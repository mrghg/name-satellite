"""Record what a running NAME process is doing, every few seconds (future method, profiling).

Where this sits: started in the background by each future-method job just before NAME
(``future.setup.add_sampler``); writes ``profile.csv`` in the run folder; read by
``future.profile_report``. Standard library only, so it runs with the job's ``python3`` without
the project environment::

    python3 -I sampler.py INPUT_FILE OUTPUT_CSV [INTERVAL_SECONDS]

It waits (up to 15 min) for a process whose command line contains a NAME executable name and
INPUT_FILE, then every INTERVAL seconds (default 5) appends one row, and stops when the process
has gone. All values come from Linux's /proc for that process:

* ``wall_s``: seconds since the first sample;
* ``cpu_s``: CPU time of the whole process (user + system, all threads), s;
* ``busy_threads``: CPU seconds used in the interval / interval length (e.g. 3.9 = four threads
  busy; 1.0 = one thread, e.g. reading met while the others wait);
* ``threads_over_half``: threads that used more than half the interval;
* ``rchar``, ``wchar``: bytes read / written through system calls (incl. the network file
  system), ``read_bytes``, ``write_bytes``: bytes the kernel fetched from or sent to storage;
* ``rss_kb``, ``hwm_kb``: resident memory now and at its peak;
* ``met_reads``: met files NAME has logged as read so far (lines "read from file" in its log),
  ``last_met_time``: the met time of the latest one (as NAME prints it).
"""

import contextlib
import os
import sys
import time
from pathlib import Path

NAME_EXECUTABLES = ("nameiii_64bit_par.exe", "nameiii_eccodes_NetCDF_64bit_par.exe")
TICKS = os.sysconf("SC_CLK_TCK")
COLUMNS = [
    "wall_s",
    "cpu_s",
    "busy_threads",
    "threads_over_half",
    "rchar",
    "wchar",
    "read_bytes",
    "write_bytes",
    "rss_kb",
    "hwm_kb",
    "met_reads",
    "last_met_time",
]


def find_process(input_file):
    """PID of the NAME process running input_file, or None (not mpirun or its launcher)."""
    for entry in os.listdir("/proc"):
        if not entry.isdigit():
            continue
        try:
            args = Path(f"/proc/{entry}/cmdline").read_bytes().split(b"\0")
        except OSError:
            continue
        words = [a.decode(errors="replace") for a in args if a]
        if not words or not any(words[0].endswith(exe) for exe in NAME_EXECUTABLES):
            continue
        if any(input_file in w for w in words[1:]):
            return int(entry)
    return None


def cpu_ticks(stat_path):
    """utime + stime (clock ticks) from a /proc .../stat file."""
    text = Path(stat_path).read_text()
    fields = text[text.rindex(")") + 2 :].split()  # skip "pid (comm)"; comm may hold spaces
    return int(fields[11]) + int(fields[12])  # utime, stime (stat fields 14 and 15)


def read_io(pid):
    values = {}
    for line in Path(f"/proc/{pid}/io").read_text().splitlines():
        key, value = line.split(":")
        values[key] = int(value)
    return values


def read_memory(pid):
    values = {}
    for line in Path(f"/proc/{pid}/status").read_text().splitlines():
        if line.startswith(("VmRSS:", "VmHWM:")):
            key, value = line.split(":")
            values[key] = int(value.split()[0])
    return values.get("VmRSS", 0), values.get("VmHWM", 0)


class LogFollower:
    """Counts met files NAME logs as read, reading only what was added since last time."""

    def __init__(self, log_path):
        self.path = Path(log_path)
        self.offset = 0
        self.count = 0
        self.last_time = ""
        self.tail = ""

    def update(self):
        if not self.path.exists():
            return
        with open(self.path, errors="replace") as f:
            f.seek(self.offset)
            new = f.read()
            self.offset = f.tell()
        text = self.tail + new
        # A log line reads "NWP met data for DD/MM/YYYY HH:MM UTC read from file"; keep the last
        # partial line for next time.
        cut = text.rfind("\n")
        complete, self.tail = (text[: cut + 1], text[cut + 1 :]) if cut >= 0 else ("", text)
        for line in complete.splitlines():
            if "read from file" in line and "met data for" in line:
                self.count += 1
                self.last_time = line.split("met data for", 1)[1].split("read from")[0].strip()


def main(input_file, out_csv, interval=5.0):
    input_file = str(Path(input_file).resolve())
    deadline = time.time() + 900
    pid = None
    while pid is None and time.time() < deadline:
        pid = find_process(input_file)
        if pid is None:
            time.sleep(1)
    if pid is None:
        print(f"sampler: no NAME process for {input_file}", file=sys.stderr)
        return 1
    log = LogFollower(Path(input_file).with_name(Path(input_file).stem + "Log.txt"))
    start = time.time()
    previous_total, previous_threads, previous_time = None, {}, start
    with open(out_csv, "w") as out:
        out.write(",".join(COLUMNS) + "\n")
        while True:
            try:
                now = time.time()
                total = cpu_ticks(f"/proc/{pid}/stat")
                threads = {}
                for tid in os.listdir(f"/proc/{pid}/task"):
                    with contextlib.suppress(OSError):  # a thread may end meanwhile
                        threads[tid] = cpu_ticks(f"/proc/{pid}/task/{tid}/stat")
                io = read_io(pid)
                rss, hwm = read_memory(pid)
            except (OSError, ValueError):
                break  # the process has ended
            log.update()
            dt = max(now - previous_time, 1e-6)
            busy = 0.0 if previous_total is None else (total - previous_total) / TICKS / dt
            used = [
                (ticks - previous_threads.get(tid, ticks)) / TICKS for tid, ticks in threads.items()
            ]
            over_half = sum(1 for seconds in used if seconds > 0.5 * dt)
            row = [f"{now - start:.1f}", f"{total / TICKS:.2f}", f"{busy:.3f}", str(over_half),
                   str(io.get("rchar", 0)), str(io.get("wchar", 0)), str(io.get("read_bytes", 0)),
                   str(io.get("write_bytes", 0)), str(rss), str(hwm), str(log.count),
                   log.last_time]  # fmt: skip
            out.write(",".join(row) + "\n")
            out.flush()
            previous_total, previous_threads, previous_time = total, threads, now
            time.sleep(interval)
    return 0


if __name__ == "__main__":
    arguments = sys.argv[1:]
    sys.exit(main(arguments[0], arguments[1], float(arguments[2]) if len(arguments) > 2 else 5.0))
