# Legacy run scripts (NAME_scripts/satellite)

Investigated 2026-11-06 in `~/code/NAME_scripts` (GitHub `ACRG-Bristol/NAME_scripts`).

## Satellite_launch.scr (bash)
- `Satellite_launch.scr [-b script] [-s setup] [-f folder] [-l LevelSplit=12] [-n nthread] start end NAME_dir Site Domain`
- Loops over days; for each `*YYYYMMDD*` CSV it calls `BackRuns_satellite.scr -l LevelSplit -d MaxAge_Days`.
- Default folder name: `${Domain}_${Site}_column`.
- MaxAge = 30 days, except near met changeovers:
  - Mk6→7, 2013-04-30..05-10: 19 days + days since 2013-04-29
  - Mk8→9, 2015-08-25..30: 24 days + days since 2015-08-24

## BackRuns_satellite.scr (ksh, 960 lines)
- Input CSV columns: `ID, timestamp, lon, d_lon, lat, d_lat, alt (Pa), d_alt (Pa)`. The level is
  the last 2 characters of the ID (`PPPLL`: 3-digit point number + 2-digit level when
  `file_per_day`).
- Particles per source: 30/s if level ≤ `LevelSplit` (script default 8; the launcher passes 12),
  otherwise 3/s. `retrieval_secs=60` → 1801 or 181 particles per source (confirmed in the 2013
  particle files).
- `MaxNumParticles=7000000`, `MinBLDepth=40`, Z0=20, SyncTime default 3 min (the 2013 test input uses 00:06).
- RunDuration = MaxAge_Hours + ceil(hours from first to last point) + 1; nIntTimes = RunDuration × 60 / Sync.
- Source window: Start = t+60 s, Stop = t (backwards run). RunStartTime = first point +10 s is used only
  for Mk selection and the met run. **It does not truncate the release:** `%StartTimeOfRun%` is
  not used in the main template, and 1801 particles per ID are observed. (PR #29 changes it to +60 s.)
- Sources: `Cuboid, Lat-Long, Z-coord Pa, dH-Metres? No, dZ-Metres? No, X, Y, Z, dX, dY, dZ`,
  `INERT-TRACER 1.0 g/s`, `Max Age infinity`, `No Reflect? Yes`.
- Requests (`SourceTermAndOutputRequests_satellite.txt`): `Req_H1Z1_<n>`, Mixing Ratio (2019) or Air
  concentration (2018), `TGridEndOfRun`, `Int` over the run, Output Format `IA`, group `Fields_grid1`.
  The 2018 run also had BL-averaged `Req_H1BL_<n>` (`Fields_gridBL`).
- Particle output: `Sets of Particle/Puff Information` with Fate Info = Yes (ACRG-modified NAME prints
  `Long Lat Ht Age(hr) Id`).
- Domain: Z Max 19000 m agl in the 2013 runs (template default 20000).
- **get_Mk_Global stops at Mk10** (bug X1). Typos: `elsemkd`, `"$Grid1_nY"mkdir`.
- A **second "met" NAME run** extracts point met at each release: one HGrid/ZGrid/TGrid per release;
  Temperature, BL depth, Pressure, Wind speed/direction. Output `Met_<X>_<Y>_<ID>_C1_X1_Y1.txt.gz`.
  It needs one H-grid per source: it fails on stock NAME (MaxHGrids=70) and even on the ACRG build
  (500) for a 109-sounding day ("too many horizontal grids", Alaska Oct 2024).
- Outputs: `RESULTDIR/FolderName/{MixR_files, Particle_files, Met/Met_<RunName>, Input_files, Fields_files_BL}`,
  `time_step.txt` (= 60).

## Met declarations
- `satellite/MetDeclarations/Use_UMG_Mk5..Mk11` (no Mk12; `site/MetDeclarations` has Mk12 and UKV).
- The Mk11 "Use_" file referencing `UMG_Mk10_L59PT*pp Whole` flow domains is **intentional** (the
  official MetDefn files reuse Mk10 domain names).
- The site script's Mk boundaries: Mk10→11 at 1651536000 (2022-05-03), Mk11→12 at 1768953600 (2026-01-21).

## PRs (as of 2026-11-06)
- #29 (open): Mk11 for satellite, boundary 2022-05-05, RunStartTime +60 s; comments out `. /etc/profile` and `ulimit`.
- #23 (closed): HiTRes satellite scripts. `OutPut_TimeGrid2` (nt=RunDuration, dt=1 h, t0=StartTimeOfRun),
  hourly `Int 01:00` requests → `MixR_hourly`. Bugs: maps every post-2017 date to Mk11; MaxNumParticles=700000;
  claims "IA2 = NetCDF" but `2` is NAME II text format.
- #7 (draft): TROPOMI array-job launchers.

## Job scripts
`job_scripts/JASMIN/satellite/*.sh` use retired partitions, modules and paths.
`python-in_progress/satellite` is an incomplete Python 2 port.
