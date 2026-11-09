"""Tests for name_satellite.lotus (Slurm job scripts)."""

import pytest

from name_satellite.lotus import write_job_script


def test_threads_tie_cpus_and_openmp(tmp_path):
    path = write_job_script(tmp_path / "job.sbatch", tmp_path, tmp_path / "input.txt", threads=4)
    text = path.read_text()
    assert "#SBATCH --cpus-per-task=4" in text
    assert "export OMP_NUM_THREADS=4" in text
    assert "#SBATCH --account=acrg" in text
    assert "#SBATCH --partition=debug" in text
    assert "shumlib_2024.11.1" in text


def test_high_qos_uses_standard_partition(tmp_path):
    text = write_job_script(tmp_path / "j", tmp_path, "in.txt", threads=8, qos="high").read_text()
    assert "#SBATCH --partition=standard" in text
    assert "#SBATCH --qos=high" in text


def test_too_many_threads_for_qos(tmp_path):
    with pytest.raises(ValueError, match="at most 1 CPU"):
        write_job_script(tmp_path / "j", tmp_path, "in.txt", threads=4, qos="standard")
    with pytest.raises(ValueError, match="at most 8 CPU"):
        write_job_script(tmp_path / "j", tmp_path, "in.txt", threads=12, qos="debug")
