"""The future method (branch ``future-method``; docs/future/scoping.md): NOT the main workflow.

NAME back-runs without side boundaries, footprints on nested output grids that follow each release
group, and the boundary term from particles' 3D end positions. Modules:

* ``grids``: the three output grids (0.1° and 0.5° boxes around each release group, global 2.5°),
  snapped to the 2.5° lines, and their fixed-size storage blocks;
* ``setup``: prepare a month of runs (NAME inputs, Slurm jobs, profiling variants);
* ``sampler``: run beside NAME in each job, recording CPU per thread, I/O and memory (stdlib only);
* ``processing``: NAME output -> one netCDF file per month in the layout of
  docs/future/file-format.md;
* ``profile_report``: turn the samples and NAME's timers into a profile of each run.

The main workflow (``prepare`` -> NAME -> ``processing``) does not use this sub-package.
"""
