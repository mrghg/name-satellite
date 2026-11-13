"""Reproducing the legacy pipeline and comparing runs: NOT needed for new runs.

The new-run workflow (``prepare`` -> NAME -> ``processing``) uses only the modules directly in
``name_satellite``. This sub-package holds what validates it (docs/workflow.md, "Validation"):

* ``legacy_files``: readers for files written by the legacy ACRG pipeline;
* ``legacy_pipeline``: the legacy satellite processing, reproduced exactly, with
  ``boundary_histograms`` and ``column`` (its steps 6b and 6c);
* ``refdata``: the registry of legacy reference data in ``data/``;
* ``compare``: measures for comparing two sets of processed footprints (totals, shapes,
  edge histograms), used against legacy and between seeds.
"""
