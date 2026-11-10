"""Several days of soundings in one NAME run (batching), and splitting the output back into days.

Where this sits (see docs/workflow.md)
--------------------------------------
    2. release plan per day  ->  **combine days**  ->  3. write NAME input  ->  4. run NAME
    5. read NAME output  ->  **split into days**  ->  6. process each day as usual

Why batch: each NAME run reads the met for its whole integration window, and a day's soundings
need the 30 days of met before them. A run per day therefore reads mostly the same met as the
day before; one run for several days reads the overlap only once. Each particle still travels
for at most the maximum travel time, so footprints are unchanged as long as the window covers
every release (``integration="all_releases"``, deviation X18).

Naming: NAME needs unique source names within a run, so the combined plan numbers the sources
consecutively (``SourceID_1 .. SourceID_N`` over all days, in day order) and keeps each source's
per-day name in ``day_source_name``. Splitting reverses this, so each day's output looks exactly
like that of a one-day run (``SourceID_1 ..`` in that day's order) and goes through the same
processing.

Limit: NAME 8.7.1 allows at most 3,000 field requests per output time
(``MaxFieldReqsPerReqInfo``, GlobalParameters.F90). Every source has one time-integrated field
request at the same time, so a batched run can hold at most 3,000 sources.

Second limit: one output group (one fields file) can hold only about 2,218 sources, because NAME
refuses output lines over 60,000 characters (MaxOutputLineLength). NAME checks this only when it
first writes, i.e. at the END of the run: the first 7-day batch (2,720 sources) failed after 8 h.
``name_input`` therefore splits the field requests into output groups of at most 2,000
(``fields_output_group``), and ``name_output.read_fields_folder`` reads them back together.
"""

import pandas as pd

MAX_SOURCES_PER_RUN = 3000  # NAME 8.7.1 MaxFieldReqsPerReqInfo


def combine_release_plans(plans_by_day):
    """Combine one-day release plans into the plan for a single NAME run.

    Args:
        plans_by_day: dict {day label (e.g. "20130101"): release plan from
            ``releases.build_release_plan``}, in the order the sources should appear.

    Returns:
        DataFrame with the same columns, plus ``day`` and ``day_source_name`` (the name in the
        one-day plan); ``source_name`` becomes ``SourceID_<n>`` numbered over the whole run.

    Raises:
        ValueError if a day's plan is not named SourceID_1..n in order, or the run would exceed
        NAME's limit of MAX_SOURCES_PER_RUN sources.
    """
    parts = []
    for day, plan in plans_by_day.items():
        expected = [f"SourceID_{n}" for n in range(1, len(plan) + 1)]
        if plan["source_name"].tolist() != expected:
            raise ValueError(f"Plan for {day} is not SourceID_1..{len(plan)} in order")
        part = plan.copy()
        part.insert(0, "day", day)
        part.insert(1, "day_source_name", part["source_name"])
        parts.append(part)
    combined = pd.concat(parts, ignore_index=True)
    if len(combined) > MAX_SOURCES_PER_RUN:
        raise ValueError(
            f"{len(combined)} sources in one run; NAME allows at most {MAX_SOURCES_PER_RUN}"
        )
    combined["source_name"] = [f"SourceID_{n}" for n in range(1, len(combined) + 1)]
    return combined


def split_batched_output(fields, particles, combined_plan, day):
    """The part of a batched run's output that belongs to one day, renamed as a one-day run.

    Args:
        fields: ``name_output.read_fields`` output of the batched run (dim ``source``).
        particles: particle DataFrame of the batched run with a 1-based ``source_number``
            (the position of the source in the input file = n in ``SourceID_<n>``).
        combined_plan: the plan from ``combine_release_plans`` used to write the run.
        day: day label, as in ``combined_plan["day"]``.

    Returns:
        (fields, particles) for that day only: sources renamed to the day's own
        ``SourceID_1..n``, and particle ``source_number`` renumbered to match.
    """
    of_day = combined_plan[combined_plan["day"] == day]
    if of_day.empty:
        raise ValueError(f"No sources for day {day} in the combined plan")

    day_fields = fields.sel(source=of_day["source_name"].to_numpy())
    day_fields = day_fields.assign_coords(source=of_day["day_source_name"].to_numpy())

    # Global source number -> number within the day (both 1-based).
    global_numbers = of_day["source_name"].str.removeprefix("SourceID_").astype(int)
    day_numbers = of_day["day_source_name"].str.removeprefix("SourceID_").astype(int)
    renumber = dict(zip(global_numbers, day_numbers, strict=True))
    day_particles = particles[particles["source_number"].isin(renumber)].copy()
    day_particles["source_number"] = day_particles["source_number"].map(renumber)
    return day_fields, day_particles.reset_index(drop=True)
