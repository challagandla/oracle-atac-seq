# Contributing

Contributions that improve correctness, reproducibility, documentation, test
coverage, portability, or scientific interpretation are welcome.

## Before changing the workflow

1. Open an issue or describe the biological and technical problem in the pull
   request.
2. Separate scientific behavior changes from formatting or dependency updates
   when possible.
3. Check whether a proposed default is valid across organisms, assemblies, and
   realistic bulk ATAC-seq designs.
4. Add primary citations for new methods or changed scientific assumptions.

Never commit real FASTQs, BAMs, results, logs, credentials, protected sample
metadata, local absolute paths, or study-specific analysis files. Tiny
synthetic fixtures belong under `.test/` and must have clear provenance.

## Development setup

Create or check the supported environments from the repository root:

```bash
bash setup.sh
bash setup.sh --check
```

Use a feature branch and keep commits focused. Do not edit generated Conda
environments under `.snakemake/`.

## Required checks

Run the checks relevant to the change:

```bash
# Shell and Python syntax.
bash -n setup.sh run.sh scripts/setup_envs.sh
conda run -n oracle-atac-runner python -m py_compile workflow/scripts/*.py tests/*.py

# Unit tests.
conda run -n oracle-atac-runner python -m pytest tests/ -q

# Resolve the test DAG through the supported runner.
bash run.sh --dry-run --configfile .test/test_config.yaml --workflow-profile none

# Whitespace and patch hygiene.
git diff --check
```

Changes to R scripts should also parse cleanly in the appropriate rule
environment. Changes to a rule, script, environment, or config key should be
covered by a focused unit test or the smallest possible synthetic DAG test.

Do not present a dry-run as end-to-end validation: it proves graph resolution,
not full tool execution or biological validity.

## Scientific review checklist

A scientific pull request should state:

- the question or failure mode being addressed;
- expected inputs and outputs;
- assembly and chromosome-name assumptions;
- effects on filtering, peak support, counting, or statistical interpretation;
- compatibility and migration consequences for existing configs;
- evidence from synthetic tests and, where appropriate, independently sourced
  benchmark data;
- new limitations or provenance requirements;
- primary references supporting the method.

Guard against silent success. Empty peak sets, missing samples, invalid
contrasts, reference mismatches, and unavailable required resources should fail
early with an actionable message.

## Documentation and public release hygiene

Update `README.md` for user-visible behavior and `TUTORIAL.md` when a change
affects setup, configuration, outputs, QC, interpretation, or limitations.
Use `bash run.sh` in public workflow examples.

Before opening a pull request, confirm that documentation links work, examples
use generic samples, and the patch contains no local paths or analysis
artifacts. Update `CITATION.cff` and `THIRD_PARTY_LICENSES.md` when authorship,
dependencies, databases, or terms change.

## Pull requests

Include a concise summary, scientific rationale, validation performed, and any
remaining limitations. CI must pass before merge. Advisory lint findings should
still be reviewed and either fixed or explained.
