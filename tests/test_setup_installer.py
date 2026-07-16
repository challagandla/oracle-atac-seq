import os
import subprocess
from pathlib import Path

import pytest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SETUP_SCRIPT = REPOSITORY_ROOT / "setup.sh"
RUN_SCRIPT = REPOSITORY_ROOT / "run.sh"
ENV_SNAKEFILE = REPOSITORY_ROOT / "workflow" / "envs.smk"


def miniforge_asset(system, architecture, version="26.3.2-2"):
    environment = os.environ.copy()
    environment["MINIFORGE_VERSION"] = version
    result = subprocess.run(
        [
            "bash",
            "-c",
            'source "$1"; miniforge_asset_name "$2" "$3"',
            "asset-test",
            str(SETUP_SCRIPT),
            system,
            architecture,
        ],
        check=True,
        capture_output=True,
        env=environment,
        text=True,
    )
    return result.stdout.strip()


@pytest.mark.parametrize(
    ("system", "architecture", "expected"),
    [
        ("Linux", "x86_64", "Miniforge3-26.3.2-2-Linux-x86_64.sh"),
        ("Linux", "aarch64", "Miniforge3-26.3.2-2-Linux-aarch64.sh"),
        ("Linux", "arm64", "Miniforge3-26.3.2-2-Linux-aarch64.sh"),
        ("Linux", "ppc64le", "Miniforge3-26.3.2-2-Linux-ppc64le.sh"),
        ("Darwin", "x86_64", "Miniforge3-26.3.2-2-MacOSX-x86_64.sh"),
        ("Darwin", "arm64", "Miniforge3-26.3.2-2-MacOSX-arm64.sh"),
        ("Darwin", "aarch64", "Miniforge3-26.3.2-2-MacOSX-arm64.sh"),
    ],
)
def test_miniforge_asset_name_is_versioned(system, architecture, expected):
    assert miniforge_asset(system, architecture) == expected


def test_miniforge_asset_name_uses_requested_version():
    assert (
        miniforge_asset("Linux", "x86_64", version="99.1-0")
        == "Miniforge3-99.1-0-Linux-x86_64.sh"
    )


@pytest.mark.parametrize(
    ("system", "architecture"),
    [("FreeBSD", "x86_64"), ("Darwin", "ppc64le"), ("Linux", "s390x")],
)
def test_miniforge_asset_name_rejects_unsupported_platforms(system, architecture):
    environment = os.environ.copy()
    environment["MINIFORGE_VERSION"] = "26.3.2-2"
    result = subprocess.run(
        [
            "bash",
            "-c",
            'source "$1"; miniforge_asset_name "$2" "$3"',
            "asset-test",
            str(SETUP_SCRIPT),
            system,
            architecture,
        ],
        check=False,
        capture_output=True,
        env=environment,
        text=True,
    )

    assert result.returncode != 0
    assert "[ERROR]" in result.stderr


@pytest.mark.parametrize("script", [SETUP_SCRIPT, RUN_SCRIPT])
def test_custom_miniforge_home_is_discovered(script, tmp_path):
    conda = tmp_path / "custom-miniforge" / "bin" / "conda"
    conda.parent.mkdir(parents=True)
    conda.write_text("#!/usr/bin/env bash\nexit 0\n")
    conda.chmod(0o755)

    environment = os.environ.copy()
    environment.pop("CONDA_EXE", None)
    environment["MINIFORGE_HOME"] = str(conda.parents[1])
    result = subprocess.run(
        [
            "bash",
            "-c",
            'source "$1"; locate_conda; printf "%s" "$CONDA_BIN"',
            "discovery-test",
            str(script),
        ],
        check=True,
        capture_output=True,
        env=environment,
        text=True,
    )

    assert result.stdout == str(conda)


def test_solver_fallback_stays_with_the_discovered_conda(tmp_path):
    conda = tmp_path / "selected" / "bin" / "conda"
    foreign_mamba = tmp_path / "foreign" / "bin" / "mamba"
    for executable in (conda, foreign_mamba):
        executable.parent.mkdir(parents=True)
        executable.write_text("#!/usr/bin/env bash\nexit 0\n")
        executable.chmod(0o755)

    environment = os.environ.copy()
    environment["PATH"] = f"{foreign_mamba.parent}:{environment['PATH']}"
    result = subprocess.run(
        [
            "bash",
            "-c",
            'source "$1"; CONDA_BIN="$2"; choose_solver; printf "%s" "$SOLVER_BIN"',
            "solver-test",
            str(SETUP_SCRIPT),
            str(conda),
        ],
        check=True,
        capture_output=True,
        env=environment,
        text=True,
    )

    assert result.stdout == str(conda)


def test_solver_prefers_mamba_from_the_discovered_conda(tmp_path):
    conda = tmp_path / "selected" / "bin" / "conda"
    sibling_mamba = conda.with_name("mamba")
    for executable in (conda, sibling_mamba):
        executable.parent.mkdir(parents=True, exist_ok=True)
        executable.write_text("#!/usr/bin/env bash\nexit 0\n")
        executable.chmod(0o755)

    result = subprocess.run(
        [
            "bash",
            "-c",
            'source "$1"; CONDA_BIN="$2"; choose_solver; printf "%s" "$SOLVER_BIN"',
            "solver-test",
            str(SETUP_SCRIPT),
            str(conda),
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    assert result.stdout == str(sibling_mamba)


def test_solver_calls_ignore_a_foreign_mamba_root_prefix(tmp_path):
    foreign_root = tmp_path / "foreign-micromamba"
    environment = os.environ.copy()
    environment["MAMBA_ROOT_PREFIX"] = str(foreign_root)

    result = subprocess.run(
        [
            "bash",
            "-c",
            'source "$1"; clean_run bash -c \'printf "%s" "${MAMBA_ROOT_PREFIX-}"\'',
            "solver-root-test",
            str(SETUP_SCRIPT),
        ],
        check=True,
        capture_output=True,
        env=environment,
        text=True,
    )

    assert result.stdout == ""


@pytest.mark.parametrize("script", [SETUP_SCRIPT, RUN_SCRIPT])
def test_runner_provides_a_writable_matplotlib_cache(script):
    source = script.read_text()
    assert 'export MPLCONFIGDIR="${MPLCONFIGDIR:-$XDG_CACHE_HOME/matplotlib}"' in source
    assert 'mkdir -p "$XDG_CACHE_HOME" "$MPLCONFIGDIR"' in source


def test_beginner_entrypoints_promote_outcome_blind_qc_review():
    setup = SETUP_SCRIPT.read_text()
    runner = RUN_SCRIPT.read_text()

    assert "qc_review" in setup
    assert "qc_review" in runner
    assert "results/qc/multiqc_report.html" not in runner


def test_runner_restarts_incomplete_outputs_by_default():
    runner = RUN_SCRIPT.read_text()

    assert '"${DEPLOYMENT_ARGS[@]}" --rerun-incomplete --printshellcmds' in runner


def test_r_environment_smoke_check_is_genome_neutral():
    env_workflow = ENV_SNAKEFILE.read_text()
    check_r = env_workflow.split("rule check_r_env:", 1)[1].split(
        "rule check_sra_env:", 1
    )[0]

    assert "org.Hs.eg.db" not in check_r
    assert "library(enrichplot)" in check_r
    assert 'export R_LIBS_USER="$CONDA_PREFIX/lib/R/library"' in check_r
