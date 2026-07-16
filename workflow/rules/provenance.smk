# =============================================================================
# provenance.smk — run record of the configured inputs and workflow state
# =============================================================================
import datetime
import glob
import hashlib
import importlib.metadata
import json
import os
import shutil
import stat
import subprocess

import yaml


REPOSITORY_ROOT = os.path.abspath(os.path.join(workflow.basedir, os.pardir))


def _json_safe(value):
    """Return a deterministic JSON-compatible copy of nested config values."""
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, os.PathLike):
        return os.fspath(value)
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def _canonical_json(value):
    return json.dumps(
        _json_safe(value),
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    )


def _sha256_bytes(value):
    return hashlib.sha256(value).hexdigest()


def _sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _installed_package_version(distribution):
    """Return a package version without depending on a CLI being on PATH."""
    try:
        return importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        return "unavailable"


def _git_text(*args):
    try:
        return subprocess.check_output(
            ["git", "-C", REPOSITORY_ROOT, *args],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def _tracked_worktree_digest():
    """Hash the current contents of Git-tracked files, never untracked files."""
    try:
        listing = subprocess.check_output(
            ["git", "-C", REPOSITORY_ROOT, "ls-files", "-z"],
            stderr=subprocess.DEVNULL,
        )
    except (OSError, subprocess.CalledProcessError):
        return 0, "unavailable"

    paths = [path for path in listing.split(b"\0") if path]
    digest = hashlib.sha256(b"oracle-atac-tracked-worktree-v1\0")
    for encoded_path in paths:
        relative_path = os.fsdecode(encoded_path)
        absolute_path = os.path.join(REPOSITORY_ROOT, relative_path)
        digest.update(encoded_path)
        digest.update(b"\0")
        try:
            metadata = os.lstat(absolute_path)
        except FileNotFoundError:
            digest.update(b"missing\0")
            continue

        executable = bool(metadata.st_mode & stat.S_IXUSR)
        digest.update(b"executable\0" if executable else b"regular\0")
        if stat.S_ISLNK(metadata.st_mode):
            content = os.fsencode(os.readlink(absolute_path))
            digest.update(b"symlink\0")
            digest.update(hashlib.sha256(content).digest())
        elif stat.S_ISREG(metadata.st_mode):
            digest.update(b"file\0")
            digest.update(bytes.fromhex(_sha256_file(absolute_path)))
        else:
            # A tracked directory is normally a submodule. Its checked-out HEAD
            # is the reproducible source identity relevant to this worktree.
            nested_head = "unavailable"
            if stat.S_ISDIR(metadata.st_mode):
                try:
                    nested_head = subprocess.check_output(
                        ["git", "-C", absolute_path, "rev-parse", "HEAD"],
                        text=True,
                        stderr=subprocess.DEVNULL,
                    ).strip()
                except (OSError, subprocess.CalledProcessError):
                    pass
            digest.update(b"other\0")
            digest.update(nested_head.encode("utf-8"))
        digest.update(b"\0")
    return len(paths), digest.hexdigest()


def _workflow_config_files():
    """Return base config followed by CLI overlays, de-duplicated in order."""
    candidates = ["config/config.yaml"]
    overwrite = getattr(workflow, "overwrite_configfiles", None)
    if overwrite:
        candidates.extend(overwrite)
    else:
        candidates.extend(getattr(workflow, "configfiles", []))

    paths = []
    seen = set()
    for candidate in candidates:
        path = os.fspath(candidate)
        absolute = os.path.abspath(path)
        if absolute in seen or not os.path.isfile(absolute):
            continue
        seen.add(absolute)
        paths.append(path)
    return paths


def _display_path(path):
    absolute = os.path.abspath(os.fspath(path))
    try:
        relative = os.path.relpath(absolute, REPOSITORY_ROOT)
    except ValueError:
        return absolute
    if relative != os.pardir and not relative.startswith(os.pardir + os.sep):
        return relative
    return absolute


def _environment_specs():
    paths = sorted(glob.glob("workflow/envs/*.yaml"))
    if os.path.isfile("environment.runner.yml"):
        paths.insert(0, "environment.runner.yml")
    return paths


PROVENANCE_CONFIG_FILES = _workflow_config_files()
PROVENANCE_ENVIRONMENT_SPECS = _environment_specs()
_tracked_file_count, _tracked_files_sha256 = _tracked_worktree_digest()
_tracked_status = _git_text("status", "--porcelain=v1", "--untracked-files=no")
PROVENANCE_SOURCE_STATE = {
    "commit": _git_text("rev-parse", "HEAD") or "unavailable",
    "dirty": None if _tracked_status is None else bool(_tracked_status),
    "tracked_file_count": _tracked_file_count,
    "tracked_files_sha256": _tracked_files_sha256,
}
PROVENANCE_SOURCE_STATE["state_sha256"] = _sha256_bytes(
    _canonical_json(PROVENANCE_SOURCE_STATE).encode("utf-8")
)

PROVENANCE_EFFECTIVE_CONFIG = _json_safe(dict(config))
PROVENANCE_PARAMETER_SNAPSHOT = {
    "effective_config": PROVENANCE_EFFECTIVE_CONFIG,
    "effective_config_sha256": _sha256_bytes(
        _canonical_json(PROVENANCE_EFFECTIVE_CONFIG).encode("utf-8")
    ),
    "selected_samples": list(SAMPLES),
    "conditions": list(conditions()),
    "resolved_paths": {
        "results_dir": RESULTS,
        "raw_dir": RAW,
        "processed_dir": PROCESSED,
        "logs_dir": LOGS,
        "reference_dir": REF,
        "genome_fasta": genome_fasta(),
        "genome_gtf": genome_gtf(),
        "blacklist_bed": blacklist_bed(),
        "footprinting_motif_db": (
            config.get("footprinting", {}).get("motif_db", "")
            if config.get("footprinting", {}).get("enabled", False) else ""
        ),
        "samples": config["samples"],
    },
    "workflow_config_files": [
        _display_path(path) for path in PROVENANCE_CONFIG_FILES
    ],
    "workflow_source_state": PROVENANCE_SOURCE_STATE,
}
PROVENANCE_PARAMETER_SNAPSHOT_JSON = _canonical_json(PROVENANCE_PARAMETER_SNAPSHOT)
PROVENANCE_PARAMETER_SNAPSHOT_SHA256 = _sha256_bytes(
    PROVENANCE_PARAMETER_SNAPSHOT_JSON.encode("utf-8")
)


def _provenance_inputs(wildcards):
    paths = [config["samples"], genome_fasta(), genome_gtf()]
    if blacklist_bed():
        paths.append(blacklist_bed())
    motif_db = config.get("footprinting", {}).get("motif_db", "")
    if config.get("footprinting", {}).get("enabled", False) and motif_db:
        paths.append(motif_db)
    paths.extend(PROVENANCE_CONFIG_FILES)
    paths.extend(PROVENANCE_ENVIRONMENT_SPECS)

    unique = []
    seen = set()
    for path in paths:
        key = os.path.abspath(os.fspath(path))
        if key not in seen:
            seen.add(key)
            unique.append(path)
    return unique


def _provenance_raw_inputs(wildcards):
    """The exact selected FASTQ bytes, ordered sample then mate."""
    return [path for sample in SAMPLES for path in raw_fastqs(sample)]


rule provenance:
    input:
        files=_provenance_inputs,
        raw=_provenance_raw_inputs,
    output:
        config=f"{RESULTS}/provenance/effective_config.yaml",
        samples=f"{RESULTS}/provenance/samples.tsv",
        manifest=f"{RESULTS}/provenance/run_manifest.json",
        envs=f"{RESULTS}/provenance/software_environments.sha256.tsv",
        raw=f"{RESULTS}/provenance/raw_inputs.sha256.tsv",
    params:
        # Snakemake's default `params` rerun trigger makes both an effective
        # config change and a tracked source-state change rebuild this record.
        effective_config_sha256=PROVENANCE_PARAMETER_SNAPSHOT[
            "effective_config_sha256"
        ],
        tracked_source_state_sha256=PROVENANCE_SOURCE_STATE["state_sha256"],
    run:
        outdir = os.path.dirname(output.manifest)
        os.makedirs(outdir, exist_ok=True)
        snapshot = json.loads(PROVENANCE_PARAMETER_SNAPSHOT_JSON)

        shutil.copyfile(config["samples"], output.samples)
        with open(output.config, "w") as handle:
            yaml.safe_dump(snapshot["effective_config"], handle, sort_keys=True)

        with open(output.envs, "w") as handle:
            handle.write("path\tsha256\n")
            for path in PROVENANCE_ENVIRONMENT_SPECS:
                handle.write(f"{_display_path(path)}\t{_sha256_file(path)}\n")

        raw_records = []
        raw_index = 0
        for sample in SAMPLES:
            for mate in ("R1", "R2"):
                path = str(input.raw[raw_index])
                raw_index += 1
                raw_records.append({
                    "sample": sample,
                    "mate": mate,
                    "path": _display_path(path),
                    "bytes": os.path.getsize(path),
                    "sha256": _sha256_file(path),
                })
        with open(output.raw, "w") as handle:
            handle.write("sample\tmate\tpath\tbytes\tsha256\n")
            for record in raw_records:
                handle.write(
                    f"{record['sample']}\t{record['mate']}\t{record['path']}\t"
                    f"{record['bytes']}\t{record['sha256']}\n"
                )

        references = {}
        reference_inputs = [genome_fasta(), genome_gtf(), blacklist_bed()]
        motif_db = config.get("footprinting", {}).get("motif_db", "")
        if config.get("footprinting", {}).get("enabled", False):
            reference_inputs.append(motif_db)
        for path in reference_inputs:
            if path and os.path.isfile(path):
                references[_display_path(path)] = _sha256_file(path)

        config_files = {}
        for path in PROVENANCE_CONFIG_FILES:
            config_files[_display_path(path)] = _sha256_file(path)

        source_state = snapshot["workflow_source_state"]
        payload = {
            "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "workflow_commit": source_state["commit"],
            "workflow_dirty": source_state["dirty"],
            "tracked_source_state_sha256": source_state["state_sha256"],
            "snakemake_version": _installed_package_version("snakemake"),
            "parameter_snapshot_sha256": PROVENANCE_PARAMETER_SNAPSHOT_SHA256,
            "parameter_snapshot": snapshot,
            "configuration_file_sha256": config_files,
            "raw_inputs": raw_records,
            "reference_sha256": references,
            "environment_spec_sha256": {
                _display_path(path): _sha256_file(path)
                for path in PROVENANCE_ENVIRONMENT_SPECS
            },
        }
        with open(output.manifest, "w") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
