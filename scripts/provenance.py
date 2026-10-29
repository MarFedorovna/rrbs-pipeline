import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def command(args):
    try:
        process = subprocess.run(args, capture_output=True, text=True, timeout=60)
        return (process.stdout + process.stderr).strip()[:4000]
    except (OSError, subprocess.TimeoutExpired) as error:
        return str(error)


root = Path.cwd()
sources = sorted(p for base in ["workflow", "scripts", "config", "profiles"]
                 for p in (root / base).rglob("*") if p.is_file() and "__pycache__" not in p.parts)
sources += [root / "environment.yml", root / "environment.linux-64.lock.txt"]
payload = {
    "completed_at_utc": datetime.now(timezone.utc).isoformat(),
    "config": snakemake.params.cfg,
    "git_commit": command(["git", "rev-parse", "HEAD"]),
    "git_status": command(["git", "status", "--porcelain"]),
    "source_sha256": {str(p.relative_to(root)): sha256(p) for p in sources},
    "input_sha256": {str(p): sha256(p) for p in snakemake.input if Path(p).is_file()},
    "software": {name: command(args) for name, args in {
        "snakemake": ["snakemake", "--version"], "bismark": ["bismark", "--version"],
        "bowtie2": ["bowtie2", "--version"], "trim_galore": ["trim_galore", "--version"],
        "cutadapt": ["cutadapt", "--version"], "samtools": ["samtools", "--version"],
        "fastqc": ["fastqc", "--version"], "multiqc": ["multiqc", "--version"]}.items()},
    "conda_explicit": command(["conda", "list", "--explicit"]),
    "conda_packages": [json.loads(p.read_text()) for p in
                       sorted((Path(os.environ.get("CONDA_PREFIX", sys.prefix)) / "conda-meta").glob("*.json"))],
}
target = Path(snakemake.output[0])
target.parent.mkdir(parents=True, exist_ok=True)
target.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
