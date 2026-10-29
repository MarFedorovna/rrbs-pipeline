"""Validation shared by the workflow and its tests. Standard library only."""
import csv
import re
from collections import Counter
from pathlib import Path


def load_samples(path):
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample", "run", "group", "treatment", "r1", "r2", "md5_r1", "md5_r2"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Missing sample sheet columns: {required - set(reader.fieldnames or [])}")
        rows = list(reader)
    if not rows:
        raise ValueError("Empty sample sheet")
    seen = set()
    for row in rows:
        sample = row["sample"]
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sample) or sample in seen:
            raise ValueError(f"Invalid or duplicate sample ID: {sample!r}")
        seen.add(sample)
        if row["treatment"] not in {"0", "1"}:
            raise ValueError("treatment must be 0 or 1")
        for mate in ["r1", "r2"]:
            if not row[mate]:
                raise ValueError(f"Missing {mate} for {sample}")
            checksum = row[f"md5_{mate}"]
            if checksum and not re.fullmatch(r"[a-fA-F0-9]{32}", checksum):
                raise ValueError(f"Invalid MD5 for {sample}/{mate}")
            if row[mate].startswith(("http://", "https://", "ftp://")) and not checksum:
                raise ValueError(f"Remote FASTQ requires MD5: {sample}/{mate}")
    counts = Counter(r["treatment"] for r in rows)
    if set(counts) != {"0", "1"} or min(counts.values()) < 2:
        raise ValueError("At least two independent biological samples in each of groups 0 and 1 are required")
    for treatment in counts:
        if len({r["group"] for r in rows if r["treatment"] == treatment}) != 1:
            raise ValueError("Each treatment code must correspond to exactly one group")
    if len({r["group"] for r in rows}) != 2:
        raise ValueError("Group labels must distinguish the two treatment codes")
    return {row["sample"]: row for row in rows}


def validate_config(config, samples):
    if type(config["non_directional"]) is not bool:
        raise ValueError("non_directional must be a YAML boolean")
    if config["assembly"] != "hg38":
        raise ValueError("This configuration and its annotation sources support hg38 only")
    if not 0 < config["qvalue"] < 1 or not 0 <= config["difference"] <= 100:
        raise ValueError("Invalid differential methylation thresholds")
    if not 0 < config["max_coverage_percentile"] <= 100:
        raise ValueError("Invalid max_coverage_percentile")
    counts = Counter(int(s["treatment"]) for s in samples.values())
    if not 2 <= config["min_per_group"] <= min(counts.values()):
        raise ValueError("min_per_group must be >=2 and fit the group sizes")
    if config["adjust"] not in {"BH", "qvalue", "SLIM"}:
        raise ValueError("Unsupported multiple-testing correction")
    if (config["overdispersion"], config["test"]) not in {("none", "Chisq"), ("MN", "F")}:
        raise ValueError("Use none/Chisq or MN/F for overdispersion/test")
    for key in ("min_coverage", "tile_size", "tile_step", "tile_min_cpg"):
        if int(config[key]) < 1:
            raise ValueError(f"{key} must be positive")
