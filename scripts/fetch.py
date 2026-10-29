#!/usr/bin/env python3
"""Atomic download/local copy with checksum validation and bounded retries."""
import argparse
import hashlib
import shutil
import time
import urllib.request
from pathlib import Path


def digest(path, algorithm="md5"):
    h = hashlib.new(algorithm)
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def fetch(source, destination, md5=""):
    target = Path(destination)
    target.parent.mkdir(parents=True, exist_ok=True)
    remote = source.startswith(("https://", "http://", "ftp://"))
    if remote and target.is_file() and md5 and digest(target) == md5.lower():
        return
    temporary = target.with_name(target.name + ".part")
    for attempt in range(3):
        try:
            if remote:
                with urllib.request.urlopen(source, timeout=120) as src, open(temporary, "wb") as out:
                    shutil.copyfileobj(src, out, 4 * 1024 * 1024)
            else:
                shutil.copyfile(source, temporary)
            if md5 and digest(temporary) != md5.lower():
                raise ValueError(f"MD5 mismatch: {source}")
            if temporary.stat().st_size == 0:
                raise ValueError(f"Empty input: {source}")
            temporary.replace(target)
            return
        except Exception:
            temporary.unlink(missing_ok=True)
            if attempt == 2:
                raise
            time.sleep(2 ** attempt)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("source")
    p.add_argument("destination")
    p.add_argument("--md5", default="")
    a = p.parse_args()
    fetch(a.source, a.destination, a.md5)
