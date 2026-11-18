"""Portable I/O, deterministic random streams, and configuration validation."""
import hashlib
from functools import lru_cache
import json
import os
from pathlib import Path
import random

TOOLS = ("minimap2", "last", "lastz", "blastn", "fastga", "bowtie")
LENGTHS = [150, 250, 500, 1000]
IDENTITIES = list(range(50, 91, 5)) + [92, 94, 96, 98, 100]


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, ensure_ascii=False,
                                     separators=(",", ":"), sort_keys=True).encode()).hexdigest()


def random_stream(master, stage, *labels):
    payload = json.dumps([master, stage, *labels], ensure_ascii=False,
                         separators=(",", ":")).encode()
    h = hashlib.sha256(payload).hexdigest()
    return random.Random(int(h[:16], 16)), h


def read_json(path):
    return json.loads(Path(path).read_text())


@lru_cache(maxsize=128)
def _cached_json(path, size, mtime_ns):
    return read_json(path)


def read_manifest(path):
    """Read a large, immutable manifest once per file version and process."""
    path = Path(path)
    stat = path.stat()
    return _cached_json(path, stat.st_size, stat.st_mtime_ns)


def write_json(path, value):
    atomic_text(path, json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def atomic_text(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".{os.getpid()}.tmp")
    tmp.write_text(text)
    tmp.replace(path)


def rc(seq):
    return seq.translate(str.maketrans("ACGTacgtNn", "TGCAtgcaNn"))[::-1]


def read_fasta(path):
    records, name, chunks = {}, None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
                if name in records:
                    raise ValueError(f"Duplicate FASTA name: {name}")
            elif line.strip():
                if name is None:
                    raise ValueError("Sequence before FASTA header")
                chunks.append(line.strip())
    if name is not None:
        records[name] = "".join(chunks)
    return records


def write_fasta(path, records):
    atomic_text(path, "".join(">" + k + "\n" + "\n".join(v[i:i+80] for i in range(0, len(v), 80))
                              + "\n" for k, v in records.items()))


def load_config(path):
    c = read_json(path)
    c["_path"] = str(Path(path).resolve())
    if c["clustered"] != {"h_percent": 50, "a_percent": 80, "w_percent": 10}:
        raise ValueError("This protocol freezes C at h50/a80/w10")
    if c["models"] != ["U", "C", "I"] or c["contexts"] != ["short", "native"]:
        raise ValueError("Expected U/C/I and exactly short/native contexts")
    n = c["loci_per_reference_length"]
    if type(n) is not int or n < 1:
        raise ValueError("Positive distinct-locus quota required")
    if c.get("profile") == "production":
        if n != 1000 or c["lengths"] != LENGTHS or c["identities"] != IDENTITIES:
            raise ValueError("Production axes are frozen at1000 distinct loci/four lengths/14 identities")
        if {r["id"] for r in c["references"]} != {"bsub", "pombe", "dmel", "zebrafish", "mouse", "human"}:
            raise ValueError("Production requires all six specified references")
        if c["enabled_tools"] != list(TOOLS):
            raise ValueError("Production requires all six declared tools")
    if any(x not in TOOLS for x in c["enabled_tools"]):
        raise ValueError("Unknown tool, including Bowtie2: original Bowtie only")
    if len({r["id"] for r in c["references"]}) != len(c["references"]):
        raise ValueError("Duplicate reference IDs")
    if len(set(c["lengths"])) != len(c["lengths"]) or any(n < 1 for n in c["lengths"]):
        raise ValueError("Distinct positive lengths required")
    if len(set(c["identities"])) != len(c["identities"]) or any(not 0 < x <= 100 for x in c["identities"]):
        raise ValueError("Distinct identity percentages required")
    if 100 not in c['identities']:
        raise ValueError('P100 source controls are mandatory')
    b = c["background"]
    if b["selected_loci_per_reference_length"] > n or b["selected_loci_per_reference_length"] % b["targets_per_sample"]:
        raise ValueError("Background quota must fit and be divisible by K")
    if not set(b["lengths"]) <= set(c["lengths"]) or not set(b["identities"]) <= set(c["identities"]):
        raise ValueError("Background axes must be a subset of the short-query bank")
    if c["batch"]["short"] < 1 or c["batch"]["native"] < 1:
        raise ValueError("Positive batch sizes required")
    return c


def public_config(c):
    return {k: v for k, v in c.items() if not k.startswith("_")}


def root(c):
    return Path(c["work_dir"]).resolve()
