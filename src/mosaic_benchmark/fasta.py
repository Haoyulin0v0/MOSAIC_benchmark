"""Seekable FASTA without loading multi-gigabase references into RAM."""
from pathlib import Path
from .common import read_json, write_json


class Fasta:
    def __init__(self, path, index_path=None):
        self.path = Path(path).resolve()
        stat = self.path.stat()
        fingerprint = [str(self.path), stat.st_size, stat.st_mtime_ns]
        if index_path and Path(index_path).exists():
            cached = read_json(index_path)
            if cached['fingerprint'] != fingerprint:
                raise RuntimeError('Indexed FASTA changed; use a new run directory')
            self.records = cached['records']
            self.handle = self.path.open('rb')
            return
        self.records = {}
        name = None
        with self.path.open("rb") as f:
            while True:
                offset = f.tell()
                line = f.readline()
                if not line:
                    break
                if line.startswith(b">"):
                    name = line[1:].split()[0].decode()
                    if name in self.records:
                        raise ValueError("Duplicate reference record: " + name)
                    self.records[name] = {"length": 0, "offset": None, "line_bases": None,
                                          "line_bytes": None, "short_line": False}
                elif line.strip():
                    if name is None:
                        raise ValueError("Missing FASTA header")
                    seq = line.rstrip(b"\r\n")
                    if b" " in seq or b"\t" in seq:
                        raise ValueError("Whitespace inside FASTA sequence")
                    r = self.records[name]
                    if r['offset'] is None:
                        r.update(offset=offset, line_bases=len(seq), line_bytes=len(line))
                    elif r['short_line']:
                        raise ValueError('Non-final short line: normalize FASTA wrapping before use')
                    elif len(seq) > r['line_bases'] or (len(seq) == r['line_bases'] and len(line) != r['line_bytes'] and line.endswith(b'\n')):
                        raise ValueError('Variable FASTA wrapping is not supported')
                    if len(seq) < r['line_bases']:
                        r['short_line'] = True
                    r["length"] += len(seq)
        if not self.records or not all(r["length"] for r in self.records.values()):
            raise ValueError("Empty FASTA record")
        self.handle = self.path.open("rb")
        if index_path:
            write_json(index_path, {'fingerprint': fingerprint, 'records': self.records})

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def close(self):
        self.handle.close()

    def fetch(self, name, start, end):
        r = self.records[name]
        if not 0 <= start <= end <= r["length"]:
            raise ValueError("FASTA interval out of bounds")
        out = []
        while start < end:
            line, delta = divmod(start, r['line_bases'])
            count = min(end-start, r['line_bases']-delta)
            self.handle.seek(r['offset'] + line*r['line_bytes'] + delta)
            out.append(self.handle.read(count).decode())
            start += count
        return "".join(out)

    def lengths(self):
        return {k: v["length"] for k, v in self.records.items()}
