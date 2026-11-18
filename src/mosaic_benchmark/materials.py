"""Checksum-pinned public inputs; downloads are explicit, never mapper side effects."""
import gzip
from pathlib import Path
import shutil
import urllib.request
from .common import read_json, sha, write_json


def fetch(manifest_path, directory):
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    receipts = []
    for entry in read_json(manifest_path)['materials']:
        archive = directory/(entry['id']+'.fa.gz')
        if not archive.exists():
            partial = archive.with_suffix('.gz.partial')
            with urllib.request.urlopen(entry['url'], timeout=120) as src, partial.open('wb') as dst:
                shutil.copyfileobj(src, dst, 4*1024*1024)
            if sha(partial) != entry['compressed_sha256']:
                raise RuntimeError('Downloaded checksum mismatch; partial retained: '+str(partial))
            partial.replace(archive)
        if sha(archive) != entry['compressed_sha256']:
            raise RuntimeError('Existing source checksum mismatch: '+str(archive))
        output = directory/(entry['id']+'.fa')
        if not output.exists():
            partial = output.with_suffix('.fa.partial')
            with gzip.open(archive, 'rb') as src, partial.open('wb') as dst:
                if not entry.get('record'):
                    shutil.copyfileobj(src, dst, 4*1024*1024)
                else:
                    selected = False
                    found = 0
                    for line in src:
                        if line.startswith(b'>'):
                            selected = line[1:].split()[0].decode() == entry['record']
                            found += selected
                        if selected:
                            dst.write(line)
                    if found != 1:
                        raise RuntimeError('Expected exactly one selected chromosome record')
            if entry.get('fasta_sha256') and sha(partial) != entry['fasta_sha256']:
                raise RuntimeError('Decompressed FASTA checksum mismatch; partial retained')
            partial.replace(output)
        if entry.get('fasta_sha256') and sha(output) != entry['fasta_sha256']:
            raise RuntimeError('Existing FASTA checksum mismatch')
        receipts.append(dict(id=entry['id'], url=entry['url'], compressed_sha256=sha(archive),
                             fasta_sha256=sha(output), record=entry.get('record')))
    write_json(directory/'material_manifest.json', dict(materials=receipts))
