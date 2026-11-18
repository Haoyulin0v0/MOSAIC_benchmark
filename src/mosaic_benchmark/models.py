"""Explicit edit budgets and generating alignment paths, not inferred gold."""
from .common import random_stream


def mutate(source, identity, model, master_seed, source_id):
    source = source.upper()
    if not source or set(source) - set("ACGT") or model not in {"U", "C", "I"}:
        raise ValueError("ACGT source and U/C/I model required")
    rng, seed_hash = random_stream(master_seed, "mutation", source_id, len(source), model, identity, 1)
    n, p = len(source), identity / 100
    gaps = round(n * min(.04, .25*(1-p))) if model == "I" else 0
    ins, dele = gaps//2, gaps-gaps//2
    substitutions = round(n-dele-p*(n+ins))
    if not 0 <= substitutions <= n-dele:
        raise ValueError("Infeasible integer edit budget")
    deleted = set()
    attempts = 0
    while len(deleted) < dele:
        attempts += 1
        if attempts > 10000:
            raise RuntimeError("Deletion placement exhausted; seed retained")
        size = 1
        while size < 20 and rng.random() > .5:
            size += 1
        size = min(size, dele-len(deleted))
        start = rng.randrange(n-size+1)
        positions = set(range(start, start+size))
        if not positions & deleted:
            deleted.update(positions)
    available = [i for i in range(n) if i not in deleted]
    hot, windows = set(), []
    if model == "C" and substitutions:
        hot_count, width = round(n*.5), max(1, round(n*.1))
        for _ in range(10000):
            if len(hot) == hot_count:
                break
            size = min(width, hot_count-len(hot))
            start = rng.randrange(n-size+1)
            hot.update(range(start, start+size))
            windows.append([start, start+size])
        if len(hot) != hot_count:
            raise RuntimeError("Hotspot placement exhausted; seed retained")
        count = round(.8*substitutions)
        changed = set(rng.sample(sorted(hot), count) + rng.sample(
            [i for i in available if i not in hot], substitutions-count))
    else:
        changed = set(rng.sample(available, substitutions))
    inserts, left, attempts = {}, ins, 0
    while left:
        attempts += 1
        if attempts > 10000:
            raise RuntimeError("Insertion placement exhausted; seed retained")
        size = 1
        while size < 20 and rng.random() > .5:
            size += 1
        size = min(size, left)
        anchor = rng.randrange(n+1)
        if anchor not in inserts:
            inserts[anchor] = "".join(rng.choices("ACGT", k=size))
            left -= size
    qa, ta = [], []
    for i in range(n+1):
        if i in inserts:
            qa.append(inserts[i]); ta.append("-"*len(inserts[i]))
        if i == n:
            break
        if i in deleted:
            qa.append("-"); ta.append(source[i])
        else:
            qa.append(rng.choice([b for b in "ACGT" if b != source[i]]) if i in changed else source[i])
            ta.append(source[i])
    qa, ta = "".join(qa), "".join(ta)
    matches = sum(a == b for a, b in zip(qa, ta))
    assert len(qa) == len(ta) and ta.replace("-", "") == source
    assert sum(a != b and a != "-" and b != "-" for a, b in zip(qa, ta)) == substitutions
    return qa.replace("-", ""), {"qaln": qa, "taln": ta, "M": matches,
        "X": substitutions, "I": ins, "D": dele, "truth_identity": matches/len(qa),
        "seed_sha256": seed_hash, "hot_windows": windows, "hot_positions": sorted(hot),
        "changed_positions": sorted(changed)}
