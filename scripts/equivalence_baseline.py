#!/usr/bin/env python
"""Gate for a change of nauty backend (or anything else under getChemEquivGroups): record
every answer chemical_equivalence gives on a fixed corpus, and diff two recordings.

    # 1. on the ATB server (read-only DB): the inputs, once
    .venv/bin/python scripts/equivalence_baseline.py extract INPUTS.json.gz [--seed S] [--sizes 10000,7000,3000]
    # 2. anywhere with no DB (helium): one recording per code version, sharded
    python scripts/equivalence_baseline.py record INPUTS.json.gz OUT_DIR --shard I/N
    # 3. compare two recordings
    python scripts/equivalence_baseline.py diff OUT_A OUT_B

The corpus: every molecule cached under core/testing/data (offline) and this package's
test molecules, plus a deterministic DB sample (`ORDER BY MD5(CONCAT(seed, id))`, so the
seed alone reproduces it) stratified by atom count, <= 50 / 51-150 / > 150, each with its
PDB (`molecule_files.pdb_normalised`) and the highest-level QM `ocoord`, read exactly as
core/testing/equivalence_sweep.py reads them.

Recorded per molecule, raw (not counts): the classes and `equivalence_details` at
`max_stereo_distance=None` and at core's DIASTEREOTOPIC_MAX_BOND_DISTANCE, or the error.
`diff` reports every molecule whose record differs in any field, and exits 1 if one does.
Which code a recording ran is in its `code` field (module paths and package versions).
"""
import argparse
import gzip
import json
import os
import sys
import time
from os.path import dirname, exists, join
from pathlib import Path

HERE = dirname(dirname(os.path.abspath(__file__)))
TEST_MOLECULES = join(HERE, 'src', 'chemical_equivalence', 'testing')
CORE_TESTING_DATA = os.environ.get('ATB_CORE_TESTING_DATA', '/home/atb/ATB/core/testing/data')
DEFAULT_SEED = 'chemical_equivalence_baseline'
STRATA = (('le50', 'atoms <= 50'), ('51to150', 'atoms BETWEEN 51 AND 150'), ('gt150', 'atoms > 150'))
#: core's DIASTEREOTOPIC_MAX_BOND_DISTANCE (core.atb.config.settings); repeated here so a
#: recording needs no core import.
DIASTEREOTOPIC_MAX_BOND_DISTANCE = 2


def highest_ocoord(qm_data):
    best = None
    for (description, level_data) in (qm_data or {}).items():
        if not isinstance(level_data, dict) or 'ocoord' not in level_data:
            continue
        try:
            level = int(str(description).replace('qm', ''))
        except ValueError:
            continue
        if best is None or level > best[0]:
            best = (level, level_data['ocoord'])
    return None if best is None else best[1]


def offline_items():
    import yaml
    loader = getattr(yaml, 'CSafeLoader', yaml.SafeLoader)
    for name in sorted(os.listdir(TEST_MOLECULES)):
        if name.endswith('.pdb'):
            yield {'source': 'test:' + name[:-4], 'pdb': open(join(TEST_MOLECULES, name)).read(), 'ocoord': None}
    if not exists(CORE_TESTING_DATA):
        return
    for molid in sorted(os.listdir(CORE_TESTING_DATA)):
        pdb_path = join(CORE_TESTING_DATA, molid, '{0}.pdb'.format(molid))
        qm_path = join(CORE_TESTING_DATA, molid, '{0}_qm_data.yml'.format(molid))
        if exists(pdb_path) and exists(qm_path):
            with open(qm_path) as handle:
                ocoord = highest_ocoord(yaml.load(handle, Loader=loader))
            yield {'source': 'core_testing:' + molid, 'pdb': open(pdb_path).read(), 'ocoord': ocoord}


def jsonable_ocoord(ocoord):
    return None if ocoord is None else {str(k): [float(x) for x in v] for (k, v) in ocoord.items()}


def extract(args) -> int:
    from atb_helpers.db import DBCon
    from atb_qm_data.interface import get_qm_data

    items = [dict(item, ocoord=jsonable_ocoord(item['ocoord'])) for item in offline_items()]
    sizes = [int(s) for s in args.sizes.split(',')]
    db = DBCon()
    cursor = db.cursor
    try:
        for ((stratum, condition), size) in zip(STRATA, sizes):
            cursor.execute('SELECT id FROM molecules WHERE error IS NULL AND {0} '
                           'ORDER BY MD5(CONCAT(%s, id)) LIMIT {1}'.format(condition, size), (args.seed,))
            molids = [row[0] for row in cursor.fetchall()]
            for start in range(0, len(molids), 500):
                chunk = molids[start:start + 500]
                cursor.execute('SELECT id, pdb_normalised FROM molecule_files WHERE id IN ({0})'.format(
                    ','.join(str(m) for m in chunk)))
                pdbs = {molid: (pdb.decode('latin-1') if isinstance(pdb, bytes) else pdb)
                        for (molid, pdb) in cursor.fetchall() if pdb is not None}
                for molid in chunk:
                    pdb = pdbs.get(molid)
                    if pdb is None:
                        continue
                    items.append({'source': 'db:{0}:{1}'.format(stratum, molid), 'pdb': pdb,
                                  'ocoord': jsonable_ocoord(highest_ocoord(get_qm_data(molid)))})
    finally:
        db.close()
    with gzip.open(args.inputs, 'wt') as fh:
        json.dump({'seed': args.seed, 'sizes': sizes, 'items': items}, fh)
    counts = {}
    for item in items:
        key = item['source'].rsplit(':', 1)[0]
        counts[key] = counts.get(key, 0) + 1
    print('wrote {0} molecules to {1}: {2}'.format(len(items), args.inputs, counts))
    return 0


def run_one(item, max_stereo_distance):
    from atb_outputs.mol_data import MolData
    from chemical_equivalence.calcChemEquivalency import getChemEquivGroups

    mol_data = MolData(item['pdb'])
    if item['ocoord']:
        ocoord = {int(k): v for (k, v) in item['ocoord'].items()}
        for (atom_id, atom) in mol_data.atoms.items():
            if atom_id in ocoord:
                atom['ocoord'] = ocoord[atom_id]
    (classes, n_runs) = getChemEquivGroups(mol_data, max_stereo_distance=max_stereo_distance)
    return {'classes': classes, 'n_runs': n_runs,
            'details': getattr(mol_data, 'equivalence_details', None)}


def canonical(value):
    """JSON round trip, so int keys and tuples compare the same before and after storage."""
    return json.loads(json.dumps(value, sort_keys=True, default=repr))


def record(args) -> int:
    import chemical_equivalence.NautyInterface as nauty_interface
    (index, count) = map(int, args.shard.split('/'))
    out_file = args.out / 'equivalence.{0}.json'.format(index)
    if out_file.exists():
        print('shard {0} already done'.format(index))
        return 0
    with gzip.open(args.inputs, 'rt') as fh:
        items = json.load(fh)['items'][index::count]
    rows = []
    started = time.time()
    for item in items:
        row = {'source': item['source']}
        for (label, distance) in (('all', None), ('near', DIASTEREOTOPIC_MAX_BOND_DISTANCE)):
            try:
                row[label] = canonical(run_one(item, distance))
            except Exception as exc:  # noqa: BLE001
                # the two backends name their failure differently; the message is kept for
                # reading but not compared
                name = 'NautyError' if type(exc).__name__ in ('DreadnautError', 'NautyError') else type(exc).__name__
                row[label] = {'error': name, 'message': str(exc)[:300]}
        rows.append(row)
    code = {'NautyInterface': nauty_interface.__file__,
            'backend': 'atb_nauty' if hasattr(nauty_interface, 'molecule_graph') else 'dreadnaut'}
    tmp = out_file.with_suffix('.tmp')
    tmp.write_text(json.dumps({'shard': args.shard, 'code': code, 'runtime_s': time.time() - started,
                               'rows': rows}))
    os.replace(tmp, out_file)
    print('shard {0}: {1} molecules, {2:.0f} s, backend {3}'.format(
        args.shard, len(rows), time.time() - started, code['backend']))
    return 0


def load(directory: Path):
    rows = {}
    codes = set()
    for path in sorted(directory.glob('equivalence.*.json')):
        data = json.loads(path.read_text())
        codes.add(data['code']['backend'])
        for row in data['rows']:
            rows[row['source']] = row
    return (rows, codes)


def diff(args) -> int:
    (a, codes_a) = load(args.a)
    (b, codes_b) = load(args.b)
    print('A: {0} molecules, backend {1}; B: {2} molecules, backend {3}'.format(
        len(a), sorted(codes_a), len(b), sorted(codes_b)))
    differences = []
    for source in sorted(set(a) | set(b)):
        if source not in a or source not in b:
            differences.append((source, 'missing from ' + ('A' if source not in a else 'B')))
            continue
        for label in ('all', 'near'):
            (x, y) = ({k: v for (k, v) in row[label].items() if k != 'message'} for row in (a[source], b[source]))
            if x != y:
                fields = sorted(k for k in set(x) | set(y) if x.get(k) != y.get(k))
                differences.append((source, '{0}: {1}'.format(label, fields)))
    errors = sum(1 for row in a.values() if 'error' in row['all'])
    print('{0} molecules compared, {1} error identically on both sides, {2} differences'.format(
        len(set(a) & set(b)), errors, len(differences)))
    for (source, what) in differences[:50]:
        print('  {0}: {1}'.format(source, what))
    return 1 if differences else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('extract')
    p.add_argument('inputs', type=Path)
    p.add_argument('--seed', default=DEFAULT_SEED)
    p.add_argument('--sizes', default='10000,7000,3000')
    p = sub.add_parser('record')
    p.add_argument('inputs', type=Path)
    p.add_argument('out', type=Path)
    p.add_argument('--shard', default='0/1')
    p = sub.add_parser('diff')
    p.add_argument('a', type=Path)
    p.add_argument('b', type=Path)
    args = parser.parse_args()
    if args.command == 'record':
        args.out.mkdir(parents=True, exist_ok=True)
    return {'extract': extract, 'record': record, 'diff': diff}[args.command](args)


if __name__ == '__main__':
    sys.exit(main())
