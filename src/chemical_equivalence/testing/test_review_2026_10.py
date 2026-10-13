"""Pins and known-defect markers from the 2026-10 code review (wave 2).

Run with:  /home/atb/ATB/.venv/bin/python -m pytest src/chemical_equivalence/testing -v
"""
import time
from os.path import dirname, join

import pytest

from atb_outputs.mol_data import MolData
from chemical_equivalence.calcChemEquivalency import getChemEquivGroups
from chemical_equivalence.config import DREADNAUT_EXECUTABLE
from chemical_equivalence import NautyInterface
from chemical_equivalence.NautyInterface import DreadnautError, _run

FIXTURES = dirname(__file__)


def _pdb_from_smiles(smiles, seed=7):
    Chem = pytest.importorskip('rdkit.Chem')
    from rdkit.Chem import AllChem
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=seed) == 0
    AllChem.MMFFOptimizeMolecule(mol)
    return Chem.MolToPDBBlock(mol)


def _methyls(smiles):
    """Methyl carbons, strict classes and MolData of a SMILES."""
    mol_data = MolData(_pdb_from_smiles(smiles))
    classes, _ = getChemEquivGroups(mol_data)
    methyls = [a for a, atom in mol_data.atoms.items()
               if atom['type'] == 'C' and sum(1 for n in atom['conn'] if mol_data.atoms[n]['type'] == 'H') == 3]
    return methyls, classes, mol_data


def test_dimethylacetamide_n_methyls_stay_merged_documented_limitation():
    # CE-C-3: an amide C-N (~1.36-1.39 A) is longer than the 0.1337 nm C=N cutoff, so it is
    # never a stereo element and the two N-methyls (diastereotopic by NMR) share a class.
    # This pins the documented limitation (README "Not detected"); changing the detector
    # must be a visible decision, so update this test together with it.
    methyls, classes, mol_data = _methyls('CC(=O)N(C)C')
    n_methyls = [m for m in methyls if any(mol_data.atoms[n]['type'] == 'N' for n in mol_data.atoms[m]['conn'])]
    assert len(n_methyls) == 2
    assert len({classes[m] for m in n_methyls}) == 1


def test_acetone_oxime_methyls_stay_merged_documented_limitation():
    # CE-C-3: C=N of an oxime is found by length, but the 2-coordinate N fails the
    # three-neighbour sp2 rule, so syn/anti methyls are not distinguished. Documented limitation.
    methyls, classes, mol_data = _methyls('CC(C)=NO')
    assert len(methyls) == 2
    assert len({classes[m] for m in methyls}) == 1


def test_run_on_malformed_input_raises():
    # CE-C-2: dreadnaut reports bad input on stderr and still exits 0; that must be an error.
    with pytest.raises(DreadnautError, match='not defined'):
        _run([DREADNAUT_EXECUTABLE], 'this is not a valid dreadnaut script ;;; %%%\n')


def test_run_raises_on_non_zero_exit_code():
    with pytest.raises(DreadnautError, match='exit code 3'):
        _run(['sh', '-c', 'exit 3'], '')


def test_run_timeout_kills_the_whole_process_group(tmp_path):
    # CE-C-2 / CE-1: a hung dreadnaut (here a shell with a background child) is killed,
    # children included, and the call raises instead of blocking forever.
    pid_file = tmp_path / 'child.pid'
    started = time.monotonic()
    with pytest.raises(DreadnautError, match='timed out'):
        _run(['sh', '-c', 'sleep 60 & echo $! > {0}; wait'.format(pid_file)], '', timeout=0.5)
    assert time.monotonic() - started < 10
    child = int(pid_file.read_text())
    deadline = time.monotonic() + 5
    while time.monotonic() < deadline:
        try:
            state = open('/proc/{0}/stat'.format(child)).read().rsplit(')', 1)[1].split()[0]
        except FileNotFoundError:
            break
        if state == 'Z':
            break
        time.sleep(0.05)
    else:
        pytest.fail('background child {0} survived the timeout'.format(child))


def test_empty_dreadnaut_output_raises(monkeypatch):
    # CE-2: no output used to come back as "" in place of the class dict.
    monkeypatch.setattr(NautyInterface, '_run', lambda *args, **kwargs: '')
    mol_data = MolData(open(join(FIXTURES, 'chlorocyclohexane.pdb')).read())
    with pytest.raises(DreadnautError, match='no output'):
        NautyInterface.calcEquivGroups(mol_data)


def test_unparseable_orbits_raise():
    # CE-C-6: the parse failure is an error carrying the text, not a print to stdout.
    with pytest.raises(DreadnautError, match='could not parse'):
        NautyInterface.nauty_equivalence_dict('cpu time = 0.00 seconds\n0 1 junk;')


def test_mirror_run_without_a_verdict_raises(monkeypatch):
    # CE-C-2: the mirror run must say "identical." or "different."; anything else used to
    # read as "no mirror image" and silently split enantiotopic atoms.
    real = NautyInterface.generate_nauty_output_from_inputstr

    def fake(text, log=None):
        return 'garbage' if '@' in text else real(text, log=log)

    monkeypatch.setattr(NautyInterface, 'generate_nauty_output_from_inputstr', fake)
    mol_data = MolData(open(join(FIXTURES, '1-chloro-1-bromopropane.pdb')).read())
    with pytest.raises(DreadnautError, match='verdict'):
        getChemEquivGroups(mol_data)
