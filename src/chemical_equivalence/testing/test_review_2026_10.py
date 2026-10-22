"""Pins and known-defect markers from the 2026-10 code review (wave 2).

Run with:  /home/atb/ATB/.venv/bin/python -m pytest src/chemical_equivalence/testing -v
"""
from os.path import dirname, join

import pytest

from atb_outputs.mol_data import MolData
from chemical_equivalence.calcChemEquivalency import getChemEquivGroups
from chemical_equivalence import NautyInterface
from chemical_equivalence.NautyInterface import DreadnautError, NautyError

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


def test_atom_ids_outside_1_to_n_raise():
    # CE-C-2 (was: dreadnaut complaining on stderr and still exiting 0). Vertices are atom
    # id - 1, so a gap in the ids is a vertex out of range, reported before nauty runs.
    mol_data = MolData(open(join(FIXTURES, 'chlorocyclohexane.pdb')).read())
    (last, atom) = max(mol_data.atoms.items())
    del mol_data.atoms[last]
    mol_data.atoms[last + 5] = dict(atom, id=last + 5)
    with pytest.raises(NautyError, match='out of range'):
        NautyInterface.calcEquivGroups(mol_data)


def test_nauty_error_is_still_a_dreadnaut_error():
    assert DreadnautError is NautyError and issubclass(NautyError, RuntimeError)


def test_oversized_graph_is_refused_before_nauty(monkeypatch):
    # The size cap replaces the dreadnaut timeout: an in-process call cannot be interrupted.
    import atb_nauty
    monkeypatch.setattr(atb_nauty, 'MAX_VERTICES', 10)
    mol_data = MolData(open(join(FIXTURES, 'chlorocyclohexane.pdb')).read())
    with pytest.raises(NautyError, match='ATB_NAUTY_MAX_VERTICES'):
        NautyInterface.calcEquivGroups(mol_data)


def test_orbit_numbering_is_by_least_member():
    assert NautyInterface.orbit_classes((0, 1, 1, 0, 4, 1)) == {0: 0, 1: 1, 2: 1, 3: 0, 4: 2, 5: 1}


@pytest.mark.parametrize(('name', 'has_mirror'), [
    ('chlorocyclohexane.pdb', True),          # Cs: the mirror maps the molecule onto itself
    ('(1R)-1-bromo-1-chloroethan-1-ol.pdb', False),
    ('(1R,3S)-1,3-dibromo-1,3-dichloropropane.pdb', True),   # meso
    ('(1S,3S)-1,3-dibromo-1,3-dichloropropane.pdb', False),
])
def test_mirror_verdict(name, has_mirror):
    getChemEquivGroups(mol_data := MolData(open(join(FIXTURES, name)).read()))
    assert mol_data.equivalence_details['has_mirror'] is has_mirror
