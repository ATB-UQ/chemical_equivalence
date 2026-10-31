"""getChemEquivGroups on a MolData built from a chemistry_data_structure Molecule3D (how
fdb_standard_pipeline and fdb_parser call it): the result is keyed by the Molecule3D atom names and
partitions the atoms as the pdb-built MolData does, including molecules with stereo centres."""
from os.path import dirname, join

import pytest

from atb_outputs.mol_data import MolData
from chemical_equivalence.calcChemEquivalency import getChemEquivGroups

input_parsers = pytest.importorskip('chemistry_data_structure.parsing.input_parsers')

FIXTURES = dirname(__file__)


def partition(classes):
    groups = {}
    for atom, equivalence_class in classes.items():
        groups.setdefault(equivalence_class, set()).add(atom)
    return {frozenset(group) for group in groups.values()}


@pytest.mark.parametrize('name', [
    '(1R)-1-bromo-1-chloroethan-1-ol',
    '(1R,3S)-1,3-dibromo-1,3-dichloropropane',
    '(1S,3S)-1,3-dibromo-1,3-dichloropropane',
])
def test_molecule3d_input_matches_pdb_input(name):
    pdb = open(join(FIXTURES, name + '.pdb')).read()
    pdb_mol_data = MolData(pdb)
    (pdb_classes, _) = getChemEquivGroups(pdb_mol_data)
    by_name = {pdb_mol_data.atoms[atom_id]['symbol']: c for (atom_id, c) in pdb_classes.items()}

    molecule = input_parsers.pdb_to_Molecule3D(pdb)
    (classes, _) = getChemEquivGroups(MolData(molecule))

    assert set(classes) == {atom.get_index('name') for atom in molecule.atom_objects}
    assert partition(classes) == partition(by_name)
