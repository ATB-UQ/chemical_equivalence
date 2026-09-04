"""Invariance properties of the stereo-aware equivalence, on molecules embedded by RDKit.
RDKit is used only to make fixtures: a SMILES becomes a 3D structure with CONECT records,
which is what the ATB receives."""
from itertools import permutations
from random import Random

import pytest

from atb_outputs.mol_data import MolData
from chemical_equivalence.calcChemEquivalency import getChemEquivGroups

rdkit = pytest.importorskip('rdkit')
from rdkit import Chem  # noqa: E402
from rdkit.Chem import AllChem  # noqa: E402


def pdb_from_smiles(smiles, seed=7, mirror=False, shuffle=None):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=seed) == 0
    AllChem.MMFFOptimizeMolecule(mol)
    conformer = mol.GetConformer()
    if mirror:
        for i in range(mol.GetNumAtoms()):
            p = conformer.GetAtomPosition(i)
            conformer.SetAtomPosition(i, (-p.x, p.y, p.z))
    order = list(range(mol.GetNumAtoms()))
    if shuffle is not None:
        Random(shuffle).shuffle(order)
        mol = Chem.RenumberAtoms(mol, order)
    # PDB serial (1-based) of the atom that was RDKit atom i before any renumbering
    serial_of_original = {original: position + 1 for (position, original) in enumerate(order)}
    return Chem.MolToPDBBlock(mol), serial_of_original


def partition(pdb, **kwargs):
    mol_data = MolData(pdb)
    classes, _ = getChemEquivGroups(mol_data, **kwargs)
    return frozenset(frozenset(a for a in classes if classes[a] == c) for c in set(classes.values())), mol_data


def in_original_labels(part, serial_of_original):
    original_of_serial = {serial: original for (original, serial) in serial_of_original.items()}
    return frozenset(frozenset(original_of_serial[a] for a in block) for block in part)


CASES = {
    'chiral_single': 'Cl[C@H](Br)CC',
    'meso': 'Cl[C@H](Br)C[C@H](Br)Cl',
    'c2': 'Cl[C@H](Br)C[C@@H](Br)Cl',
    'chlorocyclohexane': 'ClC1CCCCC1',
    'isopropyl_next_to_centre': 'C[C@H](O)C(C)C',
    'alkene_ez': 'C/C=C/CC',
    'aromatic': 'c1ccc(cc1)C(C)C',
}


@pytest.mark.parametrize('name', sorted(CASES))
def test_partition_is_invariant_to_mirror_image(name):
    pdb, _ = pdb_from_smiles(CASES[name])
    mirrored, _ = pdb_from_smiles(CASES[name], mirror=True)
    (a, _), (b, _) = partition(pdb), partition(mirrored)
    assert a == b


@pytest.mark.parametrize('name', sorted(CASES))
def test_partition_is_invariant_to_atom_relabelling(name):
    pdb, identity = pdb_from_smiles(CASES[name])
    reordered, serial_of_original = pdb_from_smiles(CASES[name], shuffle=3)
    (a, _), (b, _) = partition(pdb), partition(reordered)
    assert in_original_labels(a, identity) == in_original_labels(b, serial_of_original)


@pytest.mark.parametrize('name', sorted(CASES))
def test_strict_refines_constitutional_and_no_cutoff_equals_strict(name):
    pdb, _ = pdb_from_smiles(CASES[name])
    strict, _ = partition(pdb)
    constitutional, _ = partition(pdb, correct_symmetry=False)
    unlimited, _ = partition(pdb, max_stereo_distance=None)
    assert strict == unlimited
    for block in strict:
        assert any(block <= big for big in constitutional)


def test_meso_and_c2_diastereomers_differ_only_in_the_central_ch2():
    meso, m_meso = partition(pdb_from_smiles(CASES['meso'])[0])
    c2, m_c2 = partition(pdb_from_smiles(CASES['c2'])[0])
    assert len(meso) == len(c2) + 1
    assert m_meso.equivalence_details['has_mirror'] is True
    assert m_c2.equivalence_details['has_mirror'] is False


def test_isopropyl_methyls_next_to_a_stereocentre_split_within_two_bonds():
    pdb, _ = pdb_from_smiles(CASES['isopropyl_next_to_centre'])
    strict, mol_data = partition(pdb)
    within_2, _ = partition(pdb, max_stereo_distance=2)
    within_1, _ = partition(pdb, max_stereo_distance=1)
    # the two methyl carbons of the isopropyl group are two bonds from the stereocentre
    assert len(strict) == len(within_2)
    assert len(within_1) < len(strict)
