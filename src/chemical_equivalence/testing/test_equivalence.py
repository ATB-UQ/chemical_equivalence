"""Chemical equivalence on the fixture molecules, asserted as relations between named
atoms rather than full class dictionaries, so the tests survive renumbering.

Run with:  /home/atb/ATB/.venv/bin/python -m pytest src/chemical_equivalence/testing -v
"""
from os.path import dirname, join

import pytest

from atb_outputs.mol_data import MolData
from chemical_equivalence.calcChemEquivalency import getChemEquivGroups

FIXTURES = dirname(__file__)


def classes_for(name, **kwargs):
    mol_data = MolData(open(join(FIXTURES, name + '.pdb')).read())
    (classes, _) = getChemEquivGroups(mol_data, **kwargs)
    by_symbol = {atom['symbol']: classes[atom_id] for (atom_id, atom) in mol_data.atoms.items()}
    return by_symbol, mol_data


def n_classes(by_symbol):
    return len(set(by_symbol.values()))


def equivalent(by_symbol, *symbols):
    return len({by_symbol[s] for s in symbols}) == 1


def test_benzene_has_two_classes():
    by_symbol, _ = classes_for('benzene')
    assert n_classes(by_symbol) == 2


def test_cyclohexane_ring_flip_keeps_all_hydrogens_equivalent():
    by_symbol, mol_data = classes_for('cyclohexane')
    assert n_classes(by_symbol) == 2
    hydrogens = [atom['symbol'] for atom in mol_data.atoms.values() if atom['type'] == 'H']
    assert equivalent(by_symbol, *hydrogens)


def test_chlorocyclohexane_keeps_its_mirror_plane_but_splits_cis_trans_hydrogens():
    # C2 bears Cl1; C3/C7 are ortho, C4/C6 meta, C5 para.
    strict, _ = classes_for('chlorocyclohexane')
    assert n_classes(strict) == 12
    assert equivalent(strict, 'C3', 'C7')
    assert equivalent(strict, 'C4', 'C6')
    # the two hydrogens on an ortho carbon are cis and trans to Cl: diastereotopic
    assert not equivalent(strict, 'H9', 'H10')
    # ...but each is mirror-related to one on the other ortho carbon
    assert equivalent(strict, 'H9', 'H18') or equivalent(strict, 'H9', 'H17')
    # para CH2: still diastereotopic strictly, four bonds from C2
    assert not equivalent(strict, 'H13', 'H14')

    # distances count heavy atoms: the para CH2 is three bonds from C2, the meta two
    within_3, _ = classes_for('chlorocyclohexane', max_stereo_distance=3)
    assert n_classes(within_3) == 12

    within_2, _ = classes_for('chlorocyclohexane', max_stereo_distance=2)
    assert n_classes(within_2) == 11
    assert equivalent(within_2, 'H13', 'H14')
    assert not equivalent(within_2, 'H11', 'H12')

    within_1, _ = classes_for('chlorocyclohexane', max_stereo_distance=1)
    assert n_classes(within_1) == 10
    assert equivalent(within_1, 'H11', 'H12')
    assert not equivalent(within_1, 'H9', 'H10')


def test_gem_dichlorocyclohexane_chlorines_are_equivalent_through_the_ring_flip():
    by_symbol, mol_data = classes_for('1,1-dichlorocyclohexane')
    assert n_classes(by_symbol) == 8
    chlorines = [atom['symbol'] for atom in mol_data.atoms.values() if atom['type'] == 'CL']
    assert len(chlorines) == 2 and equivalent(by_symbol, *chlorines)


def test_meso_dibromodichloropropane_is_achiral_with_diastereotopic_central_hydrogens():
    by_symbol, mol_data = classes_for('(1R,3S)-1,3-dibromo-1,3-dichloropropane')
    assert mol_data.equivalence_details['has_mirror'] is True
    assert n_classes(by_symbol) == 7
    assert equivalent(by_symbol, 'C2', 'C5')
    assert equivalent(by_symbol, 'Br3', 'Br7')
    assert not equivalent(by_symbol, 'H9', 'H10')


def test_c2_symmetric_dibromodichloropropane_has_homotopic_central_hydrogens():
    by_symbol, mol_data = classes_for('(1S,3S)-1,3-dibromo-1,3-dichloropropane')
    assert mol_data.equivalence_details['has_mirror'] is False
    assert n_classes(by_symbol) == 6
    assert equivalent(by_symbol, 'C2', 'C5')
    assert equivalent(by_symbol, 'H9', 'H10')


def test_single_stereocentre_splits_the_adjacent_ch2_but_not_the_methyl():
    by_symbol, _ = classes_for('1-chloro-1-bromopropane')
    assert n_classes(by_symbol) == 9
    assert not equivalent(by_symbol, 'H9', 'H10')
    assert equivalent(by_symbol, 'H6', 'H7', 'H8')


def test_trans_decalin_bridgeheads_split_every_ring_ch2():
    strict, mol_data = classes_for('trans-decalin')
    assert sorted(mol_data.atoms[a]['symbol'] for a in mol_data.equivalence_details['sources']) == ['C4', 'C9']
    assert n_classes(strict) == 8
    assert equivalent(strict, 'C1', 'C2', 'C6', 'C7')
    assert not equivalent(strict, 'H11', 'H12')
    within_1, _ = classes_for('trans-decalin', max_stereo_distance=1)
    assert n_classes(within_1) == 7
    assert equivalent(within_1, 'H13', 'H14')


def test_butadiene_terminal_hydrogens_are_cis_or_trans_to_the_chain():
    by_symbol, _ = classes_for('butadiene')
    assert n_classes(by_symbol) == 5
    assert equivalent(by_symbol, 'C1', 'C4')
    assert not equivalent(by_symbol, 'H1', 'H2')
    assert equivalent(by_symbol, 'H1', 'H5') or equivalent(by_symbol, 'H1', 'H6')


def test_chloroethene_ch2_hydrogens_are_diastereotopic():
    by_symbol, _ = classes_for('chloroethene')
    assert n_classes(by_symbol) == 6


def test_constitutional_classes_ignore_stereo():
    by_symbol, _ = classes_for('chlorocyclohexane', correct_symmetry=False)
    assert n_classes(by_symbol) == 9
    by_symbol, _ = classes_for('(1R,3S)-1,3-dibromo-1,3-dichloropropane', correct_symmetry=False)
    assert n_classes(by_symbol) == 6


@pytest.mark.parametrize('name', ['glucose', 'taxol', 'CNT', 'prismane'])
def test_larger_fixtures_run(name):
    real = {'glucose': 'D-(+)-Glucose'}.get(name, name)
    strict, _ = classes_for(real)
    constitutional, _ = classes_for(real, correct_symmetry=False)
    assert n_classes(strict) >= n_classes(constitutional)
