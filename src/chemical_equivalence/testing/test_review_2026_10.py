"""Pins and known-defect markers from the 2026-10 code review (wave 2).

Run with:  /home/atb/ATB/.venv/bin/python -m pytest src/chemical_equivalence/testing -v
"""
from os.path import dirname, join

import pytest

from atb_outputs.mol_data import MolData
from chemical_equivalence.calcChemEquivalency import getChemEquivGroups
from chemical_equivalence.config import DREADNAUT_EXECUTABLE
from chemical_equivalence.NautyInterface import (
    _run, generate_nauty_output_from_inputstr, nauty_graph, partition_for_chemical_equivalence_dict,
)

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


def test_cross_molecule_partition_is_a_valid_nauty_colouring():
    # CE-C-1 :
    # comparing a molecule with itself must give dreadnaut no complaint on stderr.
    mol_data = MolData(open(join(FIXTURES, 'chlorocyclohexane.pdb')).read())
    classes, _ = getChemEquivGroups(mol_data)
    partition = partition_for_chemical_equivalence_dict(classes)
    text = '{0} c x @ {1} x ##'.format(
        nauty_graph(mol_data, nauty_node_partition=partition),
        nauty_graph(mol_data, nauty_node_partition=partition),
    )

    messages = []

    class Recorder:
        def debug(self, message):
            messages.append(message)

    generate_nauty_output_from_inputstr(text, log=Recorder())
    assert not messages, messages


test_cross_molecule_partition_is_a_valid_nauty_colouring = pytest.mark.xfail(
    strict=True,
    reason='CE-C-1: partition_for_chemical_equivalence_dict is 1-based but dreadnaut vertices are '
           '0-based ("illegal or repeated number"); only reachable via the retired Blind_RMSD / --other-pdbs path',
)(test_cross_molecule_partition_is_a_valid_nauty_colouring)


@pytest.mark.xfail(
    strict=True,
    reason='CE-C-2: _run ignores the dreadnaut return code and stderr, so malformed input is not an error',
)
def test_run_on_malformed_input_raises():
    with pytest.raises(Exception):
        _run([DREADNAUT_EXECUTABLE], 'this is not a valid dreadnaut script ;;; %%%\n')
