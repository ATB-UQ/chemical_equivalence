from atb_outputs.mol_data import MolData

from chemical_equivalence.calcChemEquivalency import getChemEquivGroups

with open("benzene.pdb") as fh:
    mol_data = MolData(fh.read())

equivalence_groups, _ = getChemEquivGroups(mol_data)
print(equivalence_groups)