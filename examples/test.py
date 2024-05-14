import pickle

from atb_outputs.mol_data import MolData
from chemistry_data_structure.parsing.input_parsers import ATB_QMData_to_Molecule3D

from chemical_equivalence.calcChemEquivalency import getChemEquivGroups


def PDB_example():
    with open("benzene.pdb") as fh:
        mol_data = MolData(fh.read())

    equivalence_groups, _ = getChemEquivGroups(mol_data)
    print(equivalence_groups)


def Molecule3D_example():
    with open("804_b3lyp_631Gd_PCM_water_hessian.pickle", "rb") as fh:
        toluene_qm_data = pickle.load(fh)
    mol3D = ATB_QMData_to_Molecule3D(toluene_qm_data, net_charge=0, name="toluene")
    mol_data = MolData(mol3D)
    equivalence_groups, _ = getChemEquivGroups(mol_data)
    print(equivalence_groups)


if __name__ == "__main__":
    Molecule3D_example()
