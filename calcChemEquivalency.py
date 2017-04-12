from argparse import ArgumentParser, Namespace
from typing import Optional, List, Dict, Tuple

from chemical_equivalence.log_helpers import print_stderr
from chemical_equivalence.NautyInterface import calcEquivGroups
from chemical_equivalence.chiral import contains_stereo_heterotopic_atoms
from chemical_equivalence.double_bond import contains_equivalence_breaking_double_bond
from chemical_equivalence.rings import contains_inversable_rings
from chemical_equivalence.helpers.types_helpers import FlavourCounter, Logger, Exception_Searching_Function, MolData
from chemical_equivalence.helpers.atoms import EQUIVALENCE_CLASS_KEY

from atb_outputs.mol_data import MolDataFailure

EXCEPTION_SEARCHING_FUNCTIONS = [
    contains_stereo_heterotopic_atoms,
    contains_equivalence_breaking_double_bond,
    contains_inversable_rings,
]

def getChemEquivGroups(molData: MolData, log: Optional[Logger] = None, correct_symmetry: bool = True, other_mol_data: Optional[MolData] = None) -> Tuple[Dict[int, int], int]:
    equivalence_dict = calcEquivGroups(molData, log)

    n_iterations = 0

    if correct_symmetry:
        # For cases with non-equivalent atoms we need to add flavour (some additional degree of freedom)
        # to distinguish them
        flavourCounter = FlavourCounter()
        while correct_chemical_equivalence_exceptions(molData, flavourCounter, log):
            n_iterations += 1
            clearEqGroupData(molData)
            equivalence_dict = calcEquivGroups(molData, log)

    return (equivalence_dict, n_iterations)

def correct_chemical_equivalence_exceptions(molData: MolData, flavourCounter: FlavourCounter, log: Logger, exception_searching_functions: List[Exception_Searching_Function] = EXCEPTION_SEARCHING_FUNCTIONS) -> bool:
    # If there is a chemical equivalence breaking groups then should_rerun = True
    should_rerun = any(
        [
            exception_searching_function(molData, flavourCounter, log)
            for exception_searching_function in exception_searching_functions
        ]
    )

    if log:
        if should_rerun:
            log.debug("Molecule contains chemical equivalence breaking groups.")
        else:
            log.debug("Molecule has NO chemical equivalence breaking groups.")

    return should_rerun

def clearEqGroupData(molData: MolData) -> None:
    molData.equivalenceGroups = {}
    for atom in list(molData.atoms.values()):
        del atom[EQUIVALENCE_CLASS_KEY]

def partial_mol_data_for_pdbstr(pdb_string: str, united_atoms: bool = True, debug: bool = False) -> MolData:
    assert pdb_string, 'Empty PDB string'

    data = MolData(pdb_string)
    getChemEquivGroups(data)
    if united_atoms:
        if debug:
            print_stderr(
                "All atoms: {0}\n".format(
                    "".join([a["type"] for a in list(data.atoms.values())]),
                ),
            )
        data.unite_atoms()
        if debug:
            print_stderr(
                "United atoms: {0}\n".format(
                    "".join([a["type"] for a in list(data.atoms.values()) if "uindex" in a]),
                ),
            )
    return data

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--pdb', type=str, help='Main PDB file.', required=True)
    parser.add_argument('--other-pdb', type=str, default=None, help='Other PDB file.')

    args = parser.parse_args()

    with open(args.pdb, "r") as fh:
        pdb_str = fh.read()

    if args.other_pdb is not None:
        with open(args.other_pdb, "r") as fh:
            other_pdb_str = fh.read()
    else:
        other_pdb_str = None

    print(
        getChemEquivGroups(
            MolData(pdb_str),
            other_mol_data=MolData(other_pdb_str) if other_pdb_str is not None else other_pdb_str,
        ),
    )
