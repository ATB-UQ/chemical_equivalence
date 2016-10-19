from optparse import OptionParser
from typing import Any, Optional, List

from chemical_equivalence.log_helpers import print_stderr
from chemical_equivalence.NautyInterface import NautyInterface
from chemical_equivalence.molData import MolData, MolDataFailure
from chemical_equivalence.chiral import contains_stereo_heterotopic_atoms
from chemical_equivalence.doubleBonds import contains_equivalence_breaking_double_bond
from chemical_equivalence.rings import contains_inversable_rings
from chemical_equivalence.helpers.types import FlavourCounter, Logger, Exception_Searching_Function

EXCEPTION_SEARCHING_FUNCTIONS = [
    contains_stereo_heterotopic_atoms,
    contains_equivalence_breaking_double_bond,
    contains_inversable_rings,
]

def getChemEquivGroups(molData: MolData, log: Optional[Logger] = None):
    # for cases with stereogenic centers we need to add flavour (some additional degree of freedom)
    # to distinguish stereoheterotopic atoms
    nautyInterface = NautyInterface(molData)

    nautyInterface.calcEquivGroups(log)

    flavourCounter = FlavourCounter()
    if chemicalEquivalenceExceptions(molData, flavourCounter, log):
        clearEqGroupData(molData)
        nautyInterface.calcEquivGroups(log)

    return molData.equivalenceGroups

def chemicalEquivalenceExceptions(molData: MolData, flavourCounter: FlavourCounter, log: Logger, exception_searching_functions: List[Exception_Searching_Function] = EXCEPTION_SEARCHING_FUNCTIONS):
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
        del atom["equivalenceGroup"]

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

def parseCommandline() -> None:
    # run CGP and symmetrization from commandline arguments
    parser = OptionParser(usage="python %prog options filename [options [filename]]")
    # required
    parser.add_option("-p", "--pdb", dest="pdb", help="PDB structure file", action="store", type="string")

    (opts, _) = parser.parse_args()

    if not opts.pdb or not opts.mtb:
        parser.print_help()
    else:
        # Try and open pdb file
        try:
            fh = open(opts.pdb, "r")
            pdbStr = fh.read()
            fh.close()

        except:
            print_stderr("ERROR: there was a problem with the PDB file.")
            parser.print_help()
            return

    data = MolData(pdbStr)

    nautyInterface = NautyInterface(data)

    # Run symmetrization
    nautyInterface.calcSym()
    print(nautyInterface.data.equivalenceGroups)
