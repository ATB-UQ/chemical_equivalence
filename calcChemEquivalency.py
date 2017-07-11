from argparse import ArgumentParser, Namespace
from typing import Optional, List, Dict, Tuple, Any
from itertools import groupby
from operator import itemgetter

from chemical_equivalence.log_helpers import print_stderr
from chemical_equivalence.NautyInterface import calcEquivGroups, nauty_graph, nauty_output, get_partition_from_nauty_output, partition_for_chemical_equivalence_dict, Partition
from chemical_equivalence.chiral import contains_stereo_heterotopic_atoms
from chemical_equivalence.double_bond import contains_equivalence_breaking_double_bond
from chemical_equivalence.rings import contains_inversable_rings
from chemical_equivalence.helpers.types_helpers import FlavourCounter, Logger, Exception_Searching_Function, MolData
from chemical_equivalence.helpers.atoms import EQUIVALENCE_CLASS_KEY

from atb_outputs.mol_data import MolDataFailure

EXCEPTION_SEARCHING_FUNCTIONS = {
    'chiral_centers': contains_stereo_heterotopic_atoms,
    'double_bonds': contains_equivalence_breaking_double_bond,
    'inversable_rings': contains_inversable_rings,
}

ALL_EXCEPTION_SEARCHING_KEYWORDS = list(EXCEPTION_SEARCHING_FUNCTIONS)

def getChemEquivGroups(
    molData: MolData,
    log: Optional[Logger] = None,
    correct_symmetry: bool = True,
    other_mol_data: Optional[MolData] = None,
    exception_searching_keywords: List[str] = ALL_EXCEPTION_SEARCHING_KEYWORDS,
) -> Tuple[Dict[int, int], int]:
    equivalence_dict = calcEquivGroups(molData, log)

    n_iterations = 0

    if correct_symmetry:
        # For cases with non-equivalent atoms we need to add flavour (some additional degree of freedom)
        # to distinguish them
        flavourCounter = FlavourCounter()
        while correct_chemical_equivalence_exceptions(molData, flavourCounter, log, exception_searching_keywords=exception_searching_keywords):
            n_iterations += 1
            clearEqGroupData(molData)
            equivalence_dict = calcEquivGroups(molData, log)

    return (equivalence_dict, n_iterations)

def atomic_equivalence_dict_to_group_equivalence_dict(atomic_equivalence_dict: Dict[int, int]) -> Dict[int, List[int]]:
    on_equivalence_group_id = itemgetter(1)
    get_atom_id = itemgetter(0)

    return {
        equivalence_group_id: [
            get_atom_id(item) for item in group
        ]
        for (equivalence_group_id, group) in
        groupby(
            sorted(
                atomic_equivalence_dict.items(),
                key=on_equivalence_group_id,
            ),
            key=on_equivalence_group_id,
        )
    }

def correct_chemical_equivalence_exceptions(
    molData: MolData,
    flavourCounter: FlavourCounter,
    log: Logger,
    exception_searching_keywords: List[str] = ALL_EXCEPTION_SEARCHING_KEYWORDS,
) -> bool:
    # If there is a chemical equivalence breaking groups then should_rerun = True
    should_rerun = any(
        [
            EXCEPTION_SEARCHING_FUNCTIONS[exception_searching_keyword](molData, flavourCounter, log)
            for exception_searching_keyword in exception_searching_keywords
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

def partial_mol_data_for_pdbstr(
    pdb_string: str,
    united_atoms: bool = True,
    debug: bool = False,
    exception_searching_keywords: List[str] = ALL_EXCEPTION_SEARCHING_KEYWORDS,
) -> MolData:
    assert pdb_string, 'Empty PDB string'

    data = MolData(pdb_string)
    getChemEquivGroups(data, exception_searching_keywords=exception_searching_keywords)
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

def get_chemical_equivalence_accross(mol_datae: List[MolData], correct_symmetry: bool) -> List[Partition]:
    chemical_equivalence_dicts = [
        getChemEquivGroups(mol_data)[0]
        for mol_data in mol_datae
    ]

    return [
        get_partition_from_nauty_output(
            nauty_output(
                '{0} c x @ {1} x ##'.format(
                    nauty_graph(
                        mol_datae[0],
                        nauty_node_partition=partition_for_chemical_equivalence_dict(chemical_equivalence_dicts[0]),
                    ),
                    nauty_graph(
                        other_mol_data,
                        nauty_node_partition=partition_for_chemical_equivalence_dict(other_chemical_equivalence_dict),
                    ),
                ),
            ),
        )
        for (other_mol_data, other_chemical_equivalence_dict) in list(zip(mol_datae, chemical_equivalence_dicts))[1:]
    ]

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--pdb', type=str, help='Main PDB file.', required=True)
    parser.add_argument('--other-pdbs', nargs='*', default=[], help='Other PDB files.')
    parser.add_argument('--index-starts-at', type=int, default=1, help='')
    parser.add_argument('--disable-symmetry', action='store_true', help='Disable symmetry correction'),

    args = parser.parse_args()

    with open(args.pdb, "r") as fh:
        pdb_str = fh.read()

    other_pdb_strs = [
        open(other_pdb, "r").read()
        for other_pdb in args.other_pdbs
    ]

    if len(other_pdb_strs) == 0:
        print(
            getChemEquivGroups(
                MolData(pdb_str),
                correct_symmetry=not args.disable_symmetry,
            ),
        )
    else:
        translate_partition = lambda partition: [(x + args.index_starts_at, y + args.index_starts_at) for (x, y) in partition]

        print(
            list(
                map(
                    translate_partition,
                    get_chemical_equivalence_accross(
                        list(map(MolData, [pdb_str] + other_pdb_strs)),
                        not args.disable_symmetry,
                    ),
                ),
            )
        )
