from argparse import ArgumentParser
from typing import Optional, List, Dict, Tuple, Any
from itertools import groupby
from operator import itemgetter

from chemical_equivalence.log_helpers import print_stderr
from chemical_equivalence.NautyInterface import (calcEquivGroups, nauty_graph, generate_nauty_output_from_inputstr,
                                                 get_partition_from_nauty_output, partition_for_chemical_equivalence_dict,
                                                 pretty_equivalence_class_dict, stereo_orbits, Partition)
from chemical_equivalence.double_bond import double_bonds as detect_double_bonds
from chemical_equivalence.stereo import (tetrahedral_orientations, cis_pairs, stereo_sources, acyclic_double_bonds,
                                         bond_distances_from, merge_distant_splits, renumber)
from chemical_equivalence.helpers.types_helpers import Logger, MolData
from chemical_equivalence.helpers.atoms import EQUIVALENCE_CLASS_KEY

from atb_outputs.mol_data import MolDataFailure

# Which stereo elements are encoded. The names are the ones the old flavour-based
# detectors went by; ring inversion needs nothing of its own any more (see stereo.py).
TETRAHEDRAL_CENTRES = 'chiral_centers'
DOUBLE_BONDS = 'double_bonds'
INVERSABLE_RINGS = 'inversable_rings'
ALL_EXCEPTION_SEARCHING_KEYWORDS = [TETRAHEDRAL_CENTRES, DOUBLE_BONDS, INVERSABLE_RINGS]


def getChemEquivGroups(
    molData: MolData,
    log: Optional[Logger] = None,
    correct_symmetry: bool = True,
    other_mol_data: Optional[MolData] = None,
    exception_searching_keywords: List[str] = ALL_EXCEPTION_SEARCHING_KEYWORDS,
    max_stereo_distance: Optional[int] = None,
) -> Tuple[Dict[int, int], int]:
    """Chemical equivalence classes, ``{atom id: class}``, written to every atom as
    ``atom['equivalenceGroup']`` as well.

    ``correct_symmetry=False`` gives the constitutional classes: nauty's orbits of the
    molecular graph. With it on, atoms that are diastereotopic -- constitutionally
    equivalent but exchanged by no rotation or reflection of the molecule -- are split,
    using the tetrahedral orientations and cis/trans relations read from the coordinates
    (``ocoord`` if present, else ``coord``); see stereo.py for the rule.

    ``max_stereo_distance`` merges such splits back for atoms further than that many
    bonds from every stereo element (the effect is negligible there, and identical
    parameters on constitutionally equivalent atoms matter more); ``None`` keeps them all.
    The distance is counted between heavy atoms: a terminal atom is as far away as the
    atom it is bonded to, so a CH2 two bonds from a stereocentre keeps its hydrogens
    split at ``max_stereo_distance=2``.

    The second value is the number of extra dreadnaut runs the stereo treatment needed
    (0 when the molecule has no stereo element), kept for the old ``n_iterations`` slot.
    ``molData.equivalence_details`` records how the answer was reached.
    """
    constitutional = calcEquivGroups(molData, log)
    if not isinstance(constitutional, dict):
        return constitutional, 0
    strict = dict(constitutional)
    details: Dict[str, Any] = {
        'constitutional': dict(constitutional),
        'sources': [],
        'split_classes': [],
        'merged_classes': [],
        'has_mirror': None,
        'max_stereo_distance': max_stereo_distance,
    }
    n_runs = 0

    if correct_symmetry:
        orientations = tetrahedral_orientations(molData, log) if TETRAHEDRAL_CENTRES in exception_searching_keywords else {}
        bonds = acyclic_double_bonds(molData, detect_double_bonds(molData, log)) if DOUBLE_BONDS in exception_searching_keywords else []
        cis = cis_pairs(molData, bonds, log)
        if orientations or cis:
            strict, has_mirror = stereo_orbits(molData, orientations, cis, log)
            strict = renumber(strict)
            n_runs = 2
            details['has_mirror'] = has_mirror
            details['split_classes'] = sorted({
                constitutional[a] for a in constitutional
                if len({strict[b] for b in constitutional if constitutional[b] == constitutional[a]}) > 1
            })
        sources = stereo_sources(molData, constitutional, orientations, bonds)
        details['sources'] = sorted(sources)
        if details['split_classes']:
            distances = bond_distances_from(molData, sources)
            effective, merged = merge_distant_splits(constitutional, strict, distances, max_stereo_distance)
            details['merged_classes'] = merged
            details['distance_to_source'] = distances
        else:
            effective = strict
    else:
        effective = strict

    details['strict'] = dict(strict)
    molData.equivalence_details = details
    for (atom_id, equivalence_class) in effective.items():
        molData.atoms[atom_id][EQUIVALENCE_CLASS_KEY] = equivalence_class

    if log and correct_symmetry:
        atoms = molData.atoms
        if details['split_classes']:
            log.info('Stereo elements at {0}; diastereotopic atoms split in {1} class(es){2}'.format(
                ' '.join(atoms[a]['symbol'] for a in details['sources']),
                len(details['split_classes']),
                ', {0} merged back beyond {1} bonds'.format(len(details['merged_classes']), max_stereo_distance)
                if details['merged_classes'] else ''))
            log.debug("Equivalence groups (stereo-aware):\n{0}".format(pretty_equivalence_class_dict(molData, effective)))
        else:
            log.debug('No diastereotopic atoms.')

    if "Molecule3D_id" in list(molData.atoms.values())[0]:
        effective = {molData.atoms[a]['Molecule3D_id']: c for (a, c) in effective.items()}
    return (effective, n_runs)


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


def partial_mol_data_for_pdbstr(
    pdb_string: str,
    united_atoms: bool = True,
    debug: bool = False,
    exception_searching_keywords: List[str] = ALL_EXCEPTION_SEARCHING_KEYWORDS,
    enforce_single_molecule: bool = True,
) -> MolData:
    assert pdb_string, 'Empty PDB string'

    data = MolData(pdb_string, enforce_single_molecule=enforce_single_molecule)
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
        getChemEquivGroups(mol_data, correct_symmetry=correct_symmetry)[0]
        for mol_data in mol_datae
    ]

    return [
        get_partition_from_nauty_output(
            generate_nauty_output_from_inputstr(
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
    parser.add_argument('--max-stereo-distance', type=int, default=None,
                        help='Merge diastereotopic splits back beyond this many bonds from a stereo element')

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
                max_stereo_distance=args.max_stereo_distance,
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
