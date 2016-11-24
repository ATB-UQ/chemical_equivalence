from glob import glob
from os.path import basename, join
from logging import basicConfig, getLogger, DEBUG, StreamHandler, Formatter
from sys import stdout
from os import devnull

from chemical_equivalence.calcChemEquivalency import getChemEquivGroups
from chemical_equivalence.helpers.types_helpers import Logger, Optional, MolData
from chemical_equivalence.helpers.atoms import EQUIVALENCE_CLASS_KEY

from atb_outputs.graph import graph_img

TESTING_DIR = 'testing'

CACHE_GRAPH_POS = True
graph_for_test = {}

def run_tests(log: Optional[Logger], correct_symmetry: bool = True) -> None:
    test_pdb_files = [filepath for filepath in sorted(glob(join(TESTING_DIR, '*.pdb')))]

    for test_pdb_file in test_pdb_files:
        mol_data = MolData(
            open(test_pdb_file).read(),
            log=log,
        )

        print('Running test for pdb file: {0}'.format(test_pdb_file), file=stdout)
        print("Rings: {0}".format(mol_data.rings), file=stdout)
        getChemEquivGroups(mol_data, log=log, correct_symmetry=correct_symmetry)
        print(list(mol_data.equivalenceGroups.values()), file=stdout)
        print("\n".join([str((atom["index"], atom[EQUIVALENCE_CLASS_KEY])) for atom in list(mol_data.atoms.values())]), file=stdout)
        print('', file=stdout)

        molecule_graph, pos = graph_for_test[test_pdb_file] if (test_pdb_file in graph_for_test and CACHE_GRAPH_POS) else (None, None)

        graph_data, molecule_graph, pos = graph_img(
            mol_data,
            return_pos=True,
            pos=pos,
            molecule_graph=None,
        )

        if test_pdb_file not in graph_for_test:
            graph_for_test[test_pdb_file] = (molecule_graph, pos)

        with open(test_pdb_file.replace('.pdb', '_' + ('not_' if not correct_symmetry else '') + 'corrected' + '.' + 'svg'), 'w' + ('b' if isinstance(graph_data, bytes) else 't')) as fh:
            fh.write(graph_data)

if __name__=="__main__":
    log = getLogger()
    log.setLevel(DEBUG)

    formatter = Formatter('[%(levelname)s] - %(message)s')
    for should_correct_symmetry in [True, False]:
        if should_correct_symmetry:
            pass
        else:
            # Do not print anything
            stdout = open(devnull, 'w')

        ch = StreamHandler(stdout)
        ch.setLevel(DEBUG)

        ch.setFormatter(formatter)
        log.addHandler(ch)

        run_tests(log if should_correct_symmetry else None, correct_symmetry=should_correct_symmetry)
