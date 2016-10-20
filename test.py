from glob import glob
from os.path import basename, join
from logging import basicConfig, getLogger, DEBUG, StreamHandler, Formatter
from sys import stdout

from chemical_equivalence.molData import MolData
from chemical_equivalence.calcChemEquivalency import getChemEquivGroups

from atb_outputs.formats import graph

TESTING_DIR = 'testing'

def run_tests(log):
    test_pdb_files = [filepath for filepath in sorted(glob(join(TESTING_DIR, '*.pdb')))]

    for test_pdb_file in test_pdb_files:
        mol_data = MolData(
            open(test_pdb_file).read(),
            log=log,
        )

        print('Running test for pdb file: {0}'.format(test_pdb_file))
        print("Rings: {0}".format(mol_data.rings))
        getChemEquivGroups(mol_data, log=log)
        print(list(mol_data.equivalenceGroups.values()))
        print("\n".join([str((atom["index"], atom["equivalenceGroup"])) for atom in list(mol_data.atoms.values())]))
        print()

        for (graph_format, graph_data) in graph(mol_data):
            with open(test_pdb_file.replace('.pdb', '.' + graph_format), 'w' + ('b' if isinstance(graph_data, bytes) else 't')) as fh:
                fh.write(graph_data)

if __name__=="__main__":
    log = getLogger()
    log.setLevel(DEBUG)

    formatter = Formatter('[%(levelname)s] - %(message)s')

    ch = StreamHandler(stdout)
    ch.setLevel(DEBUG)

    ch.setFormatter(formatter)
    log.addHandler(ch)

    run_tests(log)
