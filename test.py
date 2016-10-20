from glob import glob
from os.path import basename, join
from logging import basicConfig, getLogger, DEBUG

from chemical_equivalence.molData import MolData
from chemical_equivalence.calcChemEquivalency import getChemEquivGroups

from atb_outputs.formats import graph

TESTING_DIR = 'testing'

def run_tests():
    test_pdb_files = [filepath for filepath in glob(join(TESTING_DIR, '*.pdb'))]

    for test_pdb_file in test_pdb_files:
        mol_data = MolData(
            open(test_pdb_file).read(),
            log=getLogger(),
        )

        print("Rings: {0}".format(mol_data.rings))
        getChemEquivGroups(mol_data, log=getLogger())
        print(list(mol_data.equivalenceGroups.values()))
        print("\n".join([str((atom["index"], atom["equivalenceGroup"])) for atom in list(mol_data.atoms.values())]))

        for (graph_format, graph_data) in graph(mol_data):
            with open(test_pdb_file.replace('.pdb', '.' + graph_format), 'w' + ('b' if isinstance(graph_data, bytes) else 't')) as fh:
                fh.write(graph_data)

if __name__=="__main__":
    basicConfig(
        level=DEBUG,
        format='%(asctime)s - [%(levelname)s] - %(message)s  -->  (%(module)s.%(funcName)s: %(lineno)d)',
        datefmt='%d-%m-%Y %H:%M:%S',
    )
    run_tests()
