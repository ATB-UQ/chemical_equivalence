from os import environ
from os.path import exists

DOUBLE_BOND_LENGTH_CUTOFF = {
    frozenset(['C', 'C']): 0.1380, #nm, Source: phenix.elbow.elbow.quantum.better_bondlengths[("C", "C", 1.5)]
    frozenset(['C', 'N']): 0.1337, #nm, Source: phenix.elbow.elbow.quantum.better_bondlengths[("C", "N", 1.5)]
    frozenset(['N', 'N']): 0.1250, #nm, Source: http://www.chemikinternational.com/wp-content/uploads/2014/04/13.pdf
}

# Tier 3: chemical_equivalence-only, deployment-tunable.
NAUTY_EXECUTABLE = environ.get('ATB_NAUTY_EXECUTABLE', '/home/atb/ATB/nauty/nauty25r9/dreadnaut')

assert exists(NAUTY_EXECUTABLE), 'Could not find dreadnaut executable at: "{0}". Did you install nauty (http://users.cecs.anu.edu.au/~bdm/nauty/) ?'.format(NAUTY_EXECUTABLE)
