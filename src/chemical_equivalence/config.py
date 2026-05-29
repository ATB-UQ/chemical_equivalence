import os
from os.path import exists

DOUBLE_BOND_LENGTH_CUTOFF = {
    frozenset(['C', 'C']): 0.1380, #nm, Source: phenix.elbow.elbow.quantum.better_bondlengths[("C", "C", 1.5)]
    frozenset(['C', 'N']): 0.1337, #nm, Source: phenix.elbow.elbow.quantum.better_bondlengths[("C", "N", 1.5)]
    frozenset(['N', 'N']): 0.1250, #nm, Source: http://www.chemikinternational.com/wp-content/uploads/2014/04/13.pdf
}

# look for DEADNNAUT_EXECUTABLE environment variable
DREADNAUT_EXECUTABLE = '/usr/local/bin/dreadnaut' if os.getenv("DREADNAUT_EXECUTABLE") is None \
    else os.getenv("DREADNAUT_EXECUTABLE")

assert exists(DREADNAUT_EXECUTABLE), ('Could not find dreadnaut executable at: "{0}". '
                                      'You can find the nauty package which contains dreadnaut here: '
                                      'http://users.cecs.anu.edu.au/~bdm/nauty/.').format(DREADNAUT_EXECUTABLE)
