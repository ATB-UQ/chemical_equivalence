from os.path import exists

DOUBLE_BOND_LENGTH_CUTOFF = {
    frozenset(['C', 'C']): 0.138, #nm, Source: phenix.elbow.elbow.quantum.better_bondlengths[("C", "C", 1.5)]
    frozenset(['C', 'N']): 0.134, #nm, Source: phenix.elbow.elbow.quantum.better_bondlengths[("C", "N", 1.5)]
    frozenset(['N', 'N']): 0.125, #nm, Source: http://www.chemikinternational.com/wp-content/uploads/2014/04/13.pdf
}

NAUTY_EXECUTABLE = '/usr/local/bin/dreadnaut'

assert NAUTY_EXECUTABLE[0] == '/', 'Dreadnaut executable was not an absolute path: "{0}"'.format(NAUTY_EXECUTABLE)

assert exists(NAUTY_EXECUTABLE), 'Could not find dreadnaut executable at: "{0}". Did you install nauty (http://users.cecs.anu.edu.au/~bdm/nauty/) ?'.format(NAUTY_EXECUTABLE)
