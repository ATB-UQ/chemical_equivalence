from os.path import exists

DOUBLE_BOND_LENGTH_CUTOFF = 0.138 #nm

NAUTY_EXECUTABLE = '/usr/local/bin/dreadnaut'

assert NAUTY_EXECUTABLE[0] == '/', 'Dreadnaut executable was not an absolute path: "{0}"'.format(NAUTY_EXECUTABLE)

assert exists(NAUTY_EXECUTABLE), 'Could not find dreadnaut executable at: "{0}". Did you install nauty (http://users.cecs.anu.edu.au/~bdm/nauty/) ?'.format(NAUTY_EXECUTABLE)
