from os import environ
from os.path import exists, join

DOUBLE_BOND_LENGTH_CUTOFF = {
    frozenset(['C', 'C']): 0.1380, #nm, Source: phenix.elbow.elbow.quantum.better_bondlengths[("C", "C", 1.5)]
    frozenset(['C', 'N']): 0.1337, #nm, Source: phenix.elbow.elbow.quantum.better_bondlengths[("C", "N", 1.5)]
    frozenset(['N', 'N']): 0.1250, #nm, Source: http://www.chemikinternational.com/wp-content/uploads/2014/04/13.pdf
}

# Tier 3: chemical_equivalence-only, deployment-tunable.
#
# dreadnaut is a vendor dependency, resolved like the platform's other non-Python
# binaries (Jmol, InChI, fixnom, CyLib): a default under <ATB_ROOT>/vendor/, with an
# env override for custom installs and containers. It is built and placed there by
# ATB's scripts/fetch_vendor_deps.sh. It used to be a hard-coded path into a
# /home/atb/ATB/nauty source checkout, which tied this package to one machine's layout.
ATB_ROOT = environ.get('ATB_ROOT', '/home/atb/ATB')
NAUTY_EXECUTABLE = environ.get('ATB_NAUTY_EXECUTABLE', join(ATB_ROOT, 'vendor', 'nauty', 'dreadnaut'))

assert exists(NAUTY_EXECUTABLE), 'Could not find dreadnaut executable at: "{0}". Build it with ATB\'s scripts/fetch_vendor_deps.sh, or set ATB_NAUTY_EXECUTABLE.'.format(NAUTY_EXECUTABLE)
