from os import environ
from os.path import exists, join

DOUBLE_BOND_LENGTH_CUTOFF = {
    frozenset(['C', 'C']): 0.1380, #nm, Source: phenix.elbow.elbow.quantum.better_bondlengths[("C", "C", 1.5)]
    frozenset(['C', 'N']): 0.1337, #nm, Source: phenix.elbow.elbow.quantum.better_bondlengths[("C", "N", 1.5)]
    frozenset(['N', 'N']): 0.1250, #nm, Source: http://www.chemikinternational.com/wp-content/uploads/2014/04/13.pdf
}

# dreadnaut (from the nauty package) is resolved from the DREADNAUT_EXECUTABLE env var,
# defaulting to ATB's vendored copy — the same pattern as the platform's other non-Python
# binaries (Jmol, InChI, fixnom, CyLib), built and placed by ATB's fetch_vendor_deps.sh.
# ATB_ROOT is itself overridable, so nothing here assumes /home/atb/ATB.
#
# Outside an ATB checkout, set DREADNAUT_EXECUTABLE (e.g. to /usr/local/bin/dreadnaut).
#
# NOTE: this assert runs at IMPORT time, so a missing binary breaks `import core.atb.atb`
# outright — every topology, not one molecule. It has done exactly that once.
ATB_ROOT = environ.get('ATB_ROOT', '/home/atb/ATB')
DREADNAUT_EXECUTABLE = environ.get('DREADNAUT_EXECUTABLE', join(ATB_ROOT, 'vendor', 'nauty', 'dreadnaut'))

assert exists(DREADNAUT_EXECUTABLE), (
    'Could not find dreadnaut executable at: "{0}". Build it with ATB\'s '
    'scripts/fetch_vendor_deps.sh, or set DREADNAUT_EXECUTABLE. The nauty package that '
    'contains it is at http://users.cecs.anu.edu.au/~bdm/nauty/.'
).format(DREADNAUT_EXECUTABLE)

# Stereo perception (stereo.py). A 4-coordinate atom carries an orientation only if the
# signed volume spanned by its four neighbours, normalised by the cube of its mean bond
# length, exceeds this: a regular tetrahedron gives ~3.1, square planar gives 0.
TETRAHEDRAL_MIN_NORMALISED_VOLUME = 0.05
# Substituent pairs across a detected double bond are cis below the first dihedral (deg)
# and trans above the second; a bond with a pair in between is twisted and ignored.
CIS_MAX_DIHEDRAL = 60.0
TRANS_MIN_DIHEDRAL = 120.0
