DOUBLE_BOND_LENGTH_CUTOFF = {
    frozenset(['C', 'C']): 0.1380, #nm, Source: phenix.elbow.elbow.quantum.better_bondlengths[("C", "C", 1.5)]
    frozenset(['C', 'N']): 0.1337, #nm, Source: phenix.elbow.elbow.quantum.better_bondlengths[("C", "N", 1.5)]
    frozenset(['N', 'N']): 0.1250, #nm, Source: http://www.chemikinternational.com/wp-content/uploads/2014/04/13.pdf
}

# nauty runs in process through atb_nauty (a declared dependency), so there is no binary to
# find and no per-call timeout: a molecular graph takes well under a millisecond, and
# atb_nauty refuses graphs over ATB_NAUTY_MAX_VERTICES (default 20 000) before calling it.

# Stereo perception (stereo.py). A 4-coordinate atom carries an orientation only if the
# signed volume spanned by its four neighbours, normalised by the cube of its mean bond
# length, exceeds this: a regular tetrahedron gives ~3.1, square planar gives 0.
TETRAHEDRAL_MIN_NORMALISED_VOLUME = 0.05
# Substituent pairs across a detected double bond are cis below the first dihedral (deg)
# and trans above the second; a bond with a pair in between is twisted and ignored.
CIS_MAX_DIHEDRAL = 60.0
TRANS_MIN_DIHEDRAL = 120.0
