[![DOI](https://zenodo.org/badge/148241360.svg)](https://zenodo.org/badge/latestdoi/148241360)

# Requirements

* Python `>=3.5`

* ATB Outputs python module: `https://github.com/ATB-UQ/atb_outputs.git

* Nauty: This modules relies on the `dreadnaut` executable which is part of the `nauty` package.
Nauty can very easily be installed from [source](http://users.cecs.anu.edu.au/~bdm/nauty/) or using a package manager such as `homebrew` for Mac OS X users.

# Configuration

* `dreadnaut` is resolved from the `DREADNAUT_EXECUTABLE` environment variable, defaulting to the ATB platform's vendored copy (`<ATB_ROOT>/vendor/nauty/dreadnaut`, built by `scripts/fetch_vendor_deps.sh`). The module asserts it exists at import time.

# What "chemically equivalent" means here

Two atoms are equivalent when some rotation or reflection of the molecule, in any accessible conformer, exchanges them. That is more than the graph alone can say: nauty's orbits of the molecular graph give *constitutional* equivalence, under which the two hydrogens of a CH2 next to a stereocentre are equivalent although they are diastereotopic. `stereo.py` reads the tetrahedral orientation of every 4-coordinate atom and the cis/trans relations across acyclic double bonds from the coordinates (`ocoord` if present, else `coord`) and encodes them as gadgets on the graph, so that the surviving automorphisms are exactly the label permutations realisable by a rotation (all orientations kept) or a reflection (all flipped). Enantiotopic atoms -- meso compounds, chlorocyclohexane's two ortho carbons -- stay equivalent; diastereotopic ones split; ring inversion and bond rotation need no special treatment because tetrahedral orientation is a configurational invariant.

`getChemEquivGroups(..., max_stereo_distance=N)` merges such splits back for atoms more than `N` heavy-atom bonds from every stereo element, where the effect is negligible and identical parameters on equivalent atoms matter more; `None` keeps them all. `mol_data.equivalence_details` records the constitutional, strict and effective classes and the stereo elements found.

Not detected: axial chirality (allenes, atropisomers, spiro compounds) and pyramidal three-coordinate centres (sulfoxides, phosphines).

# Usage

A simple example is included in `examples/test.py`

The tests live in `src/chemical_equivalence/testing/` (fixture molecules as PDB files, relation-style expectations, and invariance properties on RDKit-embedded structures):

```
/home/atb/ATB/.venv/bin/python -m pytest src/chemical_equivalence/testing -v
```

 * First, a `MolData` (Molecule Data) object has to be initialised. It describes the topology and coordinate of the molecule.
Please see the documentation of the `atb_outputs` module for further description of the `MolData` object.

```
>>> from chemical_equivalence.helpers.types_helpers import MolData
>>> mol_data = MolData(open('testing/chlorocyclohexane.pdb').read())
```

* The chemical equivalence can be run in a single line of code:

```
>>> from chemical_equivalence.calcChemEquivalency import getChemEquivGroups
>>> equivalence_dict, n_iterations = getChemEquivGroups(mol_data)
```

* It returns a dictionary mapping each atom index (`atom['id']`) to an `int` per equivalence class (chlorocyclohexane: the ortho carbons C3/C7 and meta carbons C4/C6 are mirror images of each other and share a class; the two hydrogens on each CH2 are cis and trans to the chlorine and do not):

```
>>> equivalence_dict
{1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5, 7: 3, 8: 5, 9: 4, 10: 6, 11: 7, 12: 8, 13: 9, 14: 10, 15: 11, 16: 6, 17: 7, 18: 8}
```

* It can easily be mapped back to each atom by iterating over the atoms of the `MolData` object:

```
>>> list((atom['symbol'], equivalence_dict[atom['id']]) for atom in mol_data.atoms.values())
[('H8', 0), ('C2', 1), ('Cl1', 2), ('C7', 3), ('H17', 4), ('H18', 5), ('C3', 3), ('H9', 5), ('H10', 4), ('C4', 6), ('H11', 7), ('H12', 8), ('C5', 9), ('H13', 10), ('H14', 11), ('C6', 6), ('H15', 7), ('H16', 8)]
```

# Citation / Attribution

To cite this work, please use the following [Zenodo DOI](https://zenodo.org/badge/latestdoi/148241360).
