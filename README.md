[![DOI](https://zenodo.org/badge/148241360.svg)](https://zenodo.org/badge/latestdoi/148241360)

# Requirements

* Python `>=3.5`

* ATB Outputs python module: `https://github.com/ATB-UQ/atb_outputs.git

* Nauty: This modules relies on the `dreadnaut` executable which is part of the `nauty` package.
Nauty can very easily be installed from [source](http://users.cecs.anu.edu.au/~bdm/nauty/) or using a package manager such as `homebrew` for Mac OS X users.

# Configuration

* The code assumes that `dreadnaut` will be installed in `/usr/local/bin`. If installing in a different location, you can change the `NAUTY_EXECUTABLE` path in `config.py`.

# Usage

A simple example is included in `examples/test.py`

For more comprehensive examples see the `src/chemical_equivalence/test.py` script which contains code to run the chemical equivalence prediction algorithm on a series of test cases found in the `testing/` directory.

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

* It returns a dictionary mapping each atom index (`atom['id']`) to a unique `int` for each equivalence class:

```
>>> equivalence_dict
{1: 0, 2: 1, 3: 2, 4: 3, 5: 4, 6: 5, 7: 6, 8: 7, 9: 8, 10: 9, 11: 10, 12: 11, 13: 12, 14: 13, 15: 14, 16: 15, 17: 16, 18: 17}
```

* It can easily be mapped back to each atom by iterating over the atoms of the `MolData` object:

```
>>> list((atom['symbol'], equivalence_dict[atom['id']]) for atom in mol_data.atoms.values())
[('H8', 0), ('C2', 1), ('Cl1', 2), ('C7', 3), ('H17', 4), ('H18', 5), ('C3', 6), ('H9', 7), ('H10', 8), ('C4', 9), ('H11', 10), ('H12', 11), ('C5', 12), ('H13', 13), ('H14', 14), ('C6', 15), ('H15', 16), ('H16', 17)]
```

# Citation / Attribution

To cite this work, please use the following [Zenodo DOI](https://zenodo.org/badge/latestdoi/148241360).
