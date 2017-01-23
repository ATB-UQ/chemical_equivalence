# Requirements

* Python `>=3.5`

* ATB Outputs module: `https://github.com/bertrand-caron/atb_outputs`
The directory in which should be accesible by the `PYTHONPATH` variable of your shell.

* Nauty: This modules relies on the `dreadnaut` executable which is part of the `nauty` package.
Nauty can very easily be installed from [source](http://users.cecs.anu.edu.au/~bdm/nauty/) or using a package manager such as `homebrew` for Mac OS X users.
The code assumes that `dreadnaut` will be installed in `/usr/local/bin`. If installing in a different location, you can change the `NAUTY_EXECUTABLE` path in `config.py`.

# Usage

The `test.py` script contains code to run the chemical equivalence prediction algorithm on a series of test cases found in the `testing/` directory.

 * First, a `MolData` (Molecule Data) object has to be initialised. It describes the topology and coordinate of the molecule.

```
>>> from chemical_equivalence.helpers.types_helpers import MolData
>>> mol_data = MolData(open('testing/chlorocyclohexane.pdb').read())
```

* The chemical equivalence can be run in a single line of code.
It retuns a dictionary mapping each atom index (`atom['id']`) to a unique `int` for each equivalence class:

```
>>> from chemical_equivalence.calcChemEquivalency import getChemEquivGroups
>>> equivalence_dict, n_iterations = getChemEquivGroups(mol_data)
>>> equivalence_dict
>>> list((atom['symbol'], equivalence_dict[atom['id']]) for atom in mol_data.atoms.values())
```
