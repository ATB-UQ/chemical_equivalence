[![DOI](https://zenodo.org/badge/148241360.svg)](https://zenodo.org/badge/latestdoi/148241360)

# Requirements

* Python `>=3.9` (the code uses builtin generics such as `list[tuple[int, ...]]`)

* ATB Outputs python module: `https://github.com/ATB-UQ/atb_outputs.git` (provides `MolData`; declared in `setup.cfg` as a `git+ssh` dependency, so in the platform venv install it first / use `--no-deps`)

* `atb_nauty` (`https://github.com/ATB-UQ/atb_nauty.git`): nauty 2.8.8 in process, a C extension that `pip install` builds from source (it needs only a C compiler). Until 2026-10 this module ran nauty as the `dreadnaut` program from a vendored nauty 2.5r9 build; that program, its `DREADNAUT_EXECUTABLE` setting and the import-time check that it existed are gone.

# Configuration

* Nothing of its own. `ATB_NAUTY_MAX_VERTICES` (read by `atb_nauty`, default 20 000) caps the graph size nauty is given, molecule plus stereo gadgets (up to 9 extra vertices per tetrahedral centre).

# What "chemically equivalent" means here

Two atoms are equivalent when some rotation or reflection of the molecule, in any accessible conformer, exchanges them. That is more than the graph alone can say: nauty's orbits of the molecular graph give *constitutional* equivalence, under which the two hydrogens of a CH2 next to a stereocentre are equivalent although they are diastereotopic. `stereo.py` reads the tetrahedral orientation of every 4-coordinate atom and the cis/trans relations across acyclic double bonds from the coordinates (`ocoord` if present, else `coord`) and encodes them as gadgets on the graph, so that the surviving automorphisms are exactly the label permutations realisable by a rotation (all orientations kept) or a reflection (all flipped). Enantiotopic atoms -- meso compounds, chlorocyclohexane's two ortho carbons -- stay equivalent; diastereotopic ones split; ring inversion and bond rotation need no special treatment because tetrahedral orientation is a configurational invariant.

`getChemEquivGroups(..., max_stereo_distance=N)` merges such splits back for atoms more than `N` heavy-atom bonds from every stereo element, where the effect is negligible and identical parameters on equivalent atoms matter more; `None` keeps them all. `mol_data.equivalence_details` records the constitutional, strict and effective classes and the stereo elements found.

Not detected: axial chirality (allenes, atropisomers, spiro compounds), pyramidal three-coordinate centres (sulfoxides, phosphines), amide C-N bonds (too long, ~0.136-0.139 nm against the 0.1337 nm cutoff, so the two N-substituents of e.g. N,N-dimethylacetamide stay equivalent) and C=N double bonds of imines, oximes and hydrazones (the 2-coordinate N fails the three-neighbour sp2 rule, so syn/anti substituents stay equivalent). Both are pinned by tests in `testing/test_review_2026_10.py`.

# Usage

A simple example is included in `examples/test.py` (run it from `src/chemical_equivalence/testing/`: it opens `benzene.pdb` relative to the working directory). `src/chemical_equivalence/test.py` is an older timing/graph-image harness, not a pytest module.

The tests live in `src/chemical_equivalence/testing/` (fixture molecules as PDB files, relation-style expectations, and invariance properties on RDKit-embedded structures):

```
/home/atb/ATB/.venv/bin/python -m pytest src/chemical_equivalence/testing -v
```

 * First, a `MolData` (Molecule Data) object has to be initialised. It describes the topology and coordinate of the molecule.
Please see the documentation of the `atb_outputs` module for further description of the `MolData` object.

```
>>> from chemical_equivalence.helpers.types_helpers import MolData
>>> mol_data = MolData(open('src/chemical_equivalence/testing/chlorocyclohexane.pdb').read())
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

# Platform role

Status: live-production support. `core.atb.atb` calls `getChemEquivGroups` for every
topology (equivalence groups drive charge/parameter symmetrisation). Other importers:
`NMR_Interface`, the website v0.1 API (three modules import the constant
`ALL_EXCEPTION_SEARCHING_KEYWORDS`),
`core/testing/equivalence_sweep.py`, `atb_molgraph`'s benchmarks and a `PyMultifitESP` example.
`Blind_RMSD` (retired 2026-10-07) was the only caller of `partial_mol_data_for_pdbstr` and of
the cross-molecule `get_chemical_equivalence_accross` (deleted in the 2026-10 review, with the
CLI's `--other-pdbs`); `fragment_merger` does not import it.
nauty is asked two things, both invariants of the coloured graph (so the answers do not
depend on the nauty version or the vertex numbering): the automorphism orbits, and whether
the stereo graph is isomorphic to its mirror image. Each call is in process and cannot be
interrupted, so it has no timeout of its own: a molecular graph takes well under a
millisecond (the 252-atom nanotube 0.27 ms), the callers' own limits (Celery
`soft_time_limit`, the fcgid busy timeout) still apply, and `atb_nauty` refuses a graph over
`ATB_NAUTY_MAX_VERTICES` before calling nauty. Bad input (an atom id outside 1..N, for
instance) raises `atb_nauty.NautyError`, also importable as `NautyInterface.DreadnautError`.
`calcChemEquivalency.py` also has a small `__main__` argparse entry point.
`examples/test.py` (`Molecule3D_example`) additionally needs the `NXMol`
(`chemistry_data_structure`) sibling.

# Citation / Attribution

To cite this work, please use the following [Zenodo DOI](https://zenodo.org/badge/latestdoi/148241360).
