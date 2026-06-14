# test

A tester for `MCFBlock` and `MCFSolver`.

It reads a Min-Cost Flow instance (in either DIMACS or netCDF format) into an
`MCFBlock` and, from there, into an object of a `:MCFClass` chosen at compile
time, then attaches an `MCFSolver< :MCFClass >` to the `MCFBlock`. The problem
is repeatedly solved while changing costs, capacities and deficits, opening
and closing arcs, and adding and deleting arcs; the same operations are
applied to both objects and the results are compared. This mostly exercises
`MCFBlock` and `MCFSolver`, since the underlying `:MCFClass` is the same in
the two paths.

Instances live in the `data` folder. A `batch` (and `batch-l`) file is
provided that runs the tester over a set of them; all of them passing is a
good sign that no regressions have been introduced. The `makefile` builds the
executable including the `MCFBlock` module and the core SMS++ library.


## Authors

- **Antonio Frangioni**  
  Dipartimento di Informatica  
  Università di Pisa


## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html),
see the [LICENSE](../LICENSE) file for details.
