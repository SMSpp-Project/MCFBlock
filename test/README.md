# test

A tester for `MCFBlock` that needs nothing but the core SMS++ library: no
`Solver` of the MCF is attached, the `MCFBlock` being checked by itself.

It builds small Min-Cost Flow instances in memory and checks the DIMACS
(`print( 'C' )` and `load( std::istream )`) and netCDF (`serialize()` and
`deserialize()`) round trips of instances with no arcs, with unbalanced
deficits, with infinite capacities, with closed arcs and with deleted arcs;
the files are written in the working directory and removed afterwards.

Every change of the physical data (costs, capacities, deficits, closing and
opening arcs, adding and removing arcs), in its single, `Range` and `Subset`
forms, empty and unordered ones included, is checked against the abstract
representation (the flow `Variable`, the coefficients of the `Objective`, the
bound and flow conservation `Constraint`), which has to follow it, and
against the `Modification` received by a `FakeSolver` of the core registered
on the `MCFBlock`, with the `eNoMod`, `eModBlck` and `eDryRun` values of the
parameters as well; so is the other way round, a change of the abstract
representation that the `MCFBlock` has to translate into its data. The arcs
deleted before the `Constraint` are generated, and the instance whose arcs
are all dynamic and are all removed, are checked on their own.

A flow set by hand in the `Variable` is checked for its objective value and
its feasibility (flow conservation, capacities, closed and deleted arcs,
unbalanced deficits, tolerance), and so are the valid bounds on the optimal
value. Finally, a sequence of random changes, drawn with a fixed seed, is
checked step by step against a copy of the data kept by the tester.

The exit code is the number of failed checks, zero printing
`All tests passed!!`. The `makefile` builds the executable including the
`MCFBlock` module and the core SMS++ library.


## Authors

- **Donato Meoli**  
  Dipartimento di Informatica  
  Università di Pisa


## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html),
see the [LICENSE](../LICENSE) file for details.
