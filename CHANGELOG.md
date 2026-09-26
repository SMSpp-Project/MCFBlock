# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `is_direction()`, `is_direction( bool )` and `has_directions()`: a MCFBlock
  knows what a direction of its own is, i.e., a flow that conserves at every
  node with all the deficits zero, that is nonnegative on the arcs of
  infinite capacity and that is zero on all the others, an arc of finite
  capacity leaving no room to move for ever. `is_sol_feasible()` takes from
  the `MCFSolution` whether what it holds is a solution or a direction, that
  method not going through the Variable, and checks it accordingly

### Changed

- the data archive is downloaded by version: `DATA_VERSION` in CMakeLists.txt
  names the version of the Package Registry to read, and the archive and the
  marker of its extraction carry it in their name, so that a tree holding
  an older extraction (the cache of the CI, or a clone extracted before)
  downloads and extracts again instead of running on the old data;
  data/upload-dmx publishes the archive under that version

- whoever links the module keeps it: the classes of a module register
  themselves in the factory from a static initialiser, and a linker that
  drops what looks unused takes the registration away with it, so the target
  now tells whoever links it to keep the symbol that forces the module in,
  and on ELF, where naming the symbol is not enough, the library as a whole

- `chg_costs()`, `chg_ucaps()` and `chg_dfcts()` take their data as a
  `std::span< const double >`, whose length they check against the Range or
  the Subset instead of reading past the end, and are registered in the
  methods factory in that form too; the forms taking an iterator stay, and
  defer to the span ones

### Fixed

- `bound_feasible()`, `dual_feasible()`, `complementary_slackness()` and
  the DIMACS `print()` took quadratic time in the number of dynamic arcs,
  asking `is_closed()` of each arc, which walks the list of the dynamic flow
  Variable from one of its ends: they now read which arcs are closed in one
  pass [see `closed_arcs()`], and on a GOTO instance of 1024 nodes and 65536
  arcs the check of a solution goes from about 700 to about 5 million
  instructions

- on macOS a program linking the module lost the classes the module
  registers in the factories when the linker dropped the library, as it
  does under `-dead_strip_dylibs`, which conda sets: the target now asks the
  linker for the symbol that forces the module in (`-u`), which ld64,
  unlike the ELF linker, counts as a use of the library

- the step that fetches the data archive of this module says what went wrong
  when it goes wrong: the download is checked, an archive that did not arrive
  is removed instead of being left on disk for the build to take for the real
  one, and the message names the URL. A server that answers with an error page
  used to leave a file of a few bytes there, which made the next build fail
  while extracting it, with the message of `tar` and no mention of the
  download

- `load()` given no capacities or no deficits, which it documents as all
  capacities infinite and all deficits zero, read past the end of the empty
  vector

- `chg_ucap()` on a MCFBlock with no capacities returned without doing
  anything when the new capacity was finite, i.e., exactly when there was
  something to change, and `chg_dfct()` on one with no deficits sized them
  on the current number of nodes rather than on the maximum one

- `map_forward_Modification()` passed a change of the deficits on to the
  other `MCFBlock` as a change of its capacities, and, when this `MCFBlock`
  had no capacities or no deficits, read them out of the empty vector rather
  than passing on the infinite capacities and the zero deficits it has

- `map_forward_solution()` looked at the static bound Constraint where it
  had to look at the dynamic ones, setting to zero the duals of the dynamic
  bounds of the other `MCFBlock` whenever this one had no static bounds

## [0.6.1] - 2026-09-13

### Changed

- `is_sol_feasible()` and `is_sol_optimal()` say that they leave the
  Variable alone

## [0.6.0] - 2026-09-12

### Added

- accessors and setters to the flows and the potentials of `MCFSolution`, so
  that a Solver can fill it directly out of its own data structures rather
  than by writing the solution in the Variable and the Constraint of the
  `MCFBlock` and having it read back from there

- `is_sol_feasible()` and `is_sol_optimal()`, which check a Solution, i.e., a
  flow and potentials taken from the outside

### Changed

- the version of the module is the git tag of its repository, or the
  VERSION.txt of a release tarball, and the shared library carries it: its
  SONAME is major.minor while the major is 0, and it is installed with an
  RPATH relative to itself, so that an installed tree keeps working wherever
  it is moved

### Removed

- the MCFClass submodule, which is in MCFClassSolver

## [0.5.1] - 2025-12-12

### Fixed

- the order in which some conditions are checked

- avoid static vectors prone to static initialization fiasco

## [0.5.0] - 2024-02-28

### Changed

- adapted to new CMake / makefile organisation

- MCFClass is now a submodule of MCBlock rather than having to be a
  submodule of the umbrella (i.e., it is now found in ./MCFClass
  rather than in ../MCFClass)

## [0.4.4] - 2023-05-23

### Changed

- now using general RowConstraint::is_feasible()

## [0.4.3] - 2022-08-25

### Fixed

- include header <iomanip>

## [0.4.2] - 2022-06-28

Mainly a fix release after much improved testing

### Added

- added support for RelaxIV

- added un\_ModBlock and get\_objective_value()

- added MCFBlock::set\_*( Subset )

### Changed

- completely reworked arc deletion in the middle: the arc
  is no longer immediately removed from the existing flow
  conservation constraints; this only happens (if necessary)
  when a new arc is added in its place. this required several
  changes in the logic (arc being deleted in now
  C[ arc ] == NaN) that brought several changes, hopefully
  leading to less work by not changing stuff related to deleted
  arcs (which is useless)

- improved tester

- completely reworked selection of :MCFClass in MCFSolver

- adapted to new channel management

- adapted to new load/print interface

- reworked signature of get\_* and set\_* methods in MCFBlock

### Fixed

- fixed another stoopid bug in MCFBlock::open_arcs( Range )

- caught a logic error: a user may try to un-fix a variable
  corresponding to a deleted arc, which would succeed in the
  abstract representation but not in the physical one: this is
  now found out and exception is thrown

- fixed issues about open/close-ing deleted arcs: both in
  MCFBlock and MCFSolver, since an arc that is currently not
  deleted in MCFBlock can still be "seen" as deleted in a
  MCFSolver if the Modification re-adding it has not been
  processed yet

- fixed blunder in new logic for open/close_arc( Range )

- fixed minor (but significant) bug

- fixed stupid bug in set\_*( range )

- corrected blunder in map\_forward\_Modification

- added proper lock() and unlock() in MCFSolver

- fixed flaw in MCFSolution

## [0.4.1] - 2021-12-07

Minor point release to avoid the master branch to become too stale:

### Changed

- changed useabstract to hint from order

### Fixed

- fixed flaw in flow_feasible()

- fixed an issue in MCFSolution + minor changes

## [0.4.0] - 2021-02-05

### Added

- Managed vectors in Configurations.

- Added some preprocess.

### Fixed

- Using proper eps in var.is_feasible.

## [0.3.1] - 2020-09-24

### Fixed

- Workaround for default MCFSolver setting.

## [0.3.0] - 2020-09-16

### Added

- Support for MCFCplex class.

- Support for new configuration framework.

## [0.2.0] - 2020-03-06

### Fixed

- Just updated to release version.

## [0.1.2] - 2020-03-04

### Fixed

- Minor fixes in namespace use.

## [0.1.1] - 2020-02-10

### Fixed

- Minor fix in makefile support.

## [0.1.0] - 2020-02-07

### Added

- First test release.

[Unreleased]: https://gitlab.com/smspp/mcfblock/-/compare/0.6.1...develop
[0.6.1]: https://gitlab.com/smspp/mcfblock/-/compare/0.6.0...0.6.1
[0.6.0]: https://gitlab.com/smspp/mcfblock/-/compare/0.5.1...0.6.0
[0.5.1]: https://gitlab.com/smspp/mcfblock/-/compare/0.5.0...0.5.1
[0.5.0]: https://gitlab.com/smspp/mcfblock/-/compare/0.4.4...0.5.0
[0.4.4]: https://gitlab.com/smspp/mcfblock/-/compare/0.4.3...0.4.4
[0.4.3]: https://gitlab.com/smspp/mcfblock/-/compare/0.4.2...0.4.3
[0.4.2]: https://gitlab.com/smspp/mcfblock/-/compare/0.4.1...0.4.2
[0.4.1]: https://gitlab.com/smspp/mcfblock/-/compare/0.4.0...0.4.1
[0.4.0]: https://gitlab.com/smspp/mcfblock/-/compare/0.3.1...0.4.0
[0.3.1]: https://gitlab.com/smspp/mcfblock/-/compare/0.3.0...0.3.1
[0.3.0]: https://gitlab.com/smspp/mcfblock/-/compare/0.2.0...0.3.0
[0.2.0]: https://gitlab.com/smspp/mcfblock/-/compare/0.1.2...0.2.0
[0.1.2]: https://gitlab.com/smspp/mcfblock/-/compare/0.1.1...0.1.2
[0.1.1]: https://gitlab.com/smspp/mcfblock/-/compare/0.1.0...0.1.1
[0.1.0]: https://gitlab.com/smspp/mcfblock/-/tags/0.1.0
