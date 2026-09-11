# MCFBlock

`MCFBlock` is a SMS++ :Block for Linear Min-Cost Flow Problems. It is a pure
:Block depending only on the SMS++ core, and can be solved by any compatible
:Solver, for instance the `MILPSolver` family, `MCFLemonSolver`, or
`MCFSolver` (the latter provided by the separate
[MCFClassSolver](https://gitlab.com/smspp/mcfclasssolver) module, which wraps
the [MCFClass project](http://www.di.unipi.it/optimize/Software/MCF.html)).


## Getting started

These instructions will let you build MCFBlock on your system.


### Requirements

- The [SMS++ core library](https://gitlab.com/smspp/smspp) and its
  requirements.


### Build and install with CMake

Configure and build the library with:

```sh
mkdir build
cd build
cmake ..
cmake --build .
```

The library has the same configuration options of
[SMS++](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration).

Optionally, install the library in the system with:

```sh
cmake --install .
```


### Usage with CMake

After the library is built, you can use it in your CMake project with:

```cmake
find_package(MCFBlock)
target_link_libraries(<my_target> SMS++::MCFBlock)
```


### Running the tests with CMake

`MCFBlock` does not ship a tester of its own, since exercising it requires a
:Solver. The tester that loads a MCF instance into an `MCFBlock` and solves it,
comparing against direct usage of the underlying :MCFClass solver, lives in the
[MCFClassSolver](https://gitlab.com/smspp/mcfclasssolver) module.


### Build and install with makefiles

Carefully hand-crafted makefiles have also been developed for those unwilling
to use CMake. Makefiles build the executable in-source (in the same directory
tree where the code is) as opposed to out-of-source (in the copy of the
directory tree constructed in the build/ folder) and therefore it is more
convenient when having to recompile often, such as when developing/debugging
a new module, as opposed to the compile-and-forget usage envisioned by CMake.

Each executable using `MCFBlock` has to include a "main makefile" of the
module, which typically is either [makefile-c](makefile-c) including all
necessary libraries comprised the "core SMS++" one, or
[makefile-s](makefile-s) including all necessary libraries but not the "core
SMS++" one (for the common case in which this is used together with other
modules that already include them). The makefiles in turn recursively include
all the required other makefiles, hence one should only need to edit the "main
makefile" for compilation type (C++ compiler and its options) and it all
should be good to go. In case some of the external libraries are not at their
default location, it should only be necessary to create the
`../extlib/makefile-paths` out of the `extlib/makefile-default-paths-*` for
your OS `*` and edit the relevant bits (commenting out all the rest).

Check the [SMS++ installation wiki](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration#location-of-required-libraries)
for further details.


## Tools

We provide the simple `dmx2nc4` tool that converts MCF instances written in
the DIMACS standard format (DMX) into netCDF files. Optionally it hacks into
the netCDF file to change the number of static and dynamic nodes and arcs,
as well as the maximum number of nodes and arcs.

You can run the tool from the `<build-dir>/tools` directory, install it
with the library (see above), or just go in the `tools/` folder and run
`make` there (provided makefiles are properly set, see above). Run the tool
without arguments for info on its usage:

```sh
dmx2nc4
```


## Data

We provide a small sample of small-to-mid-size MCF problems in the
[data](data) folder. To get them, run

```sh
cd data
wget https://gitlab.com/smspp/mcfblock/-/package_files/215783773/download -O dmx.tgz
tar xzvf dmx.tgz
```

This creates the `data/dmx` folder containing the instances in the original
DMX text-based format. Once this is available, the netCDF versions of the
instances can be created by

```sh
mkdir nc4
cd ../tools
./batch
```

(provided the `dmx2nc4` tool is available, see above).


## Tests

`MCFBlock` is exercised through a :Solver. The
[MCFClassSolver](https://gitlab.com/smspp/mcfclasssolver) module provides a
tester that reads a MCF instance (in either DIMACS or netCDF format) into an
`MCFBlock`, solves it repeatedly under several changes in costs, capacities and
deficits, arc openings and closures and arc additions and deletions, and
compares the results against direct usage of the underlying solver. The
larger instances it uses are the netCDF files produced by the `dmx2nc4` tool
from the [data](data) folder.


## Getting help

If you need support, you want to submit bugs or propose a new feature, you can
[open a new issue](https://gitlab.com/smspp/mcfblock/-/issues/new).


## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of
conduct, and the process for submitting merge requests to us.


## Authors

### Current Lead Authors

- **Antonio Frangioni**  
  Dipartimento di Informatica  
  Università di Pisa

### Contributors


## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html) -
see the [LICENSE](LICENSE) file for details.


## Disclaimer

The code is currently provided free of charge under an open-source license.
As such, it is provided "*as is*", without any explicit or implicit warranty
that it will properly behave or it will suit your needs. The Authors of
the code cannot be considered liable, either directly or indirectly, for
any damage or loss that anybody could suffer for having used it. More
details about the non-warranty attached to this code are available in the
license description file.
