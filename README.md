# MCFBlock / MCFSolver

This project covers two conceptually different things (which may one day be
split to two different projects):

- `MCFBlock`, a SMS++ :Block for Linear Min-Cost Flow Problems

- `MCFSolver`, a SMS++ :Solver for MCFBlock based on forwarding the interface
  of (objects derived from the general abstract) [MCFClass of the
  MCFClass project](http://www.di.unipi.it/optimize/Software/MCF.html)

## Getting started

These instructions will let you build MCFBlock and MCFSolver on your system.

### Requirements

- [SMS++ core library](https://gitlab.com/smspp/smspp)
- [MCFClass](https://github.com/frangio68/Min-Cost-Flow-Class)

### Build and install

Configure and build the library with:
```sh
mkdir build
cd build
cmake ..
make
```

Optionally, install the library in the system with:
```sh
sudo make install
```

## Usage

After the library is configured and built, you can use it in your CMake project with:
```cmake
find_package(MCFBlock)
target_link_libraries(<my_target> SMS++::MCFBlock)
```

### Tools

We provide a simple tool that converts MCF instances written in the DIMACS
standard into netCDF files.
Optionally it hacks into the netCDF file to change the number of static
and dynamic nodes and arcs, as well as the maximum number of nodes and
arcs.

You can run the tool from the `<build-dir>/tools` directory or install it with
the library (see above).
Run the tool without arguments for info on its usage:

```sh
dmx2nc4
```

## Running the tests

A simple unit test will be built with the library,
To disable it, configure the library with the option `-DBUILD_TESTING=OFF`.

The test takes an instance of a MCF in DIMACS or NC4 format.
The MCF problem is then repeatedly solved with several changes in
costs/capacities/deficits, arcs openings/closures and arcs additions/deletions.
The same operations are performed on the two solvers,
and the results are compared.

## Contributing

This section is not ready yet.

## Authors

### Current Lead Authors

- **Antonio Frangioni**  
  *Operations Research Group*  
  Dipartimento di Informatica  
  Università di Pisa

### Contributors

## License

This code will provided under the LGPL license when it will be released.

## Disclaimer

The code is currently provided free of charge for academic purposes only.
As such, it is provided "*as is*", without any explicit or implicit warranty
that it will properly behave or it will suit your needs. The Authors of
the code cannot be considered liable, either directly or indirectly, for
any damage or loss that anybody could suffer for having used it. More
details about the non-warranty attached to this code are available in the
license description file.
