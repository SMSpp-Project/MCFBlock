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

## Tests

`test` directory contains a reasonably sophisticated test that takes some existing
MCF instance in the DIMACS standard format, reads it into a MCFBlock and an
:MCFClass object, create a MCFSolver using the very same :MCFClass solver,
and then repeatedly modifies both instances in the same way (possibly
passing to a "copy" R3Block of the initial MCFBlock) in order to verify that
the different Modification involved correctly do their work.

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
