
MCFBlock / MCFSolver
====================

This project covers two conceptually different things (which may one day be
splitted to two different projects):

- MCFBlock, a SMS++ :Block for Linear Min-Cost Flow Problems

- MCFSolver, a SMS++ :Solver for MCFBlock based on forwarding the interface
  of (objects derived from the general abstract) [MCFClass of the
  MCFClass project](http://www.di.unipi.it/optimize/Software/MCF.html)

The project obviously requires the ["core"
SMS++](https://gitlab.com/frangio68/sms_plus_plus); however, it also requires
(for the MCFSolver component) the [MCFClass
project](https://github.com/frangio68/Min-Cost-Flow-Class).

The arrangement of folders is assumed to be

    some root folder
     |
	 +- SMS++ (the "core" SMS++)
	 |
	 +- MCFBlock (this project)
	 |
	 +- MCFClass (the MCFClass project)

MCFBlock/test contains a reasonably sophisticated test that takes some existing
MCF instance in the DIMACS standard format, reads it into a MCFBlock and an
:MCFClass object, create a MCFSolver using the very same :MCFClass solver,
and then repeatedly modifies both instances in the same way (possibly
passing to a "copy" R3Block of the initial MCFBlock) in order to verify that
the different Modification involved correctly do their work.

License
-------

This code will provided under the LGPL license when it will be released.

Authors
=======

Lead Authors
------------

	Antonio Frangioni
	Operations Research Group
	Dipartimento di Informatica
	Universita' di Pisa
 
Contributors
------------
