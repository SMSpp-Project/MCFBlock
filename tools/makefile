##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of dmx2nc4                                                      #
#                                                                            #
#                              Antonio Frangioni                             #
#                          Dipartimento di Informatica                       #
#                              Universita' di Pisa                           #
#                                                                            #
##############################################################################

# module name
# NAME can be main, test and dmx2nc4:, the corresponding .cpp is compiled
# and the final executable named $(NAME) is constructed 
NAME = dmx2nc4

# basic directory
DIR = .

# switch to select the :MCFClass solver, supported values are:
# -DHAVE_CSCL2      for the CS2 class
# -DHAVE_CPLEX      for the MCFCplex class
# -DHAVE_MFSMX      for the MCFSimplex class
# -DHAVE_MFZIB      for the MCFZIB class
# -DHAVE_RELAX      for the RelaxIV class
# -DHAVE_CPLEX      for the MCFCplex class
# -DHAVE_SPTRE      for the SPTree class; note that SPTree cannot solve
#                   most MCF instances, except those with SPT structure
MCFCS = -DHAVE_MFSMX

# common flags
COMMON_SW = -std=c++17 $(MCFCS)

# debug switches
SW_DEBUG = -g -glldb -fno-inline $(COMMON_SW)
# production switches
SW_RELEASE = -O3 $(COMMON_SW) -DNDEBUG

# compiler
CC = clang++

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

BUILD_MODE = release

default: $(BUILD_MODE)

# debug target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

debug: SW = $(SW_DEBUG)
debug: BUILD_MODE = debug
debug: $(DIR)/$(NAME)

# release target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

release: SW = $(SW_RELEASE)
release: BUILD_MODE = release
release: $(DIR)/$(NAME)

# clean target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(DIR)/*.o $(DIR)/*~ $(DIR)/$(NAME)

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -
# if a module is not used in the current configuration, just comment out the
# corresponding include line
# each module outputs some macros to be used here:
# *OBJ is the final object(s) / library
# *LIB is the external libraries + -L< libdirs >
# *H   is the list of all include files
# *INC is the -I< include directories >

# MCFBlock
MCFBkSDR = ..
include $(MCFBkSDR)/makefile-c

# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

# object files
MOBJ =  $(MCFBkOBJ)

# libraries
MLIB =  $(MCFBkLIB)

$(DIR)/$(NAME): $(MOBJ) $(DIR)/$(NAME).o
	$(CC) -o $(DIR)/$(NAME) $(DIR)/$(NAME).o $(MOBJ) $(MLIB) $(SW)

# dependencies: every .o from its .cpp + every recursively included .h- - - -

# include directives
MINC =  $(MCFBkINC)

# includes
MH =    $(MCFBkH)

# compile command

$(DIR)/$(NAME).o: $(DIR)/$(NAME).cpp $(MH)
	$(CC) -c $*.cpp -o $@ $(MINC) $(SW)

# install/uninstall - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# check if the user is in the sudoers file
ifeq ($(shell id -u),0)
    prefix = /opt/smspp-project
else
    prefix = $(HOME)/smspp-project
endif
bindir = $(prefix)/$(BUILD_MODE)/bin

install: $(DIR)/$(NAME)
	install -s -D $(DIR)/$(NAME) $(DESTDIR)$(bindir)/$(NAME)

uninstall:
	rm -f $(DESTDIR)$(bindir)/$(NAME)

# distclean target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

distclean: clean

# phony targets - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.PHONY: debug release install uninstall clean distclean

############################ End of makefile #################################
