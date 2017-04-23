##############################################################################
################################# makefile ###################################
##############################################################################
#									     #
#   makefile of testVS10P						     #
#									     #
#   'make clean' cleans up						     #
#   'make' or 'make testVS10P' builds the module			     #
#									     #
#                                VERSION 1.00				     #
#                	        07 - 12 - 2003				     #
#									     #
##############################################################################

# basic directory
DIR = ./

# module name
NAME = $(DIR)testVS10P

# debug switches
#SW = -Wall -g -I$(DIR)

# production switches
SW = -O3 -DNDEBUG -I$(DIR)

# libreries
LIB = -lm

# C++ compiler
CC = g++

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

default: $(NAME)

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(DIR)*.o $(DIR)*~ $(NAME)

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -
# if a module is not used in the current configuration, just comment out the
# corresponding include line
# each module outputs some macros to be used here:
# *OBJ is the final object/library
# *H   is the list of all include files
# *INC is the -I<include directories>

# VrySmplP
VSPDIR = ../VrySmplP/
include $(VSPDIR)makefile

# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

# object files
OBJ =   $(VSPOBJ)\
	$(DIR)Main.o

# include files
MH  =	$(VSPH)

# inherited -I directories
MINC  =	$(VSPINC)

# inherited libreries + -L<libdirs>
MLIB  =	$(VSPLIB)\
	$(LIB)

$(NAME): $(OBJ)
	$(CC) -o $(NAME) $(OBJ) $(MLIB) $(SW)

# dependencies: every .o from its .C + every (recursively) included .h- - - -

$(DIR)Main.o: $(DIR)Main.C $(MH)
	$(CC) -c $*.C -o $@ $(MINC) $(SW)

############################# End of makefile ################################
