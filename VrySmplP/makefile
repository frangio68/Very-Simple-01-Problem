##############################################################################
################################# makefile ###################################
##############################################################################
#									     #
#   makefile of VrySmplP						     #
#									     #
#   Input:  $(CC)     = compiler command				     #
#           $(SW)     = compiler options                                     #
#           $(VSPDIR) = the directory where the VrySmplP source is	     #
#									     #
#   Output: $(VSPOBJ) = the final object(s) / library			     #
#           $(VSPLIB) = external libreries + -L<libdirs>		     #
#           $(VSPH)   = the .h files to include				     #
#           $(VSPINC) = the -I$(VrySmplP directory)			     #
#									     #
#                                VERSION 2.00				     #
#                	        20 - 05 - 2005				     #
#									     #
#			      Antonio Frangioni				     #
#   			   Operations Research Group			     #
#			  Dipartimento di Informatica			     #
#   			     Universita' di Pisa			     #
#									     #
##############################################################################

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(VSPDIR)*.o $(VSPDIR)*~

# macroes to be exported- - - - - - - - - - - - - - - - - - - - - - - - - - -

VSPOBJ = $(VSPDIR)VrySmplP.o
VSPLIB =
VSPH   = $(VSPDIR)VrySmplP.h
VSPINC = -I$(VSPDIR)

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(VSPOBJ): $(VSPDIR)VrySmplP.C $(VSPH)
	$(CC) -c $*.C -o $@ $(SW)

############################# End of makefile ################################
