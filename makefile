
#
#  Make file for two-layer diffusion model
#


#  Defining variables
#

f90comp = gfortran
CFLAGS=-I.

#
#     Compiles the fortran code
#

two_layer_diffusion: 
	$(f90comp) -o two_layer_diffusion two_layer_diffusion.f90 -I.



