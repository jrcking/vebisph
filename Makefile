# call with make const_mod=X      

# X=1 Newtonian
# X=2 Linear PTT

FC := gfortran
LD := gfortran
CFLAGS := -Wall -O3 -g -m64
FFLAGS :=-fopenmp -fbounds-check -ffpe-trap=zero -O3 -Wall -g -J./obj -I./obj -m64

ifeq ($(strip $(const_mod)),)
FFLAGS += -Dconst_mod=1
else
FFLAGS += -Dconst_mod=$(const_mod)
endif

LDFLAGS := -fopenmp -m64

SUB_DIRS := para parallel_tools base tools
SRC_DIR  := $(addprefix source/,$(SUB_DIRS))

#common_parameter needs to be first as mod file is depended upon by nearly everything.
OBJ_FILES := obj/kind_parameters.o obj/common_parameter.o obj/cell_griddin.o obj/common_2d.o
OBJ_FILES += obj/sphtools.o obj/input.o
OBJ_FILES += obj/freesurf_mod.o obj/mirror_boundaries_mod.o obj/gradient_calculation.o 
OBJ_FILES += obj/neighbour_finding_mod.o obj/nnf_models_mod.o obj/new_fick_shift.o obj/linear_solver.o
OBJ_FILES += obj/predictor_mod.o obj/correct_velocity.o 
OBJ_FILES += $(foreach sdir,$(SRC_DIR),$(patsubst $(sdir)/%.F90,obj/%.o,$(wildcard $(sdir)/*.F90)))

#HDEPS := $(OBJ_FILES:.o=.d)

vpath %.F90 $(SRC_DIR)

#-------


default: vebisph
vebisph: $(OBJ_FILES)
	$(LD) -o $@ $^ $(LDFLAGS)

obj/%.o: %.F90
	$(FC) $(FFLAGS) -c -o $@ $<

-include $(HDEPS)


cleanup:
	rm -vf ./obj/*.o
	rm -vf ./obj/*.mod
	rm -vf ./vebisph
	rm -vf I* TEMP_REC.dat
	rm -vf ./data_directory/PART*
	rm -vf ./data_directory/TIME
	rm -vf ./data_directory/TIME_OUT
	rm -vf ./data_directory/LS
	rm -vf ./paraview_files/*
	rm -vf fort*


