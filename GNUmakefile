AMREX_HOME ?= ../../amrex

CHEMISTRY_MODEL = LiDryer
#CHEMISTRY_MODEL = DRM19
#CHEMISTRY_MODEL = GRI30

CHEMISTRY_DIR = chemistry/$(CHEMISTRY_MODEL)

DEBUG	= FALSE

DIM	= 3

COMP    = gcc

USE_MPI   = FALSE
USE_OMP   = FALSE
USE_CUDA  = FALSE
USE_HIP   = FALSE
USE_SYCL  = FALSE

BL_NO_FORT = TRUE

TINY_PROFILE = FALSE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(CHEMISTRY_DIR)/Make.package
include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
