#make -f Makefile
#make clean
# Usage:
#   make
#   make clean
#
# This Makefile builds:
#   bessel_driver1
#
# Main source expected:
#   bessel_driver1.f90
#
# Required supporting sources (adjust names if needed):
#   set_rk.f90
#   Cbessel_miller.f90
#   mod_zbes.f90
#   Bessel_Cmplx.f90

# ============================================================
# Remove command
# ============================================================
RM = rm -f
# For Windows CMD, you may use instead:
 #RM = del /f

# ============================================================
# Compiler selection
# ============================================================

# ---- NAG ----
# F90 = nagfor
# F90FLAGS = -O3 -g -g90 -gline -nan -C=all -C=undefined -C=dangling -maxcontin=999
# F90LINKFLAGS = -O3 -g -g90 -gline -nan -C=all -C=undefined -C=dangling -maxcontin=999
# O = o
# EXEEXT =

# ---- Intel ifort ----
# F90 = ifort
# F90FLAGS = -O3 /WARN:ALL /WARN:DECLARATIONS /STAND:F90 /TRACEBACK /CU
# F90LINKFLAGS = -O3 /WARN:ALL /WARN:DECLARATIONS /STAND:F90 /TRACEBACK /CU
# O = obj
# EXEEXT = .exe

# ---- Intel ifx ----
# F90 = ifx
# F90FLAGS = -O3 /WARN:ALL /WARN:DECLARATIONS /STAND:F90 /TRACEBACK /CU
# F90LINKFLAGS = -O3 /WARN:ALL /WARN:DECLARATIONS /STAND:F90 /TRACEBACK /CU
# O = obj
# EXEEXT = .exe

# ---- GNU Fortran ----
F90 = gfortran
F90FLAGS = -O3
F90LINKFLAGS = -O3
O = o
EXEEXT = .exe

# For debugging with gfortran, you may use:
# F90 = gfortran
# F90FLAGS = -O0 -g -fcheck=all -Wall -Wextra -fbacktrace
# F90LINKFLAGS = -O0 -g -fcheck=all -Wall -Wextra -fbacktrace
# O = o
# EXEEXT = .exe

# ============================================================
# Target name
# ============================================================
TARGET = bessel_driver1$(EXEEXT)

# ============================================================
# Object files
# NOTE:
# If your module Bessel_Cmplx is stored in a different file,
# replace Bessel_Cmplx.$(O) below by the correct object name.
# ============================================================
OBJS = \
	set_rk.$(O) \
	Cbessel_miller.$(O) \
                  mod_zbes.$(O)\
	Bessel_Cmplx.$(O) \
	bessel_driver1.$(O)

# ============================================================
# Default target
# ============================================================
all: $(TARGET)

# ============================================================
# Compilation rule
# ============================================================
%.$(O): %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

# ============================================================
# Linking rule
# ============================================================
$(TARGET): $(OBJS)
	$(F90) $(F90LINKFLAGS) -o $(TARGET) $(OBJS)

# ============================================================
# Explicit dependencies
# ============================================================
set_rk.$(O): set_rk.f90
Cbessel_miller.$(O): Cbessel_miller.f90 set_rk.$(O)
Cbessel.$(O): Cbessel.f90 set_rk.$(O)
mod_zbes.$(O): mod_zbes.f90 set_rk.$(O)
Bessel_Cmplx.$(O): Bessel_Cmplx.f90 set_rk.$(O) Cbessel_miller.$(O) 

bessel_driver1.$(O): bessel_driver1.f90 set_rk.$(O) Cbessel_miller.$(O)  Bessel_Cmplx.$(O)


# Clean rule (Windows CMD)
.PHONY: clean
clean:
	-del /f *.o *.mod *.exe 
