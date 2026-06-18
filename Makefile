#make -f Makefile
#make clean
# Usage:
#   make
#   make clean
#
# This Makefile builds:
#   tests/bessel_driver1
#
# Main source expected:
#   tests/bessel_driver1.f90
#
# Required supporting sources (adjust names if needed):
#   src/parameters.f90
#   src/set_rk.f90
#   tests/Cbessel_miller.f90
#   tests/mod_zbes.f90
#   src/Bessel_Cmplx.f90

# ============================================================
# Platform-specific
# ============================================================
ifeq ($(OS),Windows_NT)
	RM = del /f
    EXE = .exe
else
	RM = rm -f
    EXE =
endif

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
F90FLAGS = -O3 -I src
F90LINKFLAGS = -O3
O = o
MOD = mod

# For debugging with gfortran, you may use:
# F90 = gfortran
# F90FLAGS = -O0 -g -fcheck=all -Wall -Wextra -fbacktrace
# F90LINKFLAGS = -O0 -g -fcheck=all -Wall -Wextra -fbacktrace
# O = o

# ============================================================
# Target name
# ============================================================
TARGET = tests/bessel_driver1$(EXEEXT)

# ============================================================
# Object files
# NOTE:
# If your module Bessel_Cmplx is stored in a different file,
# replace Bessel_Cmplx.$(O) below by the correct object name.
# ============================================================
OBJS = \
	src/set_rk.$(O) \
	tests/Cbessel_miller.$(O) \
    tests/mod_zbes.$(O)\
	src/Bessel_Cmplx.$(O) \
	tests/bessel_driver1.$(O)

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
src/set_rk.$(O): src/set_rk.f90
tests/Cbessel_miller.$(O): tests/Cbessel_miller.f90 src/set_rk.$(O)
src/Cbessel.$(O): src/Cbessel.f90 src/set_rk.$(O)
tests/mod_zbes.$(O): tests/mod_zbes.f90 src/set_rk.$(O)
src/Bessel_Cmplx.$(O): src/Bessel_Cmplx.f90 src/parameters.f90 src/set_rk.$(O)

tests/bessel_driver1.$(O): tests/bessel_driver1.f90 src/parameters.f90 src/set_rk.$(O) tests/mod_zbes.$(O) tests/Cbessel_miller.$(O) src/Bessel_Cmplx.$(O)

.PHONY: clean
clean:
	-$(RM) $(TARGET) $(OBJS) *.mod

check: $(TARGET)
	(cd tests; ./bessel_driver1$(EXEEXT))
