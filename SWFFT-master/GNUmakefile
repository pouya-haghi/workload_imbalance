#                  Copyright (C) 2017, UChicago Argonne, LLC
#                             All Rights Reserved
#
#            Hardware/Hybrid Cosmology Code (HACC), Version 1.0
#
#  Salman Habib, Adrian Pope, Hal Finkel, Nicholas Frontiere, Katrin Heitmann,
#       Vitali Morozov, Jeffrey Emberson, Thomas Uram, Esteban Rangel
#                         (Argonne National Laboratory)
#
#   David Daniel, Patricia Fasel, Chung-Hsing Hsu, Zarija Lukic, James Ahrens
#                       (Los Alamos National Laboratory)
#
#                                George Zagaris
#                                  (Kitware)
#
#                             OPEN SOURCE LICENSE
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#   1. Redistributions of source code must retain the above copyright notice,
#      this list of conditions and the following disclaimer. Software changes,
#      modifications, or derivative works, should be noted with comments and
#      the author and organization’s name.
#
#   2. Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in the
#      documentation and/or other materials provided with the distribution.
#
#   3. Neither the names of UChicago Argonne, LLC or the Department of Energy 
#      nor the names of its contributors may be used to endorse or promote 
#      products derived from this software without specific prior written 
#      permission.
#
#   4. The software and the end-user documentation included with the
#      redistribution, if any, must include the following acknowledgment:
#
#     "This product includes software produced by UChicago Argonne, LLC under
#      Contract No. DE-AC02-06CH11357 with the Department of Energy."
#
# *****************************************************************************
#                                DISCLAIMER
# THE SOFTWARE IS SUPPLIED "AS IS" WITHOUT WARRANTY OF ANY KIND. NEITHER THE
# UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR 
# UCHICAGO ARGONNE, LLC, NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, 
# EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE
# ACCURARY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, DATA, APPARATUS,
# PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
# PRIVATELY OWNED RIGHTS.
#
# *****************************************************************************

# output directory
DFFT_MPI_DIR ?= build

# MPI C compiler
DFFT_MPI_CC ?= mpicc

# MPI C++ compiler
DFFT_MPI_CXX ?= mpicxx

# MPI Fortran compiler
DFFT_MPI_FC ?= mpif90

# pre-processor flags
DFFT_MPI_CPPFLAGS ?= -DDFFT_TIMING=2

# C flags
DFFT_MPI_CFLAGS ?= -g -O3 -Wall -Wno-deprecated -std=gnu99

# C++ flags
DFFT_MPI_CXXFLAGS ?= -g -O3 -Wall

# Fortran flags
# -cpp seems to work with GNU and Intel
# though -fpp might be more correct for Intel
DFFT_MPI_FFLAGS ?= -g -O3 -cpp

# linker flags
DFFT_MPI_LDFLAGS ?= 

# additional Fortran linker flags
# sometimes this also needs -lmpi++, -lmpicxx, -lmpi_cxx, etc
DFFT_MPI_FLDFLAGS ?= -lstdc++

# FFTW3
DFFT_FFTW_HOME ?= /work2/07431/tg867305/stampede2/process_skew/fftw-3.3.9/build/
DFFT_FFTW_CPPFLAGS ?= -I$(DFFT_FFTW_HOME)/include
DFFT_FFTW_LDFLAGS ?= -L$(DFFT_FFTW_HOME)/lib -lfftw3 -lm



# these should not usuall require modification
DFFT_MPI_CPPFLAGS += $(DFFT_FFTW_CPPFLAGS)
DFFT_MPI_LDFLAGS += $(DFFT_FFTW_LDFLAGS)



default: nativec utilities

all: nativec utilities fortran

nativec: $(DFFT_MPI_DIR)/TestDfft

fortran: $(DFFT_MPI_DIR)/TestFDfft

utilities: $(DFFT_MPI_DIR)/CheckDecomposition

.PHONY: clean
clean: 
	rm -rf $(DFFT_MPI_DIR) *.mod



$(DFFT_MPI_DIR): 
	mkdir -p $(DFFT_MPI_DIR)

$(DFFT_MPI_DIR)/%.o: %.c | $(DFFT_MPI_DIR)
	$(DFFT_MPI_CC) $(DFFT_MPI_CFLAGS) $(DFFT_MPI_CPPFLAGS) -c -o $@ $<

$(DFFT_MPI_DIR)/%.o: %.cpp | $(DFFT_MPI_DIR)
	$(DFFT_MPI_CXX) $(DFFT_MPI_CXXFLAGS) $(DFFT_MPI_CPPFLAGS) -c -o $@ $<

$(DFFT_MPI_DIR)/%.o: %.f90 | $(DFFT_MPI_DIR)
	$(DFFT_MPI_FC) $(DFFT_MPI_FFLAGS) $(DFFT_MPI_CPPFLAGS) -c -o $@ $<



$(DFFT_MPI_DIR)/TestDfft: $(DFFT_MPI_DIR)/TestDfft.o $(DFFT_MPI_DIR)/distribution.o
	$(DFFT_MPI_CXX) $(DFFT_MPI_CXXFLAGS) -o $@ $^ $(DFFT_MPI_LDFLAGS)



$(DFFT_MPI_DIR)/CheckDecomposition: $(DFFT_MPI_DIR)/CheckDecomposition.o $(DFFT_MPI_DIR)/distribution.o
	$(DFFT_MPI_CC) $(DFFT_MPI_CFLAGS) -o $@ $^ $(DFFT_MPI_LDFLAGS)



$(DFFT_MPI_DIR)/FDfft.o: $(DFFT_MPI_DIR)/FDistribution.o

$(DFFT_MPI_DIR)/TestFDfft.o: TestFDfft.f90 $(DFFT_MPI_DIR)/FDistribution.o $(DFFT_MPI_DIR)/FDfft.o

$(DFFT_MPI_DIR)/TestFDfft: $(DFFT_MPI_DIR)/TestFDfft.o $(DFFT_MPI_DIR)/FDistribution.o $(DFFT_MPI_DIR)/FDfft.o $(DFFT_MPI_DIR)/DistributionC.o $(DFFT_MPI_DIR)/DfftC.o $(DFFT_MPI_DIR)/distribution.o
	$(DFFT_MPI_FC) $(DFFT_MPI_FFLAGS) -o $@ $^ $(DFFT_MPI_LDFLAGS) $(DFFT_MPI_FLDFLAGS)
