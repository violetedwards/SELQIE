PYTHON ?= python3
F2PY ?= $(PYTHON) -m numpy.f2py
FC ?= gfortran
FFLAGS ?= -Ofast -ffree-line-length-none -ffixed-line-length-none -Wall
LDFLAGS ?=

SRC_DIR := .
BUILD_DIR := $(SRC_DIR)

FORTRAN_SOURCES := $(SRC_DIR)/constants.f90 $(SRC_DIR)/determinants.f90 $(SRC_DIR)/hamiltonian.f90 $(SRC_DIR)/selected_ci.f90
OBJECTS := $(patsubst $(SRC_DIR)/%.f90,$(BUILD_DIR)/%.o,$(FORTRAN_SOURCES))
MODULES := $(BUILD_DIR)/*.mod

PYMODULE := selected_ci
EXTENSION := $(PYMODULE)$(shell $(PYTHON) -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX') or '.so')")

.PHONY: all f2py clean

all: f2py

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

f2py: $(OBJECTS)
	$(F2PY) -c -m $(PYMODULE) $(FORTRAN_SOURCES) only: solve : --f90flags='$(FFLAGS)' --f77flags='$(FFLAGS)'

clean:
	rm -f $(BUILD_DIR)/*.o $(BUILD_DIR)/*.mod $(BUILD_DIR)/$(EXTENSION) $(BUILD_DIR)/$(PYMODULE).so $(BUILD_DIR)/$(PYMODULE).*.so $(BUILD_DIR)/selected_ci.cpython-*.so
