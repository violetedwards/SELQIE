PYTHON ?= python3
F2PY ?= $(PYTHON) -m numpy.f2py
FC ?= gfortran
FFLAGS ?= -Ofast -ffree-line-length-none -ffixed-line-length-none -Wall
LDFLAGS ?=

SRC_DIR := .
BUILD_DIR := $(SRC_DIR)

DAVIDSON_DIR := $(SRC_DIR)/Fortran_Davidson/src
LOCAL_SOURCES := $(SRC_DIR)/constants.f90 $(SRC_DIR)/determinants.f90 $(SRC_DIR)/hamiltonian.f90 $(SRC_DIR)/selected_ci.f90
DAVIDSON_SOURCES := $(DAVIDSON_DIR)/numeric_kinds.f90 $(DAVIDSON_DIR)/lapack_wrapper.f90 \
	$(DAVIDSON_DIR)/array_utils.f90 $(DAVIDSON_DIR)/davidson.f90
FORTRAN_SOURCES := $(LOCAL_SOURCES) $(DAVIDSON_SOURCES)
OBJECTS := $(patsubst $(SRC_DIR)/%.f90,$(BUILD_DIR)/%.o,$(FORTRAN_SOURCES))
DAVIDSON_OBJECTS := $(patsubst $(SRC_DIR)/%.f90,$(BUILD_DIR)/%.o,$(DAVIDSON_SOURCES))
MODULES := $(BUILD_DIR)/*.mod

PYMODULE := selected_ci
EXTENSION := $(PYMODULE)$(shell $(PYTHON) -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX') or '.so')")

.PHONY: all f2py pyscf clean

all: f2py

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(DAVIDSON_DIR)/lapack_wrapper.o: $(DAVIDSON_DIR)/numeric_kinds.o
$(DAVIDSON_DIR)/array_utils.o: $(DAVIDSON_DIR)/lapack_wrapper.o
$(DAVIDSON_DIR)/davidson.o: $(DAVIDSON_DIR)/array_utils.o
$(BUILD_DIR)/hamiltonian.o: $(BUILD_DIR)/constants.o
$(BUILD_DIR)/selected_ci.o: $(BUILD_DIR)/constants.o $(DAVIDSON_DIR)/davidson.o $(BUILD_DIR)/hamiltonian.o

$(PYMODULE).pyf: $(PYMODULE).pyf.in
	cp $< $@

f2py: $(OBJECTS) $(PYMODULE).pyf
	$(F2PY) -c $(PYMODULE).pyf $(FORTRAN_SOURCES) --f90flags='$(FFLAGS)' --f77flags='$(FFLAGS)'

pyscf: f2py

clean:
	rm -f $(OBJECTS) $(BUILD_DIR)/*.mod $(BUILD_DIR)/*.smod $(BUILD_DIR)/$(PYMODULE).pyf $(BUILD_DIR)/$(EXTENSION) $(BUILD_DIR)/$(PYMODULE).so $(BUILD_DIR)/$(PYMODULE).*.so $(BUILD_DIR)/selected_ci.cpython-*.so
