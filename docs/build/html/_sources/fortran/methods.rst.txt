Methods
=======

This module implements the electronic structure methods available in SELQIE,
including FCI, CISD, and selected CI calculations.

sci
---

**Source:** `methods.f90:L11-L126 <https://github.com/violetedwards/SELQIE/blob/selqie-main/src/methods.f90#L11-L126>`_

``sci(n_orb, n_up, n_down, n_words, n_spin_orbs, n_basis, h, g, sel_method, conv, epsilon, cmin, energy, success, basis, c_out)``

Performs a selected configuration interaction calculation.

Starting from the Hartree-Fock reference determinant, the selected determinant
space is iteratively expanded, diagonalised, and pruned until the energy
change between successive iterations falls below the convergence threshold.

Parameters
^^^^^^^^^^

``n_orb`` : integer

     Number of spatial orbitals.

``n_up`` : integer

     Number of alpha electrons.

``n_down`` : integer

     Number of beta electrons.

``n_words`` : integer

     Number of words used to store each determinant.

``n_spin_orbs`` : integer

     Number of spin orbitals.

``h`` : real(kind=dp), dimension(:,:)

     One-electron integrals.

``g`` : real(kind=dp), dimension(:,:,:,:)

     Two-electron integrals.

``sel_method`` : character(len=*)


     Selection method used to generate new determinants.

``conv`` : real(kind=dp)

     Energy convergence threshold.

``epsilon`` : real(kind=dp)

     Selection threshold parameter.

``cmin`` : real(kind=dp)

     Minimum CI coefficient magnitude retained during pruning.

Returns
^^^^^^^

``energy`` : real(kind=dp)

     Converged SCI energy.

``success`` : logical

     Indicates whether the calculation completed successfully.

``n_basis`` : integer

     Number of determinants in the final selected space.

``basis`` : integer(kind=8), allocatable, dimension(:,:)

     Determinant basis of the final selected space.

``c_out`` : real(kind=dp), allocatable, dimension(:)

     CI coefficients corresponding to the final determinant basis.


prune
-----

**Source:** `methods.f90:L132-L156 <https://github.com/violetedwards/SELQIE/blob/selqie-main/src/methods.f90#L132-L156>`_

``prune(selected_space, n_sel_space, n_words, c, cmin)``

Removes determinants whose CI coefficient magnitude falls below a specified
threshold.

Parameters
^^^^^^^^^^

``selected_space`` : integer(kind=8), allocatable, dimension(:,:)

     Determinant basis to be pruned.

``n_sel_space`` : integer

     Number of determinants currently in the selected space.

``n_words`` : integer

     Number of words used to store each determinant.

``c`` : real(kind=dp), dimension(:)

     CI coefficient vector.

``cmin`` : real(kind=dp)

     Minimum coefficient magnitude for a determinant to be retained.

Returns
^^^^^^^

``selected_space`` : integer(kind=8), allocatable, dimension(:,:)

     Pruned determinant basis.

``n_sel_space`` : integer

     Updated number of determinants after pruning.


cisd
----

**Source:** `methods.f90:L161-L214 <https://github.com/violetedwards/SELQIE/blob/selqie-main/src/methods.f90#L161-L214>`_

``cisd(n_orb, n_up, n_down, n_words, n_spin_orbs, h_int, g_int, e_low, n_basis, success)``

Performs a configuration interaction singles and doubles (CISD)
calculation.

A determinant basis consisting of the reference determinant and all single
and double excitations is generated and diagonalised using the Davidson
solver.

Parameters
^^^^^^^^^^

``n_orb`` : integer

     Number of spatial orbitals.

``n_up`` : integer

     Number of alpha electrons.

``n_down`` : integer

     Number of beta electrons.

``n_words`` : integer

     Number of words used to store each determinant.

``n_spin_orbs`` : integer

     Number of spin orbitals.

``h_int`` : real(kind=dp), dimension(:,:)

     One-electron integrals.

``g_int`` : real(kind=dp), dimension(:,:,:,:)

     Two-electron integrals.

Returns
^^^^^^^

``e_low`` : real(kind=dp)

     Lowest-energy eigenvalue.

``n_basis`` : integer

     Number of determinants in the CISD basis.

``success`` : logical

     Indicates whether the calculation completed successfully.


fci
---

**Source:** `methods.f90:L216-L275 <https://github.com/violetedwards/SELQIE/blob/selqie-main/src/methods.f90#L216-L275>`_

``fci(n_orb, n_up, n_down, n_words, n_spin_orbs, h_int, g_int, e_low, n_basis, success)``

Performs a full configuration interaction calculation.

All possible determinants compatible with the specified numbers of alpha and
beta electrons are generated and diagonalised using the Davidson solver.

Parameters
^^^^^^^^^^

``n_orb`` : integer

     Number of spatial orbitals.

``n_up`` : integer

     Number of alpha electrons.

``n_down`` : integer

     Number of beta electrons.

``n_words`` : integer

     Number of words used to store each determinant.

``n_spin_orbs`` : integer

     Number of spin orbitals.

``h_int`` : real(kind=dp), dimension(:,:)

     One-electron integrals.

``g_int`` : real(kind=dp), dimension(:,:,:,:)

     Two-electron integrals.

Returns
^^^^^^^

``e_low`` : real(kind=dp)

     Lowest-energy eigenvalue.

``n_basis`` : integer

     Number of determinants in the FCI basis.

``success`` : logical

     Indicates whether the calculation completed successfully.


combinations
------------

**Source:** `methods.f90:L277-L306 <https://github.com/violetedwards/SELQIE/blob/selqie-main/src/methods.f90#L277-L306>`_

``combinations(n, k, combos, n_combos)``

Generates all unique combinations of ``k`` elements chosen from ``n``
elements.

This routine is used during FCI basis generation to enumerate all possible
alpha and beta occupation strings.

Parameters
^^^^^^^^^^

``n`` : integer

     Total number of available orbitals.

``k`` : integer

     Number of orbitals to occupy.

Returns
^^^^^^^

``combos`` : integer, allocatable, dimension(k,n_combos)

     Array containing all possible combinations.

``n_combos`` : integer

     Total number of combinations generated.