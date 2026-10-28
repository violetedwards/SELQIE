Selection
=========

This module contains determinant selection algorithms used by selected CI
methods.

add_selected_dets
-----------------

**Source:** `selection.f90:L9-L19 <https://github.com/violetedwards/SELQIE/blob/selqie-main/src/selection.f90#L9-L19>`_

``add_selected_dets(sel_method, selected_space, n_sel_space, n_words, old_e, old_c, n_spin_orbs, h, g, epsilon)``

Dispatch routine for determinant selection methods.

Based on the selected method, the appropriate selection algorithm is applied
to expand the current determinant space.

Parameters
^^^^^^^^^^

``sel_method`` : character(len=*)


     Selection method to apply.

``selected_space`` : integer(kind=8), allocatable, dimension(:,:)

     Current selected determinant space.

``n_sel_space`` : integer

     Number of determinants currently in the selected space.

``n_words`` : integer

     Number of words used to store each determinant.

``old_e`` : real(kind=dp)

     Current variational energy.

``old_c`` : real(kind=dp), dimension(:)

     Current CI coefficients.

``n_spin_orbs`` : integer

     Number of spin orbitals.

``h`` : real(kind=dp), dimension(:,:)

     One-electron integrals.

``g`` : real(kind=dp), dimension(:,:,:,:)

     Two-electron integrals.

``epsilon`` : real(kind=dp)

     Selection threshold.

Returns
^^^^^^^

``selected_space`` : integer(kind=8), allocatable, dimension(:,:)

     Expanded determinant space.

``n_sel_space`` : integer

     Updated size of the selected space.


hci
---

**Source:** `selection.f90:L22-L119 <https://github.com/violetedwards/SELQIE/blob/selqie-main/src/selection.f90#L22-L119>`_

``hci(selected_space, n_sel_space, n_words, old_c, n_spin_orbs, h, g, epsilon)``

Performs determinant selection using the Heat-Bath Configuration Interaction
(HCI) algorithm.

For each determinant in the current selected space, all connected single and
double excitations are generated. Candidate determinants are retained if

.. math::

   |H_{ij} c_i| > \epsilon

where :math:`H_{ij}` is the Hamiltonian coupling between determinants and
:math:`c_i` is the CI coefficient of the parent determinant.

New determinants passing this screening criterion are added to the selected
space.

Parameters
^^^^^^^^^^

``selected_space`` : integer(kind=8), allocatable, dimension(:,:)

     Current selected determinant space.

``n_sel_space`` : integer

     Number of determinants currently in the selected space.

``n_words`` : integer

     Number of words used to store each determinant.

``old_c`` : real(kind=dp), dimension(:)

     CI coefficients from the previous variational calculation.

``n_spin_orbs`` : integer

     Number of spin orbitals.

``h`` : real(kind=dp), dimension(:,:)

     One-electron integrals.

``g`` : real(kind=dp), dimension(:,:,:,:)

     Two-electron integrals.

``epsilon`` : real(kind=dp)

     HCI selection threshold.

Returns
^^^^^^^

``selected_space`` : integer(kind=8), allocatable, dimension(:,:)

     Expanded determinant space containing all selected determinants.

``n_sel_space`` : integer

     Updated number of determinants in the selected space.