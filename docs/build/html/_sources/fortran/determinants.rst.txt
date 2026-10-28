Determinants
================

This module holds all determinant operations
Determinants are stored in bitmasks, 0 for unoccupied orbitals, 1 for occupied. If the length is >64 bits then multiple 64 bit words will be used. e.g. a 104bit length mask will use 2 words for storage.

determinant_occ_to_bit
----------------------

**Source:** `determinants.f90:L8-L30 <https://github.com/violetedwards/SELQIE/blob/selqie-main/src/determinats.f90#L8-L30>`_

``determinant_occ_to_bit(occupations, n_elec, n_spin_orbs, bitstring, n_words)``

Creates a integer bitstring from an array that stores the index of all occupied orbitals.
e.g. (1,2,3,5,7) --> (1110101) --> 117

Parameters
^^^^^^^^^^

``occupations`` : integer, dimension(n_elec)

     Occupation array e.g. (1,2,3,5,7)

``n_elec`` : integer

     Number of electrons in the system

``n_spin_orbs`` : integer

     Number of spin orbitals in the basis set being used

Returns
^^^^^^^

``bitstring`` : integer, dimension(n_words)

     Determinant bitmask array e.g. (117)

``n_words`` : integer

     Number of words being used to store each determinant.


determinant_bit_to_occ
----------------------

**Source:** `determinants.f90:L32-L55 <https://github.com/violetedwards/SELQIE/blob/selqie-main/src/determinats.f90#L32-L55>`_

``determinant_bit_to_occ(bitstring, n_words, n_spin_orbs, occupations, n_elec)``

Inverse of determinant_occ_to_bit

Parameters
^^^^^^^^^^

``n_spin_orbs`` : integer

     Number of spin orbitals in the basis set being used

``bitstring`` : integer, dimension(n_words)

     Determinant bitmask array e.g. (117)

``n_words`` : integer

     Number of words being used to store each determinant.

Returns
^^^^^^^

``occupations`` : integer, dimension(n_elec)

     Occupation array e.g. (1,2,3,5,7)

``n_elec`` : integer

     Number of electrons in the system

reference_determinant
---------------------

**Source:** `determinants.f90:L57-L76 <https://github.com/violetedwards/SELQIE/blob/selqie-main/src/determinats.f90#L57-L76>`_

``reference_determinant(n_orb, n_alpha, n_beta, bitstring, n_words)``

Gets the reference (Hartree-Fock) determinant, lowest lying occupations

Parameters
^^^^^^^^^^

``n_orb`` : integer

     Number of occupied spacial orbitals

``n_alpha`` : integer

     Number of occupied alpha spin orbitals

``n_beta`` : integer

     Number of occupied beta spin orbitals

Returns
^^^^^^^

``bitstring`` : integer(kind=8), dimension(n_words)

     Reference determinant bitstring

``n_words`` : integer

     Number of words used to store reference determinant

determinant_is_occupied
-----------------------

**Source:** `determinants.f90:L78-L93 <https://github.com/violetedwards/SELQIE/blob/selqie-main/src/determinants.f90#L78-L93>`_

``determinant_is_occupied(bitstring, n_words, orbital)``

Checks whether a given spin orbital is occupied in a determinant.

Parameters
^^^^^^^^^^

``bitstring`` : integer(kind=8), dimension(n_words)

     Determinant bitmask array.

``n_words`` : integer

     Number of words used to store the determinant.

``orbital`` : integer

     Index of the spin orbital to test.

Returns
^^^^^^^

``determinant_is_occupied`` : logical

     ``.true.`` if the orbital is occupied, otherwise ``.false.``.


append_connected_det
--------------------

**Source:** `determinants.f90:L95-L117 <https://github.com/violetedwards/SELQIE/blob/selqie-main/src/determinants.f90#L95-L117>`_

``append_connected_det(candidate, connected_dets, n_connected, capacity, n_words)``

Appends a determinant to a dynamically allocated list of connected determinants,
growing the storage array when required.

Parameters
^^^^^^^^^^

``candidate`` : integer(kind=8), dimension(n_words)

     Determinant to append.

``n_words`` : integer

     Number of words used to store each determinant.

``n_connected`` : integer

     Current number of determinants in the list.

``capacity`` : integer

     Current allocated capacity of the determinant array.

Returns
^^^^^^^

``connected_dets`` : integer(kind=8), allocatable, dimension(n_words,:)

     Updated determinant list containing the appended determinant.

``n_connected`` : integer

     Updated number of determinants stored.

``capacity`` : integer

     Updated storage capacity.


get_connected_dets
------------------

**Source:** `determinants.f90:L119-L292 <https://github.com/violetedwards/SELQIE/blob/selqie-main/src/determinants.f90#L119-L292>`_

``get_connected_dets(ref_det, n_words, n_spin_orbs, order, connected_dets, n_connected)``

Generates all determinants connected to a reference determinant by spin-conserving
single or double excitations.

The routine identifies occupied and virtual orbitals and constructs all valid
excitations of the requested order.

Parameters
^^^^^^^^^^

``ref_det`` : integer(kind=8), dimension(n_words)

     Reference determinant.

``n_words`` : integer

     Number of words used to store each determinant.

``n_spin_orbs`` : integer

     Number of spin orbitals in the basis set.

``order`` : integer

     Excitation rank to generate.

     * ``1`` - Single excitations
     * ``2`` - Double excitations

Returns
^^^^^^^

``connected_dets`` : integer(kind=8), allocatable, dimension(n_words,n_connected)

     Array containing all connected determinants.

``n_connected`` : integer

     Number of connected determinants generated.


append_unique_dets
------------------

**Source:** `determinants.f90:L294-L329 <https://github.com/violetedwards/SELQIE/blob/selqie-main/src/determinants.f90#L294-L329>`_

``append_unique_dets(dets, n_dets, candidates, n_candidates, n_words)``

Adds a list of candidate determinants to an existing determinant space,
excluding determinants that are already present.

The determinant array is automatically resized if additional storage is
required.

Parameters
^^^^^^^^^^

``dets`` : integer(kind=8), allocatable, dimension(n_words,:)

     Existing determinant list.

``n_dets`` : integer

     Number of determinants currently stored.

``candidates`` : integer(kind=8), dimension(n_words,n_candidates)

     Candidate determinants to add.

``n_candidates`` : integer

     Number of candidate determinants.

``n_words`` : integer

     Number of words used to store each determinant.

Returns
^^^^^^^

``dets`` : integer(kind=8), allocatable, dimension(n_words,:)

     Updated determinant list containing all unique determinants.

``n_dets`` : integer

     Updated number of determinants in the determinant list.


    

