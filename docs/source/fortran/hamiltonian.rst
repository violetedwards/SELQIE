Hamiltonian
===========

This module contains routines for analysing determinants and constructing
Slater-Condon Hamiltonian matrix elements.

get_occupied_orbitals
---------------------

**Source:** `hamiltonian.f90:L9-L23 <https://github.com/violetedwards/SELQIE/blob/selqie-main/src/hamiltonian.f90#L9-L23>`_

``get_occupied_orbitals(det, n_words, n_spin_orbs, occ, n_occ)``

Extracts all occupied spin orbitals from a determinant bitstring.

Parameters
^^^^^^^^^^

``det`` : integer(kind=8), dimension(n_words)

     Determinant bitstring.

``n_words`` : integer

     Number of words used to store the determinant.

``n_spin_orbs`` : integer

     Number of spin orbitals in the basis set.

Returns
^^^^^^^

``occ`` : integer, dimension(n_spin_orbs)

     Array containing occupied orbital indices.

``n_occ`` : integer

     Number of occupied orbitals found.


compare_dets
------------

**Source:** `hamiltonian.f90:L25-L42 <https://github.com/violetedwards/SELQIE/blob/selqie-main/src/hamiltonian.f90#L25-L42>`_

``compare_dets(det1, n_words1, det2, n_words2, n_spin_orbs, added, removed, n_added, n_removed)``

Compares two determinants and identifies which orbitals have been added and
removed between them.

Parameters
^^^^^^^^^^

``det1`` : integer(kind=8), dimension(n_words1)

     First determinant.

``det2`` : integer(kind=8), dimension(n_words2)

     Second determinant.

``n_words1`` : integer

     Number of words used to store the first determinant.

``n_words2`` : integer

     Number of words used to store the second determinant.

``n_spin_orbs`` : integer

     Number of spin orbitals in the basis set.

Returns
^^^^^^^

``added`` : integer, dimension(n_spin_orbs)

     Orbitals occupied in ``det2`` but not ``det1``.

``removed`` : integer, dimension(n_spin_orbs)

     Orbitals occupied in ``det1`` but not ``det2``.

``n_added`` : integer

     Number of added orbitals.

``n_removed`` : integer

     Number of removed orbitals.


count_bits_before
-----------------

**Source:** `hamiltonian.f90:L44-L56 <https://github.com/violetedwards/SELQIE/blob/selqie-main/src/hamiltonian.f90#L44-L56>`_

``count_bits_before(bitstring, n_words, orbital)``

Counts the number of occupied orbitals before a specified orbital index.

Parameters
^^^^^^^^^^

``bitstring`` : integer(kind=8), dimension(n_words)

     Determinant bitstring.

``n_words`` : integer

     Number of words used to store the determinant.

``orbital`` : integer

     Orbital index.

Returns
^^^^^^^

``count_bits_before`` : integer

     Number of occupied orbitals with indices lower than ``orbital``.


set_orbital
-----------

**Source:** `hamiltonian.f90:L58-L71 <https://github.com/violetedwards/SELQIE/blob/selqie-main/src/hamiltonian.f90#L58-L71>`_

``set_orbital(bitstring, n_words, orbital)``

Returns a determinant bitstring with the specified orbital occupied.

Parameters
^^^^^^^^^^

``bitstring`` : integer(kind=8), dimension(n_words)

     Determinant bitstring.

``n_words`` : integer

     Number of words used to store the determinant.

``orbital`` : integer

     Orbital to occupy.

Returns
^^^^^^^

``new_bitstring`` : integer(kind=8), dimension(n_words)

     Modified determinant bitstring.


clear_orbital
-------------

**Source:** `hamiltonian.f90:L73-L86 <https://github.com/violetedwards/SELQIE/blob/selqie-main/src/hamiltonian.f90#L73-L86>`_

``clear_orbital(bitstring, n_words, orbital)``

Returns a determinant bitstring with the specified orbital unoccupied.

Parameters
^^^^^^^^^^

``bitstring`` : integer(kind=8), dimension(n_words)

     Determinant bitstring.

``n_words`` : integer

     Number of words used to store the determinant.

``orbital`` : integer

     Orbital to remove.

Returns
^^^^^^^

``new_bitstring`` : integer(kind=8), dimension(n_words)

     Modified determinant bitstring.


excitation_phase
----------------

**Source:** `hamiltonian.f90:L88-L111 <https://github.com/violetedwards/SELQIE/blob/selqie-main/src/hamiltonian.f90#L88-L111>`_

``excitation_phase(det, n_words, added, n_added, removed, n_removed)``

Computes the fermionic phase factor associated with an excitation between two
determinants.

Parameters
^^^^^^^^^^

``det`` : integer(kind=8), dimension(n_words)

     Reference determinant.

``n_words`` : integer

     Number of words used to store the determinant.

``added`` : integer, dimension(n_added)

     Orbitals added to the determinant.

``removed`` : integer, dimension(n_removed)

     Orbitals removed from the determinant.

``n_added`` : integer

     Number of added orbitals.

``n_removed`` : integer

     Number of removed orbitals.

Returns
^^^^^^^

``excitation_phase`` : integer

     Fermionic phase factor, ``+1`` or ``-1``.


get_matrix_element
------------------

**Source:** `hamiltonian.f90:L113-L181 <https://github.com/violetedwards/SELQIE/blob/selqie-main/src/hamiltonian.f90#L113-L181>`_

``get_matrix_element(column_det, n_words_col, row_det, n_words_row, n_spin_orbs, h, g)``

Computes a Hamiltonian matrix element between two determinants using the
Slater-Condon rules.

Parameters
^^^^^^^^^^

``column_det`` : integer(kind=8), dimension(n_words_col)

     Column determinant.

``row_det`` : integer(kind=8), dimension(n_words_row)

     Row determinant.

``n_words_col`` : integer

     Number of words used to store the column determinant.

``n_words_row`` : integer

     Number of words used to store the row determinant.

``n_spin_orbs`` : integer

     Number of spin orbitals in the basis set.

``h`` : real(kind=dp), dimension(n_spin_orbs,n_spin_orbs)

     One-electron integrals.

``g`` : real(kind=dp), dimension(n_spin_orbs,n_spin_orbs,n_spin_orbs,n_spin_orbs)

     Two-electron integrals.

Returns
^^^^^^^

``get_matrix_element`` : real(kind=dp)

     Hamiltonian matrix element between the two determinants.


build_hamiltonian
-----------------

**Source:** `hamiltonian.f90:L183-L202 <https://github.com/violetedwards/SELQIE/blob/selqie-main/src/hamiltonian.f90#L183-L202>`_

``build_hamiltonian(det_basis_set, n_det, n_words, n_spin_orbs, h, g, hamiltonian)``

Constructs the full Hamiltonian matrix for a determinant basis.

Parameters
^^^^^^^^^^

``det_basis_set`` : integer(kind=8), dimension(n_words,n_det)

     Determinant basis set.

``n_det`` : integer

     Number of determinants in the basis.

``n_words`` : integer

     Number of words used to store each determinant.

``n_spin_orbs`` : integer

     Number of spin orbitals in the basis set.

``h`` : real(kind=dp), dimension(n_spin_orbs,n_spin_orbs)

     One-electron integrals.

``g`` : real(kind=dp), dimension(n_spin_orbs,n_spin_orbs,n_spin_orbs,n_spin_orbs)

     Two-electron integrals.

Returns
^^^^^^^

``hamiltonian`` : real(kind=dp), dimension(n_det,n_det)

     Hamiltonian matrix in the determinant basis.