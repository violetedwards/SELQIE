Diagonalisers
=============

This module implements Davidson diagonalisation and
subspace expansion for selected CI calculations.


lowest_eigenvalue
-----------------

**Source:** `diagonalisers.f90:L14-L64 <https://github.com/violetedwards/SELQIE/blob/selqie-main/src/diagonalisers.f90#L14-L64>`_

``lowest_eigenvalue(h_matrix, n_basis, e_low, success, eigenvector)``

Uses LAPACK to find the lowest eigenvalue of a Hamiltonian matrix.

Parameters
^^^^^^^^^^

``h_matrix`` : real(kind=dp), dimension(n_basis, n_basis)

    Hamiltonian matrix to be diagonalised.

``n_basis`` : integer

    Number of determinants in the basis.

Returns
^^^^^^^

``e_low`` : real(kind=dp)

    Lowest eigenvalue of the Hamiltonian.

``success`` : logical

    Indicates whether diagonalisation completed successfully.

``eigenvector`` : real(kind=dp), dimension(n_basis)

    Eigenvector corresponding to ``e_low``.



build_connectivity
------------------

**Source:** `diagonalisers.f90:L66-L177 <https://github.com/violetedwards/SELQIE/blob/selqie-main/src/diagonalisers.f90#L66-L177>`_

``build_connectivity(conn,basis, n_basis, n_words, n_spin_orbs, h, g, hash)``

Builds the connectivity graphs of the hamitlonian, describing how elements interact. used in apply_hamiltonian

Parameters
^^^^^^^^^^

``basis`` : integer(kind=8), dimension(n_words,n_basis)

    The determinant basis set being used.

``n_basis`` : integer

    Number of determinants in the basis.

``n_words`` : integer

    Number of words needed to store each determinant bitstring

``n_spin_orbs`` : integer

    Number of spin orbitals in each determinant.

``h`` : real(kind=dp), dimension(n_spin_orbs,n_spin_orbs)

    One-electron integrals in spin orbital basis

``g`` : real(kind=dp), dimension(n_spin_orbs,n_spin_orbs,n_spin_orbs,n_spin_orbs)

    Two-electron integrals in spin orbital basis

``hash``: type(det_hash_t)

    Excitation lookup hash

Returns
^^^^^^^

``conn`` : type(connection_list_t), dimension(n_basis)

    Connectivity graph


davidson_lowest
---------------

**Source:** `diagonalisers.f90:L180-L371 <https://github.com/violetedwards/SELQIE/blob/selqie-main/src/diagonalisers.f90#L180-L371>`_

``davidson_lowest(basis, n_basis, n_words, n_spin_orbs, h, g, e_low, eigenvector)``

Performs davidson diagonalisation and outputs the lowest eigenvalue and eigenvector.

Parameters
^^^^^^^^^^
``basis`` : integer(kind=8), dimension(n_words,n_basis)

    The determinant basis set being used.

``n_basis`` : integer

    Number of determinants in the basis.

``n_words`` : integer

    Number of words needed to store each determinant bitstring

``n_spin_orbs`` : integer

    Number of spin orbitals in each determinant.

``h`` : real(kind=dp), dimension(n_spin_orbs,n_spin_orbs)

    One-electron integrals in spin orbital basis

``g`` : real(kind=dp), dimension(n_spin_orbs,n_spin_orbs,n_spin_orbs,n_spin_orbs)

    Two-electron integrals in spin orbital basis

Returns
^^^^^^^

``e_low`` : real(kind=dp)

    Lowest eigenvalue of the Hamiltonian.

``eigenvector`` : real(kind=dp), dimension(n_basis)

    Eigenvector corresponding to ``e_low``.


apply_hamiltonian_sparse
------------------------

**Source:** `diagonalisers.f90:L373-L410 <https://github.com/violetedwards/SELQIE/blob/selqie-main/src/diagonalisers.f90#L373-L410>`_

``apply_hamiltonian_sparse(x, y, n_basis, conn, diagH)``

applies non-explitily stored hamitlonitan to a general vector, x. Used in davidson diagonalisation.

Parameters
^^^^^^^^^^
``x`` : real(kind=dp), dimension(n_basis)

    general trial vector that the Hamiltonian is to be applied to

``n_basis`` : integer

    Number of determinants in the basis.

``conn`` : type(connection_list_t), dimension(n_basis)

    Connectivity graph

``diagH`` : real(kind=dp), dimension(n_basis)

    Diagonal part of hamitlonian

Returns
^^^^^^^
``y`` : real(kind=dp), dimension(n_basis)

    Hamiltonian applied to trial vector ``x``
