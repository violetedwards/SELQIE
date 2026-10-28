Hash
====

This module provides a hash table implementation for rapid determinant lookup
within a determinant basis.

det_hash_t
----------

**Source:** `hash.f90:L5-L9 <https://github.com/violetedwards/SELQIE/blob/selqie-main/src/hash.f90#L5-L9>`_

``type(det_hash_t)``

Hash table structure used to map determinant bitstrings to their indices within
a determinant basis.

Components
^^^^^^^^^^

``table_size`` : integer

     Number of buckets in the hash table.

``head`` : integer, allocatable, dimension(:)

     Head pointer for each bucket chain.

``next`` : integer, allocatable, dimension(:)

     Linked-list pointers used for collision resolution.


hash_det
--------

**Source:** `hash.f90:L13-L35 <https://github.com/violetedwards/SELQIE/blob/selqie-main/src/hash.f90#L13-L35>`_

``hash_det(det, n_words, table_size)``

Computes a hash value for a determinant bitstring using an FNV-style hash
function.

Parameters
^^^^^^^^^^

``det`` : integer(kind=8), dimension(n_words)

     Determinant bitstring.

``n_words`` : integer

     Number of words used to store the determinant.

``table_size`` : integer

     Number of buckets in the hash table.

Returns
^^^^^^^

``hash_det`` : integer(kind=8)

     Hash bucket index in the range ``1`` to ``table_size``.


build_det_hash
--------------

**Source:** `hash.f90:L39-L77 <https://github.com/violetedwards/SELQIE/blob/selqie-main/src/hash.f90#L39-L77>`_

``build_det_hash(hash_table, basis, n_words, n_basis)``

Constructs a hash table for a determinant basis.

Each determinant is assigned to a bucket according to its hash value and stored
using a linked-list collision scheme.

Parameters
^^^^^^^^^^

``basis`` : integer(kind=8), dimension(n_words,n_basis)

     Determinant basis set.

``n_words`` : integer

     Number of words used to store each determinant.

``n_basis`` : integer

     Number of determinants in the basis.

Returns
^^^^^^^

``hash_table`` : type(det_hash_t)

     Constructed determinant hash table.


lookup_det_index
----------------

**Source:** `hash.f90:L81-L116 <https://github.com/violetedwards/SELQIE/blob/selqie-main/src/hash.f90#L81-L116>`_

``lookup_det_index(det, basis, n_words, hash_table)``

Searches the determinant basis for a specified determinant using the hash table.

Parameters
^^^^^^^^^^

``det`` : integer(kind=8), dimension(n_words)

     Determinant to locate.

``basis`` : integer(kind=8), dimension(:,:)

     Determinant basis set.

``n_words`` : integer

     Number of words used to store each determinant.

``hash_table`` : type(det_hash_t)

     Determinant hash table.

Returns
^^^^^^^

``lookup_det_index`` : integer

     Index of the determinant within the basis set.

     Returns ``0`` if the determinant is not found.


destroy_det_hash
----------------

**Source:** `hash.f90:L120-L132 <https://github.com/violetedwards/SELQIE/blob/selqie-main/src/hash.f90#L120-L132>`_

``destroy_det_hash(hash_table)``

Deallocates memory associated with a determinant hash table.

Parameters
^^^^^^^^^^

``hash_table`` : type(det_hash_t)

     Hash table to destroy.

Returns
^^^^^^^

``hash_table`` : type(det_hash_t)

     Hash table with all storage deallocated and ``table_size`` reset to zero.