!   FILE: determinants.f90

subroutine determinant_occ_to_bit(occupations, n_spin_orbs, bitarray)
!       
!       Creates a determinant array (occupation_bitarray)
!
        use constants
        implicit none

        integer :: n_spin_orbs
        integer, allocatable :: occupations(:)
        integer, allocatable :: bitarray(:)

!f2py   intent(in) occupations, n_spin_orbs
!f2py   intent(out) bitarray
!f2py   depends(occupationsd) bitarray

        allocate(occupations(n_spin_orbs))
        allocate(bitarray(n_spin_orbs))

        

