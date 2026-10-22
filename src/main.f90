!   FILE: main.f90

module selqie_api

    use constants
    use methods

    implicit none

    integer(kind=8), allocatable, save :: stored_basis(:,:)
    real(dp), allocatable, save :: stored_coeffs(:)

contains

subroutine solve(method, sel_method, n_orb, n_up, n_down, h_int, g_int, nuc_rep, conv, epsilon, cmin, energy, n_basis)
!
!       MAIN SOLVING ROUTINE
!
        use constants
        use methods
        implicit none

        integer :: n_orb, n_up, n_down
        integer :: n_words, n_spin_orbs, n_basis
        character(len=*) :: method, sel_method
        real(kind = dp) :: h_int(:,:)
        real(kind = dp) :: g_int(:,:,:,:)
        real(kind = dp) :: nuc_rep, energy, e_low, conv, epsilon, cmin
        logical :: success
        real(kind = dp), allocatable:: c_out(:)
        integer(kind = 8), allocatable :: basis(:,:)

!f2py   intent(in) n_orb, n_up, n_down, h_int, g_int, nuc_rep, method, sel_method, conv, epsilon, cmin
!f2py   intent(out) energy, n_basis


        !get determinant words ready
        n_spin_orbs = size(h_int, 1)
        n_words = ceiling(real(n_spin_orbs) / 64.0_dp)


        n_basis = 0
        success = .false.
        if (method == "cisd") then
            call cisd(n_orb, n_up, n_down, n_words, n_spin_orbs, h_int, g_int, e_low, n_basis, success)
        else if (method == "fci") then
            call fci(n_orb, n_up, n_down, n_words, n_spin_orbs, h_int, g_int, e_low, n_basis, success)
        else if (method == "sci") then
            call sci(n_orb, n_up, n_down, n_words, n_spin_orbs, &
                n_basis, h_int, g_int, sel_method, conv, epsilon, cmin, e_low, success, basis, c_out)

            if (allocated(stored_basis)) deallocate(stored_basis)
            if (allocated(stored_coeffs)) deallocate(stored_coeffs)

            call move_alloc(basis, stored_basis)
            call move_alloc(c_out, stored_coeffs)
        end if

        !output
        if (success) then
            energy = e_low + nuc_rep
        else
            energy = 0.0_dp
        end if

end subroutine solve

end module selqie_api
!f2py END FILE main.f90