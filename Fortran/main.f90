!   FILE: main.f90

subroutine solve(method, sel_method, n_orb, n_up, n_down, h_int, g_int, nuc_rep, energy, n_basis)
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
        real(kind = dp) :: nuc_rep, energy, e_low
        logical :: success

!f2py   intent(in) n_orb, n_up, n_down, h_int, g_int, nuc_rep, method, sel_method
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
            !call sci()
        end if

        !output
        if (success) then
            energy = e_low + nuc_rep
        else
            energy = 0.0_dp
        end if

end subroutine solve
!f2py END FILE main.f90