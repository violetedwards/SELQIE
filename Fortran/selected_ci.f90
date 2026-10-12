!   FILE: selected_ci.f90

module constants
    integer, parameter :: dp = selected_real_kind(15,307)
end module constants

subroutine solve(n_orb, n_elec, spin, n_up, n_down, h_int, g_int, nuc_rep, method, conv, cmin, hci_epsilon, energy)
!   
!       MAIN SOLVING ROUTINE
!   
        use constants        
        implicit none

        integer :: n_orb, n_elec, spin, n_up, n_down
        real(kind = dp), allocatable :: h_int(:,:), g_int(:,:,:,:)
        real(kind = dp) :: nuc_rep, conv, cmin, hci_epsilon, energy
        character :: method

!f2py   intent(in) n_orb, n_elec, spin, n_up, n_down, h_int, g_int, nuc_rep, method, conv, cmin, hci_epsilon
!f2py   intent(out) energy
!f2py   depend(n_orb, n_elec, spin, n_up, n_down, h_int, g_int, nuc_rep, method, conv, cmin, hci_epsilon) energy

        allocate(h_int(2*n_orb,2*n_orb))
        allocate(g_int(2*n_orb,2*n_orb,2*n_orb,2*n_orb))


end subroutine solve
!f2py END FILE selected_ci.f90