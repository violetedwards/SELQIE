!   FILE: selected_ci.f90

subroutine solve(n_orb, n_up, n_down, h_int, g_int, nuc_rep, energy)
!
!       MAIN SOLVING ROUTINE
!
        use constants
        use hamiltonian
        implicit none

        interface
            subroutine get_connected_dets(ref_det, n_words, n_spin_orbs, order, connected_dets, n_connected)
                integer, intent(in) :: n_words, n_spin_orbs, order
                integer(kind = 8), intent(in) :: ref_det(n_words)
                integer(kind = 8), allocatable, intent(out) :: connected_dets(:,:)
                integer, intent(out) :: n_connected
            end subroutine get_connected_dets
        end interface

        integer :: n_orb, n_up, n_down
        integer :: n_words, n_spin_orbs, n_basis, info, lwork, liwork, m, il, iu, n_singles, n_doubles
        real(kind = dp) :: h_int(:,:)
        real(kind = dp) :: g_int(:,:,:,:)
        real(kind = dp) :: nuc_rep, energy
        real(kind = dp) :: vl, vu, abstol
        real(kind = dp), allocatable :: h_matrix(:,:), eval(:), z(:,:), work(:)
        integer, allocatable :: isuppz(:), iwork(:)
        integer(kind = 8), allocatable :: ref_det(:), det_basis(:,:), singles(:,:), doubles(:,:)

!f2py   intent(in) n_orb, n_up, n_down, h_int, g_int, nuc_rep
!f2py   intent(out) energy

        !get determinant words ready
        n_spin_orbs = size(h_int, 1)
        n_words = ceiling(real(n_spin_orbs) / 64.0_dp)

        !get HF reference determinant
        allocate(ref_det(n_words))
        ref_det = 0_8
        call reference_determinant(n_orb, n_up, n_down, ref_det, n_words)

        !get CISD det basis
        call get_connected_dets(ref_det, n_words, n_spin_orbs, 1, singles, n_singles)
        call get_connected_dets(ref_det, n_words, n_spin_orbs, 2, doubles, n_doubles)
        
        n_basis = 1 + n_singles + n_doubles
        allocate(det_basis(n_words, n_basis))
        det_basis(:, 1) = ref_det
        if (n_singles > 0) det_basis(:, 2:1+n_singles) = singles(:, 1:n_singles)
        if (n_doubles > 0) det_basis(:, 2+n_singles:n_basis) = doubles(:, 1:n_doubles)

        !build hamiltonian
        allocate(h_matrix(n_basis, n_basis))
        call build_hamiltonian(det_basis, n_basis, n_words, n_spin_orbs, h_int, g_int, h_matrix)

        !diagonalise
        allocate(eval(n_basis))
        allocate(z(n_basis, 1))
        allocate(isuppz(2*n_basis))
        allocate(work(1))
        allocate(iwork(1))
        lwork = -1
        liwork = -1
        vl = 0.0_dp
        vu = 0.0_dp
        abstol = 0.0_dp
        il = 1
        iu = 1
        call dsyevr('V', 'I', 'U', n_basis, h_matrix, n_basis, vl, vu, il, iu, abstol, m, &
                eval, z, n_basis, isuppz, work, lwork, iwork, liwork, info)
        lwork = max(1, int(work(1)))
        liwork = max(1, iwork(1))
        deallocate(work)
        deallocate(iwork)
        allocate(work(lwork))
        allocate(iwork(liwork))
        call dsyevr('V', 'I', 'U', n_basis, h_matrix, n_basis, vl, vu, il, iu, abstol, m, &
                eval, z, n_basis, isuppz, work, lwork, iwork, liwork, info)

        !output
        if (info == 0 .and. m == 1) then
            energy = eval(1) + nuc_rep
        else
            energy = 0.0_dp
        end if

end subroutine solve
!f2py END FILE selected_ci.f90