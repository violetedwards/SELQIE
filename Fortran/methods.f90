module methods
    use constants
    use hamiltonian
    implicit none

contains

    subroutine cisd(n_orb, n_up, n_down, n_words, n_spin_orbs, h_int, g_int, e_low, n_basis, success)
        integer, intent(in) :: n_orb, n_up, n_down, n_words, n_spin_orbs
        real(kind = dp), intent(in) :: h_int(:,:)
        real(kind = dp), intent(in) :: g_int(:,:,:,:)
        real(kind = dp), intent(out) :: e_low
        integer, intent(out) :: n_basis
        logical, intent(out) :: success

        interface
            subroutine get_connected_dets(ref_det, n_words, n_spin_orbs, order, connected_dets, n_connected)
                integer, intent(in) :: n_words, n_spin_orbs, order
                integer(kind = 8), intent(in) :: ref_det(n_words)
                integer(kind = 8), allocatable, intent(out) :: connected_dets(:,:)
                integer, intent(out) :: n_connected
            end subroutine get_connected_dets
        end interface

        integer :: n_singles, n_doubles
        real(kind = dp), allocatable :: h_matrix(:,:)
        integer(kind = 8), allocatable :: ref_det(:), det_basis(:,:), singles(:,:), doubles(:,:)

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
        call lowest_eigenvalue(h_matrix, n_basis, e_low, success)
    end subroutine cisd

    subroutine fci(n_orb, n_up, n_down, n_words, n_spin_orbs, h_int, g_int, e_low, n_basis, success)
        integer, intent(in) :: n_orb, n_up, n_down, n_words, n_spin_orbs
        real(kind = dp), intent(in) :: h_int(:,:)
        real(kind = dp), intent(in) :: g_int(:,:,:,:)
        real(kind = dp), intent(out) :: e_low
        integer, intent(out) :: n_basis
        logical, intent(out) :: success

        integer :: n_alpha_strings, n_beta_strings, ia, ib, idet, k, orb
        integer, allocatable :: alpha(:,:), beta(:,:)
        real(kind = dp), allocatable :: h_matrix(:,:)
        integer(kind = 8), allocatable :: det_basis(:,:)

        !every alpha string paired with every beta string gives each determinant once
        call combinations(n_orb, n_up, alpha, n_alpha_strings)
        call combinations(n_orb, n_down, beta, n_beta_strings)

        n_basis = n_alpha_strings * n_beta_strings
        allocate(det_basis(n_words, n_basis))
        det_basis = 0_8
        idet = 0
        do ia = 1, n_alpha_strings
            do ib = 1, n_beta_strings
                idet = idet + 1
                do k = 1, n_up
                    orb = alpha(k, ia)
                    det_basis((orb-1)/64 + 1, idet) = ibset(det_basis((orb-1)/64 + 1, idet), mod(orb-1, 64))
                end do
                do k = 1, n_down
                    orb = n_orb + beta(k, ib)
                    det_basis((orb-1)/64 + 1, idet) = ibset(det_basis((orb-1)/64 + 1, idet), mod(orb-1, 64))
                end do
            end do
        end do

        !build hamiltonian
        allocate(h_matrix(n_basis, n_basis))
        call build_hamiltonian(det_basis, n_basis, n_words, n_spin_orbs, h_int, g_int, h_matrix)

        !diagonalise
        call lowest_eigenvalue(h_matrix, n_basis, e_low, success)
    end subroutine fci

    subroutine combinations(n, k, combos, n_combos)
        integer, intent(in) :: n, k
        integer, allocatable, intent(out) :: combos(:,:)
        integer, intent(out) :: n_combos
        integer :: c(max(k, 1)), i, j, idx

        n_combos = 1
        do i = 1, k
            n_combos = n_combos * (n - k + i) / i
        end do

        allocate(combos(k, n_combos))
        do i = 1, k
            c(i) = i
        end do
        do idx = 1, n_combos
            combos(:, idx) = c(1:k)
            i = k
            do while (i >= 1)
                if (c(i) < n - k + i) exit
                i = i - 1
            end do
            if (i < 1) exit
            c(i) = c(i) + 1
            do j = i + 1, k
                c(j) = c(j-1) + 1
            end do
        end do
    end subroutine combinations

end module methods
