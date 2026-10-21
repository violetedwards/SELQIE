module methods
    use constants
    use hamiltonian
    use selection
    use diagonalisers
    use determinants
    implicit none

contains

    subroutine sci(n_orb, n_up, n_down, n_words, n_spin_orbs, n_basis, h, g, sel_method, conv, epsilon, cmin, energy, success)
        integer, intent(in) :: n_orb, n_up, n_down, n_words, n_spin_orbs
        real(kind = dp), intent(in) :: h(:,:), conv, epsilon, cmin
        real(kind = dp), intent(in) :: g(:,:,:,:)
        real(kind = dp), intent(out) :: energy
        character(len=*), intent(in) :: sel_method
        integer, intent(out) :: n_basis
        logical, intent(out) :: success

        !local variables, allocatable
        integer(kind = 8), allocatable :: selected_space(:,:), ref_det(:)
        real(kind = dp), allocatable :: sel_hamiltonian(:,:), sel_c(:), old_c(:), new_c(:)
        
        !local variables, non-allocatable
        logical :: converged
        integer :: iteration, n_sel_space
        real(kind = dp) :: sel_e, old_e, new_e
        real :: t1, t2

        !get HF reference determinant
        n_basis = 0
        allocate(ref_det(n_words))
        ref_det = 0_8
        call reference_determinant(n_orb, n_up, n_down, ref_det, n_words)

        !get initial selected space
        allocate(selected_space(n_words,1))
        selected_space(:,1) = ref_det
        n_sel_space = 1

        !main loop initial vars
        converged = .false.
        iteration = 0

        !main loop
        do while (converged .eqv. .false.)

            write(*, '(A,I0)') "Iteration ", iteration
            write(*, '(A,I0)') "Current Sel Space ", n_sel_space

            allocate(sel_c(n_sel_space))


            call davidson_lowest( &
                selected_space, &
                n_sel_space, &
                n_words, &
                n_spin_orbs, &
                h, &
                g, &
                sel_e, &
                sel_c)

            success = .true.

            !rename variables
            old_e = sel_e
            allocate(old_c(n_sel_space))
            old_c = sel_c
            deallocate(sel_c)

            !apply selection algorithm to add to selected space

            call add_selected_dets(sel_method, selected_space, n_sel_space, n_words, old_e, old_c, n_spin_orbs, h, g, epsilon)

            allocate(sel_c(n_sel_space))


            call davidson_lowest( &
                selected_space, &
                n_sel_space, &
                n_words, &
                n_spin_orbs, &
                h, &
                g, &
                sel_e, &
                sel_c)



            success = .true.

            !rename variables
            new_e = sel_e
            allocate(new_c(n_sel_space))
            new_c = sel_c
            deallocate(sel_c)
            
            !check convergence
            if (abs(old_e-new_e) < conv) then
                converged = .true.
            end if


            !prune if not converged
            if (converged .eqv. .false.) then

                call prune(selected_space, n_sel_space, n_words, new_c, cmin)

            end if

            !getting ready for next loop
            deallocate(old_c, new_c)
            iteration = iteration + 1
 
            

        end do

        !outputs
        write(*, '(A,I0)') "Energy converged"
        write(*, '(A,I0)') "----------------------------------------------------------"
        energy = new_e
        n_basis = n_sel_space



    end subroutine




    subroutine prune(selected_space, n_sel_space, n_words, c, cmin)
        integer(kind = 8), allocatable, intent(inout) :: selected_space(:,:)
        integer, intent(inout) :: n_sel_space
        integer, intent(in) :: n_words
        real(kind = dp), intent(in) :: c(:), cmin

        integer(kind = 8), allocatable :: kept(:,:)
        integer :: i, n_kept

        allocate(kept(n_words, count(abs(c) >= cmin)))
        n_kept = 0
        do i = 1, n_sel_space
            if (abs(c(i)) >= cmin) then
                n_kept = n_kept + 1
                kept(:, n_kept) = selected_space(:, i)
            end if
        end do

        call move_alloc(kept, selected_space)
        n_sel_space = n_kept

    end subroutine prune




    subroutine cisd(n_orb, n_up, n_down, n_words, n_spin_orbs, h_int, g_int, e_low, n_basis, success)
        integer, intent(in) :: n_orb, n_up, n_down, n_words, n_spin_orbs
        real(kind = dp), intent(in) :: h_int(:,:)
        real(kind = dp), intent(in) :: g_int(:,:,:,:)
        real(kind = dp), intent(out) :: e_low
        integer, intent(out) :: n_basis
        logical, intent(out) :: success

        

        integer :: n_singles, n_doubles
        real(kind = dp), allocatable :: eigvec(:)   
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

        allocate(eigvec(n_basis))

        call davidson_lowest( &
            det_basis, &
            n_basis, &
            n_words, &
            n_spin_orbs, &
            h_int, &
            g_int, &
            e_low, &
            eigvec)

        success = .true.
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
        real(kind = dp), allocatable :: eigvec(:)
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

        allocate(eigvec(n_basis))

        call davidson_lowest( &
            det_basis, &
            n_basis, &
            n_words, &
            n_spin_orbs, &
            h_int, &
            g_int, &
            e_low, &
            eigvec)

        success = .true.
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
