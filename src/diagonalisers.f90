module diagonalisers
    use constants
    use hamiltonian
    use determinants
    use hash
    implicit none

    type connection_list_t
    integer, allocatable :: idx(:)
    real(dp), allocatable :: hij(:)
    end type connection_list_t

contains
    subroutine lowest_eigenvalue(h_matrix, n_basis, e_low, success, eigenvector)
        !solve for lowest eigenvalue of explicity stored hamiltonian

        integer, intent(in) :: n_basis
        real(kind=dp), intent(inout) :: h_matrix(n_basis, n_basis)  ! destroyed by dsyevr
        real(kind=dp), intent(out) :: e_low
        logical, intent(out) :: success
        real(kind=dp), intent(out), optional :: eigenvector(n_basis)
        integer :: info, lwork, liwork, m, il, iu
        real(kind=dp) :: vl, vu, abstol
        real(kind=dp), allocatable :: eval(:), z(:,:), work(:)
        integer, allocatable :: isuppz(:), iwork(:)

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

        success = (info == 0 .and. m == 1)
        if (success) then
            e_low = eval(1)
        else
            e_low = 0.0_dp
        end if

        if (present(eigenvector)) then
            if (success) then
                eigenvector = z(:, 1)
            else
                eigenvector = 0.0_dp
            end if
        end if
    end subroutine lowest_eigenvalue

    subroutine build_connectivity( &
            conn, &
            basis, n_basis, n_words, n_spin_orbs, &
            h, g, &
            hash)

        !builds connectivity graph of the hamiltonian, used in apply_hamitlonian
        implicit none

        type(connection_list_t), allocatable, intent(out) :: conn(:)

        integer, intent(in) :: n_basis
        integer, intent(in) :: n_words
        integer, intent(in) :: n_spin_orbs

        integer(kind=8), intent(in) :: basis(n_words,n_basis)

        real(dp), intent(in) :: h(:,:)
        real(dp), intent(in) :: g(:,:,:,:)

        type(det_hash_t), intent(in) :: hash

        integer(kind=8), allocatable :: connected(:,:)

        integer, allocatable :: tmp_idx(:)
        real(dp), allocatable :: tmp_hij(:)

        integer, allocatable :: grow_idx(:)
        real(dp), allocatable :: grow_hij(:)

        integer :: i,j
        integer :: idx
        integer :: order
        integer :: n_connected
        integer :: n_keep
        integer :: capacity

        allocate(conn(n_basis))

        do i = 1,n_basis

            capacity = 128

            allocate(tmp_idx(capacity))
            allocate(tmp_hij(capacity))

            n_keep = 0

            do order = 1,2

                call get_connected_dets( &
                    basis(:,i), &
                    n_words, &
                    n_spin_orbs, &
                    order, &
                    connected, &
                    n_connected)

                do j = 1,n_connected

                    idx = lookup_det_index( &
                        connected(:,j), &
                        basis, &
                        n_words, &
                        hash)

                    if(idx /= 0) then

                        n_keep = n_keep + 1

                        if(n_keep > capacity) then

                            capacity = 2*capacity

                            allocate(grow_idx(capacity))
                            allocate(grow_hij(capacity))

                            grow_idx(1:n_keep-1) = tmp_idx(1:n_keep-1)
                            grow_hij(1:n_keep-1) = tmp_hij(1:n_keep-1)

                            call move_alloc(grow_idx,tmp_idx)
                            call move_alloc(grow_hij,tmp_hij)

                        endif

                        tmp_idx(n_keep) = idx

                        tmp_hij(n_keep) = get_matrix_element( &
                            connected(:,j), n_words, &
                            basis(:,i), n_words, &
                            n_spin_orbs, h, g)

                    endif

                end do

                if(allocated(connected)) deallocate(connected)

            end do

            allocate(conn(i)%idx(n_keep))
            allocate(conn(i)%hij(n_keep))

            conn(i)%idx(:) = tmp_idx(1:n_keep)
            conn(i)%hij(:) = tmp_hij(1:n_keep)

            deallocate(tmp_idx)
            deallocate(tmp_hij)

        end do

    end subroutine build_connectivity


    subroutine davidson_lowest( &
            basis, n_basis, n_words, n_spin_orbs, &
            h, g, e_low, eigenvector)

        !uses davidson method to get the lowest eigenvalue and eigenvector of the hamitlontain without explicitly storing it
        implicit none

        integer, intent(in) :: n_basis
        integer, intent(in) :: n_words
        integer, intent(in) :: n_spin_orbs

        integer(kind=8), intent(in) :: basis(n_words,n_basis)

        real(dp), intent(in) :: h(:,:)
        real(dp), intent(in) :: g(:,:,:,:)

        real(dp), intent(out) :: e_low
        real(dp), intent(out) :: eigenvector(n_basis)

        integer, parameter :: max_sub = 40
        integer, parameter :: max_iter = 100

        real(dp), parameter :: tol = 1.0d-8

        real(dp) :: eval(1)
        real(dp) :: resid_norm

        real(dp), allocatable :: V(:,:)
        real(dp), allocatable :: W(:,:)
        real(dp), allocatable :: Hs(:,:)
        real(dp), allocatable :: coeff(:)
        real(dp), allocatable :: q(:)
        real(dp), allocatable :: r(:)
        real(dp), allocatable :: diagH(:)

        integer :: m
        integer :: iter
        integer :: i
        integer :: total_edges
        

        logical :: success

        real(dp) :: denom
        real(dp) :: rnorm

        type(det_hash_t) :: hash_table
        type(connection_list_t), allocatable :: conn(:)

        if (n_basis == 1) then

            e_low = get_matrix_element( &
                basis(:,1), n_words, &
                basis(:,1), n_words, &
                n_spin_orbs, h, g)

            eigenvector = 0.0_dp
            eigenvector(1) = 1.0_dp

            return

        end if


        call build_det_hash( &
            hash_table, &
            basis, &
            n_words, &
            n_basis)

        call build_connectivity( &
            conn, &
            basis, n_basis, n_words, n_spin_orbs, &
            h, g, &
            hash_table)


        total_edges = 0

        do i=1,n_basis
            total_edges = total_edges + size(conn(i)%idx)
        end do

        allocate(V(n_basis,max_sub))
        allocate(W(n_basis,max_sub))
        allocate(Hs(max_sub,max_sub))
        allocate(coeff(max_sub))
        allocate(q(n_basis))
        allocate(r(n_basis))
        allocate(diagH(n_basis))
        !
        ! Davidson diagonal preconditioner
        !
        do i=1,n_basis

            diagH(i) = get_matrix_element( &
                basis(:,i), n_words, &
                basis(:,i), n_words, &
                n_spin_orbs, h, g)

        end do

        V = 0.0_dp

        V(:,1) = 0.0_dp
        V(1,1) = 1.0_dp


        m = 1
        call apply_hamiltonian_sparse( &
            V(:,m), W(:,m), &
            n_basis, &
            conn, diagH)



        do iter=1,max_iter
            

            Hs(1:m,1:m) = matmul( transpose(V(:,1:m)), W(:,1:m) )

            call lowest_eigenvalue( &
                Hs(1:m,1:m), &
                m, &
                eval(1), &
                success, &
                eigenvector=coeff(1:m))

            q = matmul(V(:,1:m), coeff(1:m))

            r = matmul(W(:,1:m), coeff(1:m)) - eval(1)*q

            resid_norm = sqrt(dot_product(r,r))
            if(resid_norm < tol) exit

            do i=1,n_basis

                denom = eval(1) - diagH(i)

                if(abs(denom) < 1.0d-2) then
                    denom = sign(1.0d-2, denom)
                endif

                r(i) = r(i) / denom

            end do

            rnorm = sqrt(dot_product(r,r))

            do i=1,m
                r = r - dot_product(V(:,i),r)*V(:,i)
            end do


            r = r / sqrt(dot_product(r,r))

            m = m + 1

            if(m > max_sub) exit

            V(:,m) = r

            call apply_hamiltonian_sparse( &
                V(:,m), W(:,m), &
                n_basis, &
                conn, diagH)


        end do

        e_low = eval(1)
        eigenvector = q

        do i = 1, n_basis

            if (allocated(conn(i)%idx)) then
                deallocate(conn(i)%idx)
            end if

            if (allocated(conn(i)%hij)) then
                deallocate(conn(i)%hij)
            end if

        end do

        if (allocated(conn)) then
            deallocate(conn)
        end if

        call destroy_det_hash(hash_table)

    end subroutine davidson_lowest

    subroutine apply_hamiltonian_sparse( &
        x, y, &
        n_basis, &
        conn, diagH)
        !applies non-explitily stored hamitlonitan to a general vector, x. Used in davidson diagonalisation

        implicit none

        integer, intent(in) :: n_basis

        real(dp), intent(in) :: x(n_basis)
        real(dp), intent(out) :: y(n_basis)

        real(dp), intent(in) :: diagH(n_basis)

        type(connection_list_t), intent(in) :: conn(n_basis)

        integer :: i
        integer :: k
        integer :: idx

        y = 0.0_dp

        do i = 1, n_basis

            y(i) = y(i) + diagH(i) * x(i)

            do k = 1, size(conn(i)%idx)

                idx = conn(i)%idx(k)

                y(idx) = y(idx) + conn(i)%hij(k) * x(i)

            end do

        end do

    end subroutine apply_hamiltonian_sparse



end module diagonalisers
