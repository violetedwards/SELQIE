!   FILE: determinants.f90

subroutine determinant_occ_to_bit(occupations, n_elec, n_spin_orbs, bitstring, n_words)
!
!       Creates a general-length integer bitstring from an occupation array.
!       The bitstring is stored as an array of 64-bit words.
!
        implicit none

        integer, intent(in) :: n_elec, n_spin_orbs, n_words
        integer, intent(in) :: occupations(n_elec)
        integer(kind=8), intent(out) :: bitstring(n_words)
        integer :: i, word_index, bit_index

        bitstring = 0_8

        do i = 1, n_elec
            if (occupations(i) >= 1 .and. occupations(i) <= n_spin_orbs) then
                word_index = (occupations(i) - 1) / 64 + 1
                bit_index = mod(occupations(i) - 1, 64)
                bitstring(word_index) = ibset(bitstring(word_index), bit_index)
            end if
        end do

end subroutine determinant_occ_to_bit

subroutine determinant_bit_to_occ(bitstring, n_words, n_spin_orbs, occupations, n_elec)
!
!       Gets occupation array from a general-length integer bitstring.
!
        implicit none

        integer, intent(in) :: n_words, n_spin_orbs, n_elec
        integer(kind=8), intent(in) :: bitstring(n_words)
        integer, intent(out) :: occupations(n_elec)
        integer :: i, index, bit_index, word_index

        occupations = 0
        index = 1

        do i = 1, n_spin_orbs
            word_index = (i - 1) / 64 + 1
            bit_index = mod(i - 1, 64)
            if (btest(bitstring(word_index), bit_index) .and. index <= n_elec) then
                occupations(index) = i
                index = index + 1
            end if
        end do

end subroutine determinant_bit_to_occ

subroutine reference_determinant(n_orb, n_alpha, n_beta, bitstring, n_words)
        implicit none

        integer, intent(in) :: n_orb, n_alpha, n_beta, n_words
        integer(kind=8), intent(out) :: bitstring(n_words)
        integer :: occupations(2*n_orb)
        integer :: i

        occupations = 0
        bitstring = 0_8

        do i = 1, n_alpha
            occupations(i) = i
        end do
        do i = 1, n_beta
            occupations(n_alpha + i) = n_orb + i
        end do

        call determinant_occ_to_bit(occupations, n_alpha + n_beta, 2*n_orb, bitstring, n_words)

end subroutine reference_determinant

logical function determinant_is_occupied(bitstring, n_words, orbital)
        implicit none
        integer, intent(in) :: n_words, orbital
        integer(kind=8), intent(in) :: bitstring(n_words)
        integer :: word_index, bit_index

        determinant_is_occupied = .false.
        if (orbital <= 0) return
        word_index = (orbital - 1) / 64 + 1
        bit_index = mod(orbital - 1, 64)
        if (word_index <= n_words) then
            determinant_is_occupied = btest(bitstring(word_index), bit_index)
        end if
end function determinant_is_occupied

subroutine append_connected_det(candidate, connected_dets, n_connected, capacity, n_words)
        implicit none
        integer, intent(in) :: n_words
        integer(kind=8), intent(in) :: candidate(n_words)
        integer(kind=8), allocatable, intent(out) :: connected_dets(:,:)
        integer, intent(inout) :: n_connected, capacity
        integer(kind=8), allocatable :: new_dets(:,:)

        if (.not. allocated(connected_dets)) then
            allocate(connected_dets(n_words, 8))
            capacity = 8
        else if (n_connected + 1 > capacity) then
            allocate(new_dets(n_words, max(capacity * 2, 8)))
            if (n_connected > 0) then
                new_dets(:, 1:n_connected) = connected_dets(:, 1:n_connected)
            end if
            call move_alloc(new_dets, connected_dets)
            capacity = size(connected_dets, 2)
        end if

        n_connected = n_connected + 1
        connected_dets(:, n_connected) = candidate
end subroutine append_connected_det

subroutine get_connected_dets(ref_det, n_words, n_spin_orbs, order, connected_dets, n_connected)
        implicit none

        integer, intent(in) :: n_words, n_spin_orbs, order
        integer(kind=8), intent(in) :: ref_det(n_words)
        integer(kind=8), allocatable, intent(out) :: connected_dets(:,:)
        integer, intent(out) :: n_connected

        logical :: determinant_is_occupied
        interface
            subroutine append_connected_det(candidate, connected_dets, n_connected, capacity, n_words)
                integer, intent(in) :: n_words
                integer(kind=8), intent(in) :: candidate(n_words)
                integer(kind=8), allocatable, intent(out) :: connected_dets(:,:)
                integer, intent(inout) :: n_connected, capacity
            end subroutine append_connected_det
        end interface

        integer :: occupied(n_spin_orbs)
        integer :: virtual(n_spin_orbs)
        integer :: n_occ, n_virt, i, j, a, b, k, capacity, max_conn
        integer :: n_orb
        integer :: removed_idx(order), added_idx(order)
        integer(kind=8) :: candidate(n_words)
        logical :: duplicate

        n_connected = 0
        capacity = 0
        n_orb = n_spin_orbs / 2

        if (order <= 0) then
            if (allocated(connected_dets)) deallocate(connected_dets)
            return
        end if

        n_occ = 0
        n_virt = 0
        do i = 1, n_spin_orbs
            if (determinant_is_occupied(ref_det, n_words, i)) then
                n_occ = n_occ + 1
                occupied(n_occ) = i
            else
                n_virt = n_virt + 1
                virtual(n_virt) = i
            end if
        end do

        if (order > min(n_occ, n_virt)) then
            if (allocated(connected_dets)) deallocate(connected_dets)
            return
        end if

        if (order == 1) then
            max_conn = max(1, n_occ * n_virt)
        else if (order == 2) then
            max_conn = max(1, (n_occ * (n_occ - 1) / 2) * (n_virt * (n_virt - 1) / 2))
        else
            max_conn = 1
        end if

        if (allocated(connected_dets)) deallocate(connected_dets)
        allocate(connected_dets(n_words, max_conn))
        connected_dets = 0_8
        capacity = max_conn

        if (order == 1) then
            do a = 1, n_occ
                removed_idx(1) = occupied(a)
                do b = 1, n_virt
                    added_idx(1) = virtual(b)
                    if ((removed_idx(1) <= n_orb) .eqv. (added_idx(1) <= n_orb)) then
                        candidate = ref_det
                        candidate((removed_idx(1)-1)/64 + 1) = ibclr(candidate((removed_idx(1)-1)/64 + 1), mod(removed_idx(1)-1,64))
                        candidate((added_idx(1)-1)/64 + 1) = ibset(candidate((added_idx(1)-1)/64 + 1), mod(added_idx(1)-1,64))

                        duplicate = .false.
                        do j = 1, n_connected
                            if (all(connected_dets(:, j) == candidate)) then
                                duplicate = .true.
                                exit
                            end if
                        end do

                        if (.not. duplicate) then
                            n_connected = n_connected + 1
                            connected_dets(:, n_connected) = candidate
                        end if
                    end if
                end do
            end do
        end if

        if (order == 2) then
            do a = 1, n_occ - 1
                do i = a + 1, n_occ
                    removed_idx(1) = occupied(a)
                    removed_idx(2) = occupied(i)
                    do b = 1, n_virt - 1
                        do j = b + 1, n_virt
                            added_idx(1) = virtual(b)
                            added_idx(2) = virtual(j)
                            if (merge(1, 0, removed_idx(1) <= n_orb) + merge(1, 0, removed_idx(2) <= n_orb) == &
                                merge(1, 0, added_idx(1) <= n_orb) + merge(1, 0, added_idx(2) <= n_orb)) then
                                candidate = ref_det
                                candidate((removed_idx(1)-1)/64 + 1) = ibclr(candidate((removed_idx(1)-1)/64 + 1), mod(removed_idx(1)-1,64))
                                candidate((added_idx(1)-1)/64 + 1) = ibset(candidate((added_idx(1)-1)/64 + 1), mod(added_idx(1)-1,64))
                                candidate((removed_idx(2)-1)/64 + 1) = ibclr(candidate((removed_idx(2)-1)/64 + 1), mod(removed_idx(2)-1,64))
                                candidate((added_idx(2)-1)/64 + 1) = ibset(candidate((added_idx(2)-1)/64 + 1), mod(added_idx(2)-1,64))

                                duplicate = .false.
                                do k = 1, n_connected
                                    if (all(connected_dets(:, k) == candidate)) then
                                        duplicate = .true.
                                        exit
                                    end if
                                end do

                                if (.not. duplicate) then
                                    n_connected = n_connected + 1
                                    connected_dets(:, n_connected) = candidate
                                end if
                            end if
                        end do
                    end do
                end do
            end do
        end if
end subroutine get_connected_dets
