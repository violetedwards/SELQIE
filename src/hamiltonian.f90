module hamiltonian
    use constants
    use determinants
    implicit none

contains


    subroutine get_occupied_orbitals(det, n_words, n_spin_orbs, occ, n_occ)
        integer, intent(in) :: n_words, n_spin_orbs
        integer(kind=8), intent(in) :: det(n_words)
        integer, intent(out) :: occ(n_spin_orbs)
        integer, intent(out) :: n_occ
        integer :: i

        n_occ = 0
        do i = 1, n_spin_orbs
            if (determinant_is_occupied(det, n_words, i)) then
                n_occ = n_occ + 1
                occ(n_occ) = i
            end if
        end do
    end subroutine get_occupied_orbitals

    subroutine compare_dets(det1, n_words1, det2, n_words2, n_spin_orbs, added, removed, n_added, n_removed)
        integer, intent(in) :: n_words1, n_words2, n_spin_orbs
        integer(kind=8), intent(in) :: det1(n_words1), det2(n_words2)
        integer, intent(out) :: added(n_spin_orbs), removed(n_spin_orbs)
        integer, intent(out) :: n_added, n_removed
        integer :: i

        n_added = 0
        n_removed = 0
        do i = 1, n_spin_orbs
            if (determinant_is_occupied(det1, n_words1, i) .and. .not. determinant_is_occupied(det2, n_words2, i)) then
                n_removed = n_removed + 1
                removed(n_removed) = i
            else if (.not. determinant_is_occupied(det1, n_words1, i) .and. determinant_is_occupied(det2, n_words2, i)) then
                n_added = n_added + 1
                added(n_added) = i
            end if
        end do
    end subroutine compare_dets

    integer function count_bits_before(bitstring, n_words, orbital)
        integer, intent(in) :: n_words, orbital
        integer(kind=8), intent(in) :: bitstring(n_words)
        integer :: i

        count_bits_before = 0
        do i = 1, orbital - 1
            if (determinant_is_occupied(bitstring, n_words, i)) then
                count_bits_before = count_bits_before + 1
            end if
        end do
    end function count_bits_before

    function set_orbital(bitstring, n_words, orbital) result(new_bitstring)
        integer, intent(in) :: n_words, orbital
        integer(kind=8), intent(in) :: bitstring(n_words)
        integer(kind=8) :: new_bitstring(n_words)
        integer :: word_index, bit_index

        new_bitstring = bitstring
        word_index = (orbital - 1) / 64 + 1
        bit_index = mod(orbital - 1, 64)
        if (word_index <= n_words) then
            new_bitstring(word_index) = ibset(new_bitstring(word_index), bit_index)
        end if
    end function set_orbital

    function clear_orbital(bitstring, n_words, orbital) result(new_bitstring)
        integer, intent(in) :: n_words, orbital
        integer(kind=8), intent(in) :: bitstring(n_words)
        integer(kind=8) :: new_bitstring(n_words)
        integer :: word_index, bit_index

        new_bitstring = bitstring
        word_index = (orbital - 1) / 64 + 1
        bit_index = mod(orbital - 1, 64)
        if (word_index <= n_words) then
            new_bitstring(word_index) = ibclr(new_bitstring(word_index), bit_index)
        end if
    end function clear_orbital

    integer function excitation_phase(det, n_words, added, n_added, removed, n_removed)
        integer, intent(in) :: n_words, n_added, n_removed
        integer(kind=8), intent(in) :: det(n_words)
        integer, intent(in) :: added(n_added), removed(n_removed)
        integer(kind=8) :: bitstring(n_words)
        integer :: i, phase, count_before

        bitstring = det
        phase = 1

        do i = 1, n_removed
            count_before = count_bits_before(bitstring, n_words, removed(i))
            if (mod(count_before, 2) /= 0) phase = -phase
            bitstring = clear_orbital(bitstring, n_words, removed(i))
        end do

        do i = n_added, 1, -1
            count_before = count_bits_before(bitstring, n_words, added(i))
            if (mod(count_before, 2) /= 0) phase = -phase
            bitstring = set_orbital(bitstring, n_words, added(i))
        end do

        excitation_phase = phase
    end function excitation_phase

    real(kind=dp) function get_matrix_element(column_det, n_words_col, row_det, n_words_row, n_spin_orbs, h, g)
        integer, intent(in) :: n_words_col, n_words_row, n_spin_orbs
        integer(kind=8), intent(in) :: column_det(n_words_col), row_det(n_words_row)
        real(kind=dp), intent(in) :: h(n_spin_orbs, n_spin_orbs)
        real(kind=dp), intent(in) :: g(n_spin_orbs, n_spin_orbs, n_spin_orbs, n_spin_orbs)
        integer :: added(n_spin_orbs), removed(n_spin_orbs)
        integer :: n_added, n_removed, i, j, alpha, beta, a, b, occ_count
        integer :: occ(n_spin_orbs)
        integer :: phase
        real(kind=dp) :: sum_term

        call compare_dets(row_det, n_words_row, column_det, n_words_col, n_spin_orbs, added, removed, n_added, n_removed)
        if (n_added > 2) then
            get_matrix_element = 0.0_dp
            return
        end if

        phase = excitation_phase(row_det, n_words_row, added, n_added, removed, n_removed)

        if (n_added == 0) then
            call get_occupied_orbitals(row_det, n_words_row, n_spin_orbs, occ, occ_count)
            sum_term = 0.0_dp
            do i = 1, occ_count
                sum_term = sum_term + h(occ(i), occ(i))
            end do
            do i = 1, occ_count
                do j = 1, occ_count
                    if (i /= j) then
                        sum_term = sum_term + 0.5_dp * (g(occ(i), occ(j), occ(i), occ(j)) - g(occ(i), occ(j), occ(j), occ(i)))
                    end if
                end do
            end do
            get_matrix_element = sum_term
            return
        end if

        if (n_added == 1) then
            alpha = added(1)
            a = removed(1)
            sum_term = h(alpha, a)
            call get_occupied_orbitals(row_det, n_words_row, n_spin_orbs, occ, occ_count)
            do i = 1, occ_count
                if (occ(i) /= a) then
                    sum_term = sum_term + (g(alpha, occ(i), a, occ(i)) - g(alpha, occ(i), occ(i), a))
                end if
            end do
            get_matrix_element = real(phase, kind=dp) * sum_term
            return
        end if

        if (n_added == 2) then
            alpha = added(1)
            beta = added(2)
            a = removed(1)
            b = removed(2)
            get_matrix_element = real(phase, kind=dp) * (g(alpha, beta, a, b) - g(alpha, beta, b, a))
            return
        end if

        get_matrix_element = 0.0_dp
    end function get_matrix_element

    subroutine build_hamiltonian(det_basis_set, n_det, n_words, n_spin_orbs, h, g, hamiltonian)
        integer, intent(in) :: n_det, n_words, n_spin_orbs
        integer(kind=8), intent(in) :: det_basis_set(n_words, n_det)
        real(kind=dp), intent(in) :: h(n_spin_orbs, n_spin_orbs)
        real(kind=dp), intent(in) :: g(n_spin_orbs, n_spin_orbs, n_spin_orbs, n_spin_orbs)
        real(kind=dp), intent(out) :: hamiltonian(n_det, n_det)
        integer :: i, j
        real(kind=dp) :: element

        hamiltonian = 0.0_dp
        do i = 1, n_det
            do j = 1, i
                element = get_matrix_element(det_basis_set(:, i), n_words, det_basis_set(:, j), n_words, n_spin_orbs, h, g)
                hamiltonian(i, j) = element
                hamiltonian(j, i) = element
            end do
        end do
    end subroutine build_hamiltonian


end module hamiltonian
