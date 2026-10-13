module hamiltonian
    use constants
    implicit none

contains

    ! Test whether a spin-orbital occupation bit is set.
    logical function determinant_is_occupied(bitstring, n_words, orbital)
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

    ! Extract occupied spin-orbital indices from a determinant.
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

    ! Find the spin orbitals removed and added when changing between determinants.
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

    ! Count occupied spin orbitals preceding an index for fermionic parity.
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

    ! Return a bitstring with the requested spin orbital occupied.
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

    ! Return a bitstring with the requested spin orbital unoccupied.
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

    ! Compute the fermionic sign for applying the specified excitation.
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

    ! Evaluate a Slater-Condon Hamiltonian element directly from two determinants.
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

        ! Apply the diagonal, single-excitation, or double-excitation formula.
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

    ! Evaluate a matrix element using precomputed excitation data and row occupations.
    real(kind=dp) function get_matrix_element_from_data(n_spin_orbs, edge_data, occ, occ_count, h, g)
        integer, intent(in) :: n_spin_orbs, edge_data(6), occ_count, occ(:)
        real(kind=dp), intent(in) :: h(n_spin_orbs, n_spin_orbs)
        real(kind=dp), intent(in) :: g(n_spin_orbs, n_spin_orbs, n_spin_orbs, n_spin_orbs)
        integer :: n_added, alpha, beta, a, b, i, j, phase
        real(kind=dp) :: sum_term

        ! Select the Slater-Condon formula from the stored excitation data.
        n_added = edge_data(1)
        phase = edge_data(6)
        if (n_added == 0) then
            sum_term = 0.0_dp
            do i = 1, occ_count
                sum_term = sum_term + h(occ(i), occ(i))
            end do
            do i = 1, occ_count
                do j = 1, occ_count
                    if (i /= j) then
                        sum_term = sum_term + 0.5_dp * &
                            (g(occ(i), occ(j), occ(i), occ(j)) - g(occ(i), occ(j), occ(j), occ(i)))
                    end if
                end do
            end do
            get_matrix_element_from_data = sum_term
        else if (n_added == 1) then
            alpha = edge_data(2)
            a = edge_data(4)
            sum_term = h(alpha, a)
            do i = 1, occ_count
                if (occ(i) /= a) then
                    sum_term = sum_term + (g(alpha, occ(i), a, occ(i)) - g(alpha, occ(i), occ(i), a))
                end if
            end do
            get_matrix_element_from_data = real(phase, kind=dp) * sum_term
        else if (n_added == 2) then
            alpha = edge_data(2)
            beta = edge_data(3)
            a = edge_data(4)
            b = edge_data(5)
            get_matrix_element_from_data = real(phase, kind=dp) * &
                (g(alpha, beta, a, b) - g(alpha, beta, b, a))
        else
            get_matrix_element_from_data = 0.0_dp
        end if
    end function get_matrix_element_from_data

    ! Apply the sparse-connectivity Hamiltonian to a block of vectors without storing its values.
    function apply_hamiltonian_matrix_free(input_vectors, det_basis_set, n_words, n_spin_orbs, h, g, &
                                           row_ptr, column_index, edge_data) result(output_vectors)
        real(kind=dp), intent(in) :: input_vectors(:,:)
        integer(kind=8), intent(in) :: det_basis_set(:,:)
        integer, intent(in) :: n_words, n_spin_orbs, row_ptr(:), column_index(:)
        integer, intent(in) :: edge_data(:,:)
        real(kind=dp), intent(in) :: h(n_spin_orbs, n_spin_orbs)
        real(kind=dp), intent(in) :: g(n_spin_orbs, n_spin_orbs, n_spin_orbs, n_spin_orbs)
        real(kind=dp) :: output_vectors(size(input_vectors, 1), size(input_vectors, 2))
        integer :: column, row, neighbor, basis_index, edge, i, occ_count
        integer :: basis_indices(size(input_vectors, 2))
        integer :: occ(n_spin_orbs)
        logical :: all_unit_vectors, is_unit_vector
        real(kind=dp) :: element

        output_vectors = 0.0_dp
        all_unit_vectors = .true.
        ! Detect unit-vector inputs used by Davidson when extracting the diagonal.
        do column = 1, size(input_vectors, 2)
            basis_index = 0
            is_unit_vector = .true.
            do i = 1, size(input_vectors, 1)
                if (input_vectors(i, column) == 1.0_dp .and. basis_index == 0) then
                    basis_index = i
                else if (input_vectors(i, column) /= 0.0_dp) then
                    is_unit_vector = .false.
                    exit
                end if
            end do

            if (is_unit_vector .and. basis_index > 0) then
                basis_indices(column) = basis_index
            else
                all_unit_vectors = .false.
            end if
        end do

        ! A unit vector needs only its corresponding Hamiltonian column.
        if (all_unit_vectors) then
            do column = 1, size(input_vectors, 2)
                basis_index = basis_indices(column)
                call get_occupied_orbitals(det_basis_set(:, basis_index), n_words, n_spin_orbs, occ, occ_count)
                do edge = row_ptr(basis_index), row_ptr(basis_index + 1) - 1
                    neighbor = column_index(edge)
                    output_vectors(neighbor, column) = get_matrix_element_from_data(n_spin_orbs, &
                        edge_data(:, edge), occ, occ_count, h, g)
                end do
            end do
        else
            ! Apply each lower-triangle edge to both rows of the symmetric product.
            do row = 1, size(input_vectors, 1)
                call get_occupied_orbitals(det_basis_set(:, row), n_words, n_spin_orbs, occ, occ_count)
                do edge = row_ptr(row), row_ptr(row + 1) - 1
                    neighbor = column_index(edge)
                    if (neighbor > row) exit
                    element = get_matrix_element_from_data(n_spin_orbs, edge_data(:, edge), occ, occ_count, h, g)
                    do column = 1, size(input_vectors, 2)
                        output_vectors(row, column) = output_vectors(row, column) + &
                            element * input_vectors(neighbor, column)
                    end do
                    if (neighbor /= row) then
                        do column = 1, size(input_vectors, 2)
                            output_vectors(neighbor, column) = output_vectors(neighbor, column) + &
                                element * input_vectors(row, column)
                        end do
                    end if
                end do
            end do
        end if
    end function apply_hamiltonian_matrix_free

    ! Build CSR connectivity and excitation metadata for determinant pairs coupled by a two-body Hamiltonian.
    subroutine build_hamiltonian_connectivity(det_basis_set, n_det, n_words, n_spin_orbs, row_ptr, column_index, edge_data)
        integer, intent(in) :: n_det, n_words, n_spin_orbs
        integer(kind=8), intent(in) :: det_basis_set(n_words, n_det)
        integer, allocatable, intent(out) :: row_ptr(:), column_index(:)
        integer, allocatable, intent(out) :: edge_data(:,:)
        integer, allocatable :: row_cursor(:)
        integer :: i, j, word, differing_bits, n_edges, n_added, n_removed, phase
        integer :: added(n_spin_orbs), removed(n_spin_orbs), edge

        ! Count connected entries in each row before allocating CSR storage.
        allocate(row_cursor(n_det))
        row_cursor = 0
        do i = 1, n_det
            do j = 1, i
                differing_bits = 0
                do word = 1, n_words
                    differing_bits = differing_bits + popcnt(ieor(det_basis_set(word, i), det_basis_set(word, j)))
                    if (differing_bits > 4) exit
                end do
                if (differing_bits <= 4) then
                    row_cursor(i) = row_cursor(i) + 1
                    if (i /= j) row_cursor(j) = row_cursor(j) + 1
                end if
            end do
        end do

        ! Convert row counts to one-based CSR offsets.
        allocate(row_ptr(n_det + 1))
        row_ptr(1) = 1
        do i = 1, n_det
            row_ptr(i + 1) = row_ptr(i) + row_cursor(i)
        end do
        n_edges = row_ptr(n_det + 1) - 1
        allocate(column_index(n_edges), edge_data(6, n_edges))
        row_cursor = row_ptr(1:n_det)

        ! Store each lower-triangle edge and its mirrored orientation.
        do i = 1, n_det
            do j = 1, i
                differing_bits = 0
                do word = 1, n_words
                    differing_bits = differing_bits + popcnt(ieor(det_basis_set(word, i), det_basis_set(word, j)))
                    if (differing_bits > 4) exit
                end do
                if (differing_bits <= 4) then
                    call compare_dets(det_basis_set(:, i), n_words, det_basis_set(:, j), n_words, &
                        n_spin_orbs, added, removed, n_added, n_removed)
                    phase = excitation_phase(det_basis_set(:, i), n_words, added, n_added, removed, n_removed)

                    edge = row_cursor(i)
                    column_index(edge) = j
                    edge_data(:, edge) = 0
                    edge_data(1, edge) = n_added
                    edge_data(6, edge) = phase
                    if (n_added >= 1) then
                        edge_data(2, edge) = added(1)
                        edge_data(4, edge) = removed(1)
                    end if
                    if (n_added == 2) then
                        edge_data(3, edge) = added(2)
                        edge_data(5, edge) = removed(2)
                    end if
                    row_cursor(i) = row_cursor(i) + 1
                    if (i /= j) then
                        edge = row_cursor(j)
                        column_index(edge) = i
                        edge_data(:, edge) = 0
                        edge_data(1, edge) = n_added
                        edge_data(6, edge) = phase
                        if (n_added >= 1) then
                            edge_data(2, edge) = removed(1)
                            edge_data(4, edge) = added(1)
                        end if
                        if (n_added == 2) then
                            edge_data(3, edge) = removed(2)
                            edge_data(5, edge) = added(2)
                        end if
                        row_cursor(j) = row_cursor(j) + 1
                    end if
                end if
            end do
        end do
    end subroutine build_hamiltonian_connectivity

    ! Assemble the full dense Hamiltonian for a determinant basis.
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
