module selection
	use constants
	use hamiltonian, only: get_matrix_element
	implicit none

contains

	subroutine hci_select(selected_space, n_selected, n_words, n_spin_orbs, h, g, epsilon, new_space, n_new)
		implicit none

		interface
			subroutine get_connected_dets(ref_det, n_words, n_spin_orbs, order, connected_dets, n_connected)
				integer, intent(in) :: n_words, n_spin_orbs, order
				integer(kind=8), intent(in) :: ref_det(n_words)
				integer(kind=8), allocatable, intent(out) :: connected_dets(:,:)
				integer, intent(out) :: n_connected
			end subroutine get_connected_dets
		end interface

		integer, intent(in) :: n_selected, n_words, n_spin_orbs
		integer(kind=8), intent(in) :: selected_space(n_words, n_selected)
		real(kind=dp), intent(in) :: h(n_spin_orbs, n_spin_orbs)
		real(kind=dp), intent(in) :: g(n_spin_orbs, n_spin_orbs, n_spin_orbs, n_spin_orbs)
		real(kind=dp), intent(in) :: epsilon
		integer(kind=8), allocatable, intent(out) :: new_space(:,:)
		integer, intent(out) :: n_new
		integer(kind=8), allocatable :: candidate_space(:,:), connected_dets(:,:), expanded(:,:)
		integer :: capacity, selected_index, order, n_connected, candidate_index, existing_index
		integer :: source_index
		logical :: duplicate, accepted

		capacity = max(8, n_selected)
		allocate(candidate_space(n_words, capacity))
		candidate_space = 0_8
		candidate_space(:, 1:n_selected) = selected_space
		n_new = n_selected

		! Examine singles and doubles connected to every currently selected determinant.
		do source_index = 1, n_selected
			do order = 1, 2
				call get_connected_dets(selected_space(:, source_index), n_words, n_spin_orbs, &
					order, connected_dets, n_connected)

				do candidate_index = 1, n_connected
					duplicate = .false.
					do existing_index = 1, n_new
						if (all(candidate_space(:, existing_index) == connected_dets(:, candidate_index))) then
							duplicate = .true.
							exit
						end if
					end do
					if (duplicate) cycle

					accepted = .false.
					do selected_index = 1, n_selected
						if (abs(get_matrix_element(connected_dets(:, candidate_index), n_words, &
							selected_space(:, selected_index), n_words, n_spin_orbs, h, g)) > epsilon) then
							accepted = .true.
							exit
						end if
					end do

					if (accepted) then
						if (n_new == capacity) then
							allocate(expanded(n_words, 2 * capacity))
							expanded(:, 1:n_new) = candidate_space(:, 1:n_new)
							call move_alloc(expanded, candidate_space)
							capacity = size(candidate_space, 2)
						end if
						n_new = n_new + 1
						candidate_space(:, n_new) = connected_dets(:, candidate_index)
					end if
				end do
			end do
		end do

		allocate(new_space(n_words, n_new))
		new_space = candidate_space(:, 1:n_new)
	end subroutine hci_select

end module selection