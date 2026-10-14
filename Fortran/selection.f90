module selection
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

        subroutine append_unique_dets(dets, n_dets, candidates, n_candidates, n_words)
            integer(kind = 8), allocatable, intent(inout) :: dets(:,:)
            integer, intent(inout) :: n_dets
            integer(kind = 8), intent(in) :: candidates(:,:)
            integer, intent(in) :: n_candidates, n_words
        end subroutine append_unique_dets
    end interface

contains

    subroutine add_selected_dets(sel_method, selected_space, n_sel_space, n_words, old_e, old_c, n_spin_orbs, h, g, epsilon)
        character(len=*), intent(in) :: sel_method
        integer(kind = 8), allocatable, intent(inout) :: selected_space(:,:)
        integer, intent(inout) :: n_sel_space
        integer, intent(in) :: n_words, n_spin_orbs
        real(kind = dp), intent(in) :: old_e, old_c(:), h(:,:), g(:,:,:,:), epsilon

        if (sel_method == "hci") then
            call hci(selected_space, n_sel_space, n_words, old_e, old_c, n_spin_orbs, h, g, epsilon)
        end if

    end subroutine add_selected_dets


    subroutine hci(selected_space, n_sel_space, n_words, old_e, old_c, n_spin_orbs, h, g, epsilon)
        integer(kind = 8), allocatable, intent(inout) :: selected_space(:,:)
        integer, intent(inout) :: n_sel_space
        integer, intent(in) :: n_words, n_spin_orbs
        real(kind = dp), intent(in) :: old_e, old_c(:), h(:,:), g(:,:,:,:), epsilon

        integer(kind = 8), allocatable :: interacting_space(:,:), connected(:,:), new_space(:,:), final_space(:,:)
        integer :: n_interacting, n_connected, n_new, order, i, j
        logical :: accept

        !interacting space: singles of every selected det, then doubles of every selected det
        allocate(interacting_space(n_words, 8))
        n_interacting = 0
        do order = 1, 2
            do i = 1, n_sel_space
                call get_connected_dets(selected_space(:, i), n_words, n_spin_orbs, order, connected, n_connected)
                if (n_connected > 0) then
                    call append_unique_dets(interacting_space, n_interacting, connected(:, 1:n_connected), n_connected, n_words)
                end if
            end do
        end do

        !selected space first, then every interacting det with a large enough element to any selected det
        allocate(new_space(n_words, max(8, n_sel_space)))
        n_new = 0
        call append_unique_dets(new_space, n_new, selected_space(:, 1:n_sel_space), n_sel_space, n_words)
        do i = 1, n_interacting
            accept = .false.
            do j = 1, n_sel_space
                if (abs(get_matrix_element(interacting_space(:, i), n_words, selected_space(:, j), n_words, &
                        n_spin_orbs, h, g)) > epsilon) then
                    accept = .true.
                    exit
                end if
            end do
            if (accept) then
                call append_unique_dets(new_space, n_new, interacting_space(:, i:i), 1, n_words)
            end if
        end do

        allocate(final_space(n_words, n_new))
        final_space = new_space(:, 1:n_new)
        call move_alloc(final_space, selected_space)
        n_sel_space = n_new

    end subroutine hci

end module selection
