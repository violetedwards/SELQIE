module selection
    use constants
    use hamiltonian
    use determinants
    implicit none

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


    subroutine hci(selected_space, n_sel_space, n_words, old_e, old_c, &
                n_spin_orbs, h, g, epsilon)

        integer(kind=8), allocatable, intent(inout) :: selected_space(:,:)
        integer, intent(inout) :: n_sel_space

        integer, intent(in) :: n_words
        integer, intent(in) :: n_spin_orbs

        real(dp), intent(in) :: old_e
        real(dp), intent(in) :: old_c(:)
        real(dp), intent(in) :: h(:,:)
        real(dp), intent(in) :: g(:,:,:,:)
        real(dp), intent(in) :: epsilon

        integer(kind=8), allocatable :: connected(:,:)
        integer(kind=8), allocatable :: new_space(:,:)
        integer(kind=8), allocatable :: final_space(:,:)

        integer :: i,j,order
        integer :: n_connected
        integer :: n_new

        real(dp) :: hij

        !--------------------------------------------------
        ! Start with current selected space
        !--------------------------------------------------

        allocate(new_space(n_words,max(8,n_sel_space)))

        n_new = 0

        call append_unique_dets( &
            new_space, n_new, &
            selected_space(:,1:n_sel_space), &
            n_sel_space, n_words)

        !--------------------------------------------------
        ! HCI screening during generation
        !--------------------------------------------------

        do i = 1, n_sel_space

            do order = 1,2

                call get_connected_dets( &
                    selected_space(:,i), &
                    n_words, &
                    n_spin_orbs, &
                    order, &
                    connected, &
                    n_connected)

                do j = 1, n_connected

                    hij = get_matrix_element( &
                        connected(:,j), &
                        n_words, &
                        selected_space(:,i), &
                        n_words, &
                        n_spin_orbs, &
                        h, g)

                    if(abs(hij*old_c(i)) > epsilon) then

                        call append_unique_dets( &
                            new_space, &
                            n_new, &
                            connected(:,j:j), &
                            1, &
                            n_words)

                    end if

                end do

                if(allocated(connected)) deallocate(connected)

            end do

        end do

        allocate(final_space(n_words,n_new))

        final_space = new_space(:,1:n_new)

        call move_alloc(final_space,selected_space)

        n_sel_space = n_new

    end subroutine hci

end module selection
