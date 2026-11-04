module selection
    use constants
    use hamiltonian
    use determinants
    implicit none

contains

    subroutine add_selected_dets(sel_method, selected_space, n_sel_space, n_words, old_e, old_c, n_spin_orbs, h, g, epsilon, sigma)

        !routes to correct selection algorithm

        character(len=*), intent(in) :: sel_method
        integer(kind = 8), allocatable, intent(inout) :: selected_space(:,:)
        integer, intent(inout) :: n_sel_space
        integer, intent(in) :: n_words, n_spin_orbs
        real(kind = dp), intent(in) :: old_e, old_c(:), h(:,:), g(:,:,:,:), epsilon, sigma

        if (sel_method == "hci") then
            call hci(selected_space, n_sel_space, n_words, old_c, n_spin_orbs, h, g, epsilon)
        else if (sel_method == "aci") then
            call aci(selected_space, n_sel_space, n_words, old_c, old_e, n_spin_orbs, h, g, sigma)
        end if

    end subroutine add_selected_dets


    subroutine hci(selected_space, n_sel_space, n_words, old_c, n_spin_orbs, h, g, epsilon)
        !Heat Bath Configuration Interaction (HCI)

        integer(kind=8), allocatable, intent(inout) :: selected_space(:,:)
        integer, intent(inout) :: n_sel_space

        integer, intent(in) :: n_words
        integer, intent(in) :: n_spin_orbs

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

        ! Start with current selected space

        allocate(new_space(n_words,max(8,n_sel_space)))

        n_new = 0

        call append_unique_dets(new_space, n_new, selected_space(:,1:n_sel_space), n_sel_space, n_words)

        ! Generate connected determinants & accept those where |matrix element * c| > epsilon

        do i = 1, n_sel_space
            do order = 1,2
                call get_connected_dets(selected_space(:,i), n_words, n_spin_orbs, order, connected, n_connected)

                do j = 1, n_connected
                    hij = get_matrix_element(connected(:,j), n_words, selected_space(:,i), n_words, n_spin_orbs, h, g)

                    if(abs(hij*old_c(i)) > epsilon) then
                        call append_unique_dets(new_space, n_new, connected(:,j:j), 1, n_words)
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

    subroutine aci(selected_space, n_sel_space, n_words, old_c, old_e, n_spin_orbs, h, g, sigma)
        !Adaptive Configuration Interaction

        integer(kind=8), allocatable, intent(inout) :: selected_space(:,:)
        integer, intent(inout) :: n_sel_space

        integer, intent(in) :: n_words
        integer, intent(in) :: n_spin_orbs

        real(dp), intent(in) :: old_c(:)
        real(dp), intent(in) :: h(:,:)
        real(dp), intent(in) :: g(:,:,:,:)
        real(dp), intent(in) :: sigma

        integer(kind=8), allocatable :: connected(:,:)
        integer(kind=8), allocatable :: interacting_space(:,:)
        integer(kind=8), allocatable :: final_space(:,:)
        integer(kind=8), allocatable :: additional_space(:,:)

        integer :: i,j,order,max_index
        integer :: n_connected, n_interacting, n_final, n_additional

        integer(kind=8) :: temp_det(n_words)
        real(dp) :: old_e, hij, ei, ep, v, delta, sum_energy_importance, temp_importance
        real(dp), allocatable :: energy_importance(:)
        real(dp), allocatable :: sorted_energy_importance(:)
        integer(kind=8), allocatable :: sorted_interacting_space(:,:)

        ! Start with current selected space

        allocate(additional_space(n_words,max(8,n_sel_space)))
        n_additional= 0

        ! Now generate interacting space with only unique elements

        allocate(interacting_space(n_words,max(8,n_sel_space)))
        n_interacting = 0
        call append_unique_dets(interacting_space, n_interacting, selected_space(:,1:n_sel_space), n_sel_space, n_words)

        do i = 1, n_sel_space
            do order=1,2
                call get_connected_dets(selected_space(:,i), n_words, n_spin_orbs, order, connected, n_connected)
                if (n_connected > 0) then
                    call append_unique_dets(interacting_space, n_interacting, connected(:,1:n_connected), n_connected, n_words)
                end if
                if (allocated(connected)) deallocate(connected)
            end do
        end do

        n_interacting = n_interacting - n_sel_space
        do i = 1, n_interacting
            interacting_space(:,i) = interacting_space(:,i+n_sel_space)
        end do

        ! Assign energy importance to each element of interacting space

        allocate(energy_importance(n_interacting))
        do i = 1, n_interacting
            ei = old_e
            ep = get_matrix_element(interacting_space(:,i), n_words, interacting_space(:,i), n_words, n_spin_orbs, h, g)
            v = 0.0_dp
            do j = 1, n_sel_space
                hij = get_matrix_element(interacting_space(:,i), n_words, selected_space(:,j), n_words, n_spin_orbs, h, g)
                v = v + hij * old_c(j)
            end do
            delta = ep - ei
            energy_importance(i) = abs(delta/2 - sqrt((delta*delta/4)+(v*v)))
        end do

        ! Sort energy importance

        !sort according to energy important sort
        allocate(sorted_energy_importance(n_interacting))
        allocate(sorted_interacting_space(n_words,n_interacting))
        sorted_energy_importance = energy_importance
        sorted_interacting_space = interacting_space(:,1:n_interacting)
        do i = 1, n_interacting - 1
            max_index = i
            do j = i + 1, n_interacting
                if (sorted_energy_importance(j) > sorted_energy_importance(max_index)) max_index = j
            end do
            if (max_index /= i) then
                temp_importance = sorted_energy_importance(i)
                sorted_energy_importance(i) = sorted_energy_importance(max_index)
                sorted_energy_importance(max_index) = temp_importance
                temp_det = sorted_interacting_space(:,i)
                sorted_interacting_space(:,i) = sorted_interacting_space(:,max_index)
                sorted_interacting_space(:,max_index) = temp_det
            end if
        end do

        !add sorted interacting space elements until remaining energy importance is <= sigma
        sum_energy_importance = sum(sorted_energy_importance)
        do i = 1, n_interacting
            if (sum_energy_importance <= sigma) exit
            call append_unique_dets(additional_space, n_additional, sorted_interacting_space(:,i:i), 1, n_words)
            sum_energy_importance = max(0.0_dp, sum_energy_importance - sorted_energy_importance(i))
        end do


        allocate(final_space(n_words,n_sel_space+n_additional))
        n_final = 0
        call append_unique_dets(final_space, n_final, selected_space(:,1:n_sel_space), n_sel_space, n_words)
        call append_unique_dets(final_space, n_final, additional_space(:,1:n_additional), n_additional, n_words)
        call move_alloc(final_space, selected_space)
        n_sel_space = n_final

        end subroutine

end module selection