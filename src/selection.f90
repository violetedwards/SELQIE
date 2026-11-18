module selection
    use constants
    use hamiltonian
    use determinants
    implicit none

contains

    subroutine add_selected_dets(sel_method, selected_space, n_sel_space, n_words, old_e, old_c, n_spin_orbs, h, g, epsilon, sigma, tdets)

        !routes to correct selection algorithm
        character(len=*), intent(in) :: sel_method
        integer(kind = 8), allocatable, intent(inout) :: selected_space(:,:)
        integer, intent(inout) :: n_sel_space
        integer, intent(in) :: n_words, n_spin_orbs, tdets
        real(kind = dp), intent(in) :: old_e, old_c(:), h(:,:), g(:,:,:,:), epsilon, sigma

        if (sel_method == "hci") then
            call hci(selected_space, n_sel_space, n_words, old_c, n_spin_orbs, h, g, epsilon)
        else if (sel_method == "aci") then
            call aci(selected_space, n_sel_space, n_words, old_c, old_e, n_spin_orbs, h, g, sigma)
        else if (sel_method == "asci") then
            call asci(selected_space, n_sel_space, n_words, old_c, old_e, n_spin_orbs, h, g, tdets)
        end if

    end subroutine add_selected_dets


    subroutine hci(selected_space, n_sel_space, n_words, old_c, n_spin_orbs, h, g, epsilon)
        !Heat Bath Configuration Interaction

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

        !start with current selected space
        allocate(new_space(n_words,max(8,n_sel_space)))

        n_new = 0

        call append_unique_dets(new_space, n_new, selected_space(:,1:n_sel_space), n_sel_space, n_words)

        !generate connected determinants & accept those where |matrix element * c| > epsilon
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

        !start with current selected space

        allocate(additional_space(n_words,max(8,n_sel_space)))
        n_additional= 0

        !now generate interacting space with only unique elements

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

        !assign energy importance to each element of interacting space

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

    subroutine asci(selected_space, n_sel_space, n_words, old_c, old_e, n_spin_orbs, h, g, tdets)
        !Adaptive Sampling Configuration Interaction

        integer(kind=8), allocatable, intent(inout) :: selected_space(:,:)
        integer, intent(inout) :: n_sel_space

        integer, intent(in) :: n_words,n_spin_orbs,tdets

        real(dp), intent(in) :: old_e
        real(dp), intent(in) :: old_c(:)
        real(dp), intent(in) :: h(:,:)
        real(dp), intent(in) :: g(:,:,:,:)

        integer(kind=8), allocatable :: connected(:,:)
        integer(kind=8), allocatable :: interacting_space(:,:)
        integer(kind=8), allocatable :: final_space(:,:)
        integer(kind=8), allocatable :: candidate_space(:,:)

        integer :: i,j,order,max_index
        integer :: n_connected, n_interacting, n_candidates, target_dets

        integer(kind=8) :: temp_det(n_words)
        real(dp) :: sum_h_c, temp_coeff, denominator
        real(dp), allocatable :: pert_amp(:), candidate_coeffs(:)

        !generate interacting space with only unique elements
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

        !get perturbed amplitudes
        allocate(pert_amp(n_interacting))
        do i = 1, n_interacting
            sum_h_c = 0.0_dp
            do j = 1, n_sel_space
                sum_h_c = sum_h_c + get_matrix_element(interacting_space(:,i),n_words, selected_space(:,j), n_words, n_spin_orbs, h, g) * old_c(j)
            end do
            denominator = get_matrix_element(interacting_space(:,i),n_words,interacting_space(:,i),n_words,n_spin_orbs,h,g)-old_e
            if (abs(denominator) > 1.0e-10_dp) then
                pert_amp(i) = sum_h_c/denominator
            else
                pert_amp(i) = 0.0_dp
            end if
        end do

        !combine old selected space and interacting space, and combine old c and perturbed amplitudes
        n_candidates = n_sel_space + n_interacting
        allocate(candidate_space(n_words, n_candidates))
        allocate(candidate_coeffs(n_candidates))

        candidate_space(:,1:n_sel_space) = selected_space(:,1:n_sel_space)
        candidate_coeffs(1:n_sel_space) = old_c(1:n_sel_space)
        if (n_interacting > 0) then
            candidate_space(:,n_sel_space+1:n_candidates) = interacting_space(:,1:n_interacting)
            candidate_coeffs(n_sel_space+1:n_candidates) = pert_amp(1:n_interacting)
        end if

        !sort amplitude list, and then sort det list by amp sorting
        do i = 1, n_candidates - 1
            max_index = i
            do j = i + 1, n_candidates
                if (abs(candidate_coeffs(j)) > abs(candidate_coeffs(max_index))) max_index = j
            end do

            if (max_index /= i) then
                temp_coeff = candidate_coeffs(i)
                candidate_coeffs(i) = candidate_coeffs(max_index)
                candidate_coeffs(max_index) = temp_coeff

                temp_det = candidate_space(:,i)
                candidate_space(:,i) = candidate_space(:,max_index)
                candidate_space(:,max_index) = temp_det
            end if
        end do

        !select largest tdets determinants
        if (tdets <= 0) then
            target_dets = n_candidates
        else
            target_dets = min(tdets,n_candidates)
        end if

        !output as new selected space
        allocate(final_space(n_words,target_dets))
        final_space = candidate_space(:,1:target_dets)
        call move_alloc(final_space, selected_space)
        n_sel_space = target_dets

        end subroutine


end module selection