module selection
    use constants
    use hamiltonian
    use determinants
    use selection_utils
    implicit none

contains

    subroutine add_selected_dets(sel_method, selected_space, n_sel_space, n_words, old_e, old_c, n_spin_orbs, h, g, epsilon, sigma, cdets, tdets)

        !routes to correct selection algorithm
        character(len=*), intent(in) :: sel_method
        integer(kind = 8), allocatable, intent(inout) :: selected_space(:,:)
        integer, intent(inout) :: n_sel_space
        integer, intent(in) :: n_words, n_spin_orbs, cdets, tdets
        real(kind = dp), intent(in) :: old_e, old_c(:), h(:,:), g(:,:,:,:), epsilon, sigma

        if (sel_method == "hci") then
            call hci(selected_space, n_sel_space, n_words, old_c, n_spin_orbs, h, g, epsilon)
        else if (sel_method == "aci") then
            call aci(selected_space, n_sel_space, n_words, old_c, old_e, n_spin_orbs, h, g, sigma)
        else if (sel_method == "asci") then
            call asci(selected_space, n_sel_space, n_words, old_c, old_e, n_spin_orbs, h, g, cdets, tdets)
        else if (sel_method == "cipsi") then
            call cipsi(selected_space, n_sel_space, n_words, old_c, old_e, n_spin_orbs, h, g, epsilon)
        else if (sel_method == "haci") then
            call haci(selected_space, n_sel_space, n_words, old_c, old_e, n_spin_orbs, h, g, sigma)
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

        integer :: i,j,order
        integer :: n_connected, n_interacting, n_final, n_additional

        real(dp) :: old_e, hij, ei, ep, v, delta, sum_energy_importance
        real(dp), allocatable :: energy_importance(:)

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

        call sort_dets_by_score(interacting_space(:,1:n_interacting), energy_importance, n_interacting, n_words, .false.)

        !add sorted interacting space elements until remaining energy importance is <= sigma
        sum_energy_importance = sum(energy_importance)
        do i = 1, n_interacting
            if (sum_energy_importance <= sigma) exit
            call append_unique_dets(additional_space, n_additional, interacting_space(:,i:i), 1, n_words)
            sum_energy_importance = max(0.0_dp, sum_energy_importance - energy_importance(i))
        end do


        allocate(final_space(n_words,n_sel_space+n_additional))
        n_final = 0
        call append_unique_dets(final_space, n_final, selected_space(:,1:n_sel_space), n_sel_space, n_words)
        call append_unique_dets(final_space, n_final, additional_space(:,1:n_additional), n_additional, n_words)
        call move_alloc(final_space, selected_space)
        n_sel_space = n_final

        end subroutine

    subroutine cipsi(selected_space, n_sel_space, n_words, old_c, old_e, n_spin_orbs, h, g, coefficient_threshold)
        !Configuration Interaction using a Perturbative Selection made Iteratively
        integer(kind=8), allocatable, intent(inout) :: selected_space(:,:)
        integer, intent(inout) :: n_sel_space
        integer, intent(in) :: n_words, n_spin_orbs
        real(dp), intent(in) :: old_c(:), old_e, h(:,:), g(:,:,:,:), coefficient_threshold

        integer(kind=8), allocatable :: connected(:,:), external_space(:,:), final_space(:,:)
        integer :: i, j, order, n_connected, n_external, n_final
        real(dp) :: coupling, denominator, first_order_coefficient, diagonal

        allocate(external_space(n_words, max(8, n_sel_space)))
        n_external = 0

        !Build the unique space of external determinants connected to the variational space.
        do i = 1, n_sel_space
            do order = 1, 2
                call get_connected_dets(selected_space(:,i), n_words, n_spin_orbs, order, connected, n_connected)
                if (n_connected > 0) then
                    call append_unique_dets(external_space, n_external, connected(:,1:n_connected), n_connected, n_words)
                end if
                if (allocated(connected)) deallocate(connected)
            end do
        end do

        !Only external determinants contribute to the first-order wavefunction correction.
        call remove_selected_dets(external_space, n_external, selected_space(:,1:n_sel_space), n_sel_space, n_words)

        allocate(final_space(n_words, n_sel_space + n_external))
        n_final = 0
        call append_unique_dets(final_space, n_final, selected_space(:,1:n_sel_space), n_sel_space, n_words)

        do i = 1, n_external
            coupling = 0.0_dp
            do j = 1, n_sel_space
                coupling = coupling + get_matrix_element(external_space(:,i), n_words, selected_space(:,j), &
                    n_words, n_spin_orbs, h, g) * old_c(j)
            end do

            diagonal = get_matrix_element(external_space(:,i), n_words, external_space(:,i), &
                n_words, n_spin_orbs, h, g)
            denominator = old_e - diagonal
            if (abs(denominator) <= 1.0e-12_dp) then
                if (abs(coupling) > tiny(1.0_dp)) then
                    first_order_coefficient = huge(1.0_dp)
                else
                    first_order_coefficient = 0.0_dp
                end if
            else
                first_order_coefficient = coupling / denominator
            end if

            if (abs(first_order_coefficient) > coefficient_threshold) then
                call append_unique_dets(final_space, n_final, external_space(:,i:i), 1, n_words)
            end if
        end do

        call move_alloc(final_space, selected_space)
        n_sel_space = n_final
    end subroutine cipsi

    subroutine asci(selected_space, n_sel_space, n_words, old_c, old_e, n_spin_orbs, h, g, cdets, tdets)
        !Adaptive Sampling Configuration Interaction

        integer(kind=8), allocatable, intent(inout) :: selected_space(:,:)
        integer, intent(inout) :: n_sel_space

        integer, intent(in) :: n_words,n_spin_orbs,cdets,tdets

        real(dp), intent(in) :: old_e
        real(dp), intent(in) :: old_c(:)
        real(dp), intent(in) :: h(:,:)
        real(dp), intent(in) :: g(:,:,:,:)

        integer(kind=8), allocatable :: connected(:,:)
        integer(kind=8), allocatable :: interacting_space(:,:)
        integer(kind=8), allocatable :: final_space(:,:)
        integer(kind=8), allocatable :: candidate_space(:,:)
        integer(kind=8), allocatable :: core_space(:,:)

        integer :: i,j,order,max_index
        integer :: n_connected, n_interacting, n_candidates, n_core, target_dets
        integer, allocatable :: core_indices(:)

        real(dp) :: sum_h_c, denominator
        real(dp), allocatable :: pert_amp(:), candidate_coeffs(:), core_coeffs(:)

        !select the cdet determinants with largest coefficients to form the core
        n_core = n_sel_space
        if (cdets > 0) n_core = min(cdets, n_sel_space)

        allocate(core_indices(n_sel_space))
        allocate(core_space(n_words,n_core))
        allocate(core_coeffs(n_core))
        do i = 1, n_sel_space
            core_indices(i) = i
        end do

        do i = 1, n_sel_space - 1
            max_index = i
            do j = i + 1, n_sel_space
                if (abs(old_c(core_indices(j))) > abs(old_c(core_indices(max_index)))) max_index = j
            end do
            if (max_index /= i) then
                j = core_indices(i)
                core_indices(i) = core_indices(max_index)
                core_indices(max_index) = j
            end if
        end do

        do i = 1, n_core
            core_space(:,i) = selected_space(:,core_indices(i))
            core_coeffs(i) = old_c(core_indices(i))
        end do

        !search only single and double excitations of the core space
        allocate(interacting_space(n_words,max(8,n_core)))
        n_interacting = 0
        do i = 1, n_core
            do order=1,2
                call get_connected_dets(core_space(:,i), n_words, n_spin_orbs, order, connected, n_connected)
                if (n_connected > 0) then
                    call append_unique_dets(interacting_space, n_interacting, connected(:,1:n_connected), n_connected, n_words)
                end if
                if (allocated(connected)) deallocate(connected)
            end do
        end do

        !core determinants already have their current CI coefficients in the target candidates
        call remove_selected_dets(interacting_space, n_interacting, core_space(:,1:n_core), n_core, n_words)

        !get perturbed amplitudes, using only core determinants in the sum
        allocate(pert_amp(n_interacting))
        do i = 1, n_interacting
            sum_h_c = 0.0_dp
            do j = 1, n_core
                sum_h_c = sum_h_c + get_matrix_element(interacting_space(:,i),n_words, core_space(:,j), n_words, n_spin_orbs, h, g) * core_coeffs(j)
            end do
            denominator = get_matrix_element(interacting_space(:,i),n_words,interacting_space(:,i),n_words,n_spin_orbs,h,g)-old_e
            if (abs(denominator) > 1.0e-10_dp) then
                pert_amp(i) = sum_h_c/denominator
            else
                pert_amp(i) = 0.0_dp
            end if
        end do

        !combine core coefficients and perturbed search amplitudes
        n_candidates = n_core + n_interacting
        allocate(candidate_space(n_words, n_candidates))
        allocate(candidate_coeffs(n_candidates))

        candidate_space(:,1:n_core) = core_space(:,1:n_core)
        candidate_coeffs(1:n_core) = core_coeffs(1:n_core)
        if (n_interacting > 0) then
            candidate_space(:,n_core+1:n_candidates) = interacting_space(:,1:n_interacting)
            candidate_coeffs(n_core+1:n_candidates) = pert_amp(1:n_interacting)
        end if

        call sort_dets_by_score(candidate_space, candidate_coeffs, n_candidates, n_words, .true.)

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

    subroutine haci(selected_space, n_sel_space, n_words, old_c, old_e, n_spin_orbs, h, g, sigma)
        !Heat Bath Adaptive Configuration Interaction

        integer(kind=8), allocatable, intent(inout) :: selected_space(:,:)
        integer, intent(inout) :: n_sel_space

        integer, intent(in) :: n_words
        integer, intent(in) :: n_spin_orbs

        real(dp), intent(in) :: old_e
        real(dp), intent(in) :: old_c(:)
        real(dp), intent(in) :: h(:,:)
        real(dp), intent(in) :: g(:,:,:,:)
        real(dp), intent(in) :: sigma

        integer(kind=8), allocatable :: connected(:,:)
        integer(kind=8), allocatable :: interacting_space(:,:)
        integer(kind=8), allocatable :: final_space(:,:)
        integer(kind=8), allocatable :: additional_space(:,:)

        integer :: i, j, order
        integer :: n_connected, n_interacting, n_final, n_additional

        real(dp) :: epsilon, hij, ep, v, delta, sum_energy_importance
        real(dp), allocatable :: energy_importance(:)

        epsilon = sigma

        !generate unique external determinants that pass the HCI heat bath threshold.
        allocate(interacting_space(n_words, max(8, n_sel_space)))
        n_interacting = 0
        do i = 1, n_sel_space
            do order = 1,2
                call get_connected_dets(selected_space(:,i), n_words, n_spin_orbs, order, connected, n_connected)
                do j = 1, n_connected
                    hij = get_matrix_element(connected(:,j), n_words, selected_space(:,i), n_words, n_spin_orbs, h, g)
                    if(abs(hij*old_c(i)) > epsilon) then
                        call append_unique_dets(interacting_space, n_interacting, connected(:,j:j), 1, n_words)
                    end if
                end do
                if(allocated(connected)) deallocate(connected)
            end do
        end do

        call remove_selected_dets(interacting_space, n_interacting, selected_space(:,1:n_sel_space), n_sel_space, n_words)

        !rank candidates by the ACI energy lowering from mixing with the current state.
        allocate(energy_importance(n_interacting))
        do i = 1, n_interacting
            ep = get_matrix_element(interacting_space(:,i), n_words, interacting_space(:,i), n_words, n_spin_orbs, h, g)
            v = 0.0_dp
            do j = 1, n_sel_space
                hij = get_matrix_element(interacting_space(:,i), n_words, selected_space(:,j), n_words, n_spin_orbs, h, g)
                v = v + hij * old_c(j)
            end do
            delta = ep - old_e
            energy_importance(i) = abs(delta/2 - sqrt((delta*delta/4)+(v*v)))
        end do

        call sort_dets_by_score(interacting_space(:,1:n_interacting), energy_importance, n_interacting, n_words, .false.)

        !add the highest-importance candidates until the omitted importance is <= sigma.
        allocate(additional_space(n_words, max(8, n_sel_space)))
        n_additional = 0
        sum_energy_importance = sum(energy_importance)
        do i = 1, n_interacting
            if (sum_energy_importance <= sigma) exit
            call append_unique_dets(additional_space, n_additional, interacting_space(:,i:i), 1, n_words)
            sum_energy_importance = max(0.0_dp, sum_energy_importance - energy_importance(i))
        end do

        allocate(final_space(n_words, n_sel_space + n_additional))
        n_final = 0
        call append_unique_dets(final_space, n_final, selected_space(:,1:n_sel_space), n_sel_space, n_words)
        call append_unique_dets(final_space, n_final, additional_space(:,1:n_additional), n_additional, n_words)
        call move_alloc(final_space, selected_space)
        n_sel_space = n_final

    end subroutine haci

end module selection