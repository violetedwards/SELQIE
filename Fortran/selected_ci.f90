!   FILE: selected_ci.f90

subroutine solve(n_orb, n_up, n_down, n_irreps, orbital_irreps, irrep_products, target_irrep, &
                 h_int, g_int, nuc_rep, energy)
!
!       MAIN SOLVING ROUTINE
!
        use constants
        use hamiltonian
        use davidson_free, only: generalized_eigensolver_free
        implicit none

        ! Declare the external determinant-excitation generator interface.
        interface
            subroutine get_connected_dets(ref_det, n_words, n_spin_orbs, order, connected_dets, n_connected)
                integer, intent(in) :: n_words, n_spin_orbs, order
                integer(kind = 8), intent(in) :: ref_det(n_words)
                integer(kind = 8), allocatable, intent(out) :: connected_dets(:,:)
                integer, intent(out) :: n_connected
            end subroutine get_connected_dets
            subroutine diagonalise_hamiltonian_from_space(n_orb, n_irreps, orbital_irreps, irrep_products, &
                    target_irrep, ref_det, n_singles, singles, n_doubles, doubles, h_int, g_int, nuc_rep, energy)
                integer, intent(in) :: n_orb, n_irreps, target_irrep, n_singles, n_doubles
                integer, intent(in) :: orbital_irreps(n_orb), irrep_products(n_irreps, n_irreps)
                integer(kind=8), intent(in) :: ref_det(:)
                integer(kind=8), allocatable, intent(in) :: singles(:,:), doubles(:,:)
                double precision, intent(in) :: h_int(2*n_orb,2*n_orb)
                double precision, intent(in) :: g_int(2*n_orb,2*n_orb,2*n_orb,2*n_orb)
                double precision, intent(in) :: nuc_rep
                double precision, intent(out) :: energy
            end subroutine diagonalise_hamiltonian_from_space
        end interface

        integer :: n_orb, n_up, n_down, n_irreps, target_irrep
        integer :: n_words, n_spin_orbs, n_singles, n_doubles
        integer :: orbital_irreps(n_orb), irrep_products(n_irreps, n_irreps)
        real(kind = dp) :: h_int(2*n_orb,2*n_orb)
        real(kind = dp) :: g_int(2*n_orb,2*n_orb,2*n_orb,2*n_orb)
        real(kind = dp) :: nuc_rep, energy
        integer(kind = 8), allocatable :: ref_det(:), singles(:,:), doubles(:,:)

!f2py   intent(in) n_orb, n_up, n_down, n_irreps, orbital_irreps, irrep_products, target_irrep, h_int, g_int, nuc_rep
!f2py   intent(out) energy

        !get determinant words ready
        n_spin_orbs = size(h_int, 1)
        n_words = ceiling(real(n_spin_orbs) / 64.0_dp)

        !get HF reference determinant
        allocate(ref_det(n_words))
        ref_det = 0_8
        call reference_determinant(n_orb, n_up, n_down, ref_det, n_words)

        !get CISD det basis
        call get_connected_dets(ref_det, n_words, n_spin_orbs, 1, singles, n_singles)
        call get_connected_dets(ref_det, n_words, n_spin_orbs, 2, doubles, n_doubles)
        
        call diagonalise_hamiltonian_from_space(n_orb, n_irreps, orbital_irreps, irrep_products, &
            target_irrep, ref_det, n_singles, singles, n_doubles, doubles, h_int, g_int, nuc_rep, energy)


            !return energy :D
        
end subroutine solve

subroutine diagonalise_hamiltonian_from_space(n_orb, n_irreps, orbital_irreps, irrep_products, &
        target_irrep, ref_det, n_singles, singles, n_doubles, doubles, h_int, g_int, nuc_rep, energy)
        use constants
        use hamiltonian
        use davidson_free, only: generalized_eigensolver_free
        implicit none

        interface
            subroutine filter_dets_by_irrep(det_basis, n_det, n_orb, n_spin_orbs, n_irreps, &
                    orbital_irreps, irrep_products, target_irrep)
                integer, intent(in) :: n_orb, n_spin_orbs, n_irreps, target_irrep
                integer, intent(in) :: orbital_irreps(n_orb), irrep_products(n_irreps, n_irreps)
                integer, intent(inout) :: n_det
                integer(kind=8), intent(inout) :: det_basis(:,:)
            end subroutine filter_dets_by_irrep
        end interface

        integer, intent(in) :: n_orb, n_irreps, target_irrep, n_singles, n_doubles
        integer, intent(in) :: orbital_irreps(n_orb), irrep_products(n_irreps, n_irreps)
        integer(kind=8), intent(in) :: ref_det(:)
        integer(kind=8), allocatable, intent(in) :: singles(:,:), doubles(:,:)
        real(kind=dp), intent(in) :: h_int(2*n_orb,2*n_orb)
        real(kind=dp), intent(in) :: g_int(2*n_orb,2*n_orb,2*n_orb,2*n_orb)
        real(kind=dp), intent(in) :: nuc_rep
        real(kind=dp), intent(out) :: energy
        integer :: n_words, n_spin_orbs, n_basis, davidson_iterations
        integer, allocatable :: row_ptr(:), column_index(:), edge_data(:,:)
        integer(kind=8), allocatable :: det_basis(:,:)
        real(kind=dp), allocatable :: eval(:), eigenvector(:,:)

        n_spin_orbs = 2 * n_orb
        n_words = size(ref_det)

        ! Assemble the reference, singles, and doubles into one determinant basis.
        n_basis = 1 + n_singles + n_doubles
        allocate(det_basis(n_words, n_basis))
        det_basis(:, 1) = ref_det
        if (n_singles > 0) det_basis(:, 2:1+n_singles) = singles(:, 1:n_singles)
        if (n_doubles > 0) det_basis(:, 2+n_singles:n_basis) = doubles(:, 1:n_doubles)

        ! Keep only determinants in the reference point-group irrep.
        call filter_dets_by_irrep(det_basis, n_basis, n_orb, n_spin_orbs, n_irreps, &
            orbital_irreps, irrep_products, target_irrep)

        ! Build sparse connectivity for matrix-free Hamiltonian products.
        call build_hamiltonian_connectivity(det_basis(:, 1:n_basis), n_basis, n_words, n_spin_orbs, &
            row_ptr, column_index, edge_data)

        allocate(eval(1), eigenvector(n_basis, 1))
        call generalized_eigensolver_free(apply_hamiltonian, eval, eigenvector, 1, 'DPR', &
            200, 1.0e-10_dp, davidson_iterations, 10, identity_operator)

        if (davidson_iterations <= 200) then
            energy = eval(1) + nuc_rep
        else
            energy = 0.0_dp
        end if

    contains

        function apply_hamiltonian(input_vectors) result(output_vectors)
            real(kind=dp), intent(in) :: input_vectors(:,:)
            real(kind=dp) :: output_vectors(size(input_vectors, 1), size(input_vectors, 2))
            output_vectors = apply_hamiltonian_matrix_free(input_vectors, det_basis, n_words, n_spin_orbs, &
                h_int, g_int, row_ptr, column_index, edge_data)
        end function apply_hamiltonian

        function identity_operator(input_vectors) result(output_vectors)
            real(kind=dp), intent(in) :: input_vectors(:,:)
            real(kind=dp) :: output_vectors(size(input_vectors, 1), size(input_vectors, 2))
            output_vectors = input_vectors
        end function identity_operator

end subroutine diagonalise_hamiltonian_from_space
!f2py END FILE selected_ci.f90