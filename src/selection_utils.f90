module selection_utils
    use constants, only: dp
    implicit none

contains

    logical function any_selected_det(det, dets, n_dets, n_words)
        integer(kind=8), intent(in) :: det(n_words), dets(n_words,n_dets)
        integer, intent(in) :: n_dets, n_words
        integer :: i

        any_selected_det = .false.
        do i = 1, n_dets
            if (all(det == dets(:,i))) then
                any_selected_det = .true.
                return
            end if
        end do
    end function any_selected_det

    subroutine remove_selected_dets(candidate_dets, n_candidates, selected_dets, n_selected, n_words)
        integer, intent(in) :: n_selected, n_words
        integer, intent(inout) :: n_candidates
        integer(kind=8), intent(inout) :: candidate_dets(:,:)
        integer(kind=8), intent(in) :: selected_dets(:,:)
        integer :: i, j

        i = 1
        do while (i <= n_candidates)
            if (any_selected_det(candidate_dets(:,i), selected_dets(:,1:n_selected), n_selected, n_words)) then
                do j = i, n_candidates - 1
                    candidate_dets(:,j) = candidate_dets(:,j+1)
                end do
                n_candidates = n_candidates - 1
            else
                i = i + 1
            end if
        end do
    end subroutine remove_selected_dets

    subroutine sort_dets_by_score(dets, scores, n_dets, n_words, sort_by_absolute_score)
        integer, intent(in) :: n_dets, n_words
        logical, intent(in) :: sort_by_absolute_score
        integer(kind=8), intent(inout) :: dets(:,:)
        real(dp), intent(inout) :: scores(:)
        integer :: i, heap_size

        !Build a min-heap, then move each minimum to the end for descending order.
        do i = n_dets / 2, 1, -1
            call sift_down(dets, scores, i, n_dets, n_words, sort_by_absolute_score)
        end do

        do heap_size = n_dets, 2, -1
            call swap_det_score(dets, scores, 1, heap_size, n_words)
            call sift_down(dets, scores, 1, heap_size - 1, n_words, sort_by_absolute_score)
        end do
    end subroutine sort_dets_by_score

    subroutine sift_down(dets, scores, root, heap_size, n_words, sort_by_absolute_score)
        integer, intent(in) :: root, heap_size, n_words
        logical, intent(in) :: sort_by_absolute_score
        integer(kind=8), intent(inout) :: dets(:,:)
        real(dp), intent(inout) :: scores(:)
        integer :: child, current_root

        current_root = root
        do while (current_root <= heap_size / 2)
            child = 2 * current_root
            if (child < heap_size) then
                if (score_less(scores(child + 1), scores(child), sort_by_absolute_score)) child = child + 1
            end if
            if (.not. score_less(scores(child), scores(current_root), sort_by_absolute_score)) exit
            call swap_det_score(dets, scores, current_root, child, n_words)
            current_root = child
        end do
    end subroutine sift_down

    pure logical function score_less(score_a, score_b, sort_by_absolute_score)
        real(dp), intent(in) :: score_a, score_b
        logical, intent(in) :: sort_by_absolute_score

        if (sort_by_absolute_score) then
            score_less = abs(score_a) < abs(score_b)
        else
            score_less = score_a < score_b
        end if
    end function score_less

    subroutine swap_det_score(dets, scores, index_a, index_b, n_words)
        integer, intent(in) :: index_a, index_b, n_words
        integer(kind=8), intent(inout) :: dets(:,:)
        real(dp), intent(inout) :: scores(:)
        integer(kind=8) :: temp_det(n_words)
        real(dp) :: temp_score

        temp_det = dets(:,index_a)
        dets(:,index_a) = dets(:,index_b)
        dets(:,index_b) = temp_det
        temp_score = scores(index_a)
        scores(index_a) = scores(index_b)
        scores(index_b) = temp_score
    end subroutine swap_det_score

end module selection_utils
