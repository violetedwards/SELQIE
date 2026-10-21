module hash

    implicit none

    type det_hash_t
        integer :: table_size
        integer, allocatable :: head(:)
        integer, allocatable :: next(:)
    end type det_hash_t

contains

    integer(kind=8) function hash_det(det,n_words,table_size)

        implicit none

        integer, intent(in) :: n_words
        integer, intent(in) :: table_size

        integer(kind=8), intent(in) :: det(n_words)

        integer(kind=8) :: h
        integer :: i

        h = 1469598103934665603_8

        do i = 1, n_words
            h = ieor(h,det(i))
            h = h * 1099511628211_8
        end do

        hash_det = modulo(h,int(table_size,8)) + 1_8

    end function hash_det



    subroutine build_det_hash(hash_table,basis,n_words,n_basis)

        implicit none

        type(det_hash_t), intent(out) :: hash_table

        integer, intent(in) :: n_words
        integer, intent(in) :: n_basis

        integer(kind=8), intent(in) :: basis(n_words,n_basis)

        integer :: i
        integer :: bucket

        hash_table%table_size = max(17,4*n_basis+1)

        allocate(hash_table%head(hash_table%table_size))
        allocate(hash_table%next(n_basis))

        hash_table%head = 0
        hash_table%next = 0

        do i = 1,n_basis

            bucket = int(hash_det( &
                basis(:,i), &
                n_words, &
                hash_table%table_size))

            hash_table%next(i) = hash_table%head(bucket)

            hash_table%head(bucket) = i

        end do

    end subroutine build_det_hash



    integer function lookup_det_index(det,basis,n_words,hash_table)

        implicit none

        integer, intent(in) :: n_words

        integer(kind=8), intent(in) :: det(n_words)

        integer(kind=8), intent(in) :: basis(:,:)

        type(det_hash_t), intent(in) :: hash_table

        integer :: bucket
        integer :: idx

        bucket = int(hash_det( &
            det, &
            n_words, &
            hash_table%table_size))

        idx = hash_table%head(bucket)

        do while(idx /= 0)

            if(all(det == basis(:,idx))) then
                lookup_det_index = idx
                return
            end if

            idx = hash_table%next(idx)

        end do

        lookup_det_index = 0

    end function lookup_det_index



    subroutine destroy_det_hash(hash_table)

        implicit none

        type(det_hash_t), intent(inout) :: hash_table

        if(allocated(hash_table%head)) deallocate(hash_table%head)

        if(allocated(hash_table%next)) deallocate(hash_table%next)

        hash_table%table_size = 0

    end subroutine destroy_det_hash

end module hash