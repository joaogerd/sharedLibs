module SortingModule
    implicit none

    interface sort
        module procedure sort_int, sort_real
    end interface sort

contains

    subroutine sort_int(arr, ascending)
        integer, intent(inout) :: arr(:)
        logical, intent(in), optional :: ascending

        integer :: i, j, temp, n
        logical :: isAscending

        n = size(arr)
        isAscending = .true.

        if (present(ascending)) then
            isAscending = ascending
        end if

        do i = 1, n - 1
            do j = 1, n - i
                if ((arr(j) > arr(j + 1)).neqv. isAscending) then
                    temp = arr(j)
                    arr(j) = arr(j + 1)
                    arr(j + 1) = temp
                end if
            end do
        end do
    end subroutine sort_int

    subroutine sort_real(arr, ascending)
        real, intent(inout) :: arr(:)
        logical, intent(in), optional :: ascending

        integer :: i, j, n
        real :: temp
        logical :: isAscending

        n = size(arr)
        isAscending = .true.

        if (present(ascending)) then
            isAscending = ascending
        end if

        do i = 1, n - 1
            do j = 1, n - i
                if ((arr(j) > arr(j + 1)).neqv. isAscending) then
                    temp = arr(j)
                    arr(j) = arr(j + 1)
                    arr(j + 1) = temp
                end if
            end do
        end do
    end subroutine sort_real

end module SortingModule

