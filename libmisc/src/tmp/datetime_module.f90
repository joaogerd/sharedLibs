module datetime_module
    use iso_fortran_env
    implicit none

    type :: datetime
        integer :: year
        integer :: month
        integer :: day
        integer :: hour
        integer :: minute
        real :: second
    end type datetime

    interface operator(-)
        module procedure datetime_difference
    end interface operator(-)

    contains

    function datetime_difference(d1, d2) result(diff)
        type(datetime), intent(in) :: d1, d2
        type(datetime) :: diff

        diff%year = d2%year - d1%year
        diff%month = d2%month - d1%month
        diff%day = d2%day - d1%day
        diff%hour = d2%hour - d1%hour
        diff%minute = d2%minute - d1%minute
        diff%second = d2%second - d1%second

        if (diff%second < 0.0) then
            diff%second = diff%second + 60.0
            diff%minute = diff%minute - 1
        endif
        if (diff%minute < 0) then
            diff%minute = diff%minute + 60
            diff%hour = diff%hour - 1
        endif
        if (diff%hour < 0) then
            diff%hour = diff%hour + 24
            diff%day = diff%day - 1
        endif
        if (diff%day < 0) then
            diff%month = diff%month - 1
            if (diff%month == 0) then
                diff%month = 12
                diff%year = diff%year - 1
            endif
            diff%day = diff%day + month_length(diff%year, diff%month)
        endif
        if (diff%month < 0) then
            diff%month = diff%month + 12
            diff%year = diff%year - 1
        endif

    end function datetime_difference

    function month_length(year, month) result(length)
        integer, intent(in) :: year, month
        integer :: length

        select case (month)
        case (1, 3, 5, 7, 8, 10, 12)
            length = 31
        case (4, 6, 9, 11)
            length = 30
        case (2)
            if (mod(year, 4) == 0 .and. (mod(year, 100) /= 0 .or. mod(year, 400) == 0)) then
                length = 29
            else
                length = 28
            endif
        end select

    end function month_length

    function total_seconds(d) result(seconds)
        type(datetime), intent(in) :: d
        real :: seconds

        seconds = real(d%year * 365 * 24 * 60 * 60 + d%month * month_length(d%year, d%month) * 24 * 60 * 60 &
                    + d%day * 24 * 60 * 60 + d%hour * 60 * 60 + d%minute * 60 + d%second, real64)

    end function total_seconds

function t2dt(atime) result(dt)
    real, intent(in) :: atime
    type(datetime) :: dt
    integer :: year
    real :: remainder

    year = int(atime)
    remainder = atime - real(year)

    dt%year = year
    dt%month = 1
    dt%day = 1
    dt%hour = 0
    dt%minute = 0
    dt%second = 0

    dt%second = dt%second + remainder * 60.0 * 60.0 * 24.0 * real(month_length(dt%year, dt%month))

    if (dt%second >= 60.0) then
        dt%second = dt%second - 60.0
        dt%minute = dt%minute + 1
    endif
    if (dt%minute >= 60) then
        dt%minute = dt%minute - 60
        dt%hour = dt%hour + 1
    endif
    if (dt%hour >= 24) then
        dt%hour = dt%hour - 24
        dt%day = dt%day + 1
    endif

    if (dt%day > month_length(dt%year, dt%month)) then
        dt%day = 1
        dt%month = dt%month + 1
        if (dt%month > 12) then
            dt%month = 1
            dt%year = dt%year + 1
        endif
    endif

end function t2dt
end module datetime_module

