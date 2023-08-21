program DateToSeconds
    implicit none
    
    integer :: year, month, day, hour, minute, second
    real :: julian_day, seconds
    
    ! Input the date
    write(*,*) "Enter the date (YYYY MM DD HH MM SS):"
    read(*,*) year, month, day, hour, minute, second
    
    ! Calculate Julian day
    julian_day = GetJulianDay(year, month, day)
    
    ! Convert Julian day to seconds
    seconds = JulianDayToSeconds(julian_day, hour, minute, second)
    
    ! Display the result
    write(*,*) "Number of seconds from Julian day:", seconds
    
    ! Convert seconds back to the date
    Call SecondsToDate(seconds, year, month, day, hour, minute, second)
    
    ! Display the reverse result
    write(*,*) "Reversed Date (YYYY MM DD HH MM SS):", year, month, day, hour, minute, second
    
contains

    function GetJulianDay(year, month, day) result(julian_day)
        implicit none
        integer, intent(in) :: year, month, day
        real :: julian_day
        
        integer :: a, y, m
        
        a = (14 - month) / 12
        y = year - a
        m = month + 12 * a - 3
        
        julian_day = day + (153 * m + 2) / 5 + 365 * y + y / 4 - 32045
        
    end function GetJulianDay

    function JulianDayToSeconds(julian_day, hour, minute, second) result(seconds)
        implicit none
        real, intent(in) :: julian_day
        integer, intent(in) :: hour, minute, second
        real :: seconds
        
        seconds = (julian_day - 2440587.5) * 24.0 * 60.0 * 60.0 + hour * 60.0 * 60.0 + minute * 60.0 + second
        
    end function JulianDayToSeconds
    
    subroutine SecondsToDate(seconds, year, month, day, hour, minute, second)
        implicit none
        real, intent(in) :: seconds
        integer, intent(out) :: year, month, day, hour, minute, second
        
        real :: julian_day, fractional_day
        integer :: a, b, c, d, e, z, alpha
        
        fractional_day = seconds / (24.0 * 60.0 * 60.0)
        
        z = int(fractional_day)
        alpha = int((fractional_day - z) * 24.0 * 60.0 * 60.0)
        
        a = z + 32044
        b = (4 * a + 3) / 146097
        c = a - (146097 * b) / 4
        d = (4 * c + 3) / 1461
        e = c - (1461 * d) / 4
        
        month = (5 * e + 2) / 153
        day = e - (153 * month + 2) / 5 + 1
        year = b * 100 + d - 4800 + (month + 2) / 12
        
        hour = alpha / (60 * 60)
        minute = mod((alpha / 60), 60)
        second = mod(alpha, 60)
        
    end subroutine SecondsToDate
    
end program DateToSeconds

