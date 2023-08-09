!------------------------------------------------------------------------------
!BOP
!
! !MODULE: dateTimeMod
!
! !DESCRIPTION: This module provides definitions for date and time manipulation.
!
! !REVISION HISTORY:
!   26 May 2023 - J. G. Z de Mattos - Initial version
!
! !SEE ALSO:
!
!EOP
!------------------------------------------------------------------------------
!BOC
module dateTimeMod
  use typeKinds
  use m_string
  implicit none
  private

  public :: dateTime
  public :: timeZone
  public :: today
  public :: timeDelta
  public :: isLeapYear
  public :: daysInYear
  public :: daysInMonth
  public :: deltaTime
  public :: gregorian
  public :: caldat
  public :: numberOfDays
  public :: dateToStr
  public :: strToDate
  public :: strptime
  public :: createDateTime
  public :: createTimeZone
  public :: ConvertToDecimal
  public :: ConvertFromDecimal
  
  TYPE timeZone
    INTEGER( kind = Long ) :: offset = 0 ! Deslocamento em relação ao UTC em minutos
    CHARACTER(LEN=10)      :: name = "UTC"    ! Nome do fuso horário, padrão é "UTC"
  contains
    procedure, public :: getName
    procedure, public :: getOffSet
  END TYPE timeZone

  type timeDelta
     integer( kind = Long ) :: year   = 0
     integer( kind = Long ) :: month  = 0
     integer( kind = Long ) :: day    = 0
     integer( kind = Long ) :: hour   = 0
     integer( kind = Long ) :: minute = 0
     integer( kind = Long ) :: second = 0
     type(timeZone)         :: timeZone
  end type timeDelta

  type dateTime
     integer( kind = Long ) :: year     = 1900
     integer( kind = Long ) :: month    = 1
     integer( kind = Long ) :: day      = 1
     integer( kind = Long ) :: hour     = 0
     integer( kind = Long ) :: minute   = 0
     integer( kind = Long ) :: second   = 0
     type(timeZone)         :: timeZone

     integer( kind = Long ) :: day_of_year
     integer( kind = Long ) :: day_of_week
   contains
     procedure, public  :: jdn => cal2jul
     procedure, public  :: jdn2 => julday2
     procedure, public  :: strftime
     procedure, public  :: yearday
     procedure, public  :: weekday
     procedure, private :: timeAdd, timePlus
     generic :: operator(+) => timeAdd, timePlus
     procedure, private :: timeSub, timeDiff
     generic :: operator(-) => timeSub, timeDiff
     procedure, private :: equal
     generic :: operator(.eq.) => equal
     procedure, private :: notEqual
     generic :: operator(.ne.) => notEqual
     procedure, private :: greaterThan
     generic :: operator(.gt.) => greaterThan
     procedure, private :: greaterEqualThan
     generic :: operator(.ge.) => greaterEqualThan
     procedure, private :: lessThan
     generic :: operator(.lt.) => lessThan
     procedure, private :: lessEqualThan
     generic :: operator(.le.) => lessEqualThan
  end type dateTime



contains

!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: timeAdd
!
! !DESCRIPTION: This function adds a time delta to a given date and time.
!
! !INTERFACE:
!
function timeAdd(a, b) result(resultTime)
  ! !INPUT PARAMETERS:
  class(dateTime), intent(in) :: a  ! Date and time
  class(timeDelta), intent(in) :: b  ! Time delta
  !
  ! !OUTPUT PARAMETERS:
  type(dateTime) :: resultTime  ! Resulting date and time
  !
  ! !REVISION HISTORY:
  !   26 May 2023 - J. G. Z de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   resultTime = timeAdd(a, b)
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  real(kind=Double) :: jdn, incr

  jdn = cal2jul(a)

  incr = (b%second + b%minute * 60.0 + b%hour * 3600.0) / 86400.0
  incr = incr + b%day + b%year * 365.0

  call jul2cal(jdn + incr, resultTime)

  resultTime%timeZone = a%timeZone

end function timeAdd
!EOC
!
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: timePlus
!
! !DESCRIPTION: This function adds two dates and times together.
!
! !INTERFACE:
!
function timePlus(a, b) result(resultTime)
  ! !INPUT PARAMETERS:
  class(dateTime), intent(in) :: a, b  ! Dates and times to be added
  !
  ! !OUTPUT PARAMETERS:
  type(dateTime) :: resultTime  ! Resulting date and time
  !
  ! !REVISION HISTORY:
  !   26 May 2023 - J. G. Z de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   resultTime = timePlus(a, b)
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  type(timeDelta) :: delta

  delta = getTimeDelta(a, b)

  resultTime = a + delta

  resultTime%timeZone = a%timeZone

end function timePlus
!
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: timeSub
!
! !DESCRIPTION: This function subtracts a time interval from a dateTime value.
!
! !INTERFACE:
!
pure function timeSub(a, b) result(resultTime)
  ! !INPUT PARAMETERS:
  class(dateTime),  intent(in) :: a  ! DateTime value
  class(timeDelta), intent(in) :: b  ! Time interval
  !
  ! !OUTPUT PARAMETERS:
  type(dateTime) :: resultTime  ! Resulting dateTime value
  !
  ! !REVISION HISTORY:
  !   26 May 2023 - J. G. Z de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   resultTime = timeSub(a, b)
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  real(kind=Double) :: jdn
  real(kind=Double) :: incr

  jdn = cal2jul(a)

  incr = (b%second + b%minute * 60.0 + b%hour * 3600.0) / 86400.0
  incr = incr + b%year * 365.0 + b%month * 30 + b%day

  call jul2cal(jdn - incr, resultTime)

  resultTime%timeZone = a%timeZone

end function timeSub
!EOC
!
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: timeDiff
!
! !DESCRIPTION: This function calculates the time difference between two dateTime values.
!
! !INTERFACE:
!
pure function timeDiff(a, b) result(resultTime)
  ! !INPUT PARAMETERS:
  class(dateTime),  intent(in) :: a  ! First dateTime value
  class(dateTime), intent(in) :: b  ! Second dateTime value
  !
  ! !OUTPUT PARAMETERS:
  type(dateTime) :: resultTime  ! Time difference as a dateTime value
  !
  ! !REVISION HISTORY:
  !   26 May 2023 - J. G. Z de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   resultTime = timeDiff(a, b)
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  type(timeDelta) :: delta

  delta = getTimeDelta(a, b)

  resultTime = a - delta

  resultTime%timeZone = a%timeZone

end function timeDiff
!EOC
!
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: equal
!
! !DESCRIPTION: This function checks if two dateTime values are equal.
!
! !INTERFACE:
!
function equal(d1, d2) result(flag)
  ! !INPUT PARAMETERS:
  class(dateTime), intent(in) :: d1  ! First dateTime value
  class(dateTime), intent(in) :: d2  ! Second dateTime value
  !
  ! !OUTPUT PARAMETERS:
  logical :: flag  ! Flag indicating if the dateTime values are equal
  !
  ! !REVISION HISTORY:
  !   26 May 2023 - J. G. Z de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   flag = equal(d1, d2)
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  flag = (cal2jul(d1) .eq. cal2jul(d2))

  return
end function equal
!EOC
!
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: notequal
!
! !DESCRIPTION: This function checks if two dateTime values are not equal.
!
! !INTERFACE:
!
function notequal(d1, d2) result(flag)
  ! !INPUT PARAMETERS:
  class(dateTime), intent(in) :: d1  ! First dateTime value
  class(dateTime), intent(in) :: d2  ! Second dateTime value
  !
  ! !OUTPUT PARAMETERS:
  logical :: flag  ! Flag indicating if the dateTime values are not equal
  !
  ! !REVISION HISTORY:
  !   26 May 2023 - J. G. Z de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   flag = notequal(d1, d2)
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  flag = (cal2jul(d1) .ne. cal2jul(d2))

  return
end function notequal
!EOC
!
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: greaterThan
!
! !DESCRIPTION: This function checks if one dateTime value is greater than another.
!
! !INTERFACE:
!
function greaterThan(d1, d2) result(flag)
  ! !INPUT PARAMETERS:
  class(dateTime), intent(in) :: d1  ! First dateTime value
  class(dateTime), intent(in) :: d2  ! Second dateTime value
  !
  ! !OUTPUT PARAMETERS:
  logical :: flag  ! Flag indicating if d1 is greater than d2
  !
  ! !REVISION HISTORY:
  !   26 May 2023 - J. G. Z de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   flag = greaterThan(d1, d2)
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  flag = (cal2jul(d1) .gt. cal2jul(d2))

  return
end function greaterThan
!EOC
!
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: greaterEqualThan
!
! !DESCRIPTION: This function checks if one dateTime value is greater than or equal to another.
!
! !INTERFACE:
!
function greaterEqualThan(d1, d2) result(flag)
  ! !INPUT PARAMETERS:
  class(dateTime), intent(in) :: d1  ! First dateTime value
  class(dateTime), intent(in) :: d2  ! Second dateTime value
  !
  ! !OUTPUT PARAMETERS:
  logical :: flag  ! Flag indicating if d1 is greater than or equal to d2
  !
  ! !REVISION HISTORY:
  !   26 May 2023 - J. G. Z de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   flag = greaterEqualThan(d1, d2)
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  flag = (cal2jul(d1) .ge. cal2jul(d2))

  return
end function greaterEqualThan
!EOC
!
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: lessThan
!
! !DESCRIPTION: This function checks if one dateTime value is less than another.
!
! !INTERFACE:
!
function lessThan(d1, d2) result(flag)
  ! !INPUT PARAMETERS:
  class(dateTime), intent(in) :: d1  ! First dateTime value
  class(dateTime), intent(in) :: d2  ! Second dateTime value
  !
  ! !OUTPUT PARAMETERS:
  logical :: flag  ! Flag indicating if d1 is less than d2
  !
  ! !REVISION HISTORY:
  !   26 May 2023 - J. G. Z de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   flag = lessThan(d1, d2)
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  flag = (cal2jul(d1) .lt. cal2jul(d2))

  return
end function lessThan
!EOC
!
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: lessEqualThan
!
! !DESCRIPTION: This function checks if one dateTime value is less than or equal to another.
!
! !INTERFACE:
!
function lessEqualThan(d1, d2) result(flag)
  ! !INPUT PARAMETERS:
  class(dateTime), intent(in) :: d1  ! First dateTime value
  class(dateTime), intent(in) :: d2  ! Second dateTime value
  !
  ! !OUTPUT PARAMETERS:
  logical :: flag  ! Flag indicating if d1 is less than or equal to d2
  !
  ! !REVISION HISTORY:
  !   26 May 2023 - J. G. Z de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   flag = lessEqualThan(d1, d2)
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  flag = (cal2jul(d1) .le. cal2jul(d2))

  return
end function lessEqualThan
!
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: today
!
! !DESCRIPTION: This function returns the current date and time as a dateTime structure.
!
! !INTERFACE:
!
function today()result(dt)
  ! !OUTPUT PARAMETERS:
  type(dateTime) :: dt  ! Current date and time
  !
  ! !REVISION HISTORY:
  !   26 May 2023 - J. G. Z de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   today = today()
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  integer(kind=Long), dimension(9) :: date

  call DATE_AND_TIME(Values=date)

  dt%year     = date(1)
  dt%month    = date(2)
  dt%day      = date(3)
  dt%hour     = date(5)
  dt%minute   = date(6)
  dt%second   = date(7)
  dt%timeZone%offset = date(4)

  dt%day_of_year = yearDay(dt)
  dt%day_of_week = weekDay(dt)

end function today
!EOC
!
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: isLeapYear
!
! !DESCRIPTION: This function determines if a year is a leap year.
!
! !INTERFACE:
!
pure elemental function isLeapYear(year) result(flag)
  ! !INPUT PARAMETERS:
  integer, intent(in) :: year  ! Year
  !
  ! !OUTPUT PARAMETERS:
  logical :: flag  ! Flag indicating if the year is a leap year
  !
  ! !REVISION HISTORY:
  !   DD Mon YYYY - J. G. Z de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   flag = isLeapYear(year)
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  flag = .false.
  if (mod(year,  4) .eq. 0 .and. &
       mod(year,100) .ne. 0 .or.  &
       mod(year,400) .eq. 0) flag = .true.

end function isLeapYear
!EOC
!
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: numberOfDays
!
! !DESCRIPTION: This function calculates the number of days between two dates, including both start and stop dates.
!
! !INTERFACE:
!
pure function numberOfDays(start, stop) result(nDays)
  ! !INPUT PARAMETERS:
  class(dateTime), intent(in) :: start, stop  ! Start and stop dates
  !
  ! !OUTPUT PARAMETERS:
  integer(kind=Long) :: nDays  ! Number of days
  !
  ! !REVISION HISTORY:
  !   26 May 2023 - J. G. Z de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   nDays = numberOfDays(start, stop)
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  integer(kind=Long) :: jdStart, jdStop

  jdStart = start%jdn()
  jdStop  = stop%jdn()

  nDays = jdStop - jdStart + 1

end function numberOfDays
!EOC
!
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: daysInYear
!
! !DESCRIPTION: This function determines the number of days in a year.
!
! !INTERFACE:
!
pure elemental function daysInYear(year) result (days)
  ! !INPUT PARAMETERS:
  integer, intent(in) :: year  ! Year
  !
  ! !OUTPUT PARAMETERS:
  integer :: days  ! Number of days in the year
  !
  ! !REVISION HISTORY:
  !   DD Mon YYYY - J. G. Z de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   days = daysInYear(year)
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  if (isLeapYear(year))then
     days = 366
  else
     days = 365
  endif
end function daysInYear
!EOC
!
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: daysInMonth
!
! !DESCRIPTION: This function determines the number of days in a month.
!
! !INTERFACE:
!
pure elemental function daysInMonth(month, year) result(days)
  ! !INPUT PARAMETERS:
  integer, intent(in) :: year   ! Year
  integer, intent(in) :: month  ! Month
  !
  ! !OUTPUT PARAMETERS:
  integer :: days  ! Number of days in the month
  !
  ! !REVISION HISTORY:
  !   DD Mon YYYY - J. G. Z de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   days = daysInMonth(month, year)
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  integer, parameter :: monthDays(*) = [31, 28, 31, 30, 31, 30, &
       31, 31, 30, 31, 30, 31]

  select case (month)
  case(:0,13:)
     days = 0
  case(2)
     if (isLeapYear(year))then
        days = 29
     else
        days = monthDays(month)
     endif
  case default
     days = monthDays(month)
  end select

end function daysInMonth
!EOC
!
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: yearDay
!
! !DESCRIPTION: This function calculates the day of the year for a given date.
!
! !INTERFACE:
!
pure function yearDay(date) result(doy)
  ! !INPUT PARAMETERS:
  class(dateTime), intent(in) :: date  ! Date and time
  !
  ! !OUTPUT PARAMETERS:
  integer :: doy  ! Day of the year
  !
  ! !REVISION HISTORY:
  !   DD Mon YYYY - J. G. Z de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   doy = yearDay(date)
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  integer :: i

  ! Initialize day of year
  doy = date%day

  ! Accumulate days of each month
  do i = 1, date%month - 1
    doy = doy + daysInMonth(i,date%year)
  end do

end function yearDay
!EOC
!
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: weekDay
!
! !DESCRIPTION: This function calculates the day of the week for a given date.
!
! !INTERFACE:
!
pure function weekDay(date) result(wday)
  ! !INPUT PARAMETERS:
  class(dateTime), intent(in) :: date  ! Date and time
  !
  ! !OUTPUT PARAMETERS:
  integer :: wday  ! Day of the week (ISO format: Monday = 1, Tuesday = 2, ..., Sunday = 7)
  !
  ! !REVISION HISTORY:
  !   DD Mon YYYY - J. G. Z de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   wday = weekDay(date)
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  integer :: century, year, month, day

  ! Extract individual components from the date
  century = date%year / 100
  year = mod(date%year, 100)
  month = date%month
  day = date%day

  ! Calculate weekday using Zeller's Congruence algorithm
  if (month < 3) then
    month = month + 12
    year = year - 1
  end if

  wday = mod(day + int(13 * (month + 1) / 5) + year + int(year / 4) + int(century / 4) - 2 * century, 7)
  wday = mod(wday + 7, 7)

  ! Convert to ISO weekday (Monday = 1, Tuesday = 2, ..., Sunday = 7)
  wday = mod(wday + 5, 7) + 1

end function weekDay
!EOC
!
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: deltaTime
!
! !DESCRIPTION: This function calculates the time difference between two dates.
!
! !INTERFACE:
!
pure elemental function deltaTime(years, months, days, hours, minutes, seconds) result(delta)
  ! !INPUT PARAMETERS:
  integer, intent(in), optional :: years, months, days, hours, minutes, seconds  ! Time components to calculate the difference
  !
  ! !OUTPUT PARAMETERS:
  type(timeDelta) :: delta  ! Time difference as a timeDelta structure
  !
  ! !REVISION HISTORY:
  !   DD Mon YYYY - J. G. Z de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   delta = deltaTime(years, months, days, hours, minutes, seconds)
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  delta%year   = merge( years  , 0, present(years))
  delta%month  = merge( months , 0, present(months))
  delta%day    = merge( days   , 0, present(days))
  delta%hour   = merge( hours  , 0, present(hours))
  delta%minute = merge( minutes, 0, present(minutes))
  delta%second = merge( seconds, 0, present(seconds))

end function deltaTime
!
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: gregorian
!
! !DESCRIPTION: This function converts a Julian day to a Gregorian date.
!
! !INTERFACE:
!
function gregorian(jd) result(dt)
  ! !INPUT PARAMETERS:
  real(kind=Double), intent(in) :: jd  ! Julian day
  !
  ! !OUTPUT PARAMETERS:
  type(dateTime) :: dt  ! Date and time in Gregorian calendar
  !
  ! !REVISION HISTORY:
  !   26 May 2023 - J. G. Z de Mattos - Initial version (Added ProTex style header)
  !
  ! !CALLING SEQUENCE:
  !   dt = gregorian(jd)
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !

  real(kind=Double) :: Z, R, G, A, B, C

  Z = floor(jd - 1721118.5)
  R = jd - 1721118.5 - Z
  G = Z - .25
  A = floor(G / 36524.25)
  B = A - floor(A / 4)
  dt%year = floor((B + G) / 365.25)
  C = B + Z - floor(365.25 * dt%year)
  dt%month = int((5 * C + 456) / 153)
  dt%day = C - int((153 * dt%month - 457) / 5) + R

  if (dt%month > 12) then
    dt%year = dt%year + 1
    dt%month = dt%month - 12
  endif

end function gregorian
!
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  cal2Jul
!
! !DESCRIPTION: This function calculates the Julian day from a given Gregorian day.
!
! !INTERFACE:
!
PURE FUNCTION cal2Jul(dt) RESULT(jd)
  !
  ! !INPUT PARAMETERS:
  !
  class(dateTime), intent(in) :: dt
  !
  ! !OUTPUT PARAMETERS:
  !
  REAL(kind=Double) :: jd  ! Julian day
  !
  ! !REVISION HISTORY: 
  !  15 Jun 2005 - J. G. de Mattos - Initial Version
  !
  ! !CALLING SEQUENCE:
  !
  !     jd = cal2Jul(dt)
  !
  ! !REMARKS:
  !   The initial Julian day in the Gregorian calendar is January 1st of the year
  !   4713 BC. This date is referred to as the "Julian Day Zero" or "Astronomical Julian Day."
  !   It is important to note that this date is used as a reference in astronomy and
  !   astronomical calculations, and it does not necessarily correspond to the start
  !   of the Gregorian calendar or other commonly used calendar systems.
  !
  !   This algorithm was adopted from Press et al.
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  integer(kind=Long) :: year
  integer(kind=Long) :: month
  integer(kind=Long) :: day
  real(kind=Double)  :: A, B, C, E, F
  real(kind=Double)  :: dayFrac
  integer(kind=LLong) :: ymd
  INTEGER(kind=Long), PARAMETER :: IGREG1 = 15 + 31 * (10 + 12 * 1582)
  INTEGER(kind=Long), PARAMETER :: IGREG2 = 04 + 31 * (10 + 12 * 1582)

  year  = dt%year
  month = dt%month
  day   = dt%day

  dayFrac = (dt%minute / 60.0 + dt%hour) / 24.0

  !-----------------------------------------------
  ! January and February are the 13th and 14th months
  ! of the previous year
  !
  if (month == 1 .or. month == 2) then
     year  = year  - 1
     month = month + 12
  endif
  !-----------------------------------------------
  
  ymd = day + 31 * (month + 12 * year)
  
  if (ymd >= IGREG1) then
     A = INT(year / 100)
     B = INT(A / 4)
     C = 2 - A + B
  endif

  if (ymd <= IGREG2) then
     C = 0
  endif

  E = INT(365.25 * (year + 4716))
  F = INT(30.6001 * (month + 1))

  jd = C + day + E + F - 1524.5 + 0.5
  jd = jd + (dt%hour / 24.0) + (dt%minute / (60.0 * 24.0)) + (dt%second / (60.0 * 60.0 * 24.0))

  RETURN

END FUNCTION cal2Jul
!
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: caldat
!
! !DESCRIPTION: This function calculates the Gregorian date from a given Julian day.
!
! !INTERFACE:
!
FUNCTION caldat(julian) RESULT(dt)
  ! !INPUT PARAMETERS:
  REAL(kind=Double), INTENT(IN) :: julian  ! Julian day
  !
  ! !OUTPUT PARAMETERS:
  TYPE(dateTime) :: dt  ! Date and time
  !
  ! !REVISION HISTORY:
  !   15 Jun 2005 - J. G. de Mattos - Initial version
  !   23 Mar 2011 - J. G. de Mattos - Modified Interface to a function call
  !   26 May 2023 - J. G. Z de Mattos - Initial version
  !
  ! !REMARKS:
  !   This algorithm was adopted from Press et al.
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  INTEGER(kind=Long), PARAMETER :: IGREG = 2299161
  INTEGER(kind=Long) :: mm, id, iyyy, ja, jd, je
  INTEGER(kind=Long) :: jalpha, jd_alpha
  INTEGER(kind=Long) :: offset

  ! Cálculo da data
  IF (julian >= IGREG) THEN
     jalpha = INT(((julian - 1867216) - 0.25_8) / 36524.25_8)
     jd_alpha = julian + 1 + jalpha - INT(0.25_8 * jalpha)
  ELSE
     jd_alpha = julian
  END IF

  ja = jd_alpha + 1524
  jd = 365 * jd_alpha + INT(0.25_8 * jd_alpha)
  je = INT((ja - jd) / 30.6001_8)
  id = ja - jd - INT(30.6001_8 * je)
  mm = je - 1
  IF (mm > 12) mm = mm - 12
  iyyy = jd_alpha - 4715
  IF (mm > 2) iyyy = iyyy - 1
  IF (iyyy <= 0) iyyy = iyyy - 1

  ! Atribui os valores calculados à estrutura dateTime
  dt%year = iyyy
  dt%month = mm
  dt%day = id
  dt%timeZone%offset = offset

END FUNCTION caldat
!
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: jul2cal
!
! !DESCRIPTION: This subroutine calculates the Gregorian date from a given Julian day.
!
! !INTERFACE:
!
pure SUBROUTINE jul2cal(jd, dt)
  ! !INPUT PARAMETERS:
  REAL(kind=Double), intent(in) :: jd  ! Julian day
  !
  ! !OUTPUT PARAMETERS:
  TYPE(dateTime), intent(out) :: dt  ! Date and time
  !
  ! !REVISION HISTORY:
  !   Date (DD Mon YYYY) - Author - Description
  !   15 Jun 2005 - J. G. de Mattos - Initial version
  !   23 Mar 2011 - J. G. de Mattos - Modified Interface to a subroutine call
  !
  ! !REMARKS:
  !   This algorithm was adopted from Press et al.
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  INTEGER(kind=Long), PARAMETER :: gregjd = 2299161
  INTEGER(kind=Long) :: year, month, day, hour, minute, second
  INTEGER(kind=Long) :: j1, j2, j3, j4, j5
  INTEGER(kind=Long) :: intgr, f, tmp
  REAL(kind=Double) :: dayfrac, frac

  ! Get the date from the Julian day number

  intgr = INT(jd)
  frac = REAL(jd - intgr, 8)

  IF (intgr >= gregjd) THEN
     tmp = INT(((intgr - 1867216) - 0.25) / 36524.25)
     j1 = intgr + 1 + tmp - INT(0.25 * tmp)
  ELSE
     j1 = intgr
  END IF

  ! Correction for half day offset
  dayfrac = frac + 0.0d0

  IF (dayfrac >= 1.0) THEN
     dayfrac = dayfrac - 1.0d0
     j1 = j1 + 1
  END IF

  j2 = j1 + 1524
  j3 = INT(6680.0 + ((j2 - 2439870) - 122.1) / 365.25)
  j4 = INT(j3 * 365.25)
  j5 = INT((j2 - j4) / 30.6001)

  day = INT(j2 - j4 - INT(j5 * 30.6001))
  month = INT(j5 - 1)

  IF (month > 12) month = month - 12

  year = INT(j3 - 4715)

  IF (month > 2) year = year - 1
  IF (year <= 0) year = year - 1

  ! Get time of day from day fraction
  hour = INT(dayfrac * 24.0d0)
  minute = INT((dayfrac * 24.0d0 - hour) * 60.0d0)
  f = NINT(((dayfrac * 24.0d0 - hour) * 60.0d0 - minute) * 60.0d0)
  second = INT(f)
  f = f - second

  IF (f > 0.5) second = second + 1

  IF (second == 60) THEN
     second = 0
     minute = minute + 1
  END IF

  IF (minute == 60) THEN
     minute = 0
     hour = hour + 1
  END IF

  IF (hour == 24) THEN
     hour = 0
     ! This could cause a bug, but probably will never happen in practice
     day = day + 1
  END IF

  IF (year < 0) THEN
     year = year * (-1)
  END IF

  dt%year = year
  dt%month = month
  dt%day = day
  dt%hour = hour
  dt%minute = minute
  dt%second = second

END SUBROUTINE jul2cal
!
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: julday2
!
! !DESCRIPTION: This function calculates the Julian day from a given date.
!
! !INTERFACE:
!
FUNCTION julday2(dt) RESULT(jd)
  ! !INPUT PARAMETERS:
  class(dateTime), intent(in) :: dt  ! Date and time
  !
  ! !OUTPUT PARAMETERS:
  REAL(kind=Double) :: jd  ! Julian day
  !
  ! !REVISION HISTORY:
  !   Date (DD Mon YYYY) - Author - Description
  !   15 Jun 2005 - J. G. de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   jd = julday2(dt)
  !
  ! !REMARKS:
  !   The initial Julian day in the Gregorian calendar is January 1st of the year
  !   4713 BC. This date is referred to as the "Julian Day Zero" or "Astronomical Julian Day."
  !   It is important to note that this date is used as a reference in astronomy and
  !   astronomical calculations, and it does not necessarily correspond to the start
  !   of the Gregorian calendar or other commonly used calendar systems.
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  INTEGER(kind=Long), PARAMETER :: IGREG = 15 + 31 * (10 + 12 * 1582)
  INTEGER(kind=Long) :: ja, jm, jy

  jy = dt%year
  if (jy == 0) print*, 'julday2: there is no year zero'
  if (jy < 0) jy = jy + 1

  if (dt%month > 2) then
     jm = dt%month + 1
  else
     jy = jy - 1
     jm = dt%month + 13
  endif

  jd = int(365.25 * jy) + int(30.6001 * jm) + dt%day + 1720995

  if (dt%day + 31 * (dt%month + 12 * dt%year) >= IGREG) then
     ja = int(0.01 * jy)
     jd = jd + 2 - ja + int(0.25 * ja)
  endif

END FUNCTION julday2
!
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: dateToStr
!
! !DESCRIPTION: This function converts a dateTime structure to a formatted string.
!
! !INTERFACE:
!
pure function dateToStr(date) result(str)
  ! !INPUT PARAMETERS:
  type(dateTime), intent(in) :: date  ! Date and time
  !
  ! !OUTPUT PARAMETERS:
  character(len=19) :: str  ! Formatted string
  !
  ! !REVISION HISTORY:
  !   26 May 2023 - J. G. Z de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   str = dateToStr(date)
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  character(len=4) :: yearStr
  character(len=2) :: monthStr, dayStr, hourStr, minuteStr, secondStr

  write(yearStr, '(I4)') date%year
  write(monthStr, '(I2)') date%month
  write(dayStr, '(I2)') date%day
  write(hourStr, '(I2)') date%hour
  write(minuteStr, '(I2)') date%minute
  write(secondStr, '(I2)') date%second

  str = trim(yearStr) // '-' // trim(monthStr) // '-' // trim(dayStr) // ' ' //&
        trim(hourStr) // ':' // trim(minuteStr) // ':' // trim(secondStr)

end function dateToStr
!
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: strToDate
!
! !DESCRIPTION: This function converts a formatted string to a dateTime structure.
!
! !INTERFACE:
!
pure function strToDate(str) result(date)
  ! !INPUT PARAMETERS:
  character(len=*), intent(in) :: str  ! Formatted string
  !
  ! !OUTPUT PARAMETERS:
  type(dateTime) :: date  ! Date and time
  !
  ! !REVISION HISTORY:
  !   26 May 2023 - J. G. Z de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   date = strToDate(str)
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  integer :: year, month, day, hour, minute, second

  read(str, *) year, month, day, hour, minute, second

  date%year     = year
  date%month    = month
  date%day      = day
  date%hour     = hour
  date%minute   = minute
  date%second   = second
  date%timeZone%offset = 0

end function strToDate
!
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: strftime
!
! !DESCRIPTION: This function formats a dateTime structure according to the specified format.
!
! !INTERFACE:
!
function strftime(date, format) result(str)
  ! !INPUT PARAMETERS:
  class(dateTime), intent(in) :: date  ! Date and time
  character(len=*), intent(in) :: format  ! Format string
  !
  ! !OUTPUT PARAMETERS:
  character(len=:), allocatable :: str  ! Formatted string
  !
  ! !REVISION HISTORY:
  !   26 May 2023 - J. G. Z de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   str = strftime(date, format)
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  integer :: i, j
  character(len=3), parameter :: monthNames(12) = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', &
                                                   'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
  character(len=3), parameter :: weekdayNames(7) = ['Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat']
  character(len=1), parameter :: weekdayDigits(7) = ['0', '1', '2', '3', '4', '5', '6']

  ! Replace format placeholders with corresponding values
  str = format
  str = replace(str, '%Y', num2str(date%year, '(I4.4)'))
  str = replace(str, '%y', num2str(mod(date%year, 100), '(I2.2)'))
  str = replace(str, '%m', num2str(date%month, '(I2.2)'))
  str = replace(str, '%d', num2str(date%day, '(I2.2)'))
  str = replace(str, '%H', num2str(date%hour, '(I2.2)'))
  str = replace(str, '%M', num2str(date%minute, '(I2.2)'))
  str = replace(str, '%S', num2str(date%second, '(I2.2)'))
  str = replace(str, '%j', num2str(date%yearDay(), '(I3.3)'))

  str = replace(str, '%b', monthNames(date%month))
  str = replace(str, '%B', adjustl(monthNames(date%month)))

  str = replace(str, '%a', weekdayNames(mod(date%weekday() + 6, 7) + 1))
  str = replace(str, '%A', adjustl(weekdayNames(mod(date%weekday() + 6, 7) + 1)))

  str = replace(str, '%w', weekdayDigits(mod(date%weekday() + 6, 7) + 1))

  str = replace(str, '%U', num2str((date%yearDay() - mod(date%weekday()+6, 7) + 7) / 7, '(I2.2)'))
  str = replace(str, '%W', num2str((date%yearDay() - mod(date%weekday()+5, 7) + 7) / 7, '(I2.2)'))

end function strftime
!
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: strptime
!
! !DESCRIPTION: This recursive function parses a formatted string and converts it to a dateTime structure.
!
! !INTERFACE:
!
recursive function strptime(datetime_string, format_string) result(date)
  USE nrtype
  USE nrutil
  USE nr
  use m_string
  !
  ! !INPUT PARAMETERS:
  character(len=*), intent(in) :: datetime_string  ! Formatted string representing date and time
  character(len=*), intent(in) :: format_string    ! Format string specifying the expected format of the datetime_string
  !
  ! !OUTPUT PARAMETERS:
  type(dateTime) :: date  ! Date and time
  !
  ! !REVISION HISTORY:
  !   26 May 2023 - J. G. Z de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   date = strptime(datetime_string, format_string)
  !
  ! !REMARKS:
  !   The function is recursive, meaning it can call itself when encountering nested formats.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !


  character(len=3), parameter :: monthNames(12) = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', &
                                                   'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

  character(len=3), parameter :: weekdayNames(7) = ['Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat']
  character(len=1), parameter :: weekdayDigits(7) = ['0', '1', '2', '3', '4', '5', '6']

  type strf
     character(len=:), allocatable :: fmt
     integer                       :: tam
  end type

  type(strf), allocatable :: f(:)
  integer,allocatable :: pos1(:), pos2(:), adj(:)
  real, allocatable, dimension(:)     :: fv, p
  integer :: i, fSize
  integer :: a, b

  type(dateTime) :: dt
  type(timeZone) :: tZ
  integer :: year, day, month
  integer :: hour, minute, second
  integer :: day_of_week, day_of_year

  ! defining format
  f = [strf( fmt='%a', tam=3 ), &
       strf( fmt='%A', tam=6 ), &
       strf( fmt='%w', tam=1 ), &
       strf( fmt='%d', tam=2 ), &
       strf( fmt='%b', tam=3 ), &
       strf( fmt='%B', tam=9 ), &
       strf( fmt='%m', tam=2 ), &
       strf( fmt='%y', tam=2 ), &
       strf( fmt='%Y', tam=4 ), &
       strf( fmt='%H', tam=2 ), &
       strf( fmt='%I', tam=2 ), &
       strf( fmt='%p', tam=2 ), &
       strf( fmt='%M', tam=2 ), &
       strf( fmt='%S', tam=2 ), &
       strf( fmt='%f', tam=6 ), &
       strf( fmt='%z', tam=5 ), &
       strf( fmt='%Z', tam=3 ), &
       strf( fmt='%j', tam=3 ), &
       strf( fmt='%U', tam=2 ), &
       strf( fmt='%W', tam=2 ), &
       strf( fmt='%c', tam=24), &
       strf( fmt='%x', tam=8 ), &
       strf( fmt='%X', tam=8 ), &
       strf( fmt='%%', tam=1 ) ]

  fSize = size(f)

  allocate(pos1(fSize))
  allocate(pos2(fSize))
  allocate(adj(fSize))

  do i = 1, fSize
     pos1(i) = index(format_string, f(i)%fmt)
  end do

  !
  ! sort
  allocate(Fv(fSize))
  allocate(p(fSize))
  p = pos1
  Fv = float([(i, i = 1, fSize)])
  Call SORT2(p, Fv)
  pos1 = int(p)

  ! ajuste da posição inicial
  adj = 0
  do i = 2, fSize
     if (pos1(i) > 0) then
        adj(i) = adj(i-1) + f(int(Fv(i-1)))%tam - 2
        pos1(i) = pos1(i) + adj(i)
     else
        adj(i) = 1
     endif
  end do

  ! reord
  p = pos1
  Call SORT2(Fv, p)
  pos1 = int(p)

  do i = 1, fSize
     pos2(i) = pos1(i) + f(i)%tam - 1
  end do

  do i = 1, fSize
     if (pos1(i) > 0) then
        a = pos1(i)
        b = pos2(i)
        select case (f(i)%fmt)
           case ('%a')
              day_of_week = maxloc(index(weekdayNames, datetime_string(a:b)), 1) - 1
           case ('%A')
              day_of_week = maxloc(index(weekdayNames, datetime_string(a:b)), 1) - 1
           case ('%w')
              day_of_week = maxloc(index(weekdayNames, datetime_string(a:b)), 1) - 1
           case ('%d')
              day = str2int(datetime_string(a:b))
           !case('%b')
           !case('%B')
           case ('%m')
              month = str2int(datetime_string(a:b))
           case ('%y')
              year = str2int(datetime_string(a:b))
           case ('%Y')
              year = str2int(datetime_string(a:b))
           case ('%H')
              hour = str2int(datetime_string(a:b))
           !case('%I')
           !case('%p')
           case ('%M')
              minute = str2int(datetime_string(a:b))
           case ('%S')
              second = str2int(datetime_string(a:b))
           !case('%f')
           !case('%z')
           case ('%Z')
              tZ%name = datetime_string(a:b)
           case ('%j')
              day_of_year = str2int(datetime_string(a:b))
           !case('%U')
           !case('%W')
           case ('%c')
              dt = strptime(datetime_string(a:b), '%d/%m/%y %H:%M:%S')
              year   = dt%year
              month  = dt%month
              day    = dt%day
              hour   = dt%hour
              minute = dt%minute
              second = dt%second
           case ('%x')
              date = strptime(datetime_string(a:b), '%d/%m/%y')
              year   = dt%year
              month  = dt%month
              day    = dt%day
           case ('%X')
              date = strptime(datetime_string(a:b), '%H:%M:%S')
              hour   = dt%hour
              minute = dt%minute
              second = dt%second
           !case('%%')
           case DEFAULT
              write(*,*) 'format not implemented yet: ', trim(f(i)%fmt)
        end select
     endif
  end do

  !
  ! Apply final values
  !

  date%year        = year
  date%month       = month
  date%day         = day
  date%hour        = hour
  date%minute      = minute
  date%second      = second
  date%timeZone    = tZ
  date%day_of_week = day_of_week
  date%day_of_year = day_of_year

end function
!
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !PURE FUNCTION: createDateTime
!
! !DESCRIPTION: This function creates a dateTime object with the provided values.
!
! !INTERFACE:
!
PURE FUNCTION createDateTime(year, month, day, hour, minute, second, tZ) RESULT(dt)
  ! !INPUT PARAMETERS:
  INTEGER(kind=Long), INTENT(IN) :: year, month, day, hour, minute, second
  TYPE(timeZone), OPTIONAL, INTENT(IN) :: tZ
  !
  ! !OUTPUT PARAMETERS:
  TYPE(dateTime) :: dt  ! Created dateTime object
  !
  ! !REVISION HISTORY:
  !   26 May 2023 - J. G. Z de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   dt = createDateTime(year, month, day, hour, minute, second, tZ)
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  dt%year   = year
  dt%month  = month
  dt%day    = day
  dt%hour   = hour
  dt%minute = minute
  dt%second = second
  dt%day_of_week = weekday(dt)
  dt%day_of_year = yearday(dt)
  IF (PRESENT(tZ)) THEN
    dt%timeZone = tZ
  END IF
  
END FUNCTION createDateTime
!
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: createTimeZone
!
! !DESCRIPTION: This function creates a timeZone object with the provided offset and name.
!
! !INTERFACE:
!
FUNCTION createTimeZone(hours, minutes, name) RESULT(tz)
  ! !INPUT PARAMETERS:
  INTEGER(kind=Long), INTENT(IN) :: hours, minutes  ! Offset from UTC in hours and minutes
  CHARACTER(*), INTENT(IN), OPTIONAL :: name  ! Time zone name (optional)
  !
  ! !OUTPUT PARAMETERS:
  TYPE(timeZone) :: tz  ! Created timeZone object
  !
  ! !REVISION HISTORY:
  !   26 May 2023 - J. G. Z de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   tz = createTimeZone(hours, minutes, name)
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  INTEGER :: offSet

  ! Assign the offset
  tz%offset = hours * 60 + minutes

  ! Check if the name is provided and assign it
  ! Calculate the offset in minutes from the provided timezone
  IF (PRESENT(name)) THEN
    offSet = calculateOffset(name)
    IF (offSet >= 0) THEN
      tz%name = 'UTC'
      tz%offSet = offSet + hours * 60 + minutes
    ELSE
      tz%name = name
    END IF
  END IF

END FUNCTION createTimeZone
!
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: calculateOffset
!
! !DESCRIPTION: This function calculates the offset in minutes from UTC based on the provided timezone name.
!
! !INTERFACE:
!
FUNCTION calculateOffset(timezone_name)result(offSet)
  ! !INPUT PARAMETERS:
  CHARACTER(len=*), INTENT(IN) :: timezone_name  ! Timezone name
  !
  ! !OUTPUT PARAMETERS:
  INTEGER(kind=Long) :: offset  ! Offset in minutes from UTC
  !
  ! !REVISION HISTORY:
  !   26 May 2023 - J. G. Z de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   offset = calculateOffset(timezone_name)
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !

  CHARACTER(len=:), ALLOCATABLE :: upper_timezone

  upper_timezone = TRANSFER(UPPERCASE(timezone_name), upper_timezone)

  SELECT CASE (TRIM(upper_timezone))
    CASE ("UTC")
      Offset = 0
    CASE ("EST")
      Offset = -5 * 60
    CASE ("PST")
      Offset = -8 * 60
    CASE ("CST")
      Offset = -6 * 60
    CASE ("IST")
      Offset = 5 * 60 + 30
    CASE ("JST")
      Offset = 9 * 60
    CASE ("BRT")
      Offset = -3 * 60
    ! Add other timezone cases as needed
    CASE DEFAULT
      WRITE(*, *) 'WARNING: UNKNOWN timeZone:', TRIM(timezone_name)
      Offset = -1
  END SELECT

END FUNCTION calculateOffset
!
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: count_digits
!
! !DESCRIPTION: This function counts the number of digits in an integer.
!
! !INTERFACE:
!
function count_digits(number) result(num_digits)
  ! !INPUT PARAMETERS:
  integer, intent(in) :: number  ! Input number
  !
  ! !OUTPUT PARAMETERS:
  integer :: num_digits  ! Number of digits in the input number
  !
  ! !REVISION HISTORY:
  !   26 May 2023 - J. G. Z de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   num_digits = count_digits(number)
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  integer :: temp_num

  temp_num = abs(number)
  num_digits = 1

  do while (temp_num >= 10)
    temp_num = temp_num / 10
    num_digits = num_digits + 1
  end do

end function count_digits
!
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: read_integer
!
! !DESCRIPTION: This function reads an integer from a string.
!
! !INTERFACE:
!
function read_integer(string) result(number)
  ! !INPUT PARAMETERS:
  character(len=*), intent(in) :: string  ! Input string
  !
  ! !OUTPUT PARAMETERS:
  integer :: number  ! Output integer
  !
  ! !REVISION HISTORY:
  !   26 May 2023 - J. G. Z de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   number = read_integer(string)
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  integer :: ierr

  read(string, *, iostat=ierr) number

  if (ierr /= 0) then
    ! Error while reading the integer
    number = -1 ! Or any other desired error value
  end if

end function read_integer
!
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: get_timezone
!
! !DESCRIPTION: This function parses a timezone string and returns the timezone information.
!
! !INTERFACE:
!
function get_timezone(tz_str) result(tz)
  ! !INPUT PARAMETERS:
  character(len=*), intent(in) :: tz_str  ! Timezone string
  !
  ! !OUTPUT PARAMETERS:
  type(timeZone) :: tz  ! Timezone information
  !
  ! !REVISION HISTORY:
  !   26 May 2023 - J. G. Z de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   tz = get_timezone(tz_str)
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  character(len=3) :: tz_sign
  integer :: tz_hours, tz_minutes

  ! Parsing the timezone string
  read(tz_str, '(*(A1))') tz_sign, tz_hours, tz_minutes

  ! Calculating the offset in minutes
  tz%offset = tz_hours * 60
  if (tz_sign == '-') tz%offset = -tz%offset
  tz%offset = tz%offset + tz_minutes

  ! Setting the timezone name
  tz%name = tz_str

end function get_timezone
!
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: str2int
!
! !DESCRIPTION: This function converts a string to an integer.
!
! !INTERFACE:
!
!function str2int(str) result(num)
!  ! !INPUT PARAMETERS:
!  character(len=*), intent(in) :: str  ! Input string
!  !
!  ! !OUTPUT PARAMETERS:
!  integer :: num  ! Output integer
!  !
!  ! !REVISION HISTORY:
!  !   26 May 2023 - J. G. Z de Mattos - Initial version
!  !
!  ! !CALLING SEQUENCE:
!  !   num = str2int(str)
!  !
!  ! !REMARKS:
!  !   Add any additional remarks here.
!  !
!  ! !SEE ALSO:
!  !
!  !EOP
!  !------------------------------------------------------------------------------
!  !BOC
!  !
!
!  read(str, *) num
!
!end function str2int
!
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: getName
!
! !DESCRIPTION: This function returns the name of a timezone.
!
! !INTERFACE:
!
function getName(tZ) result(name)
  ! !INPUT PARAMETERS:
  class(timeZone) :: tZ  ! Timezone object
  !
  ! !OUTPUT PARAMETERS:
  character(len=10) :: name  ! Timezone name
  !
  ! !REVISION HISTORY:
  !   26 May 2023 - J. G. Z de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   name = getName(tZ)
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  name = tZ%name
  return

end function getName
!
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: getOffset
!
! !DESCRIPTION: This function returns the offset of a timezone.
!
! !INTERFACE:
!
function getOffset(tZ) result(offSet)
  ! !INPUT PARAMETERS:
  class(timeZone) :: tZ  ! Timezone object
  !
  ! !OUTPUT PARAMETERS:
  integer :: offSet  ! Timezone offset
  !
  ! !REVISION HISTORY:
  !   26 May 2023 - J. G. Z de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   offSet = getOffset(tZ)
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  offSet = tZ%offset
  return

end function getOffset
!
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: ConvertToDecimal
!
! !DESCRIPTION: Converts a date and time into a decimal representation.
!
! !INTERFACE:
!
pure function ConvertToDecimal(date) result(decimal_time)
  ! !INPUT PARAMETERS:
  class(dateTime), intent(in) :: date  ! Date and time
  !
  ! !OUTPUT PARAMETERS:
  real(kind=Double) :: decimal_time  ! Decimal representation of the date and time
  !
  ! !REVISION HISTORY:
  !   31 May 2023 - J. G. Z. de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   decimal_time = ConvertToDecimal(date)
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  !

  integer :: yday, doy
  real(kind = double ) :: fhour
  real(kind = double ) :: fminute
  real(kind = double ) :: fsecond

  !
  ! Step 1: Get total number of days in that year
  !

  doy  = daysInYear(date%year)

  !
  ! The follow steps will determine the current day of the year subtracting on day 
  ! from the total number of days that have passed. The additional day accounts for 
  ! the fractional portion of the current day, which is determined by the fraction of 
  ! the hour, minute, and second.
  !

  !
  ! Step 2: Get current day of year (how many day past) and
  !         subtract by one day. 
  !

  yday = yearDay(date) - 1

  !
  ! Step 3: Divide the hour by 24
  !         To represent the hour as a fraction of a day, divide the hour value by 24. 
  !         This yields the proportion of the day that has passed. For example, if the 
  !         hour is 10, dividing by 24 gives 10/24 = 0.4167.


  fhour   = date%hour / 24.0_double

  !
  ! Step 4: Divide the minute by 1440 (24 * 60)
  !         Since there are 1440 minutes in a day (24 hours * 60 minutes), dividing the 
  !         minute value by 1440 gives the fraction of the day that has passed in terms 
  !         of minutes. For instance, if the minute is 30, dividing by 1440 gives 
  !         30/1440 = 0.0208.

  fminute = date%minute / 1440.0_double

  !
  ! Step 5: Divide the second by 86400 (24 * 60 * 60)
  !         Since there are 86400 seconds in a day (24 hours * 60 minutes * 60 seconds), 
  !         dividing the second value by 86400 gives the fraction of the day that has 
  !         passed in terms of seconds. For instance, if the second is 45, dividing by 
  !         86400 gives 45/86400 = 0.0005208333.

  fsecond = date%second / 86400.0_double

  !
  ! Step 6: sum all parts of time and divide by total number of days
  !

  decimal_time = date%year + (yday  + fhour + fminute + fsecond ) / doy


end function ConvertToDecimal

!EOC
!-----------------------------------------------------------------------------!
function ConvertFromDecimal(decimal_time)result(date)
  ! Converts a decimal representation of date and time into its components
  
  ! Inputs:
  ! decimal_time: Decimal representation of the date and time
  
  ! Outputs:
  ! year: Year component of the date
  ! month: Month component of the date
  ! day: Day component of the date
  ! hour: Hour component of the time
  ! minute: Minute component of the time
  
  real(kind=double), intent(in) :: decimal_time
  type(datetime) :: date

  real(kind=double) :: decTime, Time
  
  integer :: yday, day_accumulated, i

  ! Step 1: Extract the year component
  date%year = int(decimal_time)
  
  ! Step 2: Extract the fractional part and calculate the month

  decTime      = (decimal_time - date%year) * daysInYear(date%year)
  yday         = int(decTime)
  time         = decTime - yday
  date%hour    = int( time * 24.0_double )
  date%minute  = int( ( (time * 24.0_double) - date%hour ) * 60.0_double )
  date%second  = int( ( ( ( (time * 24.0_double) - date%hour ) * 60.0_double ) - date%minute ) * 60.0_double )

  ! Calculate the specific day
  i = 1
  day_accumulated = daysInMonth(i, date%year)
  do while (day_accumulated < yday)
      i = i + 1
      day_accumulated = day_accumulated + daysInMonth(i, date%year)
  end do

  date%month = i
  date%day   = daysInMonth(i, date%year) - (day_accumulated-yday) + 1
  
  print*,date
  
end function ConvertFromDecimal
!BOP
!
! !FUNCTION: getTimeDelta
!
! !DESCRIPTION: This function calculates the time difference between two dateTime values.
!
! !INTERFACE:
!
pure function getTimeDelta(a, b) result(delta)
  ! !INPUT PARAMETERS:
  class(dateTime), intent(in) :: a  ! Date and time A
  class(dateTime), intent(in) :: b  ! Date and time B
  !
  ! !OUTPUT PARAMETERS:
  type(timeDelta) :: delta  ! Time difference between A and B
  !
  ! !REVISION HISTORY:
  !   26 May 2023 - J. G. Z. de Mattos - Initial version
  !
  ! !CALLING SEQUENCE:
  !   delta = getTimeDelta(a, b)
  !
  ! !REMARKS:
  !   Add any additional remarks here.
  !
  ! !SEE ALSO:
  !
  !EOP
  !------------------------------------------------------------------------------
  !BOC
  !
  integer, parameter :: days_in_month(12) = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
  real(kind=Double) :: jdn_a
  real(kind=Double) :: jdn_b
  integer :: total_days, total_seconds
  integer :: remaining_days
  integer :: years, months, days
  integer :: hours, minutes, seconds

  jdn_a = cal2jul(a)
  jdn_b = cal2jul(b)

  ! Calculate the number of days between the two Julian day numbers
  total_days = abs(int(jdn_b) - int(jdn_a))

  ! Calculate the total number of seconds between the two dates
  total_seconds = (total_days) * 24 * 60 * 60 + (a%hour - b%hour) * 60 * 60 + &
                  (a%minute - b%minute) * 60 + (a%second - b%second)

  ! Calculate the number of years, months, days, hours, minutes, and seconds
  years = total_seconds / (365 * 24 * 60 * 60)
  total_seconds = total_seconds - years * (365 * 24 * 60 * 60)

  months = total_seconds / (30 * 24 * 60 * 60)
  total_seconds = total_seconds - months * (30 * 24 * 60 * 60)

  days = total_seconds / (24 * 60 * 60)
  total_seconds = total_seconds - days * (24 * 60 * 60)

  hours = total_seconds / (60 * 60)
  total_seconds = total_seconds - hours * (60 * 60)

  minutes = total_seconds / 60
  seconds = total_seconds - minutes * 60

  delta%year   = years
  delta%month  = months
  delta%day    = days
  delta%hour   = hours
  delta%minute = minutes
  delta%second = seconds

end function getTimeDelta
!EOC
!

end module dateTimeMod
