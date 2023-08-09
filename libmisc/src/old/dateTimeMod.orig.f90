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

  type(dateTime) function today()

    integer( kind = Long ) :: date(9)

    call DATE_AND_TIME(Values=date)

    today%year     = date(1)
    today%month    = date(2)
    today%day      = date(3)
    today%hour     = date(5)
    today%minute   = date(6)
    today%second   = date(7)
    today%timeZone%offset = date(4)

    today%day_of_year = yearDay(today)
    today%day_of_week = weekDay(today)

    return
  end function today

  pure elemental function isLeapYear(year) result(flag)
    integer, intent(in) :: year
    logical :: flag

    flag = .false.
    if (mod(year,  4) .eq. 0 .and. &
         mod(year,100) .ne. 0 .or.  &
         mod(year,400) .eq. 0) flag = .true.

  end function isLeapYear

  pure elemental function daysInYear(year) result (days)
    integer, intent(in) :: year
    integer :: days

    if (isLeapYear(year))then
       days = 366
    else
       days = 365
    endif
  end function daysInYear


  pure elemental function daysInMonth(month, year) result(days)
    integer, intent(in   ) :: year
    integer, intent(in   ) :: month
    integer                :: days

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
    return
  end function daysInMonth

  pure function yearDay(date) result(doy)
    class(dateTime), intent(in) :: date
    integer :: doy, i
  
    ! Initialize day of year
    doy = date%day
  
    ! Accumulate days of each month
    do i = 1, date%month - 1
      doy = doy + daysInMonth(i,date%year)
    end do
  
  end function yearDay
  
  pure function weekDay(date) result(wday)
    class(dateTime), intent(in) :: date
    integer :: wday, century, year, month, day
  
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
  
  end function weekday


  pure elemental function deltaTime(years, months, days, hours, minutes, seconds) result(delta)
      integer, intent(in), optional :: years, months, days, hours, minutes, seconds
      type(timeDelta) :: delta
  
      delta%year   = merge( years  , 0, present(years))
      delta%month  = merge( months , 0, present(months))
      delta%day    = merge( days   , 0, present(days))
      delta%hour   = merge( hours  , 0, present(hours))
      delta%minute = merge( minutes, 0, present(minutes))
      delta%second = merge( seconds, 0, present(seconds))
  
  end function deltaTime

  pure function getTimeDelta(a,b)result(delta)
     class(dateTime),  intent(in) :: a
     class(dateTime), intent(in) :: b
     type(timeDelta) :: delta
     integer, parameter :: days_in_month(12) = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
     real( kind = Double ) :: jdn_a
     real( kind = Double ) :: jdn_b
     integer :: total_days, total_seconds
     integer :: remaining_days
     integer :: years, months, days
     integer :: hours, minutes, seconds

     jdn_a = cal2jul(a)
     jdn_b = cal2jul(b)

     ! Calculate the number of days between the two Julian day numbers
     total_days = abs(int(jdn_b)-int(jdn_a))

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

  pure function numberOfDays(start, stop) result(nDays)
    class(dateTime), intent(in) :: start, stop
    integer( kind = Long ) :: nDays

    integer( kind = Long ) :: jdStart, jdStop

    jdStart = start%jdn()
    jdStop  = stop%jdn()

    nDays = jdStop - jdStart + 1

  end function numberOfDays

  function timeAdd(a, b) result(resultTime)
    class(dateTime),  intent(in) :: a
    class(timeDelta), intent(in) :: b
    type(dateTime) :: resultTime

    real( kind = Double ) :: jdn
    real( kind = Double ) :: incr

    jdn = cal2jul(a)

    incr = (b%second + b%minute * 60.0 + b%hour * 3600.0) / 86400.0
    incr = incr + b%day + b%year*365.0

    call jul2cal( jdn + incr, resultTime)

    resultTime%timeZone = a%timeZone

  end function timeAdd

  function timePlus(a, b) result(resultTime)
    class(dateTime),  intent(in) :: a
    class(dateTime), intent(in) :: b
    type(dateTime)  :: resultTime

    type(timeDelta) :: delta

    delta = getTimeDelta(a,b)

    resultTime = a + delta

    resultTime%timeZone = a%timeZone


  end function timePlus


  pure function timeSub(a, b) result(resultTime)
    class(dateTime),  intent(in) :: a
    class(timeDelta), intent(in) :: b
    type(dateTime) :: resultTime

    real( kind = Double ) :: jdn
    real( kind = Double ) :: incr

    jdn = cal2jul(a)

    incr = (b%second + b%minute * 60.0 + b%hour * 3600.0) / 86400.0
    incr = incr + b%year*365.0 + b%month*30 + b%day

    call jul2cal( jdn - incr, resultTime)

    resultTime%timeZone = a%timeZone

  end function timeSub


  pure function timeDiff(a, b) result(resultTime)
    class(dateTime),  intent(in) :: a
    class(dateTime), intent(in) :: b
    type(dateTime) :: resultTime

    type(timeDelta) :: delta

    delta = getTimeDelta(a,b)

    resultTime = a - delta

    resultTime%timeZone = a%timeZone


  end function timeDiff


  function equal(d1,d2) result(flag)
    class(dateTime), intent(in) :: d1
    class(dateTime), intent(in) :: d2
    logical :: flag

    flag = (cal2jul(d1) .eq. cal2jul(d2))

    return
  end function equal

  function notequal(d1,d2) result(flag)
    class(dateTime), intent(in) :: d1
    class(dateTime), intent(in) :: d2
    logical :: flag

    flag = (cal2jul(d1) .ne. cal2jul(d2))

    return
  end function notequal

  function greaterThan(d1,d2) result(flag)
    class(dateTime), intent(in) :: d1
    class(dateTime), intent(in) :: d2
    logical :: flag

    flag = (cal2jul(d1) .gt. cal2jul(d2))

    return
  end function greaterThan

  function greaterEqualThan(d1,d2) result(flag)
    class(dateTime), intent(in) :: d1
    class(dateTime), intent(in) :: d2
    logical :: flag

    flag = (cal2jul(d1) .ge. cal2jul(d2))

    return
  end function greaterEqualThan

  function lessThan(d1,d2) result(flag)
    class(dateTime), intent(in) :: d1
    class(dateTime), intent(in) :: d2
    logical :: flag

    flag = (cal2jul(d1) .lt. cal2jul(d2))

    return
  end function lessThan

  function lessEqualThan(d1,d2) result(flag)
    class(dateTime), intent(in) :: d1
    class(dateTime), intent(in) :: d2
    logical :: flag

    flag = (cal2jul(d1) .le. cal2jul(d2))

    return
  end function lessEqualThan


  function gregorian(jd)result(dt)
    real(kind = Double), intent(in) :: jd
    type(dateTime) :: dt

    real(kind = Double) :: Z, R, G, A, B, C

    Z = floor(jd - 1721118.5)
    R = jd - 1721118.5 - Z
    G = Z - .25
    A = floor(G / 36524.25)
    B = A - floor(A / 4)
    dt%year = floor((B+G) / 365.25)
    C = B + Z - floor(365.25 * dt%year)
    dt%month = int((5 * C + 456) / 153)
    dt%day   = C - int((153 * dt%month - 457) / 5) + R

    if (dt%month > 12)then
       dt%year  = dt%year + 1
       dt%month = dt%month - 12
    endif

    return
  end function gregorian



  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !       NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: lower_case - convert uppercase letters to lowercase.
  !
  ! !DESCRIPTION:
  !
  ! !INTERFACE:

  function lower_case(str) result(lstr)
    implicit none
    character(len=*), intent(in) :: str
    character(len=len(str))      :: lstr

    ! !REVISION HISTORY:
    !       13Aug96 - J. Guo        - (to do)
    !EOP
    !_______________________________________________________________________
    integer i
    integer,parameter :: iu2l=ichar('a')-ichar('A')

    lstr=str
    do i=1,len_trim(str)
       if(str(i:i).ge.'A'.and.str(i:i).le.'Z')   &
            lstr(i:i)=char(ichar(str(i:i))+iu2l)
    end do
  end function lower_case


  !-----------------------------------------------------------------------------!
  !
  !BOP
  !
  ! !IROUTINE:  cal2Jul
  !
  ! !DESCRIPTION: This function calculate the julian day from gregorian day
  !
  !
  ! !INTERFACE:
  !
  PURE FUNCTION cal2jul(dt) RESULT(jd)
    !
    ! !INPUT PARAMETERS:
    !
    class(dateTime), intent(in) :: dt
    !
    ! !OUTPUT PARAMETERS:
    !
    REAL(kind = Double)  :: jd
    !
    !
    ! !REVISION HISTORY: 
    !  15 Jun 2005 - J. G. de Mattos - Initial Version
    !
    !
    ! !CALLING SEQUENCE:
    !
    !     jd = dt%cal2jul()
    !
    ! !SEE ALSO:
    !
    !
    !
    ! !NEMARKS:
    !       O dia juliano inicial em data gregoriana é o dia 1º de janeiro 
    !    do ano 4713 a.C. Essa data é referenciada como o "Dia Juliano Zero" 
    !    ou "Dia Juliano Astronômico". É importante destacar que essa data é 
    !    usada como referência em astronomia e cálculos astronômicos, e não 
    !    necessariamente corresponde ao início do calendário gregoriano ou outros 
    !    sistemas de calendário comumente utilizados.
    !
    !    This algorithm was adopted from Press et al.
    !
    !EOP
    !-----------------------------------------------------------------------------!
    !BOC
    !
    integer(kind = Long) :: year
    integer(kind = Long) :: month
    integer(kind = Long) :: day

    real(kind = Double)    :: A, B, C, E, F
    real(kind = Double)    :: dayFrac
    integer(kind = LLong)  :: ymd
    INTEGER(kind = Long), PARAMETER :: IGREG1=15+31*(10+12*1582)
    INTEGER(kind = Long), PARAMETER :: IGREG2=04+31*(10+12*1582)

    year  = dt%year
    month = dt%month
    day   = dt%day

    dayFrac = (dt%minute/60.0 + dt%hour)/24.0

    !-----------------------------------------------
    ! January and February are 13th and 14th months
    ! of the previous year
    !
    if(month .eq. 1 .or. month .eq. 2)then
       year  = year  -  1
       month = month + 12
    endif
    !-----------------------------------------------
    ymd = day+31*(month+12*year)
    if(ymd>=IGREG1)then
       A = INT(year/100)
       B = INT(A/4)
       C = 2 - A + B
    endif

    if(ymd<=IGREG2)then
       c = 0
    endif

    E = INT(365.25 * (year + 4716))
    F = INT(30.6001 * (month + 1))

    jd = C + day + E + F - 1524.5 + 0.5

    jd = jd + (dt%hour/24.0)
    jd = jd + (dt%minute/(60.0*24.0))
    jd = jd + (dt%second/(60.0*60.0*24.0))

    RETURN

  END FUNCTION cal2jul
  !
  !EOC
  !-----------------------------------------------------------------------------!


  FUNCTION julday2(dt) result(jd)
    class(dateTime), intent(in) :: dt
    REAL(kind = Double) :: jd
    INTEGER(kind = Long), PARAMETER :: IGREG=15+31*(10+12*1582)
    INTEGER(kind = Long) :: ja,jm,jy
    jy=dt%year
    if (jy == 0) print*,'julday: there is no year zero'
    if (jy < 0) jy=jy+1
    if (dt%month > 2) then
       jm=dt%month+1
    else
       jy=jy-1
       jm=dt%month+13
    end if
    jd=int(365.25*jy)+int(30.6001*jm)+dt%day+1720995
    if (dt%day+31*(dt%month+12*dt%year) >= IGREG) then
       ja=int(0.01*jy)
       jd=jd+2-ja+int(0.25*ja)
    end if
  END FUNCTION julday2

  !-----------------------------------------------------------------------------!
  !BOP
  !
  ! !IROUTINE:  Jul2Cal
  !
  ! !DESCRIPTION: This function calculate the gregorian date from julian day.
  !
  !
  ! !INTERFACE:
  !
  pure SUBROUTINE jul2cal(jd, dt)
    !
    ! !INPUT PARAMETERS:
    !
    real (kind = Double), intent(in   ) :: jd
    !
    ! !OUTPUT PARAMETERS:
    !
    type(dateTime), intent(  out) :: dt
    !
    !
    ! !REVISION HISTORY: 
    !  15 Jun 2005 - J. G. de Mattos - Initial Version
    !  23 Mar 2011 - J. G. de Mattos - Modified Interface
    !                                  to a subroutine call
    !
    ! !REMARKS:
    !        This algorithm was adopted from Press et al.
    !EOP
    !-----------------------------------------------------------------------------!
    !BOC
    !
    integer (kind = Long), parameter :: gregjd = 2299161
    integer (kind = Long)            :: year
    integer (kind = Long)            :: month
    integer (kind = Long)            :: day
    integer (kind = Long)            :: hour
    integer (kind = Long)            :: minute
    integer (kind = Long)            :: second
    integer (kind = Long)            :: j1, j2, j3, j4, j5
    integer (kind = Long)            :: intgr
    integer (kind = Long)            :: f, tmp

    real (kind = Double)             :: dayfrac, frac

    !       
    ! get the date from the Julian day number
    !       
    ! jd=2453372.25

    intgr   = INT(jd)
    frac    = real(jd - intgr,8)

    IF( intgr >= gregjd )THEN              !Gregorian calendar correction
       tmp = INT( ( (intgr - 1867216) - 0.25 ) / 36524.25 )
       j1  = intgr + 1 + tmp - INT(0.25*tmp)
    ELSE
       j1 = intgr
    ENDIF
    !       correction for half day offset

    dayfrac = frac + 0.0d0
    IF( dayfrac >= 1.0 )THEN
       dayfrac = dayfrac - 1.0d0
       j1 = j1+1
    ENDIF

    j2 = j1 + 1524
    j3 = INT( 6680.0 + ( (j2 - 2439870) - 122.1 )/365.25 )
    j4 = INT(j3*365.25)
    j5 = INT( (j2 - j4)/30.6001 )

    day   = INT(j2 - j4 - INT(j5*30.6001))
    month = INT(j5 - 1)
    IF( month > 12 ) month = month- 12
    year = INT(j3 - 4715)
    IF( month > 2 )  year = year - 1
    IF( year <= 0 ) year = year - 1

    !
    ! get time of day from day fraction
    !
    hour   = INT(dayfrac * 24.0d0)
    minute = INT((dayfrac*24.0d0 - hour)*60.0d0)
    f      = NINT( ((dayfrac*24.0d0 - hour)*60.0d0 - minute)*60.0d0)
    second = INT(f)
    f      = f-second
    IF( f > 0.5 ) second = second + 1

    IF( second == 60 )THEN
       second = 0
       minute = minute + 1
    ENDIF

    IF( minute == 60 )THEN
       minute = 0
       hour   = hour + 1
    ENDIF

    IF( hour == 24 )THEN

       hour = 0
       !
       ! this could cause a bug, but probably 
       ! will never happen in practice
       !
       day = day + 1

    ENDIF

    IF( year < 0 )THEN
       year = year * (-1)
    ENDIF

    dt%year   = year
    dt%month  = month
    dt%day    = day
    dt%hour   = hour
    dt%minute = minute
    dt%second = second

  END SUBROUTINE jul2cal
  !
  !EOC
  !
  !-----------------------------------------------------------------------------!


  FUNCTION caldat(julian) result(dt)
    IMPLICIT NONE
    REAL(kind=Double), INTENT(IN) :: julian
    TYPE(dateTime) :: dt
    INTEGER(kind=Long), PARAMETER :: IGREG = 2299161
    INTEGER(kind=Long) :: mm, id, iyyy, ja, jd, je
    INTEGER(kind=Long) :: jalpha, jd_alpha
    INTEGER(kind=Long) :: offset  

   
    ! Cálculo da data
    if (julian >= IGREG) then
       jalpha = int(((julian - 1867216) - 0.25_Single) / 36524.25_Single)
       jd_alpha = julian + 1 + jalpha - int(0.25_Single * jalpha)
    else
       jd_alpha = julian
    end if
  
    ja = jd_alpha + 1524
    jd = 365 * jd_alpha + int(0.25_Single * jd_alpha)
    je = int((ja - jd) / 30.6001_Single)
    id = ja - jd - int(30.6001_Single * je)
    mm = je - 1
    if (mm > 12) mm = mm - 12
    iyyy = jd_alpha - 4715
    if (mm > 2) iyyy = iyyy - 1
    if (iyyy <= 0) iyyy = iyyy - 1
  
    ! Atribui os valores calculados à estrutura dateTime
    dt%year = iyyy
    dt%month = mm
    dt%day = id
    dt%timeZone%offset = offset
  
  END FUNCTION caldat

  
  pure function dateToStr(date) result(str)
    type(dateTime), intent(in) :: date
    character(len=19) :: str
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

  pure function strToDate(str) result(date)
    character(len=*), intent(in) :: str
    type(dateTime) :: date

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

  function strftime(date, format) result(str)
    class(dateTime), intent(in) :: date
    character(len=*), intent(in) :: format
    character(len=:), allocatable :: str
    integer :: i, j
    character(len=3), parameter :: monthNames(12) = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', &
                                                     'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    character(len=3), parameter :: weekdayNames(7) = ['Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat']
    character(len=1), parameter :: weekdayDigits(7) = ['0', '1', '2', '3', '4', '5', '6']

    ! Allocate space for str
    !allocate(character(len=len(format)) :: str)
  
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
    str = replace(str, '%A', adjustl( weekdayNames(mod(date%weekday() + 6, 7) + 1) ))
  
    str = replace(str, '%w', weekdayDigits(mod(date%weekday() + 6, 7) + 1))
  
    str = replace(str, '%U', num2str((date%yearDay() - mod(date%weekday()+6, 7) + 7) / 7, '(I2.2)'))
    str = replace(str, '%W', num2str((date%yearDay() - mod(date%weekday()+5, 7) + 7) / 7, '(I2.2)'))

!    str = replace(str, '%c', trim(strftime(date, '%a %b %d %H:%M:%S %Y')))
!    str = replace(str, '%x', trim(strftime(date, '%m/%d/%y')))
!    str = replace(str, '%X', trim(strftime(date, '%H:%M:%S')))
  
  end function strftime

    recursive function strptime(datetime_string, format_string)result(date)
        USE nrtype
        USE nrutil
        USE nr
        use m_string
        
        character(len=*), intent(in) :: datetime_string, format_string
        type(dateTime) :: date

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
           pos1(i) = index(format_string,f(i)%fmt)
        enddo

        
        !
        ! sort
        allocate(Fv(fSize))
        allocate(p(fSize))
        p = pos1
        Fv = float([(i,i=1,fSize)])        
        Call SORT2(p,Fv)
        pos1 = int(p)

        ! ajuste da posição inicial
        adj=0
        do i=2,fSize
           if (pos1(i)>0)then
              adj(i) = adj(i-1) + f(int(Fv(i-1)))%tam-2
              pos1(i)= pos1(i) + adj(i) 
           else
              adj(i) = 1
           endif
        enddo


        ! reord
        p = pos1
        Call SORT2(Fv,p)
        pos1=int(p)

        do i=1,fSize
           pos2(i) = pos1(i)+f(i)%tam-1
        enddo

        !do i=1,fSize
        !   print*,trim(f(i)%fmt),pos1(i),pos2(i)
        !enddo


        do i=1,fSize
           if(pos1(i)>0)then
              a=pos1(i)
              b=pos2(i)
              select case( f(i)%fmt )
                 case('%a')
                    day_of_week = maxloc(index(weekdayNames, datetime_string(a:b)),1)-1
                 case('%A')
                    day_of_week = maxloc(index(weekdayNames, datetime_string(a:b)),1)-1
                 case('%w')
                    day_of_week = maxloc(index(weekdayNames, datetime_string(a:b)),1)-1
                 case('%d')
                    day = str2int(datetime_string(a:b))
                 !case('%b')
                 !case('%B')
                 case('%m')
                    month = str2int(datetime_string(a:b))
                 case('%y')
                    year = str2int(datetime_string(a:b))
                 case('%Y')
                    year = str2int(datetime_string(a:b))
                 case('%H')
                    hour = str2int(datetime_string(a:b))
                 !case('%I')
                 !case('%p')
                 case('%M')
                    minute = str2int(datetime_string(a:b))
                 case('%S')
                    second = str2int(datetime_string(a:b))
                 !case('%f')
                 !case('%z')
                 case('%Z')
                    tZ%name=datetime_string(a:b)
                 case('%j')
                    day_of_year = str2int(datetime_string(a:b))
                 !case('%U')
                 !case('%W')
                 case('%c')
                    dt = strptime(datetime_string(a:b),'%d/%m/%y %H:%M:%S')
                    year   = dt%year  
                    month  = dt%month 
                    day    = dt%day   
                    hour   = dt%hour  
                    minute = dt%minute
                    second = dt%second
                 case('%x')
                    date = strptime(datetime_string(a:b),'%d/%m/%y')
                    year   = dt%year  
                    month  = dt%month 
                    day    = dt%day 
                 case('%X')
                    date = strptime(datetime_string(a:b),'%H:%M:%S')
                    hour   = dt%hour  
                    minute = dt%minute
                    second = dt%second

                 !case('%%')
                 case DEFAULT
                    write(*,*)'format not implemented yeat: ', trim(f(i)%fmt)
              end select
           endif
        !print*,f(i)%fmt,datetime_string(pos1(i):pos2(i))
        enddo

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


  PURE FUNCTION createDateTime(year, month, day, hour, minute, second, tZ) RESULT(dt)
    INTEGER( kind = Long ), INTENT(IN) :: year, month, day, hour, minute, second
    TYPE(timeZone), optional, INTENT(IN) :: tZ
    TYPE(dateTime) :: dt
  
    dt%year   = year
    dt%month  = month
    dt%day    = day
    dt%hour   = hour
    dt%minute = minute
    dt%second = second
    dt%day_of_week = weekday(dt)
    dt%day_of_year = yearday(dt)
    if(present(tZ)) dt%timeZone = tZ
  
  END FUNCTION createDateTime
  
  
  FUNCTION createTimeZone(hours, minutes, name) RESULT(tz)
    INTEGER( kind = Long ), INTENT(IN) :: hours, minutes ! Deslocamento em relação ao UTC em horas e minutos
    CHARACTER(*), INTENT(IN), OPTIONAL :: name           ! Nome do fuso horário (opcional)
    TYPE(timeZone) :: tz
    integer :: offSet

    ! Atribuir o deslocamento
    tz%offset = hours * 60 + minutes

    ! Verificar se o nome foi fornecido e atribuí-lo
    ! a partir do timezone fornecido calcular o deslocamento 
    ! em relação ao UTC em minutos
    IF (PRESENT(name)) THEN
      offSet = calculateOffset(name)
      if (offSet >= 0) then
         tz%name   = 'UTC'
         tz%offSet = offSet + hours * 60 + minutes
      else
         tz%name = name
      endif
    END IF
  END FUNCTION createTimeZone

  ! Função auxiliar para calcular o deslocamento em relação ao UTC em minutos
  INTEGER(kind=Long) FUNCTION calculateOffset(timezone_name)
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: timezone_name
  
    ! Conversão para letras maiúsculas para evitar problemas de comparação
    CHARACTER(len=:), ALLOCATABLE :: upper_timezone
    INTEGER(kind=Long) :: offset
  
    ! Converte para letras maiúsculas
    upper_timezone = TRANSFER(uppercase(timezone_name), upper_timezone)
  
    ! Define os deslocamentos para cada fuso horário
    SELECT CASE (TRIM(upper_timezone))
      CASE ("UTC")
        calculateOffset = 0
      CASE ("EST")
        calculateOffset = -5 * 60
      CASE ("PST")
        calculateOffset = -8 * 60
      CASE ("CST")
        calculateOffset = -6 * 60
      CASE ("IST")
        calculateOffset = 5 * 60 + 30
      CASE ("JST")
        calculateOffset = 9 * 60
      CASE ("BRT")
        calculateOffset = -3 * 60
      ! Adicione outros casos de fuso horário conforme necessário
      CASE DEFAULT
        write(*,*)'WARNNING : UNKNOWN timeZone:', trim(timezone_name)
        calculateOffset = -1
    END SELECT
  END FUNCTION calculateOffset

function count_digits(number) result(num_digits)
  implicit none
  integer, intent(in) :: number
  integer :: num_digits

  integer :: temp_num
  temp_num = abs(number)
  num_digits = 1

  do while (temp_num >= 10)
    temp_num = temp_num / 10
    num_digits = num_digits + 1
  end do

end function count_digits

function read_integer(string) result(number)
  implicit none
  character(len=*), intent(in) :: string
  integer :: number, ierr

  read(string, *, iostat=ierr) number

  if (ierr /= 0) then
    ! Erro ao ler o número inteiro
    number = -1 ! Ou outro valor de erro desejado
  end if

end function read_integer


function get_timezone(tz_str) result(tz)
  implicit none
  character(len=*), intent(in) :: tz_str
  type(timeZone) :: tz
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

  function str2int(str) result(num)
     character(len=*), intent(in) :: str
     integer :: num

     read(str,*)num

  end function

  function getName(tZ) result(name)
     class(timeZone) :: tZ
     character(len=10) :: name

     name = tZ%name
     return
  end function getName

  function getOffset(tZ) result(offSet)
     class(timeZone) :: tZ
     integer :: offSet

     offSet = tZ%offSet
     return
  end function getOffset

!!-----------------------------------------------------------------------------!
!!
!!BOP
!!
!! !IROUTINE: applyTimeZone
!!
!! !DESCRIPTION: This function adjusts the date and time based on the input and output time zones.
!!
!!
!! !INTERFACE:
!!
!PURE FUNCTION applyTimeZone(inputDateTime, inputTimeZone, outputTimeZone) RESULT(adjustedDateTime)
!  !
!  ! !INPUT PARAMETERS:
!  !
!  class(dateTime), intent(in) :: inputDateTime
!  class(timeZone), intent(in) :: inputTimeZone
!  class(timeZone), intent(in) :: outputTimeZone
!  !
!  ! !OUTPUT PARAMETERS:
!  !
!  class(dateTime) :: adjustedDateTime
!  !
!  !
!  ! !REVISION HISTORY: 
!  !  <Date> - <Author> - <Description of the modification>
!  !
!  !
!  ! !NOTE:
!  !       The function calculates the Julian Day for the input and output dates
!  !       using the cal2jul method defined in the dateTime type. It then calculates
!  !       the difference in time zone offsets between the input and output time zones.
!  !       The total number of seconds from the input date and time is adjusted based
!  !       on the time zone offset difference. Finally, the adjusted total number of
!  !       seconds is converted back to date and time using the convert_seconds_to_date
!  !       subroutine. The adjusted date and time are returned as the result.
!  !
!  !EOP
!  !-----------------------------------------------------------------------------!
!
!  ! Declare input and output variables
!  TYPE(dateTime), intent(in) :: inputDateTime
!  TYPE(timeZone), intent(in) :: inputTimeZone
!  TYPE(timeZone), intent(in) :: outputTimeZone
!  TYPE(dateTime) :: adjustedDateTime
!
!  ! Declare local variables
!  INTEGER( kind = Long ) :: inputJD, outputJD, offsetDiff, totalSeconds
!
!  ! Calculate Julian Day for input and output date
!  outputJD = adjustedDateTime%jdn()
!
!  ! Calculate the difference in time zone offsets
!  offsetDiff = outputTimeZone%offset - inputTimeZone%offset
!
!  ! Calculate the total number of seconds from the input date and time
!  totalSeconds = calculate_seconds_from_jdn(inputDateTime)
!
!  ! Adjust the total number of seconds based on the time zone offset difference
!  totalSeconds = totalSeconds + offsetDiff * 60
!
!  ! Convert the adjusted total number of seconds back to date and time
!  call convert_seconds_to_date(outputJD, totalSeconds, adjustedDateTime)
!
!  ! Return the adjusted date and time
!  applyTimeZone = adjustedDateTime
!end function applyTimeZone
!
!function calculate_seconds_from_jdn(date) result(total_seconds)
!    class(dateTime), intent(in) :: date
!    INTEGER, intent(in) :: julianDay, hour, minute, second
!    INTEGER :: total_seconds
!
!    ! Calculate Julian Day 
!    julianDay = date%jdn()
!
!    ! Calculate the total number of seconds from the input Julian Day and time
!    total_seconds = (julianDay - 2440587) * 24 * 60 * 60 + &
!                    hour * 60 * 60 + &
!                    minute * 60 + &
!                    second
!end function calculate_seconds_from_jdn
!
!subroutine convert_seconds_to_date(julianDay, total_seconds, dateTime)
!    INTEGER, intent(in) :: julianDay, total_seconds
!    TYPE(DateTime), intent(out) :: dateTime
!
!    ! Convert the total number of seconds to date and time
!    dateTime%year = julianDay_to_year(julianDay)
!    dateTime%month = julianDay_to_month(julianDay)
!    dateTime%day = julianDay_to_day(julianDay)
!
!    dateTime%hour = total_seconds / (60 * 60)
!    total_seconds = total_seconds - dateTime%hour * 60 * 60
!
!    dateTime%minute = total_seconds / 60
!    dateTime%second = total_seconds - dateTime%minute * 60
!end subroutine convert_seconds_to_date


end module dateTimeMod
