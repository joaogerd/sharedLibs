module dateTimeMod
  use typeKinds
  implicit none
  private

  public :: dateTime
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

  type timeDelta
     integer( kind = Long ) :: years   = 0
     integer( kind = Long ) :: months  = 0
     integer( kind = Long ) :: days    = 0
     integer( kind = Long ) :: hours   = 0
     integer( kind = Long ) :: minutes = 0
     integer( kind = Long ) :: seconds = 0
  end type timeDelta

  type dateTime
     integer( kind = Long ) :: year     = 1900
     integer( kind = Long ) :: month    = 1
     integer( kind = Long ) :: day      = 1
     integer( kind = Long ) :: hour     = 0
     integer( kind = Long ) :: minute   = 0
     integer( kind = Long ) :: second   = 0
     integer( kind = Long ) :: timeZone = 0
   contains
     procedure, public  :: jdn => cal2jul
     procedure, public  :: jdn2 => julday
     procedure, public  :: jdn3 => julday2
     procedure, private :: timeAdd
     generic :: operator(+) => timeAdd
     procedure, private :: timeSub
     generic :: operator(-) => timeSub
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
    today%timeZone = date(4)/60.0

    return
  end function today

  pure elemental function isLeapYear(year) result(flag)
    integer, intent(in) :: year
    logical :: flag

    flag = .false.
    if (mod(year,  4) == 0 .and. &
         mod(year,100) /= 0 .or.  &
         mod(year,400) == 0) flag = .true.

  end function isLeapYear

  pure elemental function daysInYear(year) result (days)
    integer, intent(in) :: year
    integer :: days

    if (isLeapYear(year)) then
        days = 366
    else
        days = 365
    endif

  end function daysInYear

  pure elemental function daysInMonth(year, month) result (days)
    integer, intent(in) :: year, month
    integer :: days

    select case (month)
      case (1, 3, 5, 7, 8, 10, 12)
         days = 31
      case (4, 6, 9, 11)
         days = 30
      case (2)
         if (isLeapYear(year)) then
            days = 29
         else
            days = 28
         endif
      case default
         days = -1
    end select

  end function daysInMonth

  pure elemental function deltaTime(years, months, days, hours, minutes, seconds) result(delta)
      integer, intent(in), optional :: years, months, days, hours, minutes, seconds
      type(timeDelta) :: delta
  
      delta%years   = merge(0, years,   present(years))
      delta%months  = merge(0, months,  present(months))
      delta%days    = merge(0, days,    present(days))
      delta%hours   = merge(0, hours,   present(hours))
      delta%minutes = merge(0, minutes, present(minutes))
      delta%seconds = merge(0, seconds, present(seconds))
  
  end function deltaTime


  pure function numberOfDays(start, stop) result(nDays)
    class(dateTime), intent(in) :: start, stop
    integer( kind = Long ) :: nDays

    integer( kind = Long ) :: jdStart, jdStop

    jdStart = start%jdn()
    jdStop  = stop%jdn()

    nDays = jdStop - jdStart + 1

  end function numberOfDays

  pure function timeAdd(a, b) result(resultTime)
    class(dateTime), intent(in) :: a, b
    type(dateTime) :: resultTime

    integer( kind = Long ) :: totalSeconds
    integer( kind = Long ) :: days, hours, minutes, seconds

    totalSeconds = a%second + b%second + &
                   a%minute * 60 + b%minute * 60 + &
                   a%hour * 3600 + b%hour * 3600 + &
                   a%day * 86400 + b%day * 86400

    days    = totalSeconds / 86400
    hours   = mod(totalSeconds, 86400) / 3600
    minutes = mod(totalSeconds, 3600) / 60
    seconds = mod(totalSeconds, 60)

    resultTime%year     = a%year + b%year
    resultTime%month    = a%month + b%month
    resultTime%day      = a%day + b%day + days
    resultTime%hour     = hours
    resultTime%minute   = minutes
    resultTime%second   = seconds
    resultTime%timeZone = a%timeZone

  end function timeAdd

  pure function timeSub(a, b) result(resultTime)
    class(dateTime), intent(in) :: a, b
    type(dateTime) :: resultTime

    integer( kind = Long ) :: totalSeconds
    integer( kind = Long ) :: days, hours, minutes, seconds

    totalSeconds = a%second - b%second + &
                   a%minute * 60 - b%minute * 60 + &
                   a%hour * 3600 - b%hour * 3600 + &
                   a%day * 86400 - b%day * 86400

    days    = totalSeconds / 86400
    hours   = mod(totalSeconds, 86400) / 3600
    minutes = mod(totalSeconds, 3600) / 60
    seconds = mod(totalSeconds, 60)

    resultTime%year     = a%year - b%year
    resultTime%month    = a%month - b%month
    resultTime%day      = a%day - b%day + days
    resultTime%hour     = hours
    resultTime%minute   = minutes
    resultTime%second   = seconds
    resultTime%timeZone = a%timeZone

  end function timeSub

!  function dateDiff(d1,d2) result(deltaT)
!    class(dateTime), intent(in) :: d1
!    class(dateTime), intent(in) :: d2
!    type(timeDelta) :: deltaT
!
!    real(kind = Double) :: days
!    real(kind = Double) :: incr
!
!    days = cal2jul(d1) - cal2jul(d2)
!    
!    deltaT = deltaTime(days)
!
!    return
!  end function dateDiff



!  pure function cal2jul(date) result(jd)
!    type(dateTime), intent(in) :: date
!    integer( kind = Long ) :: jd
!
!    jd = julday(date%year, date%month, date%day) - date%timeZone/24.0
!
!  end function cal2jul
!
!  pure function julday(year, month, day) result(jd)
!    integer( kind = Long ), intent(in) :: year, month, day
!    integer( kind = Long ) :: jd
!
!    jd = (1461 * (year + 4800 + (month - 14) / 12)) / 4 &
!       + (367 * (month - 2 - 12 * ((month - 14) / 12))) / 12 &
!       - (3 * ((year + 4900 + (month - 14) / 12) / 100)) / 4 &
!       + day - 32075
!
!  end function julday
!
!  pure function julday2(year, month, day) result(jd)
!    integer( kind = Long ), intent(in) :: year, month, day
!    integer( kind = Long ) :: jd
!
!    jd = (1461 * (year + 4800 + (month - 14) / 12)) / 4 &
!       + (367 * (month - 2 - 12 * ((month - 14) / 12))) / 12 &
!       - (3 * ((year + 4900 + (month - 14) / 12) / 100)) / 4 &
!       + day - 32075 - 0.5
!
!  end function julday2


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

  pure function gregorian(jd) result(dt)
    integer( kind = Long ), intent(in) :: jd
    type(dateTime) :: dt

    integer( kind = Long ) :: l, m, n, p, q

    l = jd + 68569
    n = (4 * l) / 146097
    l = l - (146097 * n + 3) / 4
    m = (4000 * (l + 1)) / 1461001
    l = l - (1461 * m) / 4 + 31
    p = (80 * l) / 2447
    dt%day = l - (2447 * p) / 80
    l = p / 11
    dt%month = p + 2 - (12 * l)
    dt%year = 100 * (n - 49) + m + l - 1900

  end function gregorian

  FUNCTION caldat(julian) RESULT(dt)
      REAL(kind = Double), INTENT(IN) :: julian
      TYPE(dateTime) :: dt
      INTEGER(kind = Long) :: mm, id, iyyy
      INTEGER(kind = Long) :: ja, jalpha, jb, jc, jd, je
      INTEGER(kind = Long), PARAMETER :: IGREG = 2299161
      REAL(kind = Double) :: time
  
      ! Calcula a parte inteira da data juliana
      ja = INT(julian)
  
      ! Calcula a parte fracionária que representa o tempo do dia
      time = julian - REAL(ja, kind=Double)
  
      ! Calcula a data
      if (ja >= IGREG) then
          jalpha = INT(((ja - 1867216) - 0.25_Double) / 36524.25_Double)
          ja = ja + 1 + jalpha - INT(0.25_Double * jalpha)
      end if
      jb = ja + 1524
      jc = INT(6680.0_Double + ((jb - 2439870) - 122.1_Double) / 365.25_Double)
      jd = 365 * jc + INT(0.25_Double * jc)
      je = INT((jb - jd) / 30.6001_Double)
      id = jb - jd - INT(30.6001_Double * je)
      mm = je - 1
      if (mm > 12) mm = mm - 12
      iyyy = jc - 4715
      if (mm > 2) iyyy = iyyy - 1
      if (iyyy <= 0) iyyy = iyyy - 1
  
      ! Atribui os valores à estrutura dateTime
      dt%year = iyyy
      dt%month = mm
      dt%day = id
  
      ! Extrai a hora, minutos e segundos da parte fracionária do tempo
      dt%hour = INT(time * 24.0_Double)
      dt%minute = INT((time * 24.0_Double - REAL(dt%hour, kind=Double)) * 60.0_Double)
      dt%second = INT((((time * 24.0_Double - REAL(dt%hour, kind=Double)) * 60.0_Double) &
                  - REAL(dt%minute, kind=Double)) * 60.0_Double)
  
  END FUNCTION caldat


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
  FUNCTION cal2jul(dt) RESULT(jd)
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

    ymd = (year*10000)+(month*100)+day
    if(ymd>=15821015)then
       A = INT(year/100)
       B = INT(A/4)
       C = 2 - A + B
    endif

    if(ymd<=15821004)then
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
  !
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
  SUBROUTINE jul2cal(jd, dt)
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
    date%timeZone = 0

  end function strToDate

  pure function strftime(date, format) result(str)
    type(dateTime), intent(in) :: date
    character(len=*), intent(in) :: format
    character(len=:), allocatable :: str
  
    character(len=2) :: yearStr, monthStr, dayStr, hourStr, minuteStr, secondStr
  
    ! Convert each component to string
    write(yearStr, '(I4)') date%year
    write(monthStr, '(I2)') date%month
    write(dayStr, '(I2)') date%day
    write(hourStr, '(I2)') date%hour
    write(minuteStr, '(I2)') date%minute
    write(secondStr, '(I2)') date%second
  
    ! Replace format placeholders with corresponding values
    str = format
    str = replace(str, '%Y', trim(yearStr))
    str = replace(str, '%m', trim(monthStr))
    str = replace(str, '%d', trim(dayStr))
    str = replace(str, '%H', trim(hourStr))
    str = replace(str, '%M', trim(minuteStr))
    str = replace(str, '%S', trim(secondStr))
  
  end function strftime

  pure function strptime(dateStr, format) result(date)
    character(len=*), intent(in) :: dateStr
    character(len=*), intent(in) :: format
    type(dateTime) :: date
  
    character(len=:), allocatable :: yearStr, monthStr, dayStr, hourStr, minuteStr, secondStr
    integer :: year, month, day, hour, minute, second
  
    ! Find the position of format placeholders in the given format
    integer :: yearPos, monthPos, dayPos, hourPos, minutePos, secondPos
    yearPos = index(format, '%Y')
    monthPos = index(format, '%m')
    dayPos = index(format, '%d')
    hourPos = index(format, '%H')
    minutePos = index(format, '%M')
    secondPos = index(format, '%S')
  
    ! Extract the corresponding substrings from the date string
    if (yearPos > 0) then
      yearStr = dateStr(yearPos:yearPos+3)
      read(yearStr, '(I4)') year
    end if
    if (monthPos > 0) then
      monthStr = dateStr(monthPos:monthPos+1)
      read(monthStr, '(I2)') month
    end if
    if (dayPos > 0) then
      dayStr = dateStr(dayPos:dayPos+1)
      read(dayStr, '(I2)') day
    end if
    if (hourPos > 0) then
      hourStr = dateStr(hourPos:hourPos+1)
      read(hourStr, '(I2)') hour
    end if
    if (minutePos > 0) then
      minuteStr = dateStr(minutePos:minutePos+1)
      read(minuteStr, '(I2)') minute
    end if
    if (secondPos > 0) then
      secondStr = dateStr(secondPos:secondPos+1)
      read(secondStr, '(I2)') second
    end if
  
    ! Assign the extracted values to the date object
    date%year = year
    date%month = month
    date%day = day
    date%hour = hour
    date%minute = minute
    date%second = second
  
  end function strptime
  
  function julday(d) result(jdn)
    class(dateTime), intent(in) :: d
    real(kind = Double) :: jdn

    real :: year
    real :: month

    year  = d%year
    month = d%month
    if (month < 3)then
       month = month + 12
       year  = year - 1
    endif

    jdn = d%day + int((153.0*month-457.0)/5.0) + &
         365*year + floor(year/4.0) - floor(year/100.0) + &
         floor(year/400.0) + 1721118.5

  end function julday
  
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

  end module dateTimeMod

