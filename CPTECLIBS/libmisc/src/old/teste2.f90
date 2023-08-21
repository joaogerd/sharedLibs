program datetime_parser
    implicit none
    
    character(len=100) :: datetime_string
    character(len=100) :: format_string
    integer :: year, month, day, hour, minute, second
    integer :: ierr
    
    ! Entrada da string de data e hora
    datetime_string = "23-05-2021 10:30:00"
    format_string = "%d-%m-%Y %H:%M:%S"
    
    ! Chamada para a função de parser
    call parse_datetime(datetime_string, format_string, year, month, day, hour, minute, second, ierr)
    
    ! Verificação de erros
    if (ierr /= 0) then
        write(*,*) "Erro durante o parsing da data e hora.",ierr
        stop
    endif
    
    ! Exibição dos valores obtidos
!    write(*,*) "Ano:   ", year
!    write(*,*) "Mês:   ", month
!    write(*,*) "Dia:   ", day
!    write(*,*) "Hora:  ", hour
!    write(*,*) "Minuto:", minute
!    write(*,*) "Segundo:", second


   ! Entrada da string de data e hora
    datetime_string = "05/(01)--2023)"
    format_string = "%m/(%d)--%Y)"
   
   ! Chamada para a função de parser
    call parse_datetime(datetime_string, format_string, year, month, day, hour, minute, second, ierr)
   
   ! Verificação de erros
    if (ierr /= 0) then
        write(*,*) "Erro durante o parsing da data e hora.",ierr
        stop
    endif
   
!    ! Exibição dos valores obtidos
!    write(*,*) "Ano:   ", year
!    write(*,*) "Mês:   ", month
!    write(*,*) "Dia:   ", day
!    write(*,*) "Hora:  ", hour
!    write(*,*) "Minuto:", minute
!    write(*,*) "Segundo:", second

    
contains

    subroutine parse_datetime(datetime_string, format_string, year, month, day, hour, minute, second, ierr)
        USE nrtype
        USE nrutil
        USE nr
        implicit none
        character(len=*), intent(in) :: datetime_string, format_string
        integer, intent(out) :: year, month, day, hour, minute, second
        integer, intent(out) :: ierr
        
        character(len=:), allocatable :: date_string, time_string
        character(len=:), allocatable :: year_str, month_str, day_str, hour_str, minute_str, second_str
        character(len=:), allocatable :: temp_format_string
        integer :: i, j
        character(len=2), parameter    :: fmt(*) = ['%Y','%m','%d','%H','%M','%S','%j']
        integer, parameter             :: tam(*) = [   4,   2,   2,   2,   2,   2,   1]
        integer,allocatable :: pos1(:), pos2(:), adj(:)
        integer ::  fSize
        REAL, ALLOCATABLE, DIMENSION(:)     :: Fv, p


        type strf
           character(len=:), allocatable :: fmt
           integer                       :: tam
        end type
        type(strf), allocatable :: f(:)

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

        print*,'-----AJUSTE--------'
!        do i=1,fSize
!           print*,f(int(Fv(i)))%fmt,pos1(i),adj(i), f(int(Fv(i-1)))%tam
!        enddo
!        pos1 = pos1+adj 

        ! reord
        p = pos1
        Call SORT2(Fv,p)
        pos1=int(p)

        do i=1,fSize
           pos2(i) = pos1(i)+f(i)%tam-1
        enddo


        print*,'----------Final-----------'
        do i=1,fSize
        print*,pos1(i),pos2(i)
        enddo
        print*,'--------------------------'
        do i=1,fSize
        if(pos1(i)>0)        print*,f(i)%fmt,datetime_string(pos1(i):pos2(i))
        enddo
        
    end subroutine parse_datetime

!function removeCaracteresEspeciais(string) result(stringNumeros)
!  character(len=*), intent(in) :: string
!  character(len=100) :: stringNumeros
!  character(len=100) :: tempString
!  integer :: i, j
!
!  tempString = string
!  stringNumeros = ""
!
!  do i = 1, len_trim(tempString)
!    if (tempString(i:i) >= '0' .and. tempString(i:i) <= '9') then
!      stringNumeros = trim(stringNumeros) // tempString(i:i)
!    end if
!  end do
!
!end function removeCaracteresEspeciais

    function removeCaracteresEspeciais(inputString) result(cleanedString)
        character(len=*), intent(in) :: inputString
        character(len=:), allocatable :: cleanedString
        character(len=12) :: specialCharacters
        character(len=6) :: dateFormats(6)
        integer :: i, j, n
        
        specialCharacters = ' !@#$%^&*()-_+=\|]}[{;:''"/?<>,.`~'
        dateFormats = (/'%Y', '%m', '%d', '%H', '%M', '%S'/)
        
        cleanedString = ''
        
        do i = 1, len_trim(inputString)
            if (index(specialCharacters, inputString(i:i)) == 0) then
                cleanedString = cleanedString // inputString(i:i)
            end if
        end do
        
        do i = 1, size(dateFormats)
            do j = 1, len(dateFormats(i))
                n = index(cleanedString, dateFormats(i)(j:j))
                if (n /= 0) then
                    cleanedString = cleanedString(1:n-1) // ' ' // cleanedString(n+1:)
                end if
            end do
        end do
        
    end function

end program datetime_parser

