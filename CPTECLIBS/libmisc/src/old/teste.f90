program datetime_parser
    implicit none
    
    character(len=100) :: datetime_string
    character(len=100) :: format_string
    integer :: year, month, day, hour, minute, second
    integer :: ierr
    
    ! Entrada da string de data e hora
    datetime_string = "2023-  05-   21 10:30:00"
    format_string = "%Y-  %m-   %d %H:%M:%S"
    
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


    ! Entrada da string de data e hora
    datetime_string = "05012023"
    format_string = "%m%d%Y"
    
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
        character(len=2), parameter    :: fmt(*) = ['%Y','%m','%d','%H','%M','%S']
        integer, parameter             :: tam(*) = [   4,   2,   2,   2,   2,   2]
        integer,allocatable :: pos1(:), pos2(:)
        REAL, ALLOCATABLE, DIMENSION(:)     :: Fv, p
        integer :: t1,t2,t3
        print*,format_string
        print*,datetime_string

        allocate(pos1(size(tam)))
        allocate(pos2(size(tam)))
        do i = 1, size(fmt)
           pos1(i) = index(format_string,fmt(i))
        enddo
        print*,pos1
        
        !
        ! sort
        !

        allocate(Fv(size(tam)))
        allocate(p(size(tam)))
        p = pos1
        Fv = float([(i,i=1,size(tam))])        
        Call SORT2(p,Fv)
        pos1 = int(p)

        print*,pos1
        print*,Fv
        
        ! ajuste da posição inicial
        t3=0
        do i=2,size(tam)
            t1 = pos1(i) - (pos1(i-1) + 1)
            t2 = pos1(i) - (pos1(i-1) + tam(int(Fv(i-1)))-1)
            t3 = t3 + (t1-t2)
            pos1(i) = pos1(i) + t3
        enddo
 


        do i=1,size(tam)
           pos2(i) = pos1(i)+tam(int(Fv(i)))-1
        enddo

        ! reord

        pos1 = pos1(int(Fv))
        pos2 = pos2(int(Fv))


        print*,pos1
        print*,pos2
        print*,'--------------------------'
        do i=1,size(tam)
        print*,fmt(i),datetime_string(pos1(i):pos2(i))
        enddo
        
    end subroutine parse_datetime
    
    integer function parse_value(datetime_string, format_string, ierr)
        implicit none
        character(len=*), intent(in) :: datetime_string, format_string
        integer, intent(out) :: ierr
        character(len=:), allocatable :: value_string
        integer :: pos1, pos2
        character(len=:), allocatable :: temp_string
        integer :: value
        integer :: i
        
        ierr = 0
        
        ! Encontrar a posição do valor no formato
        pos1 = index(format_string, "%")
        pos2 = index(format_string(pos1+1:), "%")
        if (pos2 == 0) then
            ierr = 1
            return
        endif
        pos2 = pos2 + pos1
        
        ! Extrair a string de valor correspondente
        value_string = datetime_string(pos1:pos2)
        
        ! Remover os caracteres de formatação
        temp_string = ""
        do i = 1, len_trim(value_string)
            if (value_string(i:i) /= "%") then
                temp_string = trim(temp_string) // value_string(i:i)
            endif
        end do
        
        ! Converter a string para um valor numérico
        read(temp_string, *) value
        
        parse_value = value
        
    end function parse_value
    
end program datetime_parser

