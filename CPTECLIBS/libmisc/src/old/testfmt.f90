program testfmt
use strtime
implicit none
character(len=50) :: timeFMT,var
integer :: ano, mes, dia
timeFMT='%d2 de %Mc de %y4 as %h2'

print*,trim(timeFMT)
call replace(timeFMT)
print*,trim(timeFMT)


ano=1993
mes=01
dia=31
write(*,'(I4.4,I2.2,I2.2)')ano,mes,dia

var='2000 02 13'
read(var,'(I4.4,1x,I2.2,1x,I2.2)')ano,mes,dia
print*,ano,mes,dia


end program
