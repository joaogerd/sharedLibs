program dateTime_example

  use dateTimeMod

  implicit none

  type(dateTime) :: dt1, dt2
  type(dateTime) :: startDateTime, currentDateTime, endDateTime
  type(timeDelta) :: timeIncrement
  type(timeZone) :: tz
  character(len=1000) :: tmp

  logical :: result

  ! Exemplo de criação de uma data e hora
  dt1 = createDateTime(2016, 2, 29, 23, 59, 0, tz)

  ! Exemplo de acesso aos campos da data e hora
  print *, "Ano: ", dt1%year
  print *, "Mês: ", dt1%month
  print *, "Dia: ", dt1%day
  print *, "Hora: ", dt1%hour
  print *, "Minuto: ", dt1%minute
  print *, "Segundo: ", dt1%second
  print *, "Fuso Horário: ", dt1%timeZone

  ! Exemplo de formatação da data e hora
!  print*,'aqio'

  print *, "Dt1 -> Data e Hora formatadas: ", dt1%strftime("%d/%m/%Y %H:%M:%S")

  print*,dt1%jdn()
  print*,'---->',ConvertToDecimal(dt1)
  dt2 = ConvertFromDecimal(ConvertToDecimal(dt1))
  print *, "Dt2 -> Data e Hora formatadas: ", dt2%strftime("%d/%m/%Y %H:%M:%S")

 ! print*, calculate_seconds_from_jdn(dt2)
  
!!  print*,'passou'
!  ! Exemplo de conversão de uma string para data e hora
!  dt2 = strptime("2023-05-21 10:30:00", "%Y-%m-%d %H:%M:%S")
!  print *, "Dt2 -> Data e Hora formatadas: ", dt2%strftime("%d/%m/%Y %H:%M:%S")
!  dt2 = strptime("21/05/23 10:30:00", "%c")
!  print *, "Dt2 -> Data e Hora formatadas: ", dt2%strftime("%d/%m/%y %H:%M:%S")
! dt2 = strptime("2023-05-21 10:30:00", "%Y-%m-%d %X")
!  print *, "Dt2 -> Data e Hora formatadas: ", dt2%strftime("%d/%m/%Y %H:%M:%S")
!  print *, "Ano: ", dt2%year
!  print *, "Mês: ", dt2%month
!  print *, "Dia: ", dt2%day
!  print *, "Hora: ", dt2%hour
!  print *, "Minuto: ", dt2%minute
!  print *, "Segundo: ", dt2%second
!  print *, "Fuso Horário: ", dt2%timeZone
!
!  ! Exemplo de comparação de datas e horas
!  result = (dt1.gt.dt2)
!  if (result) then
!    print *, "dt1 é maior que dt2"
!  else
!    print *, "dt1 é menor ou igual a dt2"
!  end if
!
  ! Definindo a data e hora inicial e final
  startDateTime = createDateTime(2023, 5, 21, 0, 0, 0)
  endDateTime = createDateTime(2023, 5, 30, 0, 0, 0)

  ! Definindo o incremento de tempo (6 horas)
  timeIncrement = deltaTime(0, 0, 0, 6, 0, 0)

  ! Loop no tempo

!  currentDateTime = startDateTime
!  do while (currentDateTime <= endDateTime )
!    print *, "Data e Hora: ", trim(currentDateTime%strftime("%d/%m/%Y -  %H:%M:%S"))
!    currentDateTime = currentDateTime + timeIncrement
!  end do
!
!   ! Definindo a data e hora atual
!  currentDateTime = createDateTime(2023, 5, 21, 12, 0, 0)
!
!  ! Definindo o timezone
!  tZ = createTimeZone(3, 0) ! UTC+3:00
!
!  ! Aplicando o timezone à data e hora atual
!!  currentDateTime = applyTimeZone(currentDateTime, tZ)
!!
!!  ! Obtendo a informação do timezone
!  print *, "Timezone Offset: ", tz%getOffset()
!  print *, "Timezone Name: ", tz%getName()
!!
!!  ! Exibindo a data e hora atual com o timezone aplicado
!!  print *, "Data e Hora: ", currentDateTime%strftime("%d/%m/%Y %H:%M:%S")
  
end program dateTime_example
