program main
    use datetime_module
    implicit none

    real :: atime = 2016.1639344262296
    type(datetime) :: result

    result = t2dt(atime)

    write(*, *) "Year: ", result%year
    write(*, *) "Month: ", result%month
    write(*, *) "Day: ", result%day
    write(*, *) "Hour: ", result%hour
    write(*, *) "Minute: ", result%minute
    write(*, *) "Second: ", result%second

end program main

