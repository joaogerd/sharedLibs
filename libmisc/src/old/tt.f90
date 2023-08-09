program tt
  type teste
    character(len=:), allocatable :: fmt
    integer          :: tam
  end type teste

  type(teste), allocatable :: forma(:)

  ! Initialize the format array using derived type constructors
  forma = [teste(fmt=char(64 + i), tam=i) | i = 1, 3]

  ! Your code logic goes here

  ! Deallocate the allocatable components when no longer needed
  deallocate(format%fmt)
  deallocate(format)
end program tt
