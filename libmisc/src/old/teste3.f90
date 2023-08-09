program teste3
  character(len=1), parameter :: a(*) = ['a','b','c']
  character(len=1), parameter :: b(*) = ['A','B','C']

  integer :: i

  do i=1,size(a)
     if ('a' == a(i) )then
        print*,'a exta em a'
     endif
     if ('A' == a(i) )then
        print*,'A exta em a'
     endif 
  enddo

  print*,index(a,'a')
  print*,index(a,'A')

end program
