program testchar
   use m_string
   character(len=60) :: iname
   character(len=:),allocatable :: oname

   iname = "Hello World"
   oname = replace(iname,'%Y','2000')
   print*,iname
   print*,oname
end program
