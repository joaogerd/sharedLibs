module strtime

   type node
      character(len=10)   :: name
      integer             :: val
      type(node), pointer :: next => null()
   end type


contains

subroutine replace(sfmt)
   CHARACTER(len=*), INTENT(INOUT) :: sfmt

   call replace_( sfmt, '%y2', '(I2.2)') 
   call replace_( sfmt, '%y4', '(I4.4)') 
   call replace_( sfmt, '%m1', '(I1.1)') 
   call replace_( sfmt, '%m2', '(I2.2)')
   call replace_( sfmt, '%mc', '(A3)')
   call replace_( sfmt, '%Mc', '(A3)')
   call replace_( sfmt, '%MC', '(A3)')
   call replace_( sfmt, '%d1', '(I1.1)') 
   call replace_( sfmt, '%d2', '(I2.2)')
   call replace_( sfmt, '%h1', '(I1.1)') 
   call replace_( sfmt, '%h2', '(I2.2)')
   call replace_( sfmt, '%h3', '(I3.3)')
   call replace_( sfmt, '%n2', '(I2.2)')


end subroutine


SUBROUTINE replace_(strg,mask,repl)

    !
    !
    IMPLICIT NONE

    ! !INPUT/OUTPUT PARAMETERS:

    CHARACTER(len=*),INTENT(INOUT)  :: strg ! String

    ! !INPUT PARAMETERS:

    CHARACTER(len=*),INTENT(IN)     :: mask ! maskout
    CHARACTER(len=*),INTENT(IN)     :: repl ! replacing string
    !
    ! !REVISION HISTORY:
    !  Joao Gerd - 20Feb2011 - Codigo Inicial
    !
    !EOP
    !---------------------------------------------------------------------!
    !BOC
    CHARACTER(len=300) ::  sub, tmp
    INTEGER :: lenRepl, lenMask, lenStrg, lenTmp
    INTEGER :: i, j

    lenStrg = LEN_TRIM(strg)
    lenRepl = LEN(repl)
    lenMask = LEN(mask)
    lenTmp  = 0

    tmp = ''
    i   = 1
    j   = 1

    do while( j .LE. lenStrg )

       sub = strg(j:j+lenMask-1)

       if( sub .EQ. mask )then

          tmp(lenTmp+1:lenTmp+lenRepl) = repl

          lenTmp = lenTmp+lenRepl
          i      = j + lenMask
          j      = i

       endif

       tmp(lenTmp+1:lenTmp+1) = strg(j:j)
       
       lenTmp = lenTmp+1
       j      = j + 1 

    enddo

    strg = tmp

    return
  END SUBROUTINE replace_

end module strtime
