module m_GrADSFiles
  use m_inpak90
  use m_stdio
  use typeKinds, only: i4 => Long, i8 => LLong
  use typeKinds, only: r4 => Single, r8 => Double
  use accessgrib, only: gribType

  implicit none
  private

  public :: GrADSFiles, GrADSVar

  !
  ! Name of this module
  !
  character(len=*), parameter :: myName = 'libGrADSFiles'
  !
  ! Undefine Value
  !
  real, parameter :: UDEF_ = -9.9e+20
  
  ! This is because grib had a problem when open
  ! many files. If we open, read and close a file,
  ! in some times grib routine retain information 
  ! about previous file and input some noise to data,
  ! or can't read corretly a file.
  ! using "exclude" we try save the number of 
  ! previous openned units, and skip it to open
  ! in a new number.
  integer, parameter :: MaxLogUnit = 254 
  integer, private   :: exclude(MaxLogUnit) = 0

  ! Size of character
  integer, parameter :: charMaxLen = 1024

  type GrADSdef
    character(len=4)  :: dName
    character(len=6)  :: dType
    integer           :: dSize
    real, pointer     :: dLevs(:)
    type(GrADSdef), pointer :: next => null()     
  end type
  
  type GrADSVar
     character(len=12)    :: name     ! Variable name
     integer              :: nlevs    ! number of levels by variable
     integer              :: pmark    ! Initial Position Marker (bytes)
     character(len=64)    :: descr    ! descrition of variable
     integer              :: units
     integer, allocatable :: param(:) ! grib parameter table
     type(GrADSVar), pointer :: next => null()
  end type GrADSVar

  type GrADSfiles

     !
     ! Common variables
     !

     character(len=1024)      :: dset    ! the filename for input
     real                     :: undef   ! missing value flag
     character(len=10)        :: dtype   ! type of file
     integer                  :: nvars   ! number of variables
     integer                  :: iosize
     integer                  :: lu    ! logical unit if already opened

     !
     ! Model definitions:
     !

     type(GrADSdef), pointer   :: dHead => null()
     type(GrADSdef), pointer   :: dTail => null()

     !
     ! Model variables:
     !

     type(GrADSVar), pointer   :: vHead => null() ! Variables
     type(GrADSVar), pointer   :: vTail => null()

     !
     !
     !

     contains
        procedure, public  :: open => open_
        procedure, private :: parseCtl_
        procedure, public  :: qctlinfo
        procedure, public  :: getDim
        procedure, public  :: getDimVec
        procedure, public  :: pVars
        procedure, public  :: getField1D
        procedure, private :: getVarPtr

!     procedure, private :: openAscii_
!     procedure, private :: reset

  end type GrADSfiles


  contains
  !----------------------------------------------------------------------!
  !                 INPE/CGCT, Data Assimilation Office                  !
  !----------------------------------------------------------------------!
  !BOP
  ! !IROUTINE: open_ - open an input GrADS "control" file for input
  !
  ! !DESCRIPTION: (to do)

  function open_(self, ctlFile)result(status)
     class(GrADSFiles), intent(inout) :: self    ! GrADSFile strutucture 
     character(len=*),  intent(in   ) :: ctlFile ! ctl file name
     integer                          :: status  ! status of open

     ! !EXAMPLES: (to do)
     ! !BUGS: (to do)
     ! !SEE ALSO: (to do)
     ! !SYSTEM ROUTINES: (to do)
     !
     ! !REVISION HISTORY:
     !    21Jan00    - Jing Guo
     !        . Added "direct" access to open_()
     !     16Jul96 - J. Guo    - modified as a Fortran 90 module.
     !    01Dec94 - Jing G.    - added zdef_gr for small values.
     !                - Original source from A. da Silva
     !    11May20 - J.G de Mattos - dapt to get all information
     !                              from vars ctl block
     !EOP
     !----------------------------------------------------------------------!
     !BOC

     character(len=*), parameter :: myname_ = myname//' :: open_'

    ! Local variables:

    character(len=64) :: str
    integer i, j, k, l, lu

    integer :: tdim, nvars
    integer :: ierr
    character(len=4) :: ext

    status = 0

    !--------------------------------------
    ! Verificando tipo de arquivo de entrada
    !--------------------------------------

    do i = len_trim(ctlFile), 1, -1
       if (ctlFile(i:i) .eq. '.') then
          j = i + 1
          exit
       endif
    end do

    ext = i90_lcase(ctlFile(j:len_trim(ctlFile)))

    select case (ext)

!       case ('txt', 'dat')
!          !
!          ! Reading ascii file, skip default read
!          !
!          lu = self%lu
!          if (lu >= 0) close (lu)
!   
!          self%lu    = i90_lua(exclude)
!          i          = minloc(exclude,1)
!          exclude(i) = self%lu
!          self%dtype = tstation
!          self%dset  = trim(ctlFile)
!   
!          call self%openAscii_(ctlFile)
!   
!          return
   
       case ('ctl')
          !
          ! reading from ctl file
          !
   
          status = self%parseCtl_(ctlFile)
          if (status /= 0)then
             call i90_perr(trim(myname_), &
                            ': Something went wrong when parsing the file: "'//trim(ctlFile)//'" , ierr =', status)
             call i90_die(myname_)
          else

          endif
       case ('netcdf')
          !
          ! reading from netCDF
          !
          
          call i90_perr(trim(myname_), &
               ': Wrong type file: "'//trim(ctlFile)//'" error, ierr =', 97)
          call i90_perr(trim(myname_),'netCDF not implemented yeat via ctl!')
          status = 97
          return

       case default
          call i90_perr(trim(myname_), &
               ': Wrong type file: "'//trim(ctlFile)//'" error, ierr =', 99)
          call i90_die(myname_)
          status = 99
          return

    end select

  end function
  !EOC
  !----------------------------------------------------------------------!
  !                 INPE/CGCT, Data Assimilation Office                  !
  !----------------------------------------------------------------------!
  !BOP
  ! !IROUTINE: parseCtl_ - parse a GrADS control (ctl) file
  !
  ! !DESCRIPTION: (to do)

  function parseCtl_(self, ctl) result(status)
    class(GrADSfiles), intent(inout) :: self
    character(len=*),  intent(in   ) :: ctl
    integer                          :: status

    character(len=*), parameter :: myname_ = myname//':: parseCtl_'

    character(len=4), parameter :: dimension(3) = ['xdef','ydef','zdef']
    character(len=4)            :: dname

    logical :: formatdefined
    integer :: ierr, val, ntokens, nvars
    integer :: i, j, k, d
    integer :: lu, ios, idx
    integer :: pmark
    integer :: currSize, lastSize
    integer :: xdef, ydef
    real    :: inip
    real    :: delta
    character(len=64) :: str, param
    character(len=64), allocatable :: tokens(:)

    !
    ! parse optional variables
    !

    status = 0

    !----------------------------------------
    ! Use m_inpak90 read the table file
    !----------------------------------------
    call i90_loadf(ctl, ierr)
    if (ierr /= 0) then

       call i90_perr(trim(myname_), ': i90_loadf('//trim(ctl)//')', ierr)
       !if (.not. present(stat)) call i90_die(myname_)
       status = ierr
       return

    endif

    !
    ! Initialize default values
    !

!    call self%reset( )

    self%iosize = 4

    !----------------------------------------
    !  Mandatory GrADS settings:
    !
    !    dset xdef ydef zdef tdef vars
    !----------------------------------------
    ! DSET

    call i90_label('DSET', ierr)
    if (ierr == 0) call i90_gtoken(self%dset, ierr)
    if (ierr /= 0) then
       call i90_perr(trim(myname_), ': DSET error with "'//trim(ctl)//'"', ierr)
       status = ierr
       return
    endif
    if (self%dset(1:1) == '^') then
       self%dset = self%dset(2:)
       i = index(ctl, '/', back=.true.)
       if (i > 0) self%dset = ctl(1:i)//self%dset
    endif

    !----------------------------------------
    !DTYPE

    call i90_label('DTYPE', ierr)
    if (ierr == 0) then
       call i90_gtoken(str, ierr)
       if (ierr /= 0)then
          call i90_perr(trim(myname_), ': DTYPE error with "'//trim(ctl)//'"', ierr)
          status = ierr
          return
       endif

       self%dtype = i90_lcase(str)

    elseif(ierr == -2)then
#ifdef DEBUG
       call i90_perr(trim(myname_), ': setting default data type, "Binary ieee"')
#endif
       self%dtype = 'binary'
    else
       call i90_perr(trim(myname_), ': DTYPE error with "'//trim(ctl)//'"', ierr)
       status = ierr
       return
    endif

    !----------------------------------------------------------------!
    ! Get dimensions Xdef, Ydef, Zdef
    select case(self%dtype)

       case('binary','grib')

          do d = 1,size(dimension)

             if (.not.associated(self%dhead))then
                allocate(self%dhead)
                self%dtail => self%dhead
             else
                allocate(self%dtail%next)
                self%dtail => self%dtail%next
             endif

             self%dtail%dName = dimension(d)
             
             call i90_label(trim(dimension(d)), ierr)
             if (ierr == 0) self%dtail%dSize = i90_gint(ierr)
             if (ierr /= 0) then
                call i90_perr(trim(myname_), ': '//TRIM(dimension(d)), ierr)
                status = ierr
                return
             endif

             if (self%dtail%dSize > 0) then
                call i90_gtoken(str, ierr)
                if (ierr /= 0) then
                   call i90_perr(trim(myname_),': '//trim(dimension(d))//' type', ierr)
                   return
                endif
         
                allocate (self%dtail%dLevs(1:self%dtail%dSize), stat=ierr)
                if (ierr /= 0) then
                   call i90_perr(trim(myname_),': allocate(levs)', ierr)
                   status = ierr
                   return
                endif
         
                self%dtail%dType = i90_lcase(str)
                select case (self%dtail%dType)
                case ('levels')
         
                   i = 1
                   do while (i <= self%dtail%dSize)
                      self%dtail%dLevs(i) = i90_gfloat(ierr)
                      if (ierr /= 0) then
                         call i90_gline(ierr)
                      else
                         i = i + 1
                      endif
                   end do
         
                   if (ierr /= 0) then
                      call i90_perr(trim(myname_), ': '//trim(dimension(d))//' level', ierr)
                      status = ierr
                      return
                   endif
         
                case ('linear')
                   inip = i90_gfloat(ierr)
                   if (ierr /= 0) then
                      call i90_perr(trim(myname_), ': '//trim(dimension(d))//' level', ierr)
                      status = ierr
                      return
                   endif
         
                   delta = i90_gfloat(ierr)
                   if (ierr /= 0) then
                      call i90_perr(trim(myname_),': '//trim(dimension(d))//' delta lev', ierr)
                      status = ierr
                      return
                   endif
         
                   do i = 1, self%dtail%dSize
                      self%dtail%dLevs(i) = inip + delta*(i - 1)
                   enddo
                
                case default
                   call i90_perr(trim(myname_),': unknown '//trim(dimension(d))//' type, "'//trim(str)//'"',99)
                   status = 99
                   return
                end select
             endif
          enddo

       case('netcdf')
          call i90_perr(trim(myname_), 'netcdf not yeat implemented via ctl!')
          status = -1
          return       
       case default
          call i90_perr(trim(myname_), &
               ': Wrong file type : "'//trim(self%dtype)//'" error, ierr =', 99)
          status = 99
          return
      
    end select
!!    !----------------------------------------------------------------!
!!    ! TDEF
!!    self%tdef%num = -1
!!    
!!    call i90_label('TDEF', ierr)
!!    if (ierr == 0) gs%tdef%num = i90_gint(ierr)
!!    if (ierr /= 0) then
!!       call i90_perr(trim(myname), ': TDEF entry error', ierr)
!!       if (.not. present(stat)) call i90_die(myname_)
!!       stat = ierr
!!       return
!!    endif
!!    !tdim = gs%tdef%num
!!
!!    !-----------------------------------------
!!    ! Optional Settings
!!    !
!!    formatdefined = .false.
!!
!!    call i90_label('OPTIONS', ierr)
!!    if (ierr == 0) then
!!       do
!!
!!          call i90_gtoken(str, ierr)
!!          if (ierr .ne. 0) exit
!!
!!          str = i90_lcase(str)
!!          select case (trim(str))
!!             !case ('pascals')
!!             !   gs%opt%pascals = .true.
!!          case ('yrev')
!!             gs%opt%yrev = .true.
!!          case ('zrev')
!!             gs%opt%zrev = .true.
!!             !case ('template')
!!             !   gs%opt%template = .true.
!!          case ('sequential')
!!             gs%opt%iacc = .false.
!!             !case ('cal365day')
!!             !   gs%opt%cal365day = .true.
!!          case ('byteswapped')
!!             gs%opt%byteswapped = .true.
!!          case ('big_endian')
!!             gs%opt%endianess = .false.
!!          case default
!!
!!             call i90_perr(trim(myname), ': unsupported option, "'//trim(str)//'"', ierr)
!!             if (.not. present(stat)) call i90_die(myname_)
!!             stat = -3
!!             return
!!
!!          end select
!!
!!       enddo
!!    endif
!!
    !----------------------------------------
    ! UNDEF

    self%undef = UDEF_


    call i90_label('UNDEF', ierr)
    if (ierr == 0) self%undef = i90_gfloat(ierr)
    if (ierr /= 0) then
       call i90_perr(trim(myname), ': UNDEF entry error', ierr)
       status = -3
       return
    endif

    !----------------------------------------
    ! VARS -- ENDVARS

    self%nvars = -1

    call i90_label('VARS', ierr)
    if (ierr == 0) self%nvars = i90_gint(ierr)
    if (ierr /= 0) then
       call i90_perr(trim(myname), ': VARS entry error', ierr)
       status = ierr
       return
    endif

!    nvars = self%nvars
!
!    allocate (self%vars(nvars), stat=ierr)
!
!    if (ierr /= 0) then
!       call i90_perr(trim(myname_), ': allocate(VARS) error', ierr)
!       status = ierr
!    endif

    !     Get variable names and labels
    !     -----------------------------
    xdef = self%getDim('xdef')
    ydef = self%getDim('ydef')

    pmark    = 0
    lastSize = 1+self%iosize
    do i = 1, self%nvars

       if (i==1)then
          allocate(self%vHead)
          self%vTail => self%vHead
       else
          allocate(self%vTail%next)
          self%vTail => self%vTail%next
       endif

       call i90_gline(ierr)
       if (ierr /= 0) then
          call i90_perr(trim(myname_), ': error to get var list', ierr)
          status = ierr
          return
       endif

       !var index
       !self%vars(i)%id = i

       !get var name
       call i90_gtoken(str, ierr)
       self%vTail%name = trim(str)

       ! get var nlevels
       k = i90_gint(ierr)
       self%vTail%nlevs = max(k, 1)

       !--------------------------------------------------!
       !get var info, if is a grid type, we should get
       !at least one parameter

       call i90_gtoken(str, ierr)
       param = trim(str)
       j = index(str, ',')
       if (j > 0) then
          do while (j > 0)
             call i90_gtoken(str, ierr)
             j = index(str, ',')
             if (j <= 0) then
                ! determine if is a number (should be a integer)
                Read (str, '(I10)', iostat=ios) val
                if (ios .eq. 0) then
                   param = trim(param)//trim(str)
                endif
             else
                param = trim(param)//trim(str)
             endif
          enddo
       endif

       !split parameters

       call split(param, ntokens, tokens, ',')
       !if a grib file will have at least 3 extra parameters
       select case(self%dtype)
          case('netcdf')
             call i90_perr(trim(myname_), 'netcdf not yeat implemented via ctl!')
             status = -1
             return  
          case default
             allocate (self%vTail%param(ntokens))
             do j = 1, ntokens
                read (tokens(j), '(I10)') self%vTail%param(j)
             enddo
       end select

       !--------------------------------------------------!
       !Get description var
       ierr = 0
       self%vTail%descr = ''
       do while (ierr == 0)
          call i90_gtoken(str, ierr)
          self%vTail%descr = trim(self%vTail%descr)//' '//trim(str)
       enddo

       !
       ! Obtain Position marker to read directly a IEEE-32 file.
       ! We assume that the integer value occupies a 32-bit (4-byte)
       ! for a single precision file and 64-bit (8-byte) for a double
       ! precion
       !

       self%vTail%pmark = pmark + lastSize
       pmark            = self%vTail%pmark
       lastSize         = (xdef*ydef*self%iosize+self%iosize*2) * self%vTail%nlevs

    enddo


    !----------------------------------------
    call i90_release(ierr)
    if (ierr /= 0) then
       call i90_perr(myname_, 'i90_release()', ierr)
       status = ierr
       return
    endif
    !----------------------------------------
!    !
!    lu = self%lu
!    if (lu >= 0) close (lu)
!    self%lu = -1
!
!    !--------------------------------------------------------
!    ! allocate the input buffer
!
!    !      allocate (gs%dbuf(gs%xdef%num, gs%ydef%num), stat=ierr)
!    !      if (ierr /= 0) then
!    !         write (stderr, '(2a,i5)') myname_, &
!    !            ': allocate(gs%dbuf) error, stat =', ierr
!    !         if (.not. present(stat)) call die(myname_)
!    !         stat = ierr
!    !         return
!    !      endif
!!    gs%lu = luavail()
!    self%lu = i90_lua(exclude)
!    i = minloc(exclude,1)
!    exclude(i) = self%lu
!
!!    self%irec = 1
!!    call opendset_(self%lu, self%dset, self%dtype, ierr)
!!    if (ierr /= 0) then
!!       call i90_perr(myname_, 'opendset_()', ierr)
!!       status = ierr
!!       return
!!    endif
!
  end function parseCtl_

  !-----------------------------------------------------------------------------!
  !             Modeling and Development Division - DMD/CPTEC/INPE              !
  !-----------------------------------------------------------------------------!
  !BOP
  !
  ! !IROUTINE: split - parse string into an array using specified delimiters
  !
  ! !DESCRIPTION: parses a string using specified delimiter characters and
  !               store tokens into an allocatable array
  !
  !
  ! !INTERFACE:
  !

  subroutine split(str, ntokens, tokens, del)

    !
    ! !INPUT PARAMETERS:
    !
    character(len=*), intent(in) :: str
    character(len=*), optional, intent(in) :: del
    !
    ! !OUTPUT PARAMETERS:
    !
    integer, intent(out) :: ntokens
    character(len=*), allocatable, intent(out) :: tokens(:)

    ! !REVISION HISTORY:
    !
    !   13 May 2020 - J. G. de Mattos -  Initial code.
    !
    !EOP
    !-----------------------------------------------------------------------------!
    !BOC

    character, parameter :: BLK = achar(32)   ! blank (space)
    character(len=1)     :: delimiter
    integer              :: i, j
    integer              :: StrLen

    ! linked list to store temporary tokens
    type token
       character(len=10)    :: tk
       type(token), pointer :: next => null()
    endtype token
    type(token), pointer :: root => null()
    type(token), pointer :: current => null()

    ! setting up delimter
    delimiter = BLK
    if (present(del)) delimiter = del

    ! get string length
    StrLen = len_trim(str)

    ! at least we has one token
    ntokens = 1

    ! find tokens using delimiter
    allocate (root)
    current => root
    j = 1
    do i = 1, StrLen

       if (str(i:i) == trim(delimiter)) then
          ntokens = ntokens + 1
          current%tk = str(j:i - 1)
          allocate (current%next)
          current => current%next
          j = i + 1
       endif

    enddo
    !get last token
    current%tk = str(j:len_trim(str))

    !copy tokens to output array
    allocate (tokens(ntokens))
    current => root
    do i = 1, ntokens
       tokens(i) = trim(current%tk)
       current => current%next
    enddo

    !
    ! deallocate temporary token list
    !
    current => root%next
    do while (associated(current))
       deallocate (root)
       root => current
       current => root%next
    enddo

  end subroutine split
  !EOC
  !-----------------------------------------------------------------------------!

  function qctlinfo(self)result(status)
     class(GrADSFiles), intent(inout) :: self
     integer                          :: status
     status = 0

     write(*,'(A7,1x,A1,A)')'dset','^',trim(adjustl(self%dset))
     write(*,'(A7,1x,ES9.3)')'undef', self%undef
     write(*,'(A7,1x,A6)')'dtype',trim(adjustl(self%dtype))
     write(*,'(A7,1x,A)')'options'
  end function

  function getDim(self, dName)result(dSize)
     class(GrADSFiles), intent(inout) :: self
     character(len=*),  intent(in   ) :: dName
     integer                          :: dSize

     type(GrADSDef), pointer :: d => null()
     character(len=4) :: iqDim

     
     iqDim = trim(i90_lcase(dName))
     dSize = -1
     d => self%dHead
     do while(associated(d))
        if (trim(iqDim) == trim(d%dName))then
           dSize = d%dSize
           exit
        endif
        d => d%next
     enddo
     return
  end function


  subroutine getDimVec(self, dName, Array, stat)
     class(GrADSFiles), intent(inout) :: self
     character(len=*),  intent(in   ) :: dName
     real, pointer,     intent(inout) :: Array(:)
     integer, optional, intent(  out) :: stat

     type(GrADSDef), pointer :: d => null()
     character(len=4) :: iqDim

     
     iqDim = trim(i90_lcase(dName))
     d => self%dHead
     do while(associated(d))
        if (trim(iqDim) == trim(d%dName))then
           Array => d%dLevs
           exit
        endif
        d => d%next
     enddo
     return
  end subroutine


  function pVars(self)result(status)
     class(GrADSFiles), intent(inout) :: self
     integer                          :: status

     type(GrADSVar), pointer :: v => null()

     status = 0

     v => self%vHead
     do while(associated(v))
        print*,trim(v%Name),v%pmark, v%nlevs, v%param
        v => v%next
     enddo

  end function

  function getField1d(self, name, klev, time, field)result(status)
     class(GrADSFiles), intent(inout) :: self
     character(len=*),  intent(in   ) :: name     ! what variable?
     integer,           intent(in   ) :: time     ! what time?
     real,              intent(in   ) :: klev     ! which level?
     real,              intent(inout) :: field(:) ! a 1-d gridded field
     integer                          :: status

     character(len=*), parameter :: myname_ = myname//':: getField1D'

     integer :: iqSize, gSize
     character(len=10) :: inqVar
     integer :: xdef, ydef

     integer(i4) :: iret
     integer(i8) :: ipos
     real, pointer :: levs(:) => null()

     type(GrADSVar), pointer :: var => null()
     type(gribType) :: grb

     status = 0

     !
     ! Sanity checks
     !

     ! check requested variable
          
     var => self%getVarPtr(name)
     if(.not.associated(var))then
       write (stderr, '(4a)') myname_, ': unknown variable "', trim(name), '"'
       status = 3
       return     
     endif
     
  
     ! Check the buffer dimensions
     iqSize = size(field)
     xdef   = self%getDim('xdef')
     ydef   = self%getDim('ydef')
     gSize  = xdef*ydef

     if (iqSize /= gSize)then
       write (stderr, '(2a)', advance="no") myname_, ': invalid arguments'
       write (stderr, '(a,2i6,a)', advance="no") ', shape(field) = (', size(field), ')'
       write (stderr, '(a,2i6,a)', advance="no") ', [xy]def = (', xdef, ydef, ')'
       write (stderr, *)
       status = 2
       return     
     endif

     ! Check the requested time

       ! time Not yeat implemented

     ! Check the requested level

     if (klev < 0 .or.  &
         klev > var%nlevs)then
        write (stderr, '(2a)', advance="no") myname_, ': invalid level request'
        write (stderr, '(a,i3,3a,i3)') ', klev =', klev, &
             ', klev("', trim(name), '") =', var%nlevs
        status = 5
        return
     endif
     
    !
    ! Get Variable
    !

    select case (self%dtype)
       case('binary')
       !--------------------------------------------------------
       ! The open statement used to access binary file in this
       ! module uses standard stream I/O,  a feature
       ! that was intruduced in fortran 2003, so we need obtain
       ! the start position of each field recorded.
       ! For single precision files we need known the size of
       ! each field plus 4-bit and for double precision files
       ! we need known the size of field plus 8-bits. Because this
       ! we have a number 4 at nrec calculus.
       !--------------------------------------------------------
       ! Compute the record number
          ipos = var%pmark + ( (((xdef*ydef*self%iosize)+self%iosize*2)) * (klev-1) )
          iret = openIEEE(unit=99,file=self%dset,status='old',action='read')
          call readIEEE(99,field,.false.,ipos,iret)
       case('grib')
          call grb%open(self%dset)
          levs => grb%getVarLevels(pds5=var%param(1))
          call grb%getField(var%param(1),levs(klev),field)
    end select


  end function

  function getVarPtr(self,name)result(ptr)
     class(GrADSFiles), intent(inout) :: self
     character(len=*),  intent(in   ) :: name     ! what variable?
     type(GrADSVar), pointer :: ptr

     character(len=10) :: inqVar


     nullify(ptr)
     inqVar =trim(adjustl(i90_lcase(name)))
     ptr => self%vHead 
     do while(associated(ptr))
       if ( trim(inqVar) == trim(ptr%name)) exit
       ptr => ptr%next
     enddo

     return
     
  end function

end module
