!-----------------------------------------------------------------------------!
!           Group on Data Assimilation Development - GDAD/CPTEC/INPE          !
!-----------------------------------------------------------------------------!
!BOI
!
! !TITLE: Documentação do pacote ``Read Grib'' \\ Version 1.0
!
! !AUTHORS: João G. Z. de Mattos
!
! !AFFILIATION: Grupo de Desenvolvimento em Assimilação de Dados, INPE/CPTEC
!
! !DATE: 20 de Julho de 2014
!
! !INTRODUCTION: Visão Geral do Pacote
!
!      O pacote ReadGrid é um modulo escrito em Fortran 90 que contém rotinas 
!      para a leitura de arquivos no formato Grib1. Faz parte deste conjunto 
!      de rotinas o arquivo ascii que contem a tabela grib utilizada pelo 
!      CPTEC/INPE. Esta tabela também é utilizada como tabela padrao caso não
!      exista uma tabela específica para o grib que está sendo lido.
!
! \subsection{Exemplo de Leitura}
!
!  A seguir é mostrado um programa teste que faz uso deste módulo para realizar
!  a leitura de um arquivo grib1 proveniente do modelo global do INPE/CPTEC.
! 
! \begin{verbatin}
!
!   program RGrib
!      use read_grib
!
!      implicit none
!
!      type(grib) :: grb
!      real, allocatable :: fld(:,:)
!
!      grb%file="GPOSCPT20130101122013010118P.fct.TQ0299L064.grb"
!
!
!      call OpenGrib(grb)
!
!      allocate(fld(grb%gds(2),grb%gds(3)))
!   
!      call ReadGrib(grb,'UVEL',1000,fld)
!      print*,minval(fld),maxval(fld)
!
!      call ReadGrib(grb,'PSLC',0000,fld)
!      print*,minval(fld),maxval(fld)
!
!      call CloseGrib(grb)
!
!   end program
!
! \end{verbatin}
!
! Este pacote faz uso de alguns módulos que são listados a seguir e devem
! ser compilados em conjunto para que o módulo ReadGrib funcione corretamente.
!
! \begin{itemize}
!  \item {\bf m\_ioutil} : módulo contendo funcoes para operacoes de I/O;
!  \item {\bf m\_die}    : módulo contendo funcoes para mensagens de saida;
!  \item {\bf m\_stdio}  : módulo com definicoes de unidades de saida;
!  \item {\bf m\_string} : módulo com rotinas para manipulacao de strings.
! \end{itemize}
!
! Os três primeiros módulos fazem parte da biblioteca {\bf MPEU} desenvolvida por 
! membros do {\it Data Assimilation Office} da {\it National Aeronautics and Space
! Administration}. Esta biblioteca está disponível para download em
! \url{http://www.nco.ncep.noaa.gov/pmb/codes//nwprod/ngac.v1.0.0/sorc/ngac_fcst.fd/chem/gocart/src/GMAO_Shared/GMAO_mpeu/}.
!
!
!
!EOI
!-----------------------------------------------------------------------------!
!           Group on Data Assimilation Development - GDAD/CPTEC/INPE          !
!-----------------------------------------------------------------------------!
!
!-----------------------------------------------------------------------------!

module accessGrib

  use m_msg, only : perr         ! modulo contendo funcoes para mensagens de saida
  use m_string, only: GetTokens, lowercase  ! modulo com rotinas para manipulacao de strings
  use m_stdio                    ! modulo contendo unidades e funcoes para operacoes de I/O
  use coord_compute, only: compute_earth_coord 



  implicit none
  private

! !ROUTINES:

!  public :: OpenGrib
!  public :: ReadGrib
!  public :: CloseGrib

! !PUBLIC TYPES:

  public :: gribType


! !REVISION HISTORY:
!   20 Jul 2014 - J. G. de Mattos - Initial Version
!   05 Jan 2016 - J. G. de Mattos - Include interface to return
!                                   1d and 2d fields
!   06 Jan 2016 - J. G. de Mattos - Bug fix
!                                   Problem to read multiple files
!                                   in same routine, like read files
!                                   to make statistical evaluation
!
!   19 May 2022 - J. G. de Mattos - transform to class 
! !BUGS:
!   Not yet
!
!EOP
!---------------------------------------------------------------------!


  type gbtab
     integer           :: parm = -1
     character(len=10) :: name = ''
     character(len=50) :: desc = ''
  end type gbtab

  type array
     real :: p
     type(array), pointer :: next => null()
  end type

  type grbDim
    character(len=4)  :: dName
    integer           :: dSize
    real, pointer     :: dLevs(:)
    type(grbDim), pointer :: next => null()     
  end type

  type :: grbvardesc
     character(len=10)    :: name
     character(len=50)    :: desc
     integer              :: nlevs
     integer              :: pds5
     integer              :: pds6
     type(array), pointer :: pds7 => null()
     type(array), pointer :: pHead => null()
     type(grbVarDesc), pointer :: next => null()
  end type grbvardesc

  type gribType
     character(len=1024)           :: File     = ''
     integer                       :: lu       = -1
     integer                       :: PDS(200) = -1
     integer                       :: GDS(200) = -1
     integer                       :: nvar     = -1
     real                          :: undef    = -1e+20
     ! Variable definition
     type(grbvardesc), pointer     :: vHead => null()
     type(grbvardesc), pointer     :: vTail => null()
     ! Dimension Definition
     type(grbDim), pointer         :: dHead => null()
     type(grbDim), pointer         :: dTail => null()
     contains
        procedure, public  :: open
        procedure, public  :: close
        procedure,private  :: getVarPtr
        procedure, public  :: getVars
        procedure, public  :: getVarLevels
        procedure, public  :: getDim
        procedure, public  :: getDimVec
        procedure, public  :: stat
        procedure          :: input1d_, input1d_gIDp,&
                              input2d_, input2d_gIDp
        generic, public    :: getField => input1d_, input1d_gIDp,&
                                          input2d_, input2d_gIDp
  end type gribType

!  interface ReadGrib
!     module procedure input1d_gip, input1d_,&
!                      input2d_gip, input2d_
!  end interface

  Character(len=*), parameter :: myname = "accessdGrib"

contains
!-----------------------------------------------------------------------------!
!           Group on Data Assimilation Development - GDAD/CPTEC/INPE          !
!-----------------------------------------------------------------------------!
!
!BOP
!
! !IROUTINE: OpenGrib - Open a Grib1 file and populate the grib data type 
!
! !DESCRIPTION:
!
! !INTERFACE:

  subroutine open( self, file,  istat )
  
    class(gribType),   intent(inout) :: self
    character(len=*),  intent(in   ) :: file
    integer, optional, intent(  out) :: istat

    Character(len=*), parameter :: myname_ = trim(myname)//"::open"

    !
    ! local
    !

    integer, parameter ::  mbuf=256*1024

    type(gbtab), allocatable :: gtb(:)

    integer   :: jpds(200),jgds(200)
    integer   :: kpds(200),kgds(200)
    real      :: kg,kf,kk
    character :: cbuf(mbuf)
    integer   :: nlen
    integer   :: nnum
    integer   :: mnum

    integer :: iret
    integer :: idx
    integer :: nrec
    integer :: i
    integer :: dsize
    real, allocatable :: rlat(:)
    real, allocatable :: rlon(:)

    type(grbvardesc), pointer :: thisVar => null()
    type(grbvardesc), pointer :: var => null()

    if(present(istat))istat=0


    self%lu   = getAvailUnit( )
    self%file = trim(file)
    call baopenr(self%lu,self%file, iret)
    if (iret.ne.0)then
       call perr(myname_,'Error to open Grib File file:'//trim(self%file),iret)
       if(present(istat))istat=iret
       return
    endif


    iret = 0
    nrec = 0
    self%nvar = 0
    do while(iret.eq.0)

       jpds = -1
       jgds = -1

       if (.not.associated(self%vHead))then

          allocate(self%vHead)
          self%vTail => self%vHead
          allocate(self%vTail%pHead)
          self%vTail%pds7 => self%vTail%pHead
          self%vTail%nlevs = 1

          call getgbmh(self%lu,0,-1,jpds,jgds,mbuf,cbuf,nlen,nnum,mnum,kg,kf,kk,kpds,kgds,iret)
          if(iret/=0)then
             call perr(myname_,'Some error to read Grib file:'//trim(self%file),iret)
             if(present(istat))istat=iret
             return
          endif

          self%gds  = kgds
          self%pds  = kpds
          self%nvar = 1
      
          call GribTable(kpds(1),gtb)

       else

          call getgbmh(self%lu,0,nrec,jpds,jgds,mbuf,cbuf,nlen,nnum,mnum,kg,kf,kk,kpds,kgds,iret)

          if (kpds(5) == self%vTail%pds5)then
             allocate(self%vTail%pds7%next)
             self%vTail%pds7 => self%vTail%pds7%next
             self%vTail%nlevs = self%vTail%nlevs + 1

          else
             allocate(self%vTail%next)
             self%vTail => self%vTail%next
             self%vTail%nlevs = 1

             allocate(self%vTail%pHead)
             self%vTail%pds7 => self%vTail%pHead
             self%nvar = self%nvar + 1
          endif

       endif

        
       self%vTail%pds5   = kpds(5)
       self%vTail%pds6   = kpds(6)
       self%vTail%pds7%p = real(kpds(7))

      
       idx = minloc(gtb%parm-kpds(5),mask=gtb%parm-kpds(5).GE.0,DIM=1)
      
       if(kpds(5).eq.gtb(idx)%parm)then
          self%vTail%name = gtb(idx)%name
          self%vTail%desc = gtb(idx)%desc
       else
          write(self%vTail%name,'(A4,I3.3)')'NONE',self%vTail%pds5
          self%vTail%desc = 'None'
       endif
          

       nrec = nrec + 1
       
    enddo

    nrec = nrec - 1

    !
    ! define dimensions
    !

    allocate(rlon(self%gds(2)))
    allocate(rlat(self%gds(3)))

    call compute_earth_coord(self%gds, self%undef, rlon, rlat, iret )

    allocate(self%dhead)
    self%dtail => self%dhead

    self%dtail%dName = 'xdef'
    self%dtail%dSize = self%gds(2)
    allocate(self%dtail%dLevs(self%dtail%dSize))
    self%dtail%dLevs = rlon
    deallocate(rlon)

    allocate(self%dtail%next)
    self%dtail => self%dtail%next
    self%dtail%dName = 'ydef'
    self%dtail%dSize = self%gds(3)
    allocate(self%dtail%dLevs(self%dtail%dSize))
    self%dtail%dLevs = rlat
    deallocate(rlat)

    allocate(self%dtail%next)
    self%dtail => self%dtail%next
    self%dtail%dName = 'zdef'
    dSize = 0
    var => self%vHead
    do while (associated(var))
       if (dSize < var%nLevs ) then
          dSize = var%nlevs
          thisVar => var
       endif
       var => var%next
    enddo
    self%dtail%dSize = dSize
    allocate(self%dtail%dLevs(self%dtail%dSize))
    thisVar%pds7 => thisVar%pHead
    i=1
    do while (associated(thisVar%pds7))
       self%dtail%dLevs(i) = thisVar%pds7%p
       i=i+1
       thisVar%pds7 => thisVar%pds7%next
    enddo
    
#ifdef VERBOSE
    write(stdout,'(A,1X,I4.3,1X,A)')' Found',nrec,'2D fields'
    write(stdout,'(A)')trim(self%file)
    var => self%vHead
    do while(associated(var))
       write(stdout,'(2x,a4,3(1x,i6.2))')var%name,var%pds5,var%pds6,var%pHead%p
       if (associated(var%pHead%next))then
          var%pds7=>var%pHead%next
          do while (associated(var%pds7))
             write(stdout,'(2x,a4,3(1x,i6.2))')var%name,var%pds5,var%pds6,var%pds7%p
             !write(stdout,'(20x,1x,i6.2)')thisGrb%pds7%p
             var%pds7 => var%pds7%next
          enddo
       endif
       var => var%next
    enddo
#endif


  end subroutine open


  subroutine input2d_gIDp(self,gIDp,hlev,vfld,stat)
       
     class(gribType),  intent(inout) :: self
     integer,          intent(in   ) :: gIDp ! grid indicator of parameter
     real,             intent(in   ) :: hlev ! which level?
     real,             intent(  out) :: vfld(:,:) ! a 2-d gridded field
     integer, optional,intent(  out) :: stat

     character(len=*), parameter :: myname_=myname//':: input2d_gIDp'
     
     integer :: i, nidx, ivar
     character(len=10) :: vname
     type(grbVarDesc), pointer :: var => null()

     ! Check/index the requested variable

     var => self%getVarPtr(pds5=gIDp)

     if(.not.associated(var)) then

        write(stderr,'(2a,i4,a)') myname_,': unknown Grid Indicator of Parameter"',gIDp,'"'
        write(stderr,'(2a)') myname_,': return undefined value to variable "'
        if(present(stat))stat = 3
        vfld = self%undef
        return
        
     endif

     vname=var%name
     
     call self%input2d_(trim(vname),hlev,vfld,stat)
     

  end subroutine

  subroutine input1d_gIDp(self,gIDp,hlev,vfld,stat)
     
     class(gribType),  intent(inout)  :: self
     integer,          intent(in   )  :: gIDp    ! grid indicator of parameter
     real,             intent(in   )  :: hlev    ! which level?
     real,             intent(  out)  :: vfld(:) ! a 1-d gridded field
     integer, optional,intent(  out)  :: stat


     character(len=*), parameter :: myname_=myname//':: input1d_gIDp'
     integer :: i, nidx, ivar
     character(len=10) :: vname
     type(grbVarDesc), pointer :: var => null()
     ! Check/index the requested variable

     var => self%getVarPtr(pds5=gIDp)

     if(.not.associated(var)) then

        write(stderr,'(2a,i4,a)') myname_,': unknown Grid Indicator of Parameter"',gIDp,'"'
        write(stderr,'(2a)') myname_,': return undefined value to variable "'
!        if(.not.present(stat)) call die(myname_)
!        stat=3
        if(present(stat))stat = 3
        vfld = self%undef
        return
        
     endif

     !
     ! Get Var Name
     !

     vname=var%name
     
     call self%input1d_(trim(vname),hlev,vfld,stat)

  end subroutine


  subroutine input2d_(self,vnam,hlev,vfld,stat)
     implicit none
     
     class(gribType),  intent(inout) :: self
     character(len=*), intent(in   ) :: vnam      ! what variable?
     real,             intent(in   ) :: hlev      ! which level?
     real,             intent(  out) :: vfld(:,:) ! a 2-d gridded field
     integer, optional,intent(  out) :: stat

! !REVISION HISTORY:
!     20 Jul 2014 - J. G. de Mattos - initial prototyping and coding
!     19 may 2022 - J. G. de Mattos - change to class
!
!
!EOP
     character(len=*), parameter :: myname_=myname//':: input2d_'

     real, allocatable      :: ftmp(:)
     integer :: ipt, jpt
     integer :: i, j, k
     integer :: iret

     if(present(stat))stat=0

     ! Check the buffer dimensions

     if(size(vfld,1) /= self%gds(2) .or. size(vfld,2) /= self%gds(3)) then

        write(stderr,'(2a)', advance="no") myname_,': invalid arguments'
        write(stderr,'(a,1i6,a)', advance="no") ', shape(vfld) = (',shape(vfld),')'
        write(stderr,'(a,2i6,a)', advance="no") ', grb%pds(2:3) = (',self%gds(2),self%gds(3),')'
        write(stderr,*)
        !if(.not.present(stat)) call die(myname_)
        if(present(stat))stat=2
        return

     endif


     ipt = size(vfld,1)
     jpt = size(vfld,2)

     allocate(ftmp(ipt*jpt))

     
     call self%input1d_(vnam,hlev,ftmp,iret)
     if(present(stat))stat=iret

     !
     ! Colocando em uma matriz 2d
     !

     k = 0
     do j=1,self%gds(3)
        do i=1,self%gds(2)
           vfld(i,j) = ftmp(i+k)
        enddo
        k = k + self%gds(2)
     enddo


  end subroutine


  subroutine input1d_(self,vnam,hlev,vfld,stat)
     
     class(gribType),  intent(inout) :: self
     character(len=*), intent(in   ) :: vnam    ! what variable?
     real,             intent(in   ) :: hlev    ! which level?
     real,dimension(:),intent(  out) :: vfld    ! a 1-d gridded field
     integer, optional,intent(  out) :: stat

! !REVISION HISTORY:
!     20 Jul 2014 - J. G. de Mattos - initial prototyping and coding
!     19 may 2022 - J. G. de Mattos - change to class
!
!EOP

     character(len=*), parameter :: myname_=myname//':: input1d_'

     !
     ! Local
     !
     integer :: i, j, k
     integer :: ierr
     integer :: ivar
     integer :: vpds
     integer :: hcount
     character(len=512) :: fmt
     type(grbVarDesc), pointer :: var => null()

     !
     ! use to Grib fields
     !

     integer :: npts
     integer :: jpds(200),jgds(200)
     integer :: kpds(200),kgds(200)
     real    :: kf,kh
     logical, allocatable :: lb(:)
     real, allocatable      :: f(:)

 
     if(present(stat))stat=0


     !
     ! Sanity checks
     !

     ! Check the buffer dimensions
     npts = self%gds(2)*self%gds(3)

     if(size(vfld) /= npts) then

        write(stderr,'(2a)') myname_,': invalid arguments'
        write(stderr,'(a,i6,a)') ', shape(vfld) = (',shape(vfld),')'
        write(stderr,'(a,i6,a)') ', grb%pds(2:3) = (',npts,')'
        write(stderr,*)
        !if(.not.present(stat)) call die(myname_)
        if(present(stat))stat=2
        return

     endif

     ! Check/index the requested variable

     var => self%getVarPtr(name=vnam)
     if(.not.associated(var)) then

        write(stderr,'(4a)') myname_,': unknown variable "',trim(vnam),'"'
        write(stderr,'(4a)') myname_,': return undefined value to variable "',trim(vnam),'"'
!        if(.not.present(stat)) call die(myname_)
!        stat=3
        if(present(stat))stat = 0
        vfld = self%undef
        return
        
     endif

     ! Check the requested level
     var%pds7=>var%pHead
     hcount = 0
     do while(associated(var%pds7))
        hcount = hcount + 1
        if(var%pds7%p == hlev) exit
        var%pds7 => var%pds7%next
     enddo

     if( .not. associated(var%pds7) )then

        write(stderr,'(2a)') myname_,': invalid level request'
        write(stderr,'( a,i4)')', Level =',int(hlev)
        write(stderr,'(3a,i3)')'grb%nlevs("',trim(vnam),'") =', hcount

        write(fmt,'(a,i5,a)')'(a,',hcount,'I4)'
        var%pds7=>var%pHead
        do while(associated(var%pds7))
           write(stderr,fmt)', grb%levs[hPa] =',var%pds7%p
           var%pds7 => var%pds7%next
        enddo

        if(present(stat))stat=5
        return

     endif

     !
     ! Get Variable
     !

     npts    = self%gds(2)*self%gds(3)
     jpds    = -1
     jgds    = -1
     kpds    = -1
     kgds    = -1

     allocate(lb(npts), f(npts), stat=ierr)

      if(ierr/=0) then
         call perr(myname_,'allocate()',ierr)
        ! if(.not.present(stat)) call die(myname_)
         if(present(stat))stat=ierr
         return
     endif

     lb = .true.

     jpds(5) = var%pds5
     jpds(6) = var%pds6
     jpds(7) = hlev
     call getgb(self%lu,0,npts,-1,jpds,jgds,kf,kh, &
                kpds,kgds,lb,f,ierr)

     if(ierr /= 0) then
        write(stderr,'(2a,i5)') myname_,    &
         ': getgb() error, ierr =',ierr
        !if(.not.present(stat)) call die(myname_)
        if(present(stat))stat=ierr
        return
     endif

     !
     ! aplicando o undef onde lb .false.
     !
     where(.not.lb) f = self%undef

     vfld = f

  end subroutine


  subroutine GribTable(center,gtb, istat)

    implicit none

    integer, intent(in)                     :: center
    type(gbtab), allocatable, intent(inout) :: gtb(:)
    integer, optional, intent(out)          :: istat

    Character(len=*), parameter :: myname_ = trim(myname)//"::GribTable"


    real                            :: x
    integer                         :: i
    integer                         :: lu
    logical                         :: iret
    integer                         :: nline
    integer                         :: nbf
    character(len=200)              :: GribTableFile
    character(len=256)              :: bf
    character(len=050), allocatable :: bf2(:)
    integer, parameter              :: len_header = 400

    if(present(istat))istat=0

    !
    ! Reading Grib Table
    !

    write(GribTableFile,'(A,I3.3,A)')'gribtab.',center,'.tab'
    inquire(file=trim(GribTableFile),exist=iret)

    if(.not.iret)then

       call perr(myname_,'GribTable '//trim(GribTableFile)//' not found! Will use gribtab.default.tab')
       GribTableFile="gribtab.default.tab"
       inquire(file=trim(GribTableFile),exist=iret)
       if(.not.iret) then
          call perr(myname_,'Default GribTable not found! Will stop...')
          if(present(istat))istat=-1
          call exit()
       endif

    endif

    lu = getAvailUnit( )
    open( unit = lu,                &
         file = trim(GribTableFile),&
         action = 'read'            &
         )

    !
    ! count lines on grib table
    !

    nline = 0
    do
       read(lu,*,end=88,err=88)
       nline = nline + 1
    enddo
88  continue
    nline = nline - 1  ! remove header from count
    rewind(lu)

    allocate(gtb(nline))
    allocate(bf2(len_header))

    !
    ! Read Table
    !

    read(lu,'(A)')bf !skip first line

    do i=1,nline

       read(lu,'(A)')bf
       call GetTokens(bf,bf2,nbf,":")

       read(bf2(1),*)x
       gtb(i)%parm = int(x)
       gtb(i)%name = trim(bf2(2))
       gtb(i)%desc = trim(bf2(3))

    enddo
    deallocate(bf2)
    close (lu)

  end subroutine GribTable
!
!  function lindex_(lsts,entr)
!    use m_chars, only : uppercase
!    implicit none
!    character(len=*),dimension(:),intent(in) :: lsts
!    character(len=*),             intent(in) :: entr
!
!    integer :: lindex_    ! the result
!
!    !--------------------------------------------------------
!    integer :: i
!
!    !--------------------------------------------------------
!    lindex_=0
!    do i=1,size(lsts)
!      if(uppercase(entr) == uppercase(lsts(i))) then
!         lindex_=i
!         return
!      endif
!    end do
!  end function lindex_
!
  function close(self)result(stat)
     class(gribType),  intent(inout) :: self
     integer                         :: stat

     Character(len=*), parameter :: myname_ = trim(myname)//"::Close"

     type(grbVarDesc), pointer :: current => null()
     type(grbVarDesc), pointer :: next => null()

     integer :: lu
     integer :: ier

     stat=0

     lu=self%lu

     call baclose(lu,ier)

     if(ier/=0) then
        call perr(myname_,'baclose()',ier)
        !if(.not.present(stat)) call die(myname_)
        stat=ier
        return
     endif

     self%File       = ''
     self%lu         = -1
     self%PDS(1:200) = -1
     self%GDS(1:200) = -1
     self%nvar       = -1
     self%undef      = -1e+20
     
     self%vTail => self%vHead%next
     do
        self%vTail%pds7 => self%vTail%pHead%next
        do
           deallocate(self%vTail%pHead)
           if (.not. associated(self%vTail%pds7)) exit
           self%vTail%pHead => self%vTail%pds7
           self%vTail%pds7 => self%vTail%pHead%next
        enddo
        deallocate(self%vHead)
        if (.not. associated(self%vTail)) exit
        self%vHead => self%vTail
        self%vTail => self%vHead%next
     enddo


  end function close
!
!  function getVarPtrN(self, name)result(ptr)
!     class(gribType),  intent(inout) :: self
!     character(len=*), intent(in   ) :: name
!     type(grbvardesc), pointer       :: ptr
!
!     character(len=20) :: inqVar
!     character(len=20) :: locVar
!
!     nullify(ptr)
!     inqVar =trim(adjustl(lowercase(name)))
!     ptr => self%vHead 
!     do while(associated(ptr))
!       locVar = trim(adjustl(lowercase(ptr%name)))
!
!       if ( trim(inqVar) == trim(locVar)) exit
!       ptr => ptr%next
!     enddo
!
!     return
!  end function

  function getVarPtr(self, pds5, name)result(ptr)
     class(gribType),            intent(inout) :: self
     character(len=*), optional, intent(in   ) :: name
     integer,          optional, intent(in   ) :: pds5
     type(grbvardesc), pointer       :: ptr

     Character(len=*), parameter :: myname_ = trim(myname)//":: getVarPtr"

     character(len=20) :: inqVar
     character(len=20) :: locVar
     logical :: ok

     !
     ! check
     !
     ok = present(name).or.present(pds5)
     if(.not.ok)then
        call perr(myname_,'Should pass VarName or PDS 5', -3)
        return
     endif
 
     nullify(ptr)

     if (present(pds5))then
        ptr => self%vHead 
        do while(associated(ptr))
          if ( pds5 == ptr%pds5) return
          ptr => ptr%next
        enddo
     endif
     
     if(present(name))then
        inqVar =trim(adjustl(lowercase(name)))
        ptr => self%vHead 
        do while(associated(ptr))
          locVar = trim(adjustl(lowercase(ptr%name)))   
          if ( trim(inqVar) == trim(locVar)) return
          ptr => ptr%next
        enddo
     endif

     return
  end function


  function getVars(self)result(vnames)
     class(gribType), intent(inout) :: self
     character(len=60),     pointer :: vnames(:)

     Character(len=*), parameter :: myname_ = trim(myname)//":: getVars"

     type(grbvardesc), pointer       :: ptr => null()
     integer :: i


     allocate(vnames(self%nvar))
     ptr => self%vHead 
     do i=1,self%nvar

       vnames(i) = trim(ptr%name)
       
       ptr => ptr%next
     enddo
     return

  end function

  function getVarLevels(self,vname,pds5)result(levels)
     class(gribType),            intent(inout) :: self
     character(len=*), optional, intent(in   ) :: vname
     integer,          optional, intent(in   ) :: pds5
     real, pointer                             :: levels(:)

     Character(len=*), parameter :: myname_ = trim(myname)//":: getVarLevels"


     type(grbvardesc), pointer  :: var => null()
     type(array),      pointer  :: lev => null()
     integer                    :: i
     logical                    :: ok

     ok=present(vname).or.present(pds5)

     if(.not.ok)then
        call perr(myname_,'Should pass VarName or PDS 5', -3)
        return
     endif

     var => self%getVarPtr(name=vname,pds5=pds5)
     allocate(levels(var%nlevs))

     lev => var%pHead
     do i=1,var%nlevs
        levels(i) = lev%p
        lev => lev%next
     enddo

     return     

  end function

  function stat(self)result(status)
     class(gribType),  intent(inout) :: self
     integer                         :: status
     Character(len=*), parameter :: myname_ = trim(myname)//":: stat"

     character(len=60), pointer :: vnames(:)
     real, pointer :: levs(:)
     integer :: i, j
     integer :: npt
     real, allocatable :: field(:)
     real :: mean, max, min, std, rms
     
     status = 0

     npt = self%gds(2)*self%gds(3)
     vnames => self%getVars()
     allocate(field(npt))
     do i = 1, self%nvar
        write(stdout,*)''
        write(stdout,*)''
        write(stdout,'(2(4x,A4),5(4x,A12))')'Name','Lev','Min','Max','MEAN','STDV','RMS'
        write(stdout,'(A96)')'------- ------- --------------- --------------- &
                              --------------- --------------- ---------------'
        levs => self%getVarLevels(vnames(i))
        do j=1,size(levs)
           call self%getField(vnames(i),levs(j),field)
           
           min  = minval(field)
           max  = maxval(field)
           mean = sum(field)/npt
           std  = sqrt( sum((field - mean)**2)/(npt-1) )
           rms  = sqrt( sum((field)**2)/(npt) )
           write(stdout,'(4x,A4,4x,I4,5(4x,G12.4))')trim(vnames(i)), levs(j), min,max,mean,std,rms
        enddo
     enddo
     
  end function

  function getDim(self, dName)result(dSize)
     class(gribType),  intent(inout) :: self
     character(len=*), intent(in   ) :: dName
     integer                         :: dSize

     Character(len=*), parameter :: myname_ = trim(myname)//"::getDim"

     type(grbDim), pointer  :: gDim => null()
     
     character(len=4) :: inqName


     dSize = -1

     select case(lowercase(dName))
        case('xdef','xdim','xpt','nlon')
           inqName = 'xdef'
        case('ydef','ydim','ypt','nlat')
           inqName = 'ydef'
        case('zdef','zdim','zpt','nlevs')
           inqName = 'zdef'
        case default
           call perr(myname_,'Unknown dim name: '//trim(dName), -3)
           return
     end select

     gDim => self%dHead
     do while(associated(gDim))
        if (gDim%dName == inqName) then
           dSize = gDim%dSize
           return
        endif
        gDim => gDim%next
     enddo

     return

  end function

 subroutine getDimVec(self, dName, Array, stat)
     class(gribType),   intent(inout) :: self
     character(len=*),  intent(in   ) :: dName
     real, pointer,     intent(inout) :: Array(:)
     integer, optional, intent(  out) :: stat

     Character(len=*), parameter :: myname_ = trim(myname)//"::getDimVec"

     type(grbDim), pointer  :: gDim => null()
     
     character(len=4) :: inqName


     if(present(stat)) stat = 0

     select case(lowercase(dName))
        case('xdef','xdim','xpt','nlon')
           inqName = 'xdef'
        case('ydef','ydim','ypt','nlat')
           inqName = 'ydef'
        case('zdef','zdim','zpt','nlevs')
           inqName = 'zdef'
        case default
           call perr(myname_,'Unknown dim name: '//trim(dName), -3)
           if(present(stat)) stat = -3
           return
     end select

     gDim => self%dHead
     do while(associated(gDim))
        if (gDim%dName == inqName) then
           Array => gDim%dLevs
           return
        endif
        gDim => gDim%next
     enddo

     return

  end subroutine



end module accessGrib
