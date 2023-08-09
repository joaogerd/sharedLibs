module accessNetcdf
   ! The standard is for COADS format,
   ! for use another format, we need a description file
   ! to inform the dimensions name. Follow an example:
   !
   !
   !
   use netcdf
   use m_string
   use m_inpak90, only: i90_lcase
   implicit none
   private

   public :: ncType
   
   character(len=4) :: Axis(0:4) = ['None','xdef','ydef','zdef','tdef']

   type def
      character(len=4)   :: name_
      character(len=60)  :: ncName_
      integer            :: ncId_
      integer            :: size_
      real, allocatable  :: array_(:)
      type(def), pointer :: next => null()
   end type

   type var
      character(len=20)    :: name_
      integer              :: id_
      integer              :: ndims
      real                 :: scaleFactor
      real                 :: addOffSet
      integer, allocatable :: dimIds(:)
      integer, allocatable :: dimSize(:)
      type(var), pointer :: next => null()
   end type

   type ncType
      private
      character(len = 1024) :: ncFile
      integer               :: ncid ! Opened Unit
      
      integer               :: ndims
      type(def), pointer    :: gridRoot => null()
      type(def), pointer    :: gridInfo => null()

      integer                :: nvars
      type(var), pointer     :: varRoot => null()
      type(var), pointer     :: varInfo => null()
      contains
         procedure, public  :: open => open_
         procedure, public  :: getDim => getDim_
         procedure, public  :: getDimVec => getDimVec_
         procedure, public  :: getField => getField_
         procedure, private :: getDimId => getDimId_
         procedure, private :: dimIsThere => dimIsThere_

   end type


   contains

   subroutine open_(self, fileName, xdim, ydim, zdim, tdim)
      class(ncType),    intent(inout) :: self
      character(len=*), intent(in   ) :: fileName
      character(len=*), optional, intent(in   ) :: xdim
      character(len=*), optional, intent(in   ) :: ydim
      character(len=*), optional, intent(in   ) :: zdim
      character(len=*), optional, intent(in   ) :: tdim

      character(len=60) :: dimName, vname
      integer :: rc
      integer :: id
      integer :: i
      integer :: ncid
      integer :: dimId
      integer :: dimSize
      type(def), pointer :: gridInfo => null()
      type(var), pointer :: varInfo => null()

      integer :: nvars
      integer :: ndims
      integer :: ngatts
      integer :: nAtts
      character(len=60) :: attName
      integer :: unlimdimid
      integer :: ndims_in, nvars_in, ngatts_in, unlimdimid_in

      ! open netCDF file

      rc = nf90_open(fileName, nf90_nowrite, ncId)
      call pdie_(rc)

      self%ncFile = trim(fileName)
      self%ncid   = ncid

      rc = nf90_inquire(ncId, ndims, self%nvars, ngatts, unlimdimid)
      self%ndims = ndims

      ! get coordinates info
      allocate(self%gridRoot)
      gridInfo => self%gridRoot
      do dimId=1,ndims

         rc = nf90_inquire_dimension(ncId, dimId, dimName, dimSize)
         
         !
         !colocar warning se a dimensao nao exitir ou for desconhecida
         !
         gridInfo%name_  = Axis(getAxis(dimName))
         gridInfo%ncName_= dimName
         call getNcDim_(ncid, dimName, gridInfo%size_, gridInfo%ncId_, rc)
         !print*,gridInfo%name_,' ',gridInfo%ncName_, gridInfo%size_, gridInfo%ncId_, dimId
         if (rc .eq. nf90_noerr)then
            allocate(gridInfo%array_(gridInfo%size_))
            call getDimField_(ncId, dimName, gridInfo%name_, gridInfo%array_,rc)
         else
            call pwrn_(rc,dimName)
         endif

         allocate(gridInfo%next)
         gridInfo => gridInfo%next

      enddo


      ! get variable info

      allocate(self%varRoot)
      varInfo => self%varRoot

      do id = 1,self%nvars

         varInfo%id_ = id

         allocate(varInfo%dimIds(NF90_MAX_VAR_DIMS))
         rc = nf90_inquire_variable(ncId, id, &
                      name   = varInfo%name_, &
                      ndims  = varInfo%ndims, &
                      dimIds = varInfo%dimIds,&
                      nAtts  = nAtts)
         varInfo%scaleFactor = 1.0
         varInfo%addOffset   = 0.0

         do i=1,nAtts
            rc = nf90_inq_attName(ncId, id, i, attName)
            if(trim(attName) .eq. 'scale_factor')then
               rc = nf90_get_att(ncId, id, attName,varInfo%scaleFactor)
            endif
            if(trim(attName) .eq. 'add_offset')then
               rc = nf90_get_att(ncId, id, attName,varInfo%addOffset)
            endif
         enddo
         allocate(varInfo%dimSize(varInfo%ndims))
         do i=1,varInfo%ndims
            rc = nf90_inquire_dimension(ncId,&
                                        dimid = varInfo%dimIds(i),&
                                        len   = varInfo%dimSize(i)&
                                        )
         enddo

         if (id .ne. nvars ) then
            allocate(varInfo%next)
            varInfo => varInfo%next
         endif
      enddo
      self%varInfo => varInfo

   end subroutine

   function getDim_(self, dimName) result(dimVal)
      class(ncType),     intent(inout) :: self
      character(len=*),  intent(in   ) :: dimName
      integer                          :: dimVal

      type(def), pointer :: gInfo => null()
      integer :: rc

      if(associated(self%gridRoot)) gInfo => self%gridRoot
      do while(associated(gInfo))
         if (trim(gInfo%name_) .eq. trim(dimName))then
            dimVal = gInfo%size_
            return
         endif
         gInfo => gInfo%next
      enddo

      call getNcDim_(self%ncId, dimName, dimVal, rc)
      if (rc .ne. nf90_noerr)then
         dimVal = -1
      endif
      !if(present(istatus)) istatus = rc

   end function

   subroutine getDimVec_(self, dimName, Array, istatus)
      class(ncType),     intent(in   ) :: self
      character(len=*),  intent(in   ) :: dimName
      real, pointer,     intent(inout) :: Array(:)
      integer, optional, intent(  out) :: istatus

      type(def), pointer :: gInfo =>  null()

      if(present(istatus)) istatus = 0
      gInfo => self%gridRoot
      do while(associated(gInfo))
         if (trim(dimName) .eq. trim(gInfo%name_))then
            Array => gInfo%array_
            return
         endif
         gInfo => gInfo%next
      end do
      if(present(istatus)) istatus = -1
   end subroutine

   function getDimId_(self,dimName,varName) result(dimId)
      class(ncType),              intent(in   ) :: self
      character(len=*),           intent(in   ) :: dimName
      character(len=*), optional, intent(in   ) :: varName
      integer                         :: dimId

      type(def), pointer :: dimInfo => null()
      type(var), pointer :: varInfo => null()
      character(len=20) :: dNameNC
      integer :: i 

      if(present(varName))then
         varInfo => self%varRoot
         do while(associated(varInfo))
            if(trim(i90_lcase(varName)) == trim(i90_lcase(varInfo%name_)))then
               do i=1,varInfo%ndims
                  dNameNC = getNcDimName_(self%ncID,varInfo%dimIds(i))
                  if(trim(dimName) == trim(Axis(getAxis(dNameNC))))then
                     dimID = I
                     return
                  endif
               enddo
            endif
            varinfo =>  varInfo%next
         enddo      
      else
         dimInfo => self%gridRoot
         do while(associated(dimInfo))
            if(trim(dimName) == trim(dimInfo%name_))then
               dimId = dimInfo%ncId_
               return
            endif
            dimInfo => dimInfo%next
         enddo
      endif
   end function

   subroutine getField_(self, varName, lev, time, varArray, istatus)
      class(ncType),     intent(inout) :: self
      character(len=*),  intent(in   ) :: varName
      real,              intent(in   ) :: lev
      integer,           intent(in   ) :: time
      real,              intent(inout) :: varArray(:)
      integer, optional, intent(  out) :: istatus

      integer :: rc
      integer :: nz, z
      integer :: iret
      type(var), pointer :: varInfo => null()
      integer, allocatable :: istart(:)
      integer, allocatable :: icount(:)
      real,    pointer :: vec(:) => null()

      varInfo => self%varRoot
      do while(associated(varInfo))
         if (trim(i90_lcase(varName)) .eq. trim(i90_lcase(varInfo%name_))) exit
         varInfo => varInfo%next
      enddo
      if(.not.associated(varInfo))then
         print*,'Trying get variable ',trim(varName)
         print*,'Not Found!!!'
         if(present(istatus)) istatus = -99
         return
      endif

      allocate(icount(varInfo%ndims))
      allocate(istart(varInfo%ndims))
      istart = 1
      icount = varInfo%dimSize
      
      if(self%dimIsThere('zdef',varName))then
         call self%getDimVec('zdef',vec,iret)
         print*,vec
         z = minloc(lev-vec,mask=(lev-vec).ge.0,DIM=1)
         istart(self%getDimId('zdef',varName)) = Z
         icount(self%getDimId('zdef',varName)) = 1
      endif

      if(self%dimIsThere('tdef',varName))then
         istart(self%getDimId('tdef',varName)) = time
         icount(self%getDimId('tdef',varName)) = 1
      endif

      rc = nf90_get_var(self%ncId, varInfo%id_, varArray, &
                        start = istart, &
                        count = icount)

      varArray = varArray * varInfo%scaleFactor + varInfo%addOffset

      if(present(istatus))then
         istatus = rc
         call pwrn_(rc,varName)
         print*,'-----------'
         print*,istart
         print*,'-----------'
         print*,icount
      else
         call pwrn_(rc,varName)
      endif

   end subroutine


   subroutine getNcDim_(ncId, dimName, dimVal, dimId, istatus)
      integer,           intent(in   ) :: ncId
      character(len=*),  intent(in   ) :: dimName
      integer,           intent(  out) :: dimId
      integer,           intent(  out) :: dimVal
      integer, optional, intent(  out) :: istatus

      integer :: rc
      integer :: id
      
      rc = nf90_inq_dimid(ncid, dimName, dimId)
      if (rc .eq. nf90_noerr)then
         rc = nf90_inquire_dimension(ncid, dimId, len=dimVal )
      endif
      if(present(istatus)) istatus = rc
!      if (rc /= nf90_noerr)then
!         print*,'Trying get ',trim(dimName),' information'
!         print*,''
!         print*,'This netcdf API use by default COARDS convention, please use a xdf file to'
!         print*,'supplement or replace any internal metadata or verify if `',trim(dimName),'`'
!         print*,'is defined inside netcdf file.'
!         call pdie_(rc, dimName)
!      else
!         rc = nf90_inquire_dimension(ncid, id, len=dimVal )
!         call pwrn_(rc,dimName)
!      endif

   end subroutine

   function getNcDimName_(ncId, dimId)result(dimName)
      integer,           intent(in   ) :: ncId
      integer,           intent(in   ) :: dimId
      character(len=20)                :: dimName

      integer :: rc
      integer :: id
      
      rc = nf90_inquire_dimension(ncid, dimId, dimName )
      call pwrn_(rc,dimName)

   end function


   subroutine getNcField_(ncId, varName, varArray, istatus)
      integer,           intent(in   ) :: ncId
      character(len=*),  intent(in   ) :: varName
      real,              intent(inout) :: varArray(:)
      integer, optional, intent(  out) :: istatus

      integer :: rc
      integer :: id

      rc = nf90_inq_varid(ncId, varName, id)
      if (rc /= nf90_noerr)then
         print*,'Trying get ',trim(varName),' information'
         print*,''
         print*,'This netcdf API use by default COARDS convention, please use a xdf file to'
         print*,'supplement or replace any internal metadata or verify if `',trim(varName),'`'
         print*,'is defined inside netcdf file.'
         call pdie_(rc, varName)
      else
         rc = nf90_get_var(ncId, id, varArray)
         call pwrn_(rc,varName)
      endif

   end subroutine

   recursive subroutine getDimField_(ncId, varName, dName, varArray, istatus)
      integer,           intent(in   ) :: ncId
      character(len=*),  intent(in   ) :: varName
      character(len=*),  intent(in   ) :: dName
      real,              intent(inout) :: varArray(:)
      integer, optional, intent(  out) :: istatus

      integer :: rc
      integer :: id
      integer :: i

      rc = nf90_inq_varid(ncId, varName, id)
      if (rc == nf90_noerr)then
         rc = nf90_get_var(ncId, id, varArray)
         call pwrn_(rc,varName)
      else
         select case(trim(dName))
            case('xdef')
               call getDimField_(ncId,'longitude', 'None', varArray, istatus)
            case('ydef')
               call getDimField_(ncId,'latitude', 'None', varArray, istatus)
            case('zdef')
               call getDimField_(ncId,'pressure', 'None', varArray, istatus)
            case default
               
               do i=1,size(varArray)
                  varArray(i) = float(i)
               enddo

         end select
      endif

   end subroutine


   function getAxis(dName)result(Axis)
      character(len=*), intent(in   ) :: dName
      integer                         :: Axis

      ! Axis : 1 - xdef
      !        2 - ydef
      !        3 - zdef
      !        4 - tdef

      character(len=11),dimension(4) :: x, y, z, t
      logical :: found(4)
      integer :: i


      !        COADS       NCEP/GFS        ECCC        IAASARS   
      x = ['longitude  ','grid_xt    ','rlon       ','west_east  ']
      y = ['latitude   ','grid_yt    ','rlat       ','south_north']
      z = ['level      ','pfull      ','pres       ','bottom_top ']
      t = ['time       ','time       ','time       ','time       ']
      do i=1,size(x)
         found=[check(dName,x(i)),&
                check(dName,y(i)),&
                check(dName,z(i)),&
                check(dName,t(i))]
         Axis = findloc(found,.true.,dim=1)
         if(Axis > 0) exit
      enddo

      return
   end function

!   subroutine getNcVarInfo_(ncId)
!      integer, intent(in   ) :: ncId
!
!      call 
!   end subroutine

   subroutine pdie_(istatus, mess)
      integer, intent (in) :: istatus
      character(len=*), optional, intent(in) :: mess

      call pwrn_(istatus, mess)
      if (istatus /= nf90_noerr)then
         print*,'Mandatory information necessary'
         print*,'will stop ...'
         stop
      endif

   end subroutine pdie_    

    subroutine pwrn_(istatus, mess)
      integer, intent (in) :: istatus
      character(len=*), optional, intent(in) :: mess
      character(len=60) :: msc

      msc = ''
      if (present(mess)) msc = mess
      
      if (istatus /= nf90_noerr) then
         write(*,*)
         write(*,'(3A)') trim(adjustl(nf90_strerror(istatus))),' : ',trim(msc)
         write(*,*)
      end if

   end subroutine pwrn_ 

   function dimIsThere_(self,dName,vName)result(answer)
      class(ncType),              intent(inout) :: self
      character(len=*),           intent(in   ) :: dName
      character(len=*), optional, intent(in   ) :: vName
      logical                         :: answer
      
      type(def), pointer :: dimInfo => null()
      type(var), pointer :: varInfo => null()
      character(len=20) :: dNameNC
      integer :: i 

      answer = .false.
      if(present(vName))then
         varInfo => self%varRoot
         do while(associated(varInfo))
            if(trim(vName) == trim(varInfo%name_))then
               do i=1,varInfo%ndims
                  dNameNC = getNcDimName_(self%ncID,varInfo%dimIds(i))

                  if(trim(dName) == Axis(getAxis(dNameNC)))then
                     answer = .true.
                     return
                  endif

               enddo
            endif
            varinfo =>  varInfo%next
         enddo
      else
         dimInfo => self%gridRoot
         do while(associated(dimInfo))
            if(trim(dName) == trim(dimInfo%name_))then
               answer = .true.
               return
            endif
            dimInfo => dimInfo%next
         enddo
      endif

   end function

end module
