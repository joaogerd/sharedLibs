module fileAccess
   use m_GrADSFiles
   use accessGrib
   use accessNetCDF
   use m_inpak90

   implicit none
   private

   public :: ax
   
   type ax
      private
      character(len=512)    :: fileName

      class(*), allocatable :: data
      contains
      procedure, public :: open => open_
      procedure, public :: close => close_
      procedure, public :: getDim => getDim_
      procedure, public :: getDimVec => getDimVec_
      procedure         :: getField1D_, getField2D_
      generic,   public :: getField => getField1D_, getField2D_
   endtype


   contains

   function open_(self,fileName)result(iret)
      class(ax),        intent(inout) :: self
      character(len=*), intent(in   ) :: fileName
      integer                         :: iret

      character(len=4) :: ext
      integer :: i, j

      !--------------------------------------
      ! Verificando tipo de arquivo de entrada
      !--------------------------------------
  
      do i = len_trim(fileName), 1, -1
         if (fileName(i:i) .eq. '.') then
            j = i + 1
            exit
         endif
      end do
  
      ext = i90_lcase(fileName(j:len_trim(fileName)))

      select case(ext)
         case('ctl')
            allocate(GrADSFiles::self%data)
         case('grb','grib')
         print*,'GRib File'
            allocate(gribType::self%data)
         case('nc')
            allocate(ncType::self%data)
         case default
            write(*,*)'Wrong data file type! Abort ... ', trim(ext)
            write(*,*)'Should be ctl, grb (or grib), or nc file type!'
            iret = -3
            return
      end select

      select type (ptr => self%data)
         type is (GrADSFiles)
            iret = ptr%open(fileName)
         type is (gribType)
            call ptr%open(fileName, iret)
         type is (ncType)
            call ptr%open(fileName)
      end select

   end function

   function getDim_(self,dName)result(dSize)
      class(ax),        intent(inout) :: self
      character(len=*), intent(in   ) :: dName
      integer                         :: dSize

      select type(ptr => self%data)
         type is (GrADSFiles)
            dSize = ptr%getDim(dName)
         type is (gribType)
            dSize = ptr%getDim(dName)
         type is (ncType)
            dSize = ptr%getDim(dName)
      end select
   endfunction

   subroutine getDimVec_(self,dName, dVec)
      class(ax),        intent(inout) :: self
      character(len=*), intent(in   ) :: dName
      real,             intent(inout) :: dVec(:)

      real, pointer :: array(:) => null()
      integer :: npt
      
      select type(ptr => self%data)
         type is (GrADSFiles)
            call ptr%getDimVec(dName, array)
         type is (gribType)
            call ptr%getDimVec(dName, array)
         type is (ncType)
            call ptr%getDimVec(dName, array)
      end select

      dVec = array

   end subroutine


   function close_(self)result(iret)
      class(ax),        intent(inout) :: self
      integer                         :: iret
      
      deallocate(self%data, stat=iret)
   end function

   subroutine getField1D_(self,vName,vLev,fld, stat)
      class(ax),        intent(inout) :: self
      character(len=*), intent(in   ) :: vName
      real,             intent(in   ) :: vLev
      real,             intent(  out) :: fld(:)
      integer,          intent(  out) :: stat

      stat = 0

      select type(ptr => self%data)
         type is (GrADSFiles)
            stat = ptr%getField1D(vName, vLev, 1, fld)
         type is (gribType)
            call ptr%getField(vName,vLev, fld, stat)
         type is (ncType)
            call ptr%getField(vName, real(vLev,4), 1, fld, stat)
      end select
   end subroutine

   subroutine getField2D_(self,vName,vLev,fld, stat)
      class(ax),        intent(inout) :: self
      character(len=*), intent(in   ) :: vName
      real,             intent(in   ) :: vLev
      real,             intent(  out) :: fld(:,:)
      integer,          intent(  out) :: stat

      real, allocatable :: tmp(:)
      integer :: i, j, k
      integer :: x, y

      x = size(fld,1)
      y = size(fld,2)

      allocate(tmp(x*y))

      stat = 0

      call self%getField1D_(vName,vLev,tmp, stat)
      
      k=0
      do j=1,y
         do i=1,x
               fld(i,j) = tmp(i+k)
         enddo
         k=k+1
      enddo

   end subroutine

end module
