MODULE MiscMod

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: GetLongitudes
   PUBLIC :: GetImaxJmax
   PUBLIC :: GetGaussianLatitudes
   PUBLIC :: lterp

   !
   !precisao dos dados do modelo BAM
   integer, public, parameter :: I4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers
   integer, public, parameter :: I8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
   integer, public, parameter :: R4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
   integer, public, parameter :: R8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
   integer, public, parameter :: strlen = 1024

   !Logical Units 
   integer, parameter :: stderr = 0 ! Error Unit
   integer, parameter :: stdinp = 5 ! Input Unit
   integer, parameter :: stdout = 6 ! Output Unit


   ! Constants
   real(kind=r8), public, parameter :: rd     = 45_r8/ATAN(1.0_r8) ! convert to radian
   real(kind=r8), public, parameter :: EMRad   = 6.37E6_r8         ! Earth Mean Radius (m)
   real(kind=r8), public, parameter :: EMRad1  = 1.0_r8/EMRad      ! 1/EMRad (1/m)
   real(kind=r8), public, parameter :: EMRad12 = EMRad1*EMRad1     ! EMRad1**2 (1/m2)
   real(kind=r8), public, parameter :: EMRad2  = EMRad*EMRad       ! EMRad**2 (m2)

!
! Undef value
!

   integer, parameter :: iUdef = -9999
   real,    parameter :: Udef  = -9.99E9

   integer :: iMax
   integer :: jMax
!   integer :: kMax
!   integer :: Mend
   integer :: Mend1
   integer :: Mend2
   integer :: Mend3
   integer :: MnWv0
   integer :: MnWv1
   integer :: MnWv2
   integer :: MnWv3
   integer :: jMaxHf
   integer :: iMx


CONTAINS


SUBROUTINE InitParameters ( mend )

  IMPLICIT NONE

  integer, intent(in) :: mend

  
  CALL GetImaxJmax (Mend, Imax, Jmax)

  Mend1 = Mend+1
  Mend2 = Mend+2
  Mend3 = Mend+3
  Mnwv2 = Mend1*Mend2
  Mnwv0 = Mnwv2/2
  Mnwv3 = Mnwv2+2*Mend1
  Mnwv1 = Mnwv3/2

  Imx    = Imax+2
  Jmaxhf = Jmax/2

END SUBROUTINE InitParameters


SUBROUTINE GetImaxJmax (Mend, Imax, Jmax)

  IMPLICIT NONE

  INTEGER, INTENT (IN) :: Mend
  INTEGER, INTENT (OUT) :: Imax, Jmax

  INTEGER :: Nx, Nm, N2m, N3m, N5m, &
             n2, n3, n5, j, n, Check, Jfft

  INTEGER, SAVE :: Lfft=40000

  INTEGER, DIMENSION (:), ALLOCATABLE, SAVE :: Ifft

  N2m=CEILING(LOG(REAL(Lfft,r8))/LOG(2.0_r8))
  N3m=CEILING(LOG(REAL(Lfft,r8))/LOG(3.0_r8))
  N5m=CEILING(LOG(REAL(Lfft,r8))/LOG(5.0_r8))
  Nx=N2m*(N3m+1)*(N5m+1)

  ALLOCATE (Ifft (Nx))
  Ifft=0

  n=0
  DO n2=1,N2m
     Jfft=(2**n2)
     IF (Jfft > Lfft) EXIT
     DO n3=0,N3m
        Jfft=(2**n2)*(3**n3)
        IF (Jfft > Lfft) EXIT
        DO n5=0,N5m
           Jfft=(2**n2)*(3**n3)*(5**n5)
           IF (Jfft > Lfft) EXIT
           n=n+1
           Ifft(n)=Jfft
        END DO
     END DO
  END DO
  Nm=n

  n=0
  DO 
     Check=0
     n=n+1
     DO j=1,Nm-1
        IF (Ifft(j) > Ifft(j+1)) THEN
           Jfft=Ifft(j)
           Ifft(j)=Ifft(j+1)
           Ifft(j+1)=Jfft
           Check=1
        END IF
     END DO
     IF (Check == 0) EXIT
  END DO

  Imax=3*Mend+1
  DO n=1,Nm
     IF (Ifft(n) >= Imax) THEN
        Imax=Ifft(n)
        EXIT
     END IF
  END DO
  Jmax=Imax/2
  IF (MOD(Jmax, 2) /= 0) Jmax=Jmax+1

  DEALLOCATE (Ifft)

END SUBROUTINE GetImaxJmax

SUBROUTINE GetGaussianLatitudes (Jmax, Lat)

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: Jmax ! Number of Gaussian Latitudes

   REAL (KIND=r8), DIMENSION (Jmax), INTENT(OUT) :: &
                   Lat ! Gaussian Latitudes In Degree

   INTEGER :: j

   REAL (KIND=r8) :: eps, rd, dCoLatRadz, CoLatRad, dCoLatRad, p2, p1

   eps=1.0e-12_r8
   rd=45.0_r8/ATAN(1.0_r8)
   dCoLatRadz=((180.0_r8/REAL(Jmax,r8))/rd)/10.0_r8
   CoLatRad=0.0_r8
   DO j=1,Jmax/2
      dCoLatRad=dCoLatRadz
      DO WHILE (dCoLatRad > eps)
         CALL LegendrePolynomial (Jmax,CoLatRad,p2)
         DO
            p1=p2
            CoLatRad=CoLatRad+dCoLatRad
            CALL LegendrePolynomial (Jmax,CoLatRad,p2)
            IF (SIGN(1.0_r8,p1) /= SIGN(1.0_r8,p2)) EXIT
         END DO
         CoLatRad=CoLatRad-dCoLatRad
         dCoLatRad=dCoLatRad*0.25_r8
      END DO
      Lat(j)=90.0_r8-CoLatRad*rd
      Lat(Jmax-j+1)=-Lat(j)
      CoLatRad=CoLatRad+dCoLatRadz
   END DO

END SUBROUTINE GetGaussianLatitudes


SUBROUTINE LegendrePolynomial (N, CoLatRad, LegPol)

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: N ! Order of the Ordinary Legendre Function

   REAL (KIND=r8), INTENT(IN) :: CoLatRad ! Colatitude (In Radians)

   REAL (KIND=r8), INTENT(OUT) :: LegPol ! Value of The Ordinary Legendre Function

   INTEGER :: i

   REAL (KIND=r8) :: x, y1, y2, g, y3

   x=COS(CoLatRad)
   y1=1.0_r8
   y2=X
   DO i=2,N
      g=x*y2
      y3=g-y1+g-(g-y1)/REAL(i,r8)
      y1=y2
      y2=y3
   END DO
   LegPol=y3

END SUBROUTINE LegendrePolynomial

SUBROUTINE GetLongitudes (Imax, Lon0, Lon)

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: Imax ! Number of Longitudes

   REAL (KIND=r8), INTENT(In) :: Lon0 ! First Longitude In Degree

   REAL (KIND=r8), DIMENSION (Imax), INTENT(OUT) :: &
                   Lon ! Longitudes In Degree

   INTEGER :: i

   REAL (KIND=r8) :: dx

   dx=360.0_r8/REAL(Imax,r8)
   DO i=1,Imax
      Lon(i)=Lon0+REAL(i-1,r8)*dx
   END DO

END SUBROUTINE GetLongitudes


subroutine lterp(iField, iMend, oMend, oField)

   use coord_compute,   only: compute_grid_coord

   real(r8), intent(in   ) :: iField(:,:)
   integer,  intent(in   ) :: iMend
   integer,  intent(in   ) :: oMend
   real(r8), intent(inout) :: oField(:,:)

   ! Local variables

   real    :: y1, y2, y3, y4
   real    :: t, u, d
   integer :: iIMax, iJMax
   integer :: oIMax, oJMax
   integer :: i, i1, ii
   integer :: j, j1, jj
   integer :: iret
   
   real,              dimension(200)     :: gDesci             ! Grid description parameters of input field
   real,              dimension(200)     :: gDesco             ! Grid description parameters of output field
   real(r8), allocatable, dimension(:)   :: olat               ! latitudes in degrees of output field
   real(r8), allocatable, dimension(:)   :: olon               ! longitudes in degrees of output field
   real(r8), allocatable, dimension(:)   :: ilat               ! latitudes in degrees of input field
   real(r8), allocatable, dimension(:)   :: ilon               ! longitudes in degrees of input field

   real,     allocatable, dimension(:)   :: rlat               ! latitudes in degrees of input field
   real,     allocatable, dimension(:)   :: rlon               ! longitudes in degrees of input field

   real,     allocatable, dimension(:)   :: xpts
   real,     allocatable, dimension(:)   :: ypts
   real(r8),     allocatable, dimension(:,:) :: ya
   real(r8),     allocatable, dimension(:)   :: x1a
   real(r8),     allocatable, dimension(:)   :: x2a

!   real, parameter    :: udef  = 9.999E20

   !
   ! Get input info
   !
   call GetImaxJmax(iMend, iIMax, iJMax)
   if(size(iField,1).ne.iImax.or.size(iField,2).ne.iJMax)then
      write(stdout,*)'iField error:'
      write(stdout,*)'isize in  :',size(iField,1),size(iField,2)
      write(stdout,*)'isize jcap:',iIMax,iJMax
      stop
   endif

   allocate(ilon(iIMax))
   call GetLongitudes (iIMax, 0.0_r8, ilon)
   allocate(ilat(iJMax))
   call GetGaussianLatitudes(iJMax,ilat)
   gDesci    = 0
   gDesci(1) = 4
   gDesci(2) = iIMax
   gDesci(3) = iJMax
   gDesci(4) = maxval(ilat)
   gDesci(5) = minval(ilon)
   gDesci(6) = 128
   gDesci(7) = minval(ilat)
   gDesci(8) = maxval(ilon)
   gDesci(9) = abs(ilon(2)-ilon(1))
   gDesci(10)= iJMax/2
   gDesci(11)= 64
   gDesci(20)= 0

   DeAllocate(ilat)
   DeAllocate(ilon)

   !
   ! Get output info
   !
   call GetImaxJmax(oMend, oIMax, oJMax)
   if(size(oField,1).ne.oImax.or.size(oField,2).ne.oJMax)then
      write(stdout,*)'oField error:'
      write(stdout,*)'osize in  :',size(oField,1),size(oField,2)
      write(stdout,*)'osize jcap:',oIMax,oJMax
      stop
   endif

   allocate(olon(oIMax))
   call GetLongitudes (oIMax, 0.0_r8, olon)
   allocate(olat(oJMax))
   call GetGaussianLatitudes(oJMax,olat)

   Allocate(xpts(oIMax))
   Allocate(ypts(oJMax))

   call compute_grid_coord(gDesci, real(olon,4), real(olat,4), udef, xpts, ypts, iret)

!   write(stdout,*)'xpts:',minval(xpts),maxval(xpts)
   ypts = ypts + 1 ! <- to use with a ghost zone

   !
   ! Create a ghost zone at north and south poles
   !

   Allocate(ya(iIMax,iJMax+2))
   d               = 1/float(iIMax)
   ya(:,        1) = sum(iField(:,1))*d
   ya(:,2:iJMax+1) = iField
   ya(:,  iJMax+2) = sum(iField(:,iJMax))*d

   !
   ! Create input Grid 
   !

   Allocate(x1a(iIMax+1))
   x1a(1:iIMax+1) = (/(i,i=1,iIMax+1)/) !  <- Global field last point + 1 = 1st point
   Allocate(x2a(iJMax+2))
   x2a(1:iJMax+2) =(/(i,i=1,iJMax+2)/)

   do jj=1,oJMax
      do ii=1,oIMax

         if(xpts(ii).eq.udef .or. ypts(jj).eq. udef)then
            write(stdout,'(A1,x,2I6)')'wrong x,y points at:',ii,jj
            oField(ii,jj) = udef
            cycle
         endif

         i  = int(xpts(ii))
         i1 = mod(i,iIMax) + 1 ! <- Global field last point + 1 = 1st point
         j  = int(ypts(jj))
         j1 = j + 1


         y1 = ya( i, j)
         y2 = ya(i1, j)
         y3 = ya(i1,j1)
         y4 = ya( i,j1)

         t  = (real(xpts(ii),r8) - x1a(i))/(x1a(i+1)-x1a(i))
         u  = (real(ypts(jj),r8) - x2a(j))/(x2a(j+1)-x2a(j))

         oField(ii,jj) = (1-t)*(1-u)*y1 + t*(1-u)*y2 + t*u*y3 + (1-t)*u*y4

      enddo
   enddo

   DeAllocate(xpts)
   DeAllocate(ypts)
   DeAllocate(ya)
   DeAllocate(x1a)
   DeAllocate(x2a)

end subroutine

END MODULE MiscMod
