program cspline

use Precision
IMPLICIT                           NONE
integer(IntKi) ::                   i,j,k,numL,unitIN,unitOUT,ErrStat
real(ReKi) ::                      dummy,ClInterp
real(ReKi), allocatable ::         AoA(:), Cl(:), Cd(:),Coef(:,:)
character(len=36) ::               Fn
character(len=1) ::                dumchar
character(len=64) ::                    ErrMSG

! read the polar data
unitIN=11
unitOUT=33
numL=62
k=0
fn='testPolar.txt'
allocate(AoA(numL) &
        ,Cl(numL) &
        ,Cd(numL));

open(unit=unitIN, file=fn,status='old',action='read');
do i=1,numL
  if(i==1)then
    read(unitIN,*) dumchar
  else
    k=k+1
  read(unitIN,*) AoA(k), Cl(k),Cd(k)
  endif
enddo
close(unitIN)
!check input data
do i=1,numL-1
  print*,AoA(i),Cl(i),Cd(i)
enddo

!get coefs
allocate(Coef(numL-2,4))
call CubicSplineInit(numL-1,AoA,Cl,Coef,ErrStat,ErrMSG)
!check coefs
do j=1,numL-2
  print*,Coef(j,1),Coef(j,2),Coef(j,3),Coef(j,4)
enddo

!get interpolated value
do i=1,numL-1
   dummy=AoA(i)-0.3*i+0.01*sqrt(real(i,8))
   call CubicSplineInterps(dummy, numL-1, AoA, Cl, Coef,ClInterp, ErrStat, ErrMSG )
!   ClInterp=CubicSplineInterp(X, AryLen, XAry, YAry, Coef, ErrStat, ErrMsg )
   print*,'DummyAOA = ',dummy, 'Clinterp = ',ClInterp
enddo




deallocate(AoA,Cl,Cd)
deallocate(Coef)

end program

!=======================================================================
!> This routine calculates the parameters needed to compute a irregularly-spaced natural cubic spline.
!! Natural cubic splines are used in that the curvature at the end points is zero.
!! This routine does not require that the XAry be regularly spaced.
   SUBROUTINE CubicSplineInit ( AryLen, XAry, YAry, Coef, ErrStat, ErrMsg )
     use Precision
      ! Argument declarations:

   INTEGER, INTENT(IN)          :: AryLen                                     !< Length of the array

   REAL(ReKi), INTENT(OUT)      :: Coef  (AryLen-1,0:3)                       !< The coefficients for the cubic polynomials
   REAL(ReKi), INTENT(IN)       :: XAry  (AryLen)                             !< Input array of x values
   REAL(ReKi), INTENT(IN)       :: YAry  (AryLen)                             !< Input array of y values

   INTEGER(IntKi), INTENT(OUT)  :: ErrStat                                    !< Error status

   CHARACTER(*), INTENT(OUT)    :: ErrMsg                                     !< Error message


      ! Local declarations.

   REAL(ReKi), ALLOCATABLE      :: DelX  (:)                                  ! The distances between the randomly spaced points.
   REAL(ReKi), ALLOCATABLE      :: Slope (:)                                  ! The AryLen-1 length array of slopes between points.
   REAL(ReKi), ALLOCATABLE      :: U     (:)                                  ! An AryLen-1 length array used in the Gaussian elimination.
   REAL(ReKi), ALLOCATABLE      :: V     (:)                                  ! An AryLen-1 length array used in the Gaussian elimination.
   REAL(ReKi)                   :: ZHi                                        ! A parameter used to calculate the polynomial coefficients.
   REAL(ReKi)                   :: ZLo                                        ! A parameter used to calculate the polynomial coefficients.

   INTEGER(IntKi)               :: ErrStatLcL                                 ! Local error status.
   INTEGER                      :: I                                          ! The index into the arrays.
   CHARACTER(*), PARAMETER      :: RoutineName = 'CubicSplineInit'

   ErrStat = ErrID_None
   ErrMsg  = ""

      ! Allocate the various intermediate arrays.

   ALLOCATE ( DelX( AryLen - 1 ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, RoutineName//':Error allocating memory for the DelX array.' )
      RETURN
   ENDIF

   ALLOCATE ( Slope( AryLen - 1 ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, RoutineName//':Error allocating memory for the Slope array.' )
      RETURN
   ENDIF

   ALLOCATE ( U( AryLen - 1 ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, RoutineName//':Error allocating memory for the U array.' )
      RETURN
   ENDIF

   ALLOCATE ( V( AryLen - 1 ), STAT=ErrStatLcL )
   IF ( ErrStatLcL /= 0 )  THEN
      CALL ExitThisRoutine ( ErrID_Fatal, RoutineName//':Error allocating memory for the V array.' )
      RETURN
   ENDIF


      ! Compute the distance between XAry values and the slopes between points.

   DO I=1,AryLen-1
      DelX (I) =   XAry(I+1) - XAry(I)
      Slope(I) = ( YAry(I+1) - YAry(I) )/DelX(I)
   END DO ! I


      ! Use Gaussian elimination to solve the tri-diagonal matrix.

   U(1) = 2.0_ReKi*( DelX (2) + DelX (1) )
   V(1) = 6.0_ReKi*( Slope(2) - Slope(1) )

   DO I=2,AryLen-1
      U(I) = 2.0_ReKi*( DelX(I-1) + DelX(I)    ) - DelX(I-1)*DelX(I-1)/U(I-1)
      V(I) = 6.0_ReKi*( Slope(I)  - Slope(I-1) ) - DelX(I-1)*   V(I-1)/U(I-1)
   END DO ! I


      ! Determine the coefficients of the polynomials.

   Coef(:,0) = YAry(:)

   ZHi = 0.0_ReKi

   DO I=AryLen-1,1,-1
      ZLo       = ( V(I) - DelX(I)*ZHi )/U(I)
      Coef(I,1) = Slope(I) - DelX(I)*( ZHi/6.0_ReKi + ZLo/3.0_ReKi )
      Coef(I,2) = 0.5_ReKi*ZLo
      Coef(I,3) = ( ZHi - ZLo )/( 6.0_ReKi*DelX(I) )
      ZHi       = ZLo
   END DO ! I



   CALL ExitThisRoutine ( ErrID_None, 'No Problemo' )

   RETURN

   !=======================================================================
   CONTAINS
   !=======================================================================
      SUBROUTINE ExitThisRoutine ( ErrID, Msg )

         ! This subroutine cleans up the parent routine before exiting.


            ! Argument declarations.

         INTEGER(IntKi), INTENT(IN)       :: ErrID                            ! The error identifier (ErrLev)

         CHARACTER(*),   INTENT(IN)       :: Msg                              ! The error message (ErrMsg)


            ! Local declarations.

         LOGICAL                          :: IsOpen                           ! A flag that indicates if the input unit is still open.


            ! Set error status/message

         ErrStat = ErrID
         ErrMsg  = Msg


            ! Deallocate the Words array if it had been allocated.

         IF ( ALLOCATED( DelX  ) )  DEALLOCATE( DelX  )
         IF ( ALLOCATED( Slope ) )  DEALLOCATE( Slope )
         IF ( ALLOCATED( U     ) )  DEALLOCATE( U     )
         IF ( ALLOCATED( V     ) )  DEALLOCATE( V     )


         RETURN

      END SUBROUTINE ExitThisRoutine ! ( ErrID, Msg )

   END SUBROUTINE CubicSplineInit ! ( AryLen, XAry, YAry, YAry, Coef, ErrStat, ErrMsg )

!=======================================================================
!> This routine interpolates a pair of arrays using cubic splines to find the function value at X.
!! One must call cubicsplineinit first to compute the coefficients of the cubics.
!! This routine does not require that the XAry be regularly spaced.
   FUNCTION CubicSplineInterp ( X, AryLen, XAry, YAry, Coef, ErrStat, ErrMsg )
     use Precision

      ! Function declaration.

   REAL(ReKi)                   :: CubicSplineInterp                          !  This function


      ! Argument declarations:

   INTEGER, INTENT(IN)          :: AryLen                                     !< Length of the array

   REAL(ReKi), INTENT(IN)       :: Coef  (AryLen-1,0:3)                       !< The coefficients for the cubic polynomials
   REAL(ReKi), INTENT(IN)       :: X                                          !< The value we are trying to interpolate for
   REAL(ReKi), INTENT(IN)       :: XAry (AryLen)                              !< Input array of regularly spaced x values
   REAL(ReKi), INTENT(IN)       :: YAry (AryLen)                              !< Input array of y values

   INTEGER(IntKi), INTENT(OUT)  :: ErrStat                                    !< Error status

   CHARACTER(*), INTENT(OUT)    :: ErrMsg                                     !< Error message


      ! Local declarations.

   REAL(ReKi)                   :: XOff                                       ! The distance from X to XAry(ILo).

!   INTEGER(IntKi)               :: ErrStatLcL                                 ! Local error status.
   INTEGER                      :: ILo                                        ! The index into the array for which X is just above or equal to XAry(ILo).


   ErrStat = ErrID_None
   ErrMsg  = ""

      ! See if X is within the range of XAry.  Return the end point if it is not.

   IF ( X <= XAry(1) )  THEN
      CubicSplineInterp = YAry(1)
      RETURN
   ELSEIF ( X >= XAry(AryLen) )  THEN
      CubicSplineInterp = YAry(AryLen)
      RETURN
   ENDIF ! ( X <= XAry(1) )


      ! We are somewhere inside XAry.  Find the segment that bounds X using binary search.

   CALL LocateBin( X, XAry, ILo, AryLen )

   XOff = X - XAry(ILo)

   CubicSplineInterp = Coef(ILo,0) + XOff*( Coef(ILo,1) + XOff*( Coef(ILo,2) + XOff*Coef(ILo,3) ) )


   RETURN
   END FUNCTION CubicSplineInterp ! ( X, AryLen, XAry, YAry, Coef, ErrStat, ErrMsg )


!=======================================================================
!> This routine interpolates a pair of arrays using cubic splines to find the function value at X.
!! One must call cubicsplineinit first to compute the coefficients of the cubics.
!! This routine does not require that the XAry be regularly spaced.
   subroutine CubicSplineInterps ( X, AryLen, XAry, YAry, Coef, InterpdVal,ErrStat, ErrMsg )
     use Precision

      ! Function declaration.

      ! Argument declarations:

   INTEGER, INTENT(IN)          :: AryLen                                     !< Length of the array

   REAL(ReKi), INTENT(IN)       :: Coef  (AryLen-1,0:3)                       !< The coefficients for the cubic polynomials
   REAL(ReKi), INTENT(IN)       :: X                                          !< The value we are trying to interpolate for
   REAL(ReKi), INTENT(IN)       :: XAry (AryLen)                              !< Input array of regularly spaced x values
   REAL(ReKi), INTENT(IN)       :: YAry (AryLen)                              !< Input array of y values
   REAL(ReKi), intent(out)                  :: InterpdVal                         !  This function

   INTEGER(IntKi), INTENT(OUT)  :: ErrStat                                    !< Error status

   CHARACTER(*), INTENT(OUT)    :: ErrMsg                                     !< Error message


      ! Local declarations.

   REAL(ReKi)                   :: XOff                                       ! The distance from X to XAry(ILo).

!   INTEGER(IntKi)               :: ErrStatLcL                                 ! Local error status.
   INTEGER                      :: ILo                                        ! The index into the array for which X is just above or equal to XAry(ILo).


   ErrStat = ErrID_None
   ErrMsg  = ""

      ! See if X is within the range of XAry.  Return the end point if it is not.

   IF ( X <= XAry(1) )  THEN
      CubicSplineInterp = YAry(1)
      RETURN
   ELSEIF ( X >= XAry(AryLen) )  THEN
      CubicSplineInterp = YAry(AryLen)
      RETURN
   ENDIF ! ( X <= XAry(1) )


      ! We are somewhere inside XAry.  Find the segment that bounds X using binary search.

   CALL LocateBin( X, XAry, ILo, AryLen )

   XOff = X - XAry(ILo)

   InterpdVal = Coef(ILo,0) + XOff*( Coef(ILo,1) + XOff*( Coef(ILo,2) + XOff*Coef(ILo,3) ) )


   RETURN
   END subroutine CubicSplineInterps ! ( X, AryLen, XAry, YAry, Coef, ErrStat, ErrMsg )


!=======================================================================
!> This subroutine finds the lower-bound index of an input x-value located in an array.
!! On return, Ind has a value such that
!!           XAry(Ind) <= XVal < XAry(Ind+1), with the exceptions that
!!             Ind = 0 when XVal < XAry(1), and
!!          Ind = AryLen when XAry(AryLen) <= XVal.
!!
!! It uses a binary interpolation scheme that takes about log(AryLen)/log(2) steps to converge.
!! If the index doesn't change much between calls, LocateStp() (nwtc_num::locatestp) may be a better option.
   SUBROUTINE LocateBin( XVal, XAry, Ind, AryLen )
     use Precision

      ! Argument declarations.

   INTEGER, INTENT(IN)          :: AryLen                                          !< Length of the array.
   INTEGER, INTENT(OUT)         :: Ind                                             !< Final (low) index into the array.

   REAL(ReKi), INTENT(IN)       :: XAry    (AryLen)                                !< Array of X values to be interpolated.
   REAL(ReKi), INTENT(IN)       :: XVal                                            !< X value to be interpolated.


      ! Local declarations.

   INTEGER                      :: IHi                                             ! The high index into the arrays.
   INTEGER                      :: IMid                                            ! The mid-point index between IHi and Ind.



      ! Let's check the limits first.

   IF ( XVal < XAry(1) )  THEN
      Ind = 0
   ELSE IF ( XVal >= XAry(AryLen) )  THEN
      Ind = AryLen
   ELSE
         ! Let's interpolate!

      Ind  = 1
      IHi  = AryLen

      DO WHILE ( IHi-Ind > 1 )

         IMid = ( IHi + Ind )/2

         IF ( XVal >= XAry(IMid) ) THEN
            Ind = IMid
         ELSE
            IHi = IMid
         END IF

      END DO

   END IF

   RETURN
   END SUBROUTINE LocateBin
!=======================================================================


