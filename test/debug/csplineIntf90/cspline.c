  int LocateBin( double XVal, double XAry[], int AryLen)
  {
  /*
  INTEGER, INTENT(IN)          :: AryLen                                          !< Length of the array.
  INTEGER, INTENT(OUT)         :: Ind                                             !< Final (low) index into the array.

  REAL(ReKi), INTENT(IN)       :: XAry (AryLen)                                   !< Array of X values to be interpolated.
  REAL(ReKi), INTENT(IN)       :: XVal                                            !< X value to be interpolated.
  */

  /* Local declarations. */
  int  Ind;
  int  IHi;  // The high index into the arrays.
  int  IMid; // The mid-point index between IHi and Ind.

  // Let's check the limits first.

  if ( XVal < XAry[0] ){
     Ind = 0
  }
  else if(XVal>=XAry[AryLen]){
     Ind = AryLen
  }
  else{
     // Let's interpolate!
  Ind  = 0
  IHi  = AryLen
  do while(IHi-Ind>1){
     IMid = ( IHi + Ind )/2
     if( XVal >= XAry[IMid]){
        Ind = IMid
     }
     else{
        IHi = IMid
     }
  }
  }
  return Ind;
  }


  double CubicSplineInterps (double X, int AryLen, double XAry[], double YAry[], double Coef[][], double InterpdVal)
  {
   /*
   INTEGER, INTENT(IN)          :: AryLen                                     !< Length of the array

   REAL(ReKi), INTENT(IN)       :: Coef  (AryLen-1,0:3)                       !< The coefficients for the cubic polynomials
   REAL(ReKi), INTENT(IN)       :: X                                          !< The value we are trying to interpolate for
   REAL(ReKi), INTENT(IN)       :: XAry (AryLen)                              !< Input array of regularly spaced x values
   REAL(ReKi), INTENT(IN)       :: YAry (AryLen)                              !< Input array of y values
   REAL(ReKi), intent(out)                  :: InterpdVal                         !  This function

   INTEGER(IntKi), INTENT(OUT)  :: ErrStat                                    !< Error status

   CHARACTER(*), INTENT(OUT)    :: ErrMsg                                     !< Error message
   */

   //! Local declarations.
/*
   REAL(ReKi)                   :: XOff                                       ! The distance from X to XAry(ILo).

!   INTEGER(IntKi)               :: ErrStatLcL                                 ! Local error status.
   INTEGER                      :: ILo                                        ! The index into the array for which X is just above or equal to XAry(ILo).
*/
  double XOff;
  int    ErrStatLcL;
  int    ILo;

  int ErrStat = 0
  char *ErrMsg = ""

  // ! See if X is within the range of XAry.  Return the end point if it is not.

  if(X<=XAry[0]){
    InterpVal=YAry[0]
     return InterpVal;
  }
  else if(X>=XAry[AryLen-1]]{
     InterpVal=YAry[AryLen-1]
     return InterpVal;
  } // ( X <= XAry(1) )
  //! We are somewhere inside XAry.  Find the segment that bounds X using binary search.

  ILo=LocateBin( X, XAry,AryLen )

  XOff=X-XAry[ILo]

  InterpdVal=Coef[ILo][0] + XOff*(Coef[ILo][1] + XOff*(Coef[ILo][2] + XOff*Coef[ILo][3]));

  return InterpVal;
}



   double **CubicSplineInit(int AryLen, double XAry[], double YAry[])
   /*
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
   */
   //! Allocate the various intermediate arrays
   double **Coef=NULL; Coef=dmatrix(0,Arylen-1,0,4);
   double *DelX=NULL; DelX=dvector(0,AryLen-1);
   double *Slope=NULL; Slope=dvector(0,AryLen-1);
   double *U=NULL; U=dvector(0,AryLen-1);
   double *V=NULL; V=dvector(0,AryLen-1);
   double ZHi,ZLo;

   for(int i=0;i<AryLen-2;i++){
     DelX[i]=XAry[i+1]-XAry[i];
     Slope[i]=(YAry[i+1]-YAry[i])/DelX[i];
   }
   
   //! Use Gaussian elimination to solve the tri-diagonal matrix.

   U[0]=2.0*(DelX[1] + DelX[0]);
   V[0]=6.0*(Slope[1] - Slope[0]);
   
   //I=2,AryLen-1
   for(int i=1;i<AryLen-2;i++){
     U[i] = 2.0*(DelX[i-1]+ DelX[i]) - DelX[i-1]*DelX[i-1]/U[i-1];
     V[i] = 6.0*(Slope[i] - Slope[i-1]) - DelX[i-1]*V[i-1]/U[I-1];
   }
   //! Determine the coefficients of the polynomials.

   for(int i=0;i<AryLen;i++){
     Coef[i][0]=YAry[i];
     ZHi=0.0;
   }

//   DO I=AryLen-1,1,-1
   for(int i=AryLen-1;i>=0;i--){
      ZLo       = ( V[i] - DelX[i]*ZHi )/U[i];
      Coef[i][1]= Slope[i]-DelX[i]*(ZHi/6.0+ZLo/3.0);
      Coef[i][2]= 0.5*ZLo;
      Coef[I][3]= ( ZHi - ZLo )/( 6.0*DelX[I]);
      ZHi=ZLo
   }

/*
   CALL ExitThisRoutine ( ErrID_None, 'No Problemo' )
*/

  return Ceof;
}


