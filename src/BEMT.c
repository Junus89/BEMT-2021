#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>

#include "BEMT.h"

/* utility functions */
  void BEMT_error(char error_text[],\
                  const char* routine)
  /*-------------------------------------------------------------------
  Purpose:
          reports error for the BEMT Solver implementation.

  Written by:
             Yunusi Fuerkaiti
             email: y.fuerkaiti@tudelft.nl
  date:
       12.06.2020
  -------------------------------------------------------------------*/
  {
    fprintf(stderr,"                                    \n");
    fprintf(stderr,"BEMT run-time error in routine:\t%s\n",\
      routine);
    fprintf(stderr,"%s\n",error_text);
    fprintf(stderr,"Now exiting to system ...\n");
    exit(1);
  }

  void BEMT_warning(char error_text[],\
                    const char* routine)
  /*-------------------------------------------------------------------
  Purpose:
          reports warning for the BEMT Solver implementation.

  Written by:
             Yunusi Fuerkaiti
             email: y.fuerkaiti@tudelft.nl
  date:
       12.06.2020
  -------------------------------------------------------------------*/
  {
    fprintf(stderr,"                                    \n");
    fprintf(stderr,"BEMT run-time warning in routine:\t%s\n",\
      routine);
    fprintf(stderr,"%s\n",error_text);
    fprintf(stderr,"Please check your routine and run again! \n");
    //exit(1);
  }

/* memory unitilty functions */
  double *dvector(long nl,long nh)
  /*-------------------------------------------------------------------
  Purpose:
          allocates an array of double type.

  Written by:
             Yunusi Fuerkaiti
             email: y.fuerkaiti@tudelft.nl
  date:
       12.06.2020
  -------------------------------------------------------------------*/
  {
    double *outputAr=NULL;
    const char* thisroutine="double *dvector(...)";
    outputAr = (double *)malloc((nh+1)*sizeof(double));
    if(outputAr == NULL)
    {
      BEMT_error("No enough memory space.",thisroutine);
    }
    return outputAr;
  }


  void free_dvector(double *inputAr,\
                    long nl,\
                    long nh)
  {
      free(inputAr);
      return;
  }

  /* complex dvector */
  double complex *dcvector(long nl,long nh)
  /*-------------------------------------------------------------------
  Purpose:
          allocates an array of double complex type.

  Written by:
             Yunusi Fuerkaiti
             email: y.fuerkaiti@tudelft.nl
  date:
       12.06.2020
  -------------------------------------------------------------------*/
  {
    const char* thisroutine="complex *dvector(...)";

    double complex *outputAr=NULL;

    outputAr = (double complex*)malloc((nh+1)*sizeof(double));
    if (outputAr == NULL)
    {
      BEMT_error("No enough memory space.",thisroutine);
    }
    return outputAr;
  }


  void free_dcvector(double complex *inputAr,\
                     long nl,\
                     long nh)
  {
    free(inputAr);
    return;
  }
  /*-----------------------------2D----------------------------------*/

  double **dmatrix(long nrl,\
                   long nrh,\
                   long ncl,\
                   long nch)
  /*-------------------------------------------------------------------
  Purpose:
          allocates a matrix of double type.

  Written by:
             Yunusi Fuerkaiti
             email: y.fuerkaiti@tudelft.nl
  date:
       12.06.2020
  -------------------------------------------------------------------*/
  {
    const char* thisroutine="double **dmatrix(...)";

    int i;
    double **outputAr=NULL;
    outputAr = (double **)malloc((nrh+1)*sizeof(double *));
    if (outputAr == NULL)
    {
      BEMT_error("No enough memory space.",thisroutine);
    }

    for (i = 0; i < nrh+1; i++){
      outputAr[i] = (double *)malloc((nch+1)*sizeof(double));
      if (outputAr[i] == NULL){
      BEMT_error("No enough memory space.",thisroutine);
      }
    }
    return outputAr;
  }
  /* */
  void free_dmatrix(double **inputAr,\
                    long nrl,\
                    long nrh,\
                    long ncl,\
                    long nch)
  /*      
   */
  {
      int i;
      for (i = 0; i <= nrh; i++) free(inputAr[i]);
      free(inputAr);              
      return;
  }

/* */
/* complex dmatrix */
/* complex dmatrix */
  double complex **dcmatrix(long nrl,\
                            long nrh,\
                            long ncl,\
                            long nch)
  /*-------------------------------------------------------------------
  Purpose:
          allocates a matrix of double type.

  Written by:
             Yunusi Fuerkaiti
             email: y.fuerkaiti@tudelft.nl
  date:
       12.06.2020
  -------------------------------------------------------------------*/
  {
    const char* thisroutine="double complex **dcmatrix(...)";
    int i;
    double complex **outputAr=NULL;
    outputAr = (double complex**)malloc((nrh+1)*sizeof(double complex*));
    if (outputAr == NULL)
    {
    BEMT_error("No enough memory space.",thisroutine);
    }

    for (i = 0; i < nrh+1; i++){
      outputAr[i] = (double complex*)malloc((nch+1)*sizeof(double complex));
      if (outputAr[i] == NULL){
        BEMT_error("No enough memory space.",thisroutine);
      }
    }
    return outputAr;
  }
/* */

  /* */
  void free_dcmatrix(double complex **inputAr,\
                     long nrl,\
                     long nrh,\
                     long ncl,\
                     long nch)
  /*      
   */
  {
      int i;
      for (i = 0; i <= nrh; i++)  free(inputAr[i]);
      free(inputAr);              
      return;
  }

  /* */
  d2_t *d2_tvector(long nl,long nh)
  /*-------------------------------------------------------------------
  Purpose:
          allocate a vector of d2_t type.

  Written by:
             Yunusi Fuerkaiti
             email: y.fuerkaiti@tudelft.nl
  date:
       12.06.2020
  -------------------------------------------------------------------*/
  {
    const char* thisroutine="double complex ***dc3dtensor(...)";
    d2_t *outputAr=NULL;

    outputAr=(d2_t *)malloc((nh+1)*sizeof(d2_t));
    if (outputAr==NULL){
      BEMT_error("No enough memory space.",thisroutine);
    }
    return outputAr;
  }


  void free_d2_tvector(d2_t *inputAr,\
                       long nl,\
                       long nh)
  {
      free(inputAr);
      return;
  }

  d2_t **d2_tmatrix(long nrl,\
                    long nrh,\
                    long ncl,\
                    long nch)
    /*-------------------------------------------------------------------
    Purpose:
            allocate a matrix of d2_t type.

    Written by:
               Yunusi Fuerkaiti
               email: y.fuerkaiti@tudelft.nl
    date:
         12.06.2020
    -------------------------------------------------------------------*/
  {
    int i;
    const char* thisroutine="d2_t **d2_tmatrix(...)";
    d2_t **outputAr=NULL;
    outputAr=(d2_t **)malloc((nrh+1)*sizeof(d2_t *));
    if (outputAr==NULL)
    {
      BEMT_error("No enough memory space.",thisroutine);
    }

    for (i=0;i<nrh+1;i++){
      outputAr[i]=(d2_t *)malloc((nch+1)*sizeof(d2_t));
      if (outputAr[i] == NULL){
        BEMT_error("No enough memory space.",thisroutine);
      }
    }
    return outputAr;
  }

  /* */
  void free_d2_tmatrix(d2_t **inputAr,\
                       long nrl,\
                       long nrh,\
                       long ncl,\
                       long nch)
  /*      
   */
  {
    int i;
    for (i=0;i<=nrh;i++)  free(inputAr[i]);
    free(inputAr);
    return;
  }


  double ***d3dtensor(long nrl,
                      long nrh,\
                      long ncl,\
                      long nch,\
                      long ndl,\
                      long ndh)
  /*-------------------------------------------------------------------
  Purpose:
          allocates a tensor(3D vector) of double type.

  Written by:
             Yunusi Fuerkaiti
             email: y.fuerkaiti@tudelft.nl
  date:
       12.06.2020
  -------------------------------------------------------------------*/
  {
    const char* thisroutine="double ***d3dtensor(...)";
    int i,j;
    double ***outputAr=NULL;
    outputAr = (double ***)malloc((nrh+1)*sizeof(double **));
    if (outputAr == NULL)
    {
        BEMT_error("No enough memory space.",thisroutine);
        exit(1);
    }

    for (i = 0; i < nrh+1; i++){
      outputAr[i] = (double **)malloc((nch+1)*sizeof(double*));
      if (outputAr[i] == NULL)
              {
      BEMT_error("No enough memory space.",thisroutine);
              }
      }

    for (i = 0; i<nrh+1;i++){
      for (j = 0;j<nch+1;j++){
        outputAr[i][j]=(double *)malloc((ndh+1)*sizeof(double));
        if (outputAr[i][j]==NULL){
          BEMT_error("No enough memory space.",thisroutine);
        }
      }
    }
    return outputAr;
  }

  void free_d3dtensor(double ***inputAr,
                      long nrl,\
                      long nrh,\
                      long ncl,\
                      long nch,\
                      long ndl,\
                      long ndh)
  /*      
   */
  {
    int i,j;

    for (i = 0; i < nrh+1; i++) {
      for (j = 0; j < nch+1; j++){
        free(inputAr[i][j]);
      }
    }
    for (i = 0; i < nrh+1; i++) free(inputAr[i]);
    free(inputAr);
    return;
  }






/*File reader functions */
/*    */
  d3_t *freader3col(char* fn,int* n)
  /*-------------------------------------------------------------------
  Purpose:
             reads a data file with 3 columns and returns the data
             and number of rows. 
  Written by:
             Yunusi Fuerkaiti
             email: y.fuerkaiti@tudelft.nl
  date:
       12.06.2020
  -------------------------------------------------------------------*/
  {
  const char* thisroutine="d3_t *freader3col(...)";

  d3_t *data=NULL;
  d3_t data_tmp;
  int index,ii=0;
  const int ncol=3;
  for(int i=0;i<ncol;i++)data_tmp.x[i]=0.0;

  /* */
  FILE *file;
  if((file=fopen(fn,"r"))==NULL){
        BEMT_error("Error on opening input file.\n", thisroutine);
  }
  file=fopen(fn,"r");

  while(fscanf(file,"%lf %lf %lf",&data_tmp.x[0], \
        &data_tmp.x[1], &data_tmp.x[2])!=EOF){
    if(data==NULL){
      data=malloc(sizeof(d3_t));
      *data=data_tmp;
    }
    if(data!=NULL){
      ii++;
      data=realloc(data,sizeof(d3_t)*ii);
      index=ii-1;
      *(data+index)=data_tmp;
    }
  }
  *n=ii;
  fclose(file);
  return data;
  }
/* */

  d4_t *freader4col(char* fn,int* n)
  /*-------------------------------------------------------------------
  Purpose:
             reads a data file with 6 columns and returns the data
             and number of rows. 
  Written by:
             Yunusi Fuerkaiti
             email: y.fuerkaiti@tudelft.nl
  Date:
       12.06.2020
  Modification: 
                It reads and discards the first two line that contains 
  names of each column.
  Date of modification:
  26.02.2021
  -------------------------------------------------------------------*/
  {
  const char* thisroutine="freader4col";

  d4_t *data=NULL;
  d4_t data_tmp;
  int index,ii=0;
  const int ncol=4;
  for(int i=0;i<ncol;i++) data_tmp.x[i]=0.0;
  /* */
  /* */
  FILE *file;
  if((file=fopen(fn,"r"))==NULL){
        BEMT_error("Error on opening input file.\n", thisroutine);
  }
  file=fopen(fn,"r");
  /* skip the first line*/
  fscanf(file,"%*[^\n]");//read and discard the first line
  fscanf(file,"%*[^\n]");//read and discard the second line
  while(fscanf(file,"%lf %lf %lf %lf",&data_tmp.x[0], \
        &data_tmp.x[1], &data_tmp.x[2],&data_tmp.x[3])!=EOF){
    if(data==NULL){
      data=malloc(sizeof(d4_t));
      *data=data_tmp;
    }
    if(data!=NULL){
      ii++;
      data=realloc(data,sizeof(d4_t)*ii);
      index=ii-1;
      *(data+index)=data_tmp;
    }
  }
  *n=ii;
  fclose(file);
  return data;
  }


  d5_t *freader5col(char* fn,int* n)
  /*-------------------------------------------------------------------
  Purpose:
             reads a data file with 6 columns and returns the data
             and number of rows. 
  Written by:
             Yunusi Fuerkaiti
             email: y.fuerkaiti@tudelft.nl
  Date:
       12.06.2020
  Modification: 
                It reads and discards the first line that contains 
  names of each column.
  Date of modification:
  26.02.2021
  -------------------------------------------------------------------*/
  {
  const char* thisroutine="d5_t *freader5col(...)";

  d5_t *data=NULL;
  d5_t data_tmp;
  int index,ii=0;
  const int ncol=5;
  for(int i=0;i<ncol;i++) data_tmp.x[i]=0.0;
  /* */
  /* */
  FILE *file;
  if((file=fopen(fn,"r"))==NULL){
        BEMT_error("Error on opening input file.\n", thisroutine);
  }
  file=fopen(fn,"r");
  /* skip the first line*/
  fscanf(file,"%*[^\n]");//read and discard the first line
  while(fscanf(file,"%lf %lf %lf %lf %lf",&data_tmp.x[0], \
        &data_tmp.x[1], &data_tmp.x[2],&data_tmp.x[3], \
        &data_tmp.x[4])!=EOF){
    if(data==NULL){
      data=malloc(sizeof(d5_t));
      *data=data_tmp;
    }
    if(data!=NULL){
      ii++;
      data=realloc(data,sizeof(d5_t)*ii);
      index=ii-1;
      *(data+index)=data_tmp;
    }
  }
  *n=ii;
  fclose(file);
  return data;
  }



  d4_t *polarDataReader(char* fn,int* n)
  /*-------------------------------------------------------------------
  Purpose:
             reads a polar data file with 4 columns and returns the data
             and number of rows. 
  Written by:
             Yunusi Fuerkaiti
             email: y.fuerkaiti@tudelft.nl
  date:
       21.02.2021
  -------------------------------------------------------------------*/
  {
  const char* thisroutine="d4_t *polarDataReader(char* fn,int* n)";

  d4_t *data=NULL;
  d4_t data_tmp;
  int index,ii=0;
  const int ncol=4;
  for(int i=0;i<ncol;i++)data_tmp.x[i]=0.0;

  /* */
  FILE *file;
  if((file=fopen(fn,"r"))==NULL){
        BEMT_error("Error on opening input file.\n", thisroutine);
  }
  file=fopen(fn,"r");
  /* skipping the first two columns that start with "#" */
  fscanf(file,"%*[^\n]");
//  fscanf(file,"%*[^\n]");
  while(fscanf(file,"%lf %lf %lf %lf",&data_tmp.x[0], \
        &data_tmp.x[1], &data_tmp.x[2],&data_tmp.x[3])!=EOF){
    if(data==NULL){
      data=malloc(sizeof(d4_t));
      *data=data_tmp;
    }
    if(data!=NULL){
      ii++;
      data=realloc(data,sizeof(d4_t)*ii);
      index=ii-1;
      *(data+index)=data_tmp;
    }
  }
  *n=ii;
  fclose(file);
  return data;
  }
/* */
  /*

  */
  Airfoil_t *AirfoilDBreader(char* fn,int* n)
  /*-------------------------------------------------------------------
  Purpose:
             reads airfoils data long a blade section. The first 4 
             columns are: STATION/R  CHORD/R  TWIST (DEG)  SWEEP/R 
             and the last column is airfoil PROFILE file name.
  Written by:
             Yunusi Fuerkaiti
             email: y.fuerkaiti@tudelft.nl
  date:
       19.12.2020
  -------------------------------------------------------------------*/
  {
  const char* thisroutine="Airfoil_t *AirFoilDBreader(...)";

  Airfoil_t *data=NULL;
  Airfoil_t data_tmp;
  int index,ii=0;
  const int ncol=5;
  for(int i=0;i<ncol-1;i++){
    data_tmp.x[i]=0.0;
    // data_tmp.FoilFN=" ";
  }
  /* */
  /* */
  FILE *file;
  if((file=fopen(fn,"r"))==NULL){
        BEMT_error("Error on opening input file.\n", thisroutine);
  }
  file=fopen(fn,"r");

  while(fscanf(file,"%lf %lf %lf %lf %s",&data_tmp.x[0], \
        &data_tmp.x[1],&data_tmp.x[2],&data_tmp.x[3], \
        data_tmp.FoilFN)!=EOF){
    if(data==NULL){
      data=malloc(sizeof(Airfoil_t));
      *data=data_tmp;
    }
    if(data!=NULL){
      ii++;
      data=realloc(data,sizeof(Airfoil_t)*ii);
      index=ii-1;
      *(data+index)=data_tmp;
    }
  }
  *n=ii;
  fclose(file);
  return data;
  }

/*--------------------------------------------------------------------*/
/*                                                                    */
/*                         Math utility functions                     */
/*                                                                    */
/*--------------------------------------------------------------------*/
/*--------------- Barycentric interpolation functions ----------------*/
  d2_t blii1d(d2_t x,\
              d2_t y,
              double xi)
/***********************************************************************
! Barycentric Linear Interpolation 1D
! Written by YF
! returns the interpolated value and its first derivatives:yi and yxi;
******************************************************************/
{
  double x1,x2;
  d2_t yi;
  double y1,y2;

  x1=x.x[0];
  x2=x.x[1];
  
  y1=y.x[0];
  y2=y.x[1];

  yi.x[1]=(y2-y1)/(x2-x1);
  yi.x[0]=y1+yi.x[1]*(xi-x1);
  
  return yi;
}

  int brcket(int n,\
             double x[],\
             double xi)
/* **********************************************************************
c      BRaCKET subroutine 
c      Written by YF
********************************************************************** */
{
  int i,ia,im,ib;
  
  ia = 0;
  ib = n-2;
  printf("%lf %lf\n",xi,x[n-1]);

/*If xi is outside the interval [x(1) x(n)] just output zero and
      return to the calling program:*/
  if((xi<x[0])||(xi>x[n-1])){
    i=0;
    printf("outside--brkt");
  return i;
  }
/*If xi is inside the interval [x(1) x(n)] just divide it by half
 until ib - ia = 1:*/
  im=(ia+ib)/2;
  if((ib-ia)>1){
    if((xi>x[ia]) && (xi<x[im])){
      ib=im;}
    else{
      ia=im;
    }
    im=(ia+ib)/2;
  }
  i=im;
  return i;
}

    double bcui1d(d4_t x,\
                d4_t f,\
                double xi)
/*****************************************************************
! Barycentric Cubic Interpolation 1D 
! Written by Yunusi Fuerkaiti
! returns the interpolated value for given xi and its first and
! second derivates to x[0],x[1],x[2], respectively.
! 
!*****************************************************************/
{
  int i=0;
  d3_t a,px,sx,qx;
  d3_t fi;
  double x1,x2,x3,x4;
  x1=x.x[0];
  x2=x.x[1];
  x3=x.x[2];
  x4=x.x[3];

  px.x[0]= ( x2 - x1 )*( x2 - x3 )*( x2 - x4 );
  px.x[1]= ( x3 - x1 )*( x3 - x2 )*( x3 - x4 );
  px.x[2]= ( x4 - x1 )*( x4 - x2 )*( x4 - x3 );
  
  for(i=0;i<3;i++) a.x[i]=(f.x[i+1]-f.x[0])/px.x[i];
  px.x[0]=(xi-x1)*( xi - x3 )*( xi - x4 );
  px.x[1]=(xi-x1)*( xi - x2 )*( xi - x4 );
  px.x[2]=(xi-x1)*( xi - x2 )*( xi - x3 );
  
  qx.x[0]=(xi-x1)*(xi-x3)+(xi-x1)*(xi-x4)+(xi-x3)*(xi-x4);
  qx.x[1]=(xi-x1)*(xi-x2)+(xi-x1)*(xi-x4)+(xi-x2)*(xi-x4);
  qx.x[2]=(xi-x1)*(xi-x2)+(xi-x1)*(xi-x3)+(xi-x2)*(xi-x3);
  
  sx.x[0]=2.0*( 3.0*xi - x1 - x3 - x4 );
  sx.x[1]=2.0*( 3.0*xi - x1 - x2 - x4 );
  sx.x[2]=2.0*( 3.0*xi - x1 - x2 - x3 );
  
  fi.x[0]=f.x[0];//interpolated value
  fi.x[1]=0.0;//first derivative
  fi.x[2]=0.0;//second derivative
  
  for(i=0;i<3;i++){
    fi.x[0]=fi.x[0]+a.x[i]*px.x[i];
    fi.x[1]=fi.x[1]+a.x[i]*qx.x[i];
    fi.x[2]=fi.x[2]+a.x[i]*sx.x[i];}
  return fi.x[0];
}

  double interp1D(int nx,\
                double tabx[],\
                double tabc[],\
                double xi)
  {
  /*
  */
  double ci;
  d2_t c2;
  d4_t xc,fc;
  d2_t xl,fl;
  double cxxxi=0;
  int i;
  if(xi<tabx[0]){
    xl.x[0]=tabx[0];
    xl.x[1]=tabx[1];
    fl.x[0]=tabc[0];
    fl.x[1]=tabc[1];
    c2=blii1d(xl,fl,xi);
    ci=c2.x[0];
    if(DEBUG==0) printf("test1 %lf %lf\n",xi,tabx[0]);
  }
  else if(xi>tabx[nx-2]){
    xl.x[0]=tabx[nx-2];
    xl.x[1]=tabx[nx];
    fl.x[0]=tabc[nx-2];
    fl.x[1]=tabc[nx];
    c2=blii1d(xl,fl,xi);
    ci=c2.x[0];
    if(DEBUG==0) printf("test2 %lf %lf \n",xi,tabx[nx-2]);

  }
  else{
    i=brcket(nx,tabx,xi);
    xc.x[0]=tabx[i-1];
    xc.x[1]=tabx[i];
    xc.x[2]=tabx[i+1];
    xc.x[3]=tabx[i+2];
    
    fc.x[0]=tabc[i-1];
    fc.x[1]=tabc[i];
    fc.x[2]=tabc[i+1];
    fc.x[3]=tabc[i+2];
    ci=bcui1d(xc,fc,xi);
    if(DEBUG==0) printf("test3\n");

    }
    return ci;
}

  void polint(double xa[],\
              double ya[],\
              int n,\
              double x,\
              double *y)
  {
  char *thisroutine="void polint(...)";

  double *dy=NULL;
  double TINY=1.0e-25;
  int i,m,ns=0;
  double den,dif,dift,ho,hp,w;
  double *c=NULL,*d=NULL;

  dif=fabs(x-xa[0]);
  c=dvector(0,n);
  d=dvector(0,n);
  for (i=0;i<n;i++){
          if ( (dift=fabs(x-xa[i])) < dif) {
                  ns=i;
                  dif=dift;
          }
          c[i]=ya[i];
          d[i]=ya[i];
  }
  *y=ya[ns--];
  for (m=1;m<n;m++) {
          for (i=0;i<=n-m;i++) {
                  ho=xa[i]-x;
                  hp=xa[i+m]-x;
                  w=c[i+1]-d[i];
                  if((den=ho-hp) == 0.0) den=den+TINY;//BEMT_warning("Error in routine %s\n",thisroutine);
                  den=w/den;
                  d[i]=hp*den;
                  c[i]=ho*den;
          }
          // *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
          *y += (2*(ns+1) < (n-m) ? c[ns+1] : d[ns--]);

  }
  free_dvector(d,0,n);
  free_dvector(c,0,n);
  }


  void ratint(double xa[],\
              double ya[],\
              int n,\
              double x,\
              double *y)
{
  char *thisroutine="void ratint(...)";
  double *dy=NULL;
  double TINY=1.0e-25;
  int m,i,ns=0;
  double w,t,hh,h,dd,*c=NULL,*d=NULL;

  c=dvector(0,n);
  d=dvector(0,n);
  hh=fabs(x-xa[0]);
  for (i=0;i<n;i++) {
    h=fabs(x-xa[i]);
    if (h == 0.0) {
            *y=ya[i];
            *dy=0.0;
            free_dvector(d,0,n);free_dvector(c,0,n);
    } else if (h < hh) {
            ns=i;
            hh=h;
    }
    c[i]=ya[i];
    d[i]=ya[i]+TINY;
  }
  *y=ya[ns--];
  for (m=1;m<n;m++) {
          for (i=0;i<=n-m;i++) {
                  w=c[i+1]-d[i];
                  h=xa[i+m]-x;
                  t=(xa[i]-x)*d[i]/h;
                  dd=t-c[i+1];
                  if (dd == 0.0) dd=dd+TINY;// APnoise_warning("Error in routine %s",thisroutine);
                  dd=w/dd;
                  d[i]=c[i+1]*dd;
                  c[i]=t*dd;
          }
          *y += (2*(ns+1) < (n-m) ? c[ns+1] : d[ns--]);
  if(DEBUG==0) printf("y=%lf\n",*y);
  }
  if(x==xa[n-1]) *y=ya[n-1];
  free_dvector(d,0,n);free_dvector(c,0,n);
  }

  void spline(double x[], double y[], int n, double yp1, double ypn, double y2[])
  {
          int i,k;
          double p,qn,sig,un,*u;
  
          u=dvector(0,n-1);
          if (yp1 > 0.99e30)
                  y2[0]=u[0]=0.0;
          else {
                  y2[0] = -0.5;
                  u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
          }
          for (i=1;i<n-1;i++) {
                  sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
                  p=sig*y2[i-1]+2.0;
                  y2[i]=(sig-1.0)/p;
                  u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
                  u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
          }
          if (ypn > 0.99e30)
                  qn=un=0.0;
          else {
                  qn=0.5;
                  un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
          }
          y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
          for (k=n-2;k>=0;k--){
                  y2[k]=y2[k]*y2[k+1]+u[k];}
          free_dvector(u,0,n-1);
  }


  void splint(double xa[], double ya[], double y2a[], int n, double x, double *yi)
  {
    const char* thisroutine="splint";
    int klo,khi,k;
    double h,b,a;

    klo=0;
    khi=n;
    while (khi-klo > 1) {
            k=(khi+klo) >> 1;
            if (xa[k] > x) khi=k;
            else klo=k;
    }
    h=xa[khi]-xa[klo];
    if (h == 0.0) BEMT_error("Bad xa input to routine %s\n",thisroutine);
    a=(xa[khi]-x)/h;
    b=(x-xa[klo])/h;
    *yi=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
  }

  void cspline(double x[], double y[], int n, double yp1, double ypn, double xi, double *yi)
  {
  /*
  */
  const char* thisroutine="cspline";
  int code=0;
  if(DEBUG==0){
    printf(" yp1  = %lf , ypN = %lf\n",yp1, ypn);
  }
  double y2[n];
  spline(x,y,n,yp1,ypn,y2);
  if(DEBUG==0){
    for(int i=0;i<n;i++) printf(" y2[%d] = %lf\n",i,y2[i]);
  }
  /* Call splint for interpolations */ 
  splint(x,y,y2,n,xi,yi); 
  if(code!=0) BEMT_error(" error in this routine %s\n",thisroutine);
  return;
  }


/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
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
     Ind = 0;
  }
  else if(XVal>=XAry[AryLen-1]){
     Ind = AryLen;
  }
  else{
     // Let's interpolate!
  Ind  = 0;
  IHi  = AryLen;
  do {
     IMid = ( IHi + Ind )/2;
     if( XVal >= XAry[IMid]){
        Ind = IMid;
     }
     else{
       IHi = IMid;
     }
    }while(IHi-Ind>1);
  }
  return Ind;
  }


  double CubicSplineInterps (double X, int AryLen, double XAry[], double YAry[], double **Coef, double InterpVal)
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

  int ErrStat = 0;
  char *ErrMsg = "";

  // ! See if X is within the range of XAry.  Return the end point if it is not.

  if(X<=XAry[0]){
    InterpVal=YAry[0];
     return InterpVal;
  }
  else if(X>=XAry[AryLen-1]){
     InterpVal=YAry[AryLen-1];
     return InterpVal;
  } // ( X <= XAry(1) )
  //! We are somewhere inside XAry.  Find the segment that bounds X using binary search.

  ILo=LocateBin( X, XAry,AryLen );

  XOff=X-XAry[ILo];

  InterpVal=Coef[ILo][0] + XOff*(Coef[ILo][1] + XOff*(Coef[ILo][2] + XOff*Coef[ILo][3]));

  return InterpVal;
}



   double **CubicSplineInit(int AryLen, double *XAry, double *YAry)
   {
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

   double **Coef=NULL; Coef=dmatrix(0,AryLen-1,0,4);
   double *DelX=NULL; DelX=dvector(0,AryLen-1);
   double *Slope=NULL; Slope=dvector(0,AryLen-1);
   double *U=NULL; U=dvector(0,AryLen-1);
   double *V=NULL; V=dvector(0,AryLen-1);
   double ZHi,ZLo;

   for(int i=0;i<AryLen-1;i++){
     DelX[i]=XAry[i+1]-XAry[i];
     Slope[i]=(YAry[i+1]-YAry[i])/DelX[i];
     printf(" Slope = %lf\n",Slope[i]);

   }
   
   //! Use Gaussian elimination to solve the tri-diagonal matrix.

   U[0]=2.0*(DelX[1] + DelX[0]);
   V[0]=6.0*(Slope[1] - Slope[0]);
   
   //I=2,AryLen-1
   for(int i=1;i<AryLen-2;i++){
     U[i] = 2.0*(DelX[i-1]+ DelX[i]) - DelX[i-1]*DelX[i-1]/U[i-1];
     V[i] = 6.0*(Slope[i] - Slope[i-1]) - DelX[i-1]*V[i-1]/U[i-1];
   }
   //! Determine the coefficients of the polynomials.

   for(int i=0;i<AryLen-1;i++){
     Coef[i][0]=YAry[i];
  }
   ZHi=0.0;

//   DO I=AryLen-1,1,-1
   for(int i=AryLen-2;i>0;--i){
      ZLo       = ( V[i] - DelX[i]*ZHi )/U[i];
      Coef[i][1]= Slope[i]-DelX[i]*(ZHi/6.0+ZLo/3.0);
      Coef[i][2]= 0.5*ZLo;
      Coef[i][3]= ( ZHi - ZLo )/( 6.0*DelX[i]);
      ZHi=ZLo;
   }

/*
   CALL ExitThisRoutine ( ErrID_None, 'No Problemo' )
*/

  return Coef;
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/







  double *diff(double *vec,\
               int N)
  /* */
  {
  double *diffvec=NULL; diffvec=dvector(0,N);
  int i;
  for(i=0;i<N;i++){
    diffvec[i]=vec[i+1]-vec[i];
  }
  return diffvec;
  }

  double *linspace(double start,\
                   double end,\
                   int N)
  /* */
  {
  double *vec=NULL;vec=dvector(0,N);
  for(int i=0;i<N;i++) vec[i]=start+(end-start)/(N-1)*i;
  return vec;
  }

  minValLoc findMinValLoc(double *arr,\
                          int numArr)
  /* */
  {
  minValLoc mVL;
  int c, location=0;

  for (c=0;c<numArr;c++)
    if (arr[c]<arr[location])
      location=c;
  mVL.minValIdx=location;
  mVL.minVal=arr[location];
  return mVL;
  }

char* cutoffstr(const char* str,\
                int from,\
                int to)
/* */
{
    if (from >= to)
        return  NULL;

    char* cut = calloc(sizeof(char), (to - from) + 1);
    char* begin = cut;
    if (!cut)
        return  NULL;
    const char* fromit = str+from;
    const char* toit = str+to;
    (void)toit;
    memcpy(cut, fromit, to);
    return begin;
}

  int strFind(Airfoil_t *profile,\
              int  N,\
              int  idx,\
              char pattern[128])
  /*
  strcmp function is a build in function that returns 0 if both strings are identical,
  for more information, please google it "strcmp()")
  */
  {
  const char* thisroutine="int strFind(...)";
  int strIdx=-1;
  char* first4chars;
  //for(int k=0;k<N;k++){
    first4chars=cutoffstr(profile[idx].FoilFN,0,4);
    if(strcmp(first4chars,pattern)==0){
      strIdx=idx;}
    if(strcmp(first4chars,pattern)!=0){
      //BEMT_warning("String can not be found.\n", thisroutine);
    }
  //}
  //printf("%s \n",first4chars);
  //printf("true %s %s\n",first4chars,pattern);
  return strIdx;
  }

  /* Computation routine functions */
  geometryOut_t geometry(BladeData Blade,\
                         Airfoil_t *BladeData,\
                         double R_hub,\
                         double h_hub,\
                         int n,\
                         double step)
  /* 
  */
  {
  geometryOut_t BladeGeom;
  int N;N=Blade.NumBladeSections;
  BladeGeom.Mu_blade_interp=NULL;
  BladeGeom.Ml_blade_interp=NULL;
  BladeGeom.r_interp=NULL;
  BladeGeom.Z_grid_int=NULL;
  BladeGeom.Line_u_x_interp=NULL;
  BladeGeom.Y_u_int=NULL;
  BladeGeom.Line_l_x_interp=NULL;
  BladeGeom.Y_l_int=NULL;
  double *m=NULL;m=dvector(0,N);
  double *p=NULL;p=dvector(0,N);
  double *nc=NULL;nc=dvector(0,N);
  char *ncstr=NULL,*mstr=NULL,*pstr=NULL;
  double *deltaGeom=NULL;deltaGeom=dvector(0,N);
  double *beta=NULL;beta=dvector(0,n);
  double *x=NULL;x=dvector(0,n);
  double *yc=NULL;yc=dvector(0,n);
  double *dyc=NULL;dyc=dvector(0,n);
  double *theta=NULL;theta=dvector(0,n);
  double *yt=NULL;yt=dvector(0,n);
  d2_t *upper=NULL;upper=d2_tvector(0,n);
  d2_t *lower=NULL;lower=d2_tvector(0,n);
  d2_t *upper_f=NULL;upper_f=d2_tvector(0,n);
  d2_t *lower_f=NULL;lower_f=d2_tvector(0,n);
  double *z_f=NULL;z_f=dvector(0,n);
  d2_t **NACAprofiles=NULL;NACAprofiles=d2_tmatrix(0,N,0,n);
  d2_t **rotorAirfoils=NULL;rotorAirfoils=d2_tmatrix(0,N,0,n);
  double **Rot=NULL;Rot=dmatrix(0,2,0,2);
  double **Rot_u=NULL;Rot_u=dmatrix(0,2,0,n);
  double **Rot_l=NULL;Rot_l=dmatrix(0,2,0,n);
  double **Xu_rot=NULL;Xu_rot=dmatrix(0,2,0,n);
  double **Xl_rot=NULL;Xl_rot=dmatrix(0,2,0,n);

  double ***Mu_blade=NULL;Mu_blade=d3dtensor(0,n,0,3,0,N);
  double ***Ml_blade=NULL;Ml_blade=d3dtensor(0,n,0,3,0,N);

  /* INITIALIZATION */
  for(int k=0;k<N;k++) deltaGeom[k]=BladeData[k].x[2]*PI/180.0;
  beta=linspace(0,PI,n);
  for(int i=0;i<n;i++) x[i]=(1-cos(beta[i]))/2.0;

  d5_t a;
  a.x[0]=0.2969;
  a.x[1]=-0.126;
  a.x[2]=-0.3516;
  a.x[3]=0.2843;
  a.x[4]=-0.1015;// Open TE
  //a4=-0.1036;//   Closed TE
  int strIdx=-1;



  /* */
  for(int j=0;j<N;j++){
    strIdx=strFind(BladeData,N,j,"NACA");
    printf("strIdx=%d\n",strIdx);
    /* NACA Airfoils */
    //if(N-strIdx==1){
    if(strIdx!=-1){
      rmdir("NACAprofiles",0777);
      mkdir("NACAprofiles",0777);
      mstr=&BladeData[j].FoilFN[4];
      pstr=&BladeData[j].FoilFN[5];
      ncstr=cutoffstr(BladeData[j].FoilFN,6,7);
      mstr=cutoffstr(mstr,0,1);
      pstr=cutoffstr(pstr,0,1);
      m[j]=atof(mstr)/100.0;
      p[j]=atof(pstr)/10.0;
      nc[j]=strtof(ncstr,NULL)/100.0;
      printf("%s\n",mstr);
      printf("%s\n",pstr);
      printf("%s\n",ncstr);
      printf("%lf, %lf, %lf\n",m[j],p[j],nc[j]);
      
      /* Mean Chamber line coordinates and gradients */
      for(int i=0;i<n;i++){
        if(x[i]<p[j]){
          yc[i]=(m[j]/pow(p[j],2))*(2*p[j]*x[i]-pow(x[i],2));
          dyc[i]=2*m[j]/(pow(p[j],2))*(p[j]-x[i]);
          theta[i]=atan(dyc[i]);
        }
        else{
          yc[i]=(m[j]/pow((1-p[j]),2))*(1-2*p[j]+2*p[j]*x[i]-pow(x[i],2));
          dyc[i]=2*m[j]/(pow((1-p[j]),2))*(p[j]-x[i]);
          theta[i]=atan(dyc[i]);
        }
        /* Thickness distribution and radius of curvature */
        yt[i]=nc[j]/0.2*(a.x[0]*pow(x[i],0.5)\
             +a.x[1]*x[i]+a.x[2]*pow(x[i],2)\
             +a.x[3]*pow(x[i],3)\
             +a.x[4]*pow(x[i],4));
        /*Upper and lower surface coordinates*/
        upper[i].x[0]=x[i]-yt[i]*sin(theta[i]);
        upper[i].x[1]=yc[i]+yt[i]*cos(theta[i]);
        lower[i].x[0]=x[i]+yt[i]*sin(theta[i]);
        lower[i].x[1]=yc[i]-yt[i]*cos(theta[i]);
      }
      if(lower[0].x[0] != upper[0].x[0] || lower[0].x[1] != upper[0].x[1]){
        lower[0].x[0]=upper[0].x[0];
        lower[0].x[1]=upper[0].x[1];
      }
      FILE *fNACA[j];
      char fnNACA[256];
      sprintf(fnNACA, "NACAprofiles/NACA%s%s%s_%d.txt",mstr,pstr,ncstr,j+1);
/*
      if(j>0 && j<10){sprintf(fnNACA, "NACAprofiles/NACA%s%s%s_00%d.txt",mstr,pstr,ncstr,j);}
      if(j>9 && j<100){sprintf(fnNACA, "NACAprofiles/NACA%s%s%s_0%d.txt",mstr,pstr,ncstr,j);}
*/
      printf("%s ->----------test\n",fnNACA);

      fNACA[j] = fopen(fnNACA,"w");

      fprintf(fNACA[j],"NACA%s%s%s\n",mstr,pstr,ncstr);
      for(int k=n-1;k>0;k--){
        fprintf(fNACA[j],"%12.10f\t %12.10f\n",upper[k].x[0],upper[k].x[1]);
      }
      for(int k=0;k<n-1;k++){
        fprintf(fNACA[j],"%12.10f\t %12.10f\n",lower[k].x[0],lower[k].x[1]);
      }
      fclose(fNACA[j]);
    }
    /* If Airfoil from file >>>> There is a bug, needs to be resolved! does not work now!*/
    double xx=0.0,yy=0.0;
    int nl=0;
    if(strIdx==-1){
      rmdir("Rotor_airfoils",0777);
      mkdir("Rotor_airfoils",0777);
      printf("rotor-profiles need to be read!\n");
      int nl=0;
      int jj=j+1;
      FILE *files[j];
      char filename[256];
      if(jj>0 && jj<10){sprintf(filename, "rotor_profiles/rotor_profile_00%d.txt",jj);}
      if(jj>9 && jj<100){sprintf(filename, "rotor_profiles/rotor_profile_0%d.txt",jj);}
      printf("%s\n",filename);
      files[j] = fopen(filename,"r");
      fscanf(files[j], "%*[^\n]"); /* Read and discard a line */
      while(fscanf(files[j],"%lf %lf",&xx,&yy)!=EOF){
        nl++;
        printf("xu = %lf, yu = %lf, %d\n",xx,yy,nl);
        rotorAirfoils[j][nl].x[0]=xx;
        rotorAirfoils[j][nl].x[1]=yy;
      }
      fclose(files[j]);
      printf("xx = %lf, yy = %lf, %d\n",rotorAirfoils[j][nl].x[0],rotorAirfoils[j][nl].x[1],nl);
      printf("nl=%d\n",nl);
      d3_t upp,low;for(int i=0;i<3;i++) upp.x[i]=0.0,low.x[i]=0.0;
      int nll;nll=(nl+1)/2;
      printf("nll = %d\n",nll);
      double tblLx[nll],tblLy[nll];
      double tblUx[nll],tblUy[nll];
      for(int i=1;i<nl;i++){
        if(i<nll){
          tblLx[i]=rotorAirfoils[j][i].x[0];
          tblLy[i]=rotorAirfoils[j][i].x[1];
          printf("tblLy=%lf\n",tblLy[i]);
        }
        else{
          tblUx[i]=rotorAirfoils[j][i].x[0];
          tblUy[i]=rotorAirfoils[j][i].x[1];
          printf("tblUy=%lf\n",tblUy[i]);

        }
      }
      for(int k=0;k<n;k++){
      d3_t upp,low;for(int i=0;i<3;i++) upp.x[i]=0.0,low.x[i]=0.0;
        upper[k].x[0]=x[k];
        lower[k].x[0]=x[k];
//        upp=interp1D(nll,tblUx,tblUy,x[k]);
 //       low=interp1D(nll,tblLx,tblLy,x[k]);
        upper[k].x[1]=upp.x[0];
        lower[k].x[1]=low.x[0];
      }
      FILE *fRotor[j];
      char fnRotor[256];
      sprintf(fnRotor, "Rotor_airfoils/rotorAirfoils_%d.txt",j+1);

      fRotor[j] = fopen(fnRotor,"w");

      fprintf(fRotor[j],"Rotor_airfoils_%d\n",j+1);
      for(int k=n-1;k>0;k--){
        fprintf(fRotor[j],"%12.10f\t %12.10f\n",upper[k].x[0],upper[k].x[1]);
      /*for(int k=1;k<=nl;k++){
      //fprintf(fRotor[j],"%12.10f\t %12.10f\n",rotorAirfoils[j][k].x[0],rotorAirfoils[j][k].x[1]);
*/
      }
      for(int k=0;k<n-1;k++){
        fprintf(fRotor[j],"%12.10f\t %12.10f\n",lower[k].x[0],lower[k].x[1]);
      }
      fclose(fRotor[j]);
    }
    /*Scaling airfoil (with the real chord in mm)*/
    for(int i=0;i<n;i++){
      upper[i].x[0]=upper[i].x[0]*Blade.c[j]*1000;
      upper[i].x[1]=upper[i].x[1]*Blade.c[j]*1000;
      lower[i].x[0]=lower[i].x[0]*Blade.c[j]*1000;
      lower[i].x[1]=lower[i].x[1]*Blade.c[j]*1000;
    }
    /*Rotating airfoil respect 0 = c/4 */
    for(int i=0;i<n;i++){
      upper[i].x[0]= upper[i].x[0]-Blade.c[j]*1000/4;
      upper[i].x[1]= upper[i].x[1];
      lower[i].x[0]= lower[i].x[0]-Blade.c[j]*1000/4;
      lower[i].x[1]= lower[i].x[1];
      printf("xl = %lf\n",lower[i].x[0]);
    }
    Rot[0][0]=cos(deltaGeom[j]);
    Rot[1][0]=-sin(deltaGeom[j]);
    Rot[0][1]=sin(deltaGeom[j]);
    Rot[1][1]=cos(deltaGeom[j]);
    printf("Rot = %lf %lf\n",Rot[0][0],Rot[0][1]);
    printf("Rot = %lf %lf\n",Rot[1][0],Rot[1][1]);
    /* */
    for(int k=0;k<n;k++){
      Rot_u[0][k]=upper[k].x[0];
      Rot_u[1][k]=upper[k].x[1];
      Rot_l[0][k]=lower[k].x[0];
      Rot_l[1][k]=lower[k].x[1];
      /* */
      Xu_rot[0][k]=Rot[0][0]*Rot_u[0][k]+Rot[0][1]*Rot_u[1][k];
      Xu_rot[1][k]=Rot[1][0]*Rot_u[1][k]+Rot[1][1]*Rot_u[1][k];
      Xl_rot[0][k]=Rot[0][0]*Rot_l[0][k]+Rot[0][1]*Rot_l[1][k];
      Xl_rot[1][k]=Rot[1][0]*Rot_l[1][k]+Rot[1][1]*Rot_l[1][k];
      
    }
    printf("Xu_rot = %lf, %lf \n",Xu_rot[0][0],Xl_rot[0][n-1]);
    printf("Xu_rot = %lf, %lf \n",Xu_rot[1][0],Xl_rot[1][n-1]);
    /* Applying sweep */
    for(int i=0;i<n;i++){
      upper_f[i].x[0]=(Xu_rot[0][i]-Xu_rot[0][0])-Blade.sweep[j]*1000.0;
      upper_f[i].x[1]=Xu_rot[1][i]-(Xu_rot[1][n-1]+(Xu_rot[1][n-1]-Xl_rot[1][n-1])/2.0);
      lower_f[i].x[0]=(Xl_rot[0][i]-Xl_rot[0][0])-Blade.sweep[j]*1000.0;
      lower_f[i].x[1]=Xl_rot[1][i]-(Xl_rot[1][n-1]+(Xu_rot[1][n-1]-Xl_rot[1][n-1])/2.0);
      z_f[i]=-1.0*i/i*Blade.r[j]*1000;
    }
    for(int i=0;i<n;i++){
          Mu_blade[i][0][j]=upper_f[i].x[0];
          Mu_blade[i][1][j]=upper_f[i].x[1];
          Mu_blade[i][2][j]=z_f[i];

          Ml_blade[i][0][j]=lower_f[i].x[0];
          Ml_blade[i][1][j]=lower_f[i].x[1];
          Ml_blade[i][2][j]=z_f[i];
    printf("Mu %lf\n",Mu_blade[i][0][j]);
    printf("Ml %lf\n",Ml_blade[i][0][j]);
    }


  }






/* Freeing mems */
  free_dvector(deltaGeom,0,N);
  free_dvector(beta,0,n);
  free_dvector(x,0,n);
  free_dvector(yc,0,n);
  free_dvector(dyc,0,n);
  free_dvector(theta,0,n);
  free_dvector(m,0,N);
  free_dvector(p,0,N);
  free_dvector(nc,0,N);
  free_dvector(z_f,0,n);
  free_dmatrix(Rot,0,2,0,2);
  free_dmatrix(Rot_u,0,2,0,n);
  free_dmatrix(Rot_l,0,2,0,n);
  free_dmatrix(Xu_rot,0,2,0,n);
  free_dmatrix(Xl_rot,0,2,0,n);


  
  free_d2_tvector(upper,0,n);
  free_d2_tvector(lower,0,n);
  free_d2_tvector(upper_f,0,n);
  free_d2_tvector(lower_f,0,n);
  free_d2_tmatrix(NACAprofiles,0,N,0,n);
  free_d2_tmatrix(rotorAirfoils,0,N,0,n);

  free_d3dtensor(Mu_blade,0,n,0,3,0,N);
  free_d3dtensor(Ml_blade,0,n,0,3,0,N);

  return BladeGeom;
  }
/*===================================================================*/
/*                                                                   */
/*                          BEMT ROUTINES IN USE                     */
/*                                                                   */
/*===================================================================*/

  double ainduction(double CT)
  {
  /*
  This function calculates the induction factor 'a' as a function of thrust coefficient CT 
  including Glauert's correction
  */
  double a = 0.0;
  double CT1=1.816;
  double CT2=2*sqrt(CT1)-CT1;
  if(CT>=CT2){
    a=1 + (CT-CT1)/(4*(sqrt(CT1)-1));
  }
  else if(CT<CT2){
    a = 0.5-0.5*sqrt(1-CT);
  }
  return a;
 }

  PrndlCorr PrandtlTipRootCorrection(RotorData rotor, 
                                     double r_R, 
                                     double rootradius_R, 
                                     double tipradius_R, 
                                     double axial_induction)
  {
  /*
  This function calcualte steh combined tip and root Prandtl correction
  at agiven radial position 'r_R' (non-dimensioned by rotor radius), 
  given a rotor data that contains root and tip radius 
  (also non-dimensioned), a tip speed ratio
  TSR, the number lf blades NBlades and the axial induction factor
  The correction factor F can be written as:
                                                                       
    F=2/PI*arccos(exp(-f)), where     f=B/2*(R-r)/(r*sin(inflowAngle)) 
                                                                       
  
  */
  const char* thisroutine="PrandtlTipRootCorrection";
  double TINY=1.0e-25;
  if(DEBUG==0){
    printf("\n\n Debug info from routine %s\n\n",thisroutine);
    printf("TSR = %lf\n",rotor.TSR);
    printf("omega = %lf\n",rotor.omega);
    printf("rootR = %lf, tipR = %lf\n",rootradius_R,tipradius_R);
    printf(" 0/0 = %lf\n",(0.0/0.0));
    
  }
  
  PrndlCorr Pcorr; double Froot=0.0, Ftip=0.0, temp1=0.0, temp2=0.0;
  temp1=-rotor.NB/2*(tipradius_R-r_R)/r_R*sqrt(1+pow((rotor.TSR*r_R),2)/(pow((1-axial_induction),2)));
  Ftip=2/PI*acos(exp(temp1));
  temp2=rotor.NB/2*(rootradius_R-r_R)/r_R*sqrt(1+pow((rotor.TSR*r_R),2)/(pow((1-axial_induction),2)));
  Froot=2/PI*acos(exp(temp2));

  Pcorr.PrandtlRoot=Froot;
  Pcorr.PrandtlTip=Ftip;
  Pcorr.PrandtlRxT=Froot*Ftip;

  return Pcorr;
  }

  SkewedWakeCorr WakeCorr(RotorData rotor, 
                          double Vx,\
                          double Vy,\
                          double azimuth,\
                          double r_R,\
                          double a)
  {
  /* 
  */
  char *thisroutine="WakeCorr";
  SkewedWakeCorr SWCorr;
  SWCorr.skewAngle=rotor.yawAngle*PI/180;
  double yawCorr=0.0;
  double sineOfAz=sin(azimuth*PI/180);
  if(fabs(sineOfAz>0.005)){
    SWCorr.skewAngle=(0.6*a+1.0)*rotor.yawAngle*PI/180;
    yawCorr=(15.0*PI/32.0+tan(SWCorr.skewAngle/2.0)*r_R*sineOfAz);
    SWCorr.a=a*(1.0+yawCorr);
    if(DEBUG==1)printf(" axial induction: before %lf updated %lf from routine %s\n",a,SWCorr.a,thisroutine);
  }else{
    if(DEBUG==0)printf(" got a as it is\n");
    SWCorr.a=a;
    SWCorr.skewAngle=rotor.yawAngle;
  }
  return SWCorr;
  }

  bladeForce loadBladeElement(double v_n,\
                              double v_t,\
                              double chord,\
                              double twist,\
                              int numPolar,\
                              polarData polar)
  {
  /* computes the load in a blade element
  */
  int code=0;
  const char* thisroutine="loadBladeElement";
  if(DEBUG==0){
    printf("\n\n Debug info from routine %s\n\n",thisroutine);

  }
  bladeForce bladeF;
  /*reading and getting the polar data */
  double vmag2 = v_n*v_n + v_t*v_t;
  double inflowangle = atan2(v_n,v_t);
  double alpha =inflowangle*180.0/PI-twist;
  double Cl=0.0, Cd=0.0;
  /* interpolate Cl and Cd at alpha */
  cspline(polar.AoA, polar.Cl, numPolar, (polar.Cl[1]-polar.Cl[0]) \
          ,(polar.Cl[numPolar-1]-polar.Cl[numPolar-2]), alpha, &Cl);
  cspline(polar.AoA, polar.Cd, numPolar, (polar.Cd[1]-polar.Cd[0]) \
          ,(polar.Cd[numPolar-1]-polar.Cd[numPolar-2]), alpha, &Cd);

  if(DEBUG==0){
  printf(" alpha = %lf interpCL = %lf\n",alpha,Cl);
  printf(" alpha = %lf interCD = %lf\n",alpha,Cd);
  }
  double lift = 0.5*vmag2*Cl*chord;
  double drag = 0.5*vmag2*Cd*chord;
  bladeF.F_n=lift*cos(inflowangle)+drag*sin(inflowangle);
  bladeF.F_t=lift*sin(inflowangle)-drag*cos(inflowangle);
  bladeF.gamma=0.5*sqrt(vmag2)*Cl*chord;
  bladeF.AoA=alpha;
  bladeF.Cl=Cl;
  bladeF.Cd=Cd;
  
  return bladeF;
  }

  streamTube solveStreamTube(RotorData rotor,\
                             double r1_R,\
                             double r2_R,\
                             double rootR,\
                             double tipR,\
                             double omega,\
                             double chord,\
                             double twist,\
                             int numPolar,\
                             polarData polar)
{
/*

*/
  streamTube sTube;
  PrndlCorr Pcorr;
  SkewedWakeCorr YawCorr;
  bladeForce BForce;
  double anew=0.0;

  double AnularRingArea = PI*(pow(r2_R*rotor.R,2)-pow(r1_R*rotor.R,2)); /* area streamtube */
  double r_R = (r1_R+r2_R)/2; /* centroide */
  /* initiatlize variables */
  double a = 0.30; /* axial induction */
  double aline = 0.0; /* tangential induction factor */
  
    
  int Niter = 1000;
  /* */
  double Erroriterations =0.9e-10; /* error limit for iteration rpocess, in absolute value of induction */
  double Vrotor=0.0, Vtan=0.0, load3Daxial=0.0, CT=0.0, TSR=0.0;
  for(int i=0;i<Niter;i++){
     /* calculate the velocity and loads at blade element */
    Vrotor = rotor.Vinf*(1-a); /* axial velocity at rotor */
    //Vrotor=rotor.Vinf*(cos(rotor.yawAngle*PI/180)-a);
    //printf(" vrotor = %lf\n",Vrotor);
    Vtan = (1+aline)*omega*r_R*rotor.R; /* tangential velocity at rotor */
    //Vtan=(cos(rotor.yawAngle*PI/180)+aline)*omega*r_R*rotor.R;
    //printf(" Vtan = %lf\n",Vtan);
    /* calculate loads in blade segment in 2D (N/m) */
    BForce=loadBladeElement(Vrotor,Vtan,chord,twist,numPolar,polar);
    //printf("Fn = %lf Ft = %lf gamma = %lf\n",BForce.F_n, BForce.F_t,BForce.gamma);
    load3Daxial =BForce.F_n*rotor.R*(r2_R-r1_R)*rotor.NB; /* 3D force in axial direction */
    /*Calculate new estimate of axial and azimuthal induction */
    /* calculate thrust coefficient at the streamtube */
    //printf(" load3D = %lf\n",load3Daxial);
    CT = load3Daxial/(0.5*AnularRingArea*rotor.Vinf*rotor.Vinf);
    //printf("CT = %lf\n",CT);
    anew=ainduction(CT);
    /* correct new axial induction with Prandtl's correction */
    Pcorr=PrandtlTipRootCorrection(rotor,\
                                   r_R,\
                                   rootR,\
                                   tipR,\
                                   anew);
    if(Pcorr.PrandtlRxT<0.0001) Pcorr.PrandtlRxT=0.0001;/* avoid devide by zero */
    anew=anew/Pcorr.PrandtlRxT;/*correct estimate of axial induction */
    /* apply yaw correction */
    YawCorr=WakeCorr(rotor, 
                     Vrotor,\
                     Vtan,\
                     rotor.yawAngle,\
                     r_R,\
                     anew);
    anew=YawCorr.a;
    a=0.75*a+0.25*anew; /* for improving convergence, weigh current and previous iteration of axial induction */
    /* calculate azimuthal/tangential induction*/
    aline = BForce.F_t*rotor.NB/(2*PI*rotor.Vinf*(1-a)*omega*2*pow((r_R*rotor.R),2));
    aline =aline/Pcorr.PrandtlRxT; /* correct estimate of azimuthal induction with Prandtl's correction */
    /* test convergence of solution, by checking convergence of axial induction*/
    if (fabs(a-anew) < Erroriterations){
    if(DEBUG==0){
      printf("iterations %d    a-anew = %lf\n",i,fabs(a-anew));
      printf(" a\t %lf ap\t %lf r_R\t %lf F_n\t %lf F_t\t %lf gamma\t %lf\n",a\
             ,aline,r_R,BForce.F_n,BForce.F_t,BForce.gamma);
      }
       break;
    }
    sTube.a=a;if(a==0.0/0.0) sTube.a=0.0;
    sTube.a_p=aline;if(aline==0.0/0.0) sTube.a_p=0.0;
    sTube.r_R=r_R;
    sTube.F_n=BForce.F_n;
    sTube.F_t=BForce.F_t;
    sTube.gamma=BForce.gamma;
    sTube.AoA=BForce.AoA;
    sTube.Cl=BForce.Cl;
    sTube.Cd=BForce.Cd;
  }

  return sTube;
}


  RotorData getBladeGeom(RotorData rotor,caseData caseIN)
  {
  /*
  */
  FILE *fbg;
  char fnm[256];strcpy(fnm,caseIN.caseName);

  strcat(fnm,"_BladeGeom.txt");
  fbg=fopen(fnm,"w");

  double tipN=rotor.R/rotor.R;
  rotor.BGeom.r_R=NULL; rotor.BGeom.r_R=linspace(rotor.hubR,tipN,rotor.BGeom.NumStations);

  rotor.BGeom.c_R=dvector(0,rotor.BGeom.NumStations);
  rotor.BGeom.twist=dvector(0,rotor.BGeom.NumStations);
  
  if(DEBUG==0){
    for(int i=0;i<rotor.BGeom.NumStations;i++){
      printf("Numstations =  %d, r_R = %lf\n",i,rotor.BGeom.r_R[i]);
    }
  }
  fprintf(fbg,"%8s %8s %8s \n","r/R", "c/R", "twist [deg]");
  for(int i=0;i<rotor.BGeom.NumStations;i++){
    rotor.BGeom.c_R[i]=3*(1-rotor.BGeom.r_R[i])+1; /* [meters] */
    rotor.BGeom.twist[i]=14*(1-rotor.BGeom.r_R[i])+rotor.pitchAngle; /* [degrees] */
    fprintf(fbg,"%8.4f %8.4f %8.4f\n",rotor.BGeom.r_R[i],rotor.BGeom.c_R[i],rotor.BGeom.twist[i]);
  }
  fclose(fbg);
  return rotor;
  }

  /* */

  double get_cpu_time(clock_t start,\
                      clock_t end)
  {
  /*
  */
  double cpu_time_used=0.0;
  cpu_time_used=((double) (end-start))/CLOCKS_PER_SEC;
  printf("\n\n");
  printf(" Computation completed sucessfully! \n");
  printf(" CPU time %lf seconds\n", cpu_time_used);
  printf(" -----------------------------------------------\n");
  return cpu_time_used;
}

