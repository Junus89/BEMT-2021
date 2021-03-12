#pragma once
#ifndef BEMT_H_
#include<complex.h>
#include<math.h>
#include<string.h>
#include<time.h>
#define BEMT_H_

  #define PI      3.14159265358979323846
  #define TWOPI   (2.0*PI)
  #define DEBUG   1

/* User defined data types */
  typedef struct{
    double x[2];
  }d2_t;

  typedef struct{
   double minVal;
   int    minValIdx;
  }minValLoc;

  typedef struct{
    double complex x[2];
  }dc2_t;

  typedef struct{
    double x[3];
  }d3_t;

  typedef struct{
    double complex x[3];
  }dc3_t;

  typedef struct{
    double x[4];
  }d4_t;

  typedef struct{
    double x[5];
  }d5_t;

  typedef struct{
    double x[4];
    char   FoilFN[128];
  }Airfoil_t;


  typedef struct{
    double complex x[5];
  }dc5_t;

  typedef struct{
    double x[6];
  }d6_t;

  typedef struct{
    double complex x[6];
  }dc6_t;

  typedef struct{
    char caseName[256];
    int  ifCheckGeomSTL;
    int  ifSrfcPressureGen;
    int  ifFWHcomp;
    int  ifEPNLcomp;
    int  ifBBnoiseComp;
    int  ifNearFieldNoise;
    int  ifSTLgen;
    int  ifPolarsComp;
    int  ifBEMTcomp;
    int  ifTEBLparam;
    int  ifInstallXFOIL;
  }INCheck;
  
  typedef struct{
    double FlightSpeed;
    double AC_PitchAngle;
    int    ifUnsteadyLoad;
    double relativeRandomForce;
  }FlightData;
  
  typedef struct{
    double Pressure; //[Pa]
    double Temp; //[degK]
    double DynViscos;// [Pa*s]
    double RelativeHumi; // for atmospheric sound attenuation [%]
    double Z;
    double rho;
  }ATMData;

  typedef struct{
    d3_t Coords;
    int  NumProp;
    int  ifUncorr; //Uncorreltaed left/right propellers
  }ACData;

  typedef struct{
    d3_t   ReynoldsRange;
    d3_t   AlphaRange;
    int    XfoilMaxNumIter;
    int    NumPHalfFoil;
    double BLtransitPsuc;
    double BLtransitPpre;
  }ADynData;

  typedef struct{
    d3_t   Coords;
    double Diameter;
    double FirstBladeAngle2y;
    double RPM;
    double RotationSign;
    int    NumBlades;
    double Hub2TipRatio;
    double HubHoleDiameter;
    double HubHeight;
    int NumPHalfFoilSTL;
    double FoilInterpStep;
    int NumBladeSections;
    double n_j;
    double omega;
    double HubHoleRadius;
    double J;
  }PropellerData;
  
  typedef struct{
    char BladeDataFN[256];
    char PropellerJetFlowFN[256];
    int  ifBladeJetInteract;
    int  NumBladeSections;
    double *r;
    double *c;
    double *sweep;
    double R;
  }BladeData;

  typedef struct{
    int NumRotorRevol;
    int NumTimeStepPerRevol;
  }TimeIntegData;
  
  typedef struct{
    double Radius;
    int    NumParallel;
    int    NumMeridian;
    char   MicsFN[256];
    int    MicsFileFormat;
  }HemisphereData;

  typedef struct{
    char  TrajectoryFN[256];
    double FlightAltitude;
    d2_t   GrndMeshSize;
    d2_t   GrndMeshCellSize;
    double ReceiverHeight;
    double PNLTimeWindow;
    int    ifGroundReflection;
    int    EnhancedAurelization;
  }FootprintData;

  typedef struct{
    int VisualFormat;
  }OutputData;

  typedef struct{
    char MpirunPath[256];
    int NumProcs;
  }SystemInfo;

/* */
  typedef struct{
    double ***Mu_blade_interp;
    double ***Ml_blade_interp;
    double *r_interp;
    double **Z_grid_int;
    double **Line_u_x_interp;
    double *Y_u_int; 
    double **Line_l_x_interp; 
    double *Y_l_int;
  }geometryOut_t;

  typedef struct{
    double T0;    /* Temperature at MSL [K] */
    double rho;   /* Density [kg/m^3]*/
    double mu;    /* Dry air dynamic viscosity [kg/m-s] */
    double nu;    /* Dry air knematic viscosity [m^2/s] */
    double R;     /* Air specific constant [J/kg-K] */
    double gamma; /* Air specific heat ratio [-] */
  }air;

  typedef struct{
    int    NumStations;
    double *r_R;
    double *c_R;
    double *twist;
  }BladeGeom;


  typedef struct{
    char caseName[256];
    int runID;
  }caseData;


  typedef struct{
    //double r_R;   /* Non-dimensional radial distance [-] */
    //double c_R;   /* chord distribution [-] */
    //double theta; /* geometric pitch angle [rad] */
    double R;     /* Propeller Diameter [m] */
    double hubR;  /* hub radius [m] */
    double tipR;  /* tip radius [m] */
    int    NB;    /* Number of blades */
    int    RPM;   /* rotor RPM [-] */
    double omega; /* rotor rotational speed [rad/s] */
    double TSR;   /* Tip Speed Ratio (Lambda) */
    double Vinf;  /* Free stream velocity */
    double pitchAngle; /* rotor pitch angle */
    double yawAngle;   /* rotor yaw angle */
    char polarDataFN[256]; /* File name of the polar data file */
    BladeGeom BGeom;
    double CT;
    double CP;
    caseData caseIN;
  }RotorData;

  typedef struct{
    double F_n;
    double F_t;
    double gamma;
    double AoA;
    double Cl;
    double Cd;
  }bladeForce;

  typedef struct{
    double PrandtlRxT;
    double PrandtlTip;
    double PrandtlRoot;
  }PrndlCorr;

  typedef struct{
    double a;
    double a_p;
    double skewAngle;
    double wakeAngle;
  }SkewedWakeCorr;

  typedef struct{
    double *AoA;
    double *Cl;
    double *Cd;
    double *Cm;
  }polarData;

  typedef struct{
    double a;
    double a_p;
    double r_R;
//    bladeForce *bF; 
    double F_n;
    double F_t;
    double gamma;
    double AoA;
    double Cl;
    double Cd;
  }streamTube;

  typedef struct{
    double *a;
    double *a_p;
    double *r_R;
//    bladeForce *bF; 
    double *F_n;
    double *F_t;
    double *gamma;
  }BEMTresults;

  enum clLinearRange{
    foundmin=1,\
    foundmax=2,\
    finished=3
  };

  void BEMT_error(char error_text[],\
                  const char* routine);


  void BEMT_warning(char error_text[],\
                    const char* routine);
  /*------------------------------------------------------------------*/
  /* Functions/Routines for memory allocation for 1D-4D arrays*/

  /*--------------------------------1D--------------------------------*/
  double *dvector(long nl,long nh);

  void free_dvector(double *inputAr,\
                    long nl,long nh);

  /* complex dvector */
  double complex *dcvector(long nl,long nh);

  void free_dcvector(double complex *inputAr,\
                     long nl,long nh);

  /*--------------------------------2D--------------------------------*/

  double **dmatrix(long nrl,long nrh, \
                   long ncl,long nch);
  /* */
  void free_dmatrix(double **inputAr,\
                    long nrl,long nrh,\
                    long ncl,long nch);

  d2_t **d2_tmatrix(long nrl,\
                    long nrh,\
                    long ncl,\
                    long nch);
  /* */
  void free_d2_tmatrix(d2_t **inputAr,\
                       long nrl,\
                       long nrh,\
                       long ncl,\
                       long nch);
/*------------------------------3D------------------------------------*/

  double ***d3dtensor(long nrl,long nrh,\
                      long ncl,long nch,\
                      long ndl,long ndh);

  void free_d3dtensor(double ***inputAr,\
                      long nrl,long nrh,\
                      long ncl,long nch,\
                      long ndl,long ndh);

  /* */
  /* complex dmatrix */
  double complex **dcmatrix(long nrl,long nrh,\
                            long ncl,long nch);

  void free_dcmatrix(double complex **inputAr,\
                     long nrl,long nrh,\
                     long ncl,long nch);


  /* Mem functions for user defined data types */
  d2_t *d2_tvector(long nl,long nh);
  void free_d2_tvector(d2_t *inputAr,long nl,long nh);

  d3_t *d3_tvector(long nl,long nh);
  void free_d3_tvector(d3_t *inputAr,long nl,long nh);
/*--------------------------- File reader utility functions ----------*/
  d3_t *freader3col(char* fn,int* n);

  d4_t *freader4col(char* fn,int* n);

  d5_t *freader5col(char* fn,int* n);

  d4_t *polarDataReader(char* fn,int* n);

  Airfoil_t *AirfoilDBreader(char* fn,int* n);


/* --------------------------- Math utility functions ----------------*/
  d2_t blii1d(d2_t x,\
              d2_t y,
              double xi);

  double bcui1d(d4_t x,\
              d4_t f,\
              double xi);

  d2_t blii1d(d2_t x,\
              d2_t y,
              double xi);

  int brcket(int n,\
             double x[],\
             double xi);

  double interp1D(int nx,\
                double tabx[],\
                double tabc[],\
                double xi);

  void ratint(double xa[],\
              double ya[],\
              int n,\
              double x,\
              double *y);


  void polint(double xa[],\
              double ya[],\
              int n,\
              double x,\
              double *y);

  void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
  void splint(double xa[], double ya[], double y2a[], int n, double xi, double *yi);
  void cspline(double x[], double y[], int n, double yp1, double ypn, double xi, double *yi);

  double **CubicSplineInit(int AryLen, double *XAry, double *YAry);
  double CubicSplineInterps (double X, int AryLen, double XAry[], double YAry[], double **Coef, double InterpdVal);
  int LocateBin( double XVal, double XAry[], int AryLen);


  double *diff(double *vec,\
               int N);

  double *linspace(double start,\
                   double end,\
                   int N);

  minValLoc findMinValLoc(double *arr,\
                          int numArr);

  char* cutoffstr(const char* str, int from , int to);

  int strFind(Airfoil_t *profile,\
              int  N,\
              int  idx,\
              char pattern[128]);

  /* Computation routine functions */
  geometryOut_t geometry(BladeData Blade,\
                         Airfoil_t *BladeData,\
                         double R_hub,\
                         double h_hub,\
                         int n,\
                         double step);
  /* */
  double ainduction(double CT);

  PrndlCorr PrandtlTipRootCorrection(RotorData rotor, 
                                     double r_R, 
                                     double rootradius_R, 
                                     double tipradius_R, 
                                     double axial_induction);

  SkewedWakeCorr WakeCorr(RotorData rotor, 
                          double Vx,\
                          double Vy,\
                          double azimuth,\
                          double r_R,\
                          double a;

  RotorData getBladeGeom(RotorData rotor, caseData caseIN);

  bladeForce loadBladeElement(double v_n,\
                              double v_t,\
                              double chord,\
                              double twist,\
                              int numPolar,\
                              polarData polar);

  streamTube solveStreamTube(RotorData rotor,\
                             double r1_R,\
                             double r2_R,\
                             double rootR,\
                             double tipR,\
                             double omega,\
                             double chord,\
                             double twist,\
                             int numPolar,\
                             polarData polar);

/* */
  double get_cpu_time(clock_t start,\
                      clock_t end);


#endif
