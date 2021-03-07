#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>

#include "BEMT.h"

int main(int argc, char **argv)
{
  /*cpu time variables */
  clock_t start=0.0;
  clock_t end=0.0;
  double cpu_time;
  start=clock();
  /*reading input file */
  const char* thisroutine="int main(...)";
  /* input data declerations */
  INCheck         IN;
  FlightData      FL;
  ATMData         atmNG;
  ATMData         atmNP;
  ACData          AC;
  ADynData        AD;
  PropellerData   Prop;
  BladeData       Blade;
  TimeIntegData   TI;
  HemisphereData  Hemisphere;
  FootprintData   Footprint;
  OutputData      Output;
  SystemInfo      SysInfo;
  char fnSrcGeom[128],skipLine[256];
  FILE *fp=NULL;
  if(argc>=2){
    fp=fopen(argv[1],"r");
  }
  else{
    BEMT_error("Can no read the input file.\n\
Please check your input file!\n",thisroutine);
  }
  while(1==fscanf(fp,"%s%*[^\n] \
                      %d%*[^\n] \
                      %d%*[^\n] \
                      %d%*[^\n] \
                      %d%*[^\n] \
                      %d%*[^\n] \
                      %d%*[^\n] \
                      %d%*[^\n] \
                      %d%*[^\n] \
                      %d%*[^\n] \
                      %d%*[^\n] \
                      %d%*[^\n] \
                      %s%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %d%*[^\n] \
                      %lf%*[^\n] \
                      %s%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %s%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %s%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %d%*[^\n] \
                      %d%*[^\n] \
                      %s%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %d%*[^\n] \
                      %d%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %s%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %d%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %d%*[^\n] \
                      %lf%*[^\n] \
                      %d%*[^\n] \
                      %s%*[^\n] \
                      %s%*[^\n] \
                      %s%*[^\n] \
                      %d%*[^\n] \
                      %s%*[^\n] \
                      %d%*[^\n] \
                      %d%*[^\n] \
                      %s%*[^\n] \
                      %lf%*[^\n] \
                      %d%*[^\n] \
                      %d%*[^\n] \
                      %s%*[^\n] \
                      %d%*[^\n] \
                      %s%*[^\n] \
                      %s%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %d%*[^\n] \
                      %d%*[^\n] \
                      %s%*[^\n] \
                      %d%*[^\n] \
                      %s%*[^\n] \
                      %s%*[^\n] \
                      %d%*[^\n] \
                      %s%*[^\n]\n", \
                      IN.caseName,\
                      &IN.ifCheckGeomSTL,\
                      &IN.ifSrfcPressureGen,\
                      &IN.ifFWHcomp,\
                      &IN.ifEPNLcomp,\
                      &IN.ifBBnoiseComp,\
                      &IN.ifNearFieldNoise,\
                      &IN.ifSTLgen,\
                      &IN.ifPolarsComp,\
                      &IN.ifBEMTcomp,\
                      &IN.ifTEBLparam,\
                      &IN.ifInstallXFOIL,\
                      skipLine,\
                      &FL.FlightSpeed,\
                      &FL.AC_PitchAngle,\
                      &FL.ifUnsteadyLoad,\
                      &FL.relativeRandomForce,\
                      skipLine,\
                      &atmNG.Pressure,\
                      &atmNG.Temp,\
                      &atmNG.DynViscos,\
                      &atmNG.RelativeHumi,\
                      skipLine,\
                      &atmNP.Pressure,\
                      &atmNP.Temp,\
                      &atmNP.DynViscos,\
                      &atmNP.RelativeHumi,\
                      skipLine,\
                      &AC.Coords.x[0],\
                      &AC.Coords.x[1],\
                      &AC.Coords.x[2],\
                      &AC.NumProp,\
                      &AC.ifUncorr,\
                      skipLine,\
                      &AD.ReynoldsRange.x[0],\
                      &AD.ReynoldsRange.x[1],\
                      &AD.ReynoldsRange.x[2],\
                      &AD.AlphaRange.x[0],\
                      &AD.AlphaRange.x[0],\
                      &AD.AlphaRange.x[0],\
                      &AD.XfoilMaxNumIter,\
                      &AD.NumPHalfFoil,\
                      &AD.BLtransitPsuc,\
                      &AD.BLtransitPpre,\
                      skipLine,\
                      &Prop.Coords.x[0],\
                      &Prop.Coords.x[1],\
                      &Prop.Coords.x[2],\
                      &Prop.Diameter,\
                      &Prop.FirstBladeAngle2y,\
                      &Prop.RPM,\
                      &Prop.RotationSign,\
                      &Prop.NumBlades,\
                      &Prop.Hub2TipRatio,\
                      &Prop.HubHoleDiameter,\
                      &Prop.HubHeight,\
                      &Prop.NumPHalfFoilSTL,\
                      &Prop.FoilInterpStep,\
                      &Prop.NumBladeSections,\
                      skipLine,\
                      Blade.BladeDataFN,\
                      Blade.PropellerJetFlowFN,\
                      &Blade.ifBladeJetInteract,\
                      skipLine,\
                      &TI.NumRotorRevol,\
                      &TI.NumTimeStepPerRevol,\
                      skipLine,\
                      &Hemisphere.Radius,\
                      &Hemisphere.NumParallel,\
                      &Hemisphere.NumMeridian,\
                      Hemisphere.MicsFN,\
                      &Hemisphere.MicsFileFormat,\
                      skipLine,\
                      Footprint.TrajectoryFN,\
                      &Footprint.FlightAltitude,\
                      &Footprint.GrndMeshSize.x[0],\
                      &Footprint.GrndMeshSize.x[1],\
                      &Footprint.GrndMeshCellSize.x[0],\
                      &Footprint.GrndMeshCellSize.x[1],\
                      &Footprint.ReceiverHeight,\
                      &Footprint.PNLTimeWindow,\
                      &Footprint.ifGroundReflection,\
                      &Footprint.EnhancedAurelization,\
                      skipLine,\
                      &Output.VisualFormat,\
                      skipLine,\
                      SysInfo.MpirunPath,\
                      &SysInfo.NumProcs,\
                      skipLine)){}
  fclose(fp);
  printf("---------------------------------------------------\n");
  printf("!                                                 !\n");
  printf("!       Blade Element Momentum Theory Solver      !\n");
  printf("!                                                 !\n");
  printf("!                      BEMT                       !\n");
  printf("!                                                 !\n");
  printf("!          Written by: Yunusi Fuerkaiti           !\n");
  printf("!                                                 !\n");
  printf("!              y.fuerkaiti@tudelft.nl             !\n");
  printf("!                                                 !\n");
  printf("---------------------------------------------------\n");
  printf("---------------------------------------------------\n");
  printf("!                                                 !\n");
  printf("!  %s is running ...                             !\n",IN.caseName);
  printf("!                                                 !\n");
  printf("---------------------------------------------------\n");
  printf("!                                                 !\n");
/* Initialization */
  printf("%s \n",skipLine);
  printf("%lf \n",FL.FlightSpeed);
  printf("%lf \n",atmNG.Pressure);
  printf("%lf \n",atmNP.Pressure);
  for(int i=0;i<3;i++) printf("Aircraft center coords: %lf\n",AC.Coords.x[i]);
  for(int i=0;i<3;i++) printf("Reynolds number range: %lf\n",AD.ReynoldsRange.x[i]);
  for(int i=0;i<3;i++) printf("Alpha range: %lf\n",AD.AlphaRange.x[i]);
  printf("%d \n",Prop.NumBlades);
  printf("%lf \n",Prop.RPM);
  printf("%s\n",Blade.BladeDataFN);
/*reading blade geometry */
  Airfoil_t *BladeData=NULL;
  BladeData=AirfoilDBreader(Blade.BladeDataFN,&Blade.NumBladeSections);
  printf("Number of blade sections %d\n",Blade.NumBladeSections);
  for(int i=0;i<Blade.NumBladeSections;i++){
    printf(" %lf %lf %lf %lf %s \n",BladeData[i].x[0],BladeData[i].x[1],BladeData[i].x[2],BladeData[i].x[3],BladeData[i].FoilFN);}
  printf("%d\n",TI.NumRotorRevol);
  printf("%d\n",TI.NumTimeStepPerRevol);
  printf("%lf\n",Hemisphere.Radius);
  printf("%s\n",Hemisphere.MicsFN);
  printf("%s\n",Footprint.TrajectoryFN);
  printf("%lf\n",Footprint.FlightAltitude);
  printf("%d\n",Footprint.ifGroundReflection);
  printf("%d\n",Output.VisualFormat);
  printf("%s\n",SysInfo.MpirunPath);
  printf("%d\n",SysInfo.NumProcs);
/* ---------End of Input data check---------- */
/* Additional Input data processing and Initializations */
/* Atm data processing */
  atmNG.Z=(288.15 - atmNG.Temp)/0.0065;
  atmNG.Rho=1.22*pow((1 - 0.0000226*atmNG.Z),4.256);
  printf("Z=%lf \n",atmNG.Z);
  printf("Rho=%lf \n",atmNG.Rho);
/* Propeller data processing */
  Prop.n_j=Prop.RPM/60.0;
  printf("n_j=%lf\n",Prop.n_j);
  Prop.omega=2*PI*Prop.n_j;
  printf("n_j=%lf\n",Prop.omega);
  Prop.HubHoleRadius=Prop.HubHoleDiameter/2.0;
  printf("HubHoleR=%lf\n",Prop.HubHoleRadius);
  Prop.FoilInterpStep=Prop.FoilInterpStep*1000.0;
  printf("Foil interp step=%lf\n",Prop.FoilInterpStep);
  Prop.J=FL.FlightSpeed/(Prop.n_j*Prop.Diameter);
  printf("J=%lf\n",Prop.J);
/* Blade data processing */
  Blade.r=NULL;Blade.r=dvector(0,Blade.NumBladeSections);
  Blade.c=NULL;Blade.c=dvector(0,Blade.NumBladeSections);
  Blade.sweep=NULL;Blade.sweep=dvector(0,Blade.NumBladeSections);
  Blade.R=0.0;Blade.R=Prop.Diameter/2.0;
  printf("%lf\n",Blade.R);
  for(int k=0;k<Blade.NumBladeSections;k++){
  Blade.r[k]=BladeData[k].x[0]*Blade.R;
  Blade.c[k]=BladeData[k].x[1]*Blade.R;
  Blade.sweep[k]=BladeData[k].x[3]*Blade.R;
  printf("r[%d] = %lf, %lf, %lf\n",k,Blade.r[k],Blade.c[k],Blade.sweep[k]);
  }

  double RefBladeAspectRatio=0.75;
  /* r/R=BladeData.x[0] */
  double *r75_dist=NULL; r75_dist=dvector(0,Blade.NumBladeSections);
  for(int j=0;j<Blade.NumBladeSections;j++){
    r75_dist[j]=fabs(RefBladeAspectRatio-BladeData[j].x[0]);
    //printf("r75_dist[%d] = %lf\n",j,r75_dist[j]);
  }
  minValLoc mVLAR; mVLAR=findMinValLoc(r75_dist,Blade.NumBladeSections);
  printf("MinVal=%lf, minValLoc=%d\n",mVLAR.minVal,mVLAR.minValIdx);
  double AR=0; AR=Blade.R/Blade.c[mVLAR.minValIdx];
  printf("AR=%lf\n",AR);
/* Warning */
  if(IN.ifBEMTcomp==1 && IN.ifTEBLparam==0){
    BEMT_error("In order to compute BL parameters the load computation is necessary.\n\
    Please set line 10 of the input file as 0 and run again.'\n",thisroutine);}
  if(IN.ifBEMTcomp==1 && (IN.ifFWHcomp==0 &&IN.ifBBnoiseComp==0)){
    BEMT_error("In order to perform noise analysis the load computation is necessary.\n\
    Please set line 10 of the input file as 0 and run again.'\n",thisroutine);}
  if(IN.ifTEBLparam==1 && (IN.ifFWHcomp==0 &&IN.ifBBnoiseComp==0)){
    BEMT_error("In order to perform noise analysis the parameters computation is necessary.\n\
    Please set line 11 of the input file as 0 and run again.'\n",thisroutine);}
  if((IN.ifFWHcomp==0 && IN.ifBBnoiseComp==1) && (IN.ifFWHcomp==1 &&IN.ifBBnoiseComp==0)){
    BEMT_error("Broadband and tonal noise modules work together.\n\
    Please set both line 4 and 6 of the input file as 0 and run again.'\n",thisroutine);}
  printf("%d\n",IN.ifBEMTcomp);
/* Generation of the Blade surface ->worked well*/
/* checking linspace function */
/*
  double *vec=NULL;vec=dvector(0,5);
  vec=linspace(0,PI,5);
  for(int k=0;k<5;k++) printf("%lf\n",vec[k]);
  */
/*checking strFind function */
  //int strIdx=0;strIdx=strFind(BladeData,Blade.NumBladeSections,"NACA");
  //printf("%d\n",strIdx);
  geometryOut_t BladeGeom;
  geometry(Blade,\
           BladeData,\
           Prop.HubHoleRadius,\
           Prop.HubHeight,\
           Prop.NumPHalfFoilSTL,\
           Prop.FoilInterpStep);/* later reduce it with just two input blade and prop or just prop that includes blade data */


/* test interpocu1*/
  d4_t x,y;
  for(int i=0;i<4;i++){
    x.x[i]=i*0.1;
    y.x[i]=x.x[i]*cos(x.x[i]);
    printf("%lf %lf \n",x.x[i],y.x[i]);}
  d3_t fi; double xi=0.29;
  fi=bcui1d(x,y,xi);
  printf("xi=%lf\n",xi);
  for(int i=0;i<3;i++) printf("fi = %lf\n",fi.x[i]);
  d2_t x1,y1,f1;
  for(int i=0;i<2;i++){
    x1.x[i]=i*0.1+0.3;
    y1.x[i]=x1.x[i]*cos(x1.x[i]);
    printf("%lf %lf \n",x1.x[i],y1.x[i]);}
  f1=blii1d(x1,y1,xi);
  printf("x1 = %lf, %lf\n",x1.x[0],x1.x[1]);
  printf("y1 = %lf, %lf\n",y1.x[0],y1.x[1]);
  printf("f1 = %lf, %lf\n",f1.x[0],f1.x[1]);
  /*test brcket */















/*Freeing memories */
/*1D */
  free_dvector(r75_dist,0,Blade.NumBladeSections);
  free_dvector(Blade.r,0,Blade.NumBladeSections);
  free_dvector(Blade.c,0,Blade.NumBladeSections);
  free_dvector(Blade.sweep,0,Blade.NumBladeSections);

  free(BladeData);
  end=clock();
  cpu_time=get_cpu_time(start,end);

  return 0;

}
