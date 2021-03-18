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
  const char* thisroutine="main";
  /* input data declerations */
  caseData        caseIN;
  RotorData       rotor;
  BladeData       blade;
  ATMData         atm;
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
                      %s%*[^\n] \
                      %lf%*[^\n] \
                      %d%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %d%*[^\n] \
                      %s%*[^\n] \
                      %s%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %s%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n] \
                      %lf%*[^\n]\n", \
                      caseIN.caseName,\
                      &caseIN.runID,\
                      skipLine,\
                      &rotor.R,\
                      &rotor.NB,\
                      &rotor.hubR,\
                      &rotor.pitchAngle,\
                      &rotor.yawAngle,\
                      &rotor.BGeom.NumStations,\
                      rotor.polarDataFN,\
                      skipLine,\
                      &rotor.Vinf,\
                      &rotor.TSR,\
                      &rotor.RPM,\
                      skipLine,\
                      &atm.Pressure,\
                      &atm.Temp,\
                      &atm.rho,\
                      &atm.DynViscos,\
                      &atm.RelativeHumi)){};
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
  printf("!  %s is running ...                              !\n",caseIN.caseName);
  printf("!                                                 !\n");
  printf("---------------------------------------------------\n");
/* */
  /*reading and getting the polar data */
  rotor.caseIN.runID=caseIN.runID;
  
  int numPolarP;
  d4_t *polar=NULL;polar=polarDataReader(rotor.polarDataFN,&numPolarP);
  polarData Polar;
  Polar.AoA=NULL; Polar.AoA=dvector(0,numPolarP);
  Polar.Cl=NULL; Polar.Cl=dvector(0,numPolarP);
  Polar.Cd=NULL; Polar.Cd=dvector(0,numPolarP);
  for(int i=0;i<numPolarP;i++){
     Polar.AoA[i]=polar[i].x[0];
     Polar.Cl[i]=polar[i].x[1];
     Polar.Cd[i]=polar[i].x[2];
  }
  if(DEBUG==0){
    printf(" Number of blades %d\n",rotor.NB);
    printf(" Free-stream velocity %lf\n",rotor.Vinf);
    printf(" Number of blade sections %d\n", rotor.BGeom.NumStations);
    printf(" Name of the polar data %s \n", rotor.polarDataFN);
    for(int i=0;i<numPolarP;i++){
      printf(" %d AoA = %lf Cl = %lf Cd = %lf Cm = %lf\n",i,polar[i].x[0],polar[i].x[1],polar[i].x[2],polar[i].x[3]);
    }
    //printf("\n\nCheck Prandtle Corrections\n\n");
    printf("\n\n Check Polar interpolation \n\n");
    double Cl=0.0, alpha=0.0, **Coef=NULL; Coef=dmatrix(0,numPolarP,0,4);
    Coef=CubicSplineInit(numPolarP, Polar.AoA, Polar.Cl);
    for(int i=0;i<numPolarP;i++){
      alpha=0.01*i+1.4-0.5*i;
      //printf("%d %lf %lf %lf %lf\n",i,Coef[i][0],Coef[i][1],Coef[i][2],Coef[i][3]);
      cspline(Polar.AoA, Polar.Cl, numPolarP, (Polar.Cl[1]-Polar.Cl[0]) \
              ,(Polar.Cl[numPolarP-1]-Polar.Cl[numPolarP-2]), alpha, &Cl);

     //printf(" alpha = %lf interpCL = %lf\n",alpha,Cl);
     //printf(" alpha = %lf interCD = %lf\n",alpha,Cd);
    }
  }
  /* compute and generate the blade geometry */
  rotor=getBladeGeom(rotor,caseIN);

  enum WTorProp WTP=caseIN.runID;
  switch(WTP){
    case(Windturbine):
      rotor.omega=rotor.Vinf*rotor.TSR/rotor.R;
      break;
    case(Propeller):
      rotor.omega=rotor.RPM/(30.0/PI);
      rotor.J=rotor.Vinf/(rotor.RPM/60*rotor.R*2);
      printf("rotor RPM %lf rotor omega %lf J %lf\n",rotor.RPM,rotor.omega, rotor.J);
      break;
    default:
      BEMT_error("Please specify run ID: 0->WT; 1-> Prop\n",thisroutine);
  }
  //printf("%lf\n",rotor.BGeom.r_R[0]);

  /* compute stream tube */
  streamTube sTube;
if(DEBUG==1){
  double chord=0.0, twist=0.0,r_R=0.0;
  double CT_r=0.0, CP_r=0.0,CQ_r=0.0;
  FILE *fm;
  char fnm[256];strcpy(fnm,caseIN.caseName);

  strcat(fnm,"_Restuls.txt");
  fm=fopen(fnm,"w");
  
  fprintf(fm,"%8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s\n","r/R"\
    ,"a","a'","Fn [N]","Ft[N]","CT","CQ","CP","Phi","AoA","Cl","Cd");
  for(int i=0;i<rotor.BGeom.NumStations-1;i++){
  //for(int i=0;i<1;i++){
  CT_r=0.0, CQ_r=0.0,CP_r=0.0;
    r_R=(rotor.BGeom.r_R[i]+rotor.BGeom.r_R[i+1])/2;
    ratint(rotor.BGeom.r_R,rotor.BGeom.c_R,rotor.BGeom.NumStations,r_R,&chord);
    ratint(rotor.BGeom.r_R,rotor.BGeom.twist,rotor.BGeom.NumStations,r_R,&twist);
    if(DEBUG==0) printf(" r/R = %lf chod = %lf twist = %lf\n",r_R, chord, twist);

    sTube=solveStreamTubeInWake(rotor,\
                          rotor.BGeom.r_R[i],\
                          rotor.BGeom.r_R[i+1],\
                          rotor.hubR,\
                          1.0,\
                          rotor.omega,\
                          chord,\
                          twist,\
                          numPolarP,\
                          Polar);

    if(sTube.a==-1/0.0) sTube.a=0.0;
    if(sTube.a_p!=sTube.a_p) sTube.a_p=0.0;
    if(sTube.F_n!=sTube.F_n) sTube.F_n=0.0;
    if(sTube.F_t!=sTube.F_t) sTube.F_t=0.0;
    if(sTube.gamma!=sTube.gamma) sTube.gamma=0.0;
    if(sTube.AoA!=sTube.AoA) sTube.AoA=0.0;
    if(sTube.Cl!=sTube.Cl) sTube.Cl=0.0;
    if(sTube.Cd!=sTube.Cd) sTube.Cd=0.0;

    double dr=(rotor.BGeom.r_R[i+1]-rotor.BGeom.r_R[i])*rotor.R;

    //CT_r=dr*sTube.F_n*rotor.NB/(0.5*pow(rotor.Vinf,2)*PI*pow(rotor.R,2));
    CT_r=dr*sTube.F_n*rotor.NB/(pow((rotor.RPM/60),2)*pow(2*rotor.R,4));
    CQ_r=dr*sTube.F_t*r_R*rotor.R*rotor.NB/(pow((rotor.RPM/60),2)*pow(2*rotor.R,5));
    CP_r=2*PI*CQ_r;

    rotor.CT+=CT_r; rotor.CQ+=CQ_r,rotor.CP+=CP_r;
    rotor.Thrust+=dr*sTube.F_n*rotor.NB;
    rotor.Torque+=dr*sTube.F_t*r_R*rotor.R*rotor.NB;
    rotor.Power=rotor.Torque*2*PI;
    fprintf(fm,"%8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f %8.4f  %8.4f\n",sTube.r_R\
             ,sTube.a,sTube.a_p,sTube.F_n,sTube.F_t,CT_r,CQ_r,CP_r, sTube.phi,sTube.AoA, sTube.Cl, sTube.Cd);

  }
  fclose(fm);
  if(DEBUG==0){
    printf("Omega = %lf\n",rotor.omega);
  }
    printf("\n\n--------------------------------------------------\n");
    printf("Thrust [N] = %lf,\t Torque [N.m] = %lf,\t Power [N.m/s]= %lf\n",rotor.Thrust,rotor.Torque,rotor.Power);
    printf("CT = %lf,\t CQ = %lf,\t CP = %lf\n",rotor.CT, rotor.CQ, rotor.CP);
    double eta_prop=(rotor.CT/rotor.CP)*rotor.J;
    printf("Propeller efficiency = %lf\n",eta_prop);

  }
  if(DEBUG==0){
    BEMT(rotor, atm, Polar, numPolarP,caseIN);
  }



/* dealloc */
//  free(polar);
   free_dvector(Polar.AoA,0,numPolarP);
   free_dvector(Polar.Cl,0,numPolarP);
   free_dvector(Polar.Cd,0,numPolarP);
/* */
  end=clock();
  cpu_time=get_cpu_time(start,end);

  return 0;

}
