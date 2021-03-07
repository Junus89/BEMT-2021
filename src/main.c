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
                      %s%*[^\n] \
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
                      skipLine,\
                      &atm.Pressure,\
                      &atm.Temp,\
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
/*      ratint(Polar.AoA,\
             Polar.Cl,\
             numPolarP,\
             alpha,
             &Cl);
*/
    
      cspline(Polar.AoA, Polar.Cl, numPolarP, (Polar.Cl[1]-Polar.Cl[0]) \
              ,(Polar.Cl[numPolarP-1]-Polar.Cl[numPolarP-2]), alpha, &Cl);

     printf(" alpha = %lf interpCL = %lf\n",alpha,Cl);
      //printf(" alpha = %lf interCD = %lf\n",alpha,Cd);
    }
  }
  /* compute and generate the blade geometry */
  rotor=getBladeGeom(rotor,caseIN);
  rotor.omega=rotor.Vinf*rotor.TSR/rotor.R;
  //printf("%lf\n",rotor.BGeom.r_R[0]);

  /* compute stream tube */
  streamTube sTube;

  double chord=0.0, twist=0.0,r_R=0.0;
  FILE *fm;
  char fnm[256];strcpy(fnm,caseIN.caseName);

  strcat(fnm,"_Restuls.txt");
  fm=fopen(fnm,"w");
  
  fprintf(fm,"%8s %8s %8s %8s %8s %8s %8s %8s %8s\n","r/R", "a", "a'", "Fn [N]", "Ft[N]", "Gamma", "AoA", "Cl", "Cd");
  for(int i=0;i<rotor.BGeom.NumStations-1;i++){
    r_R=(rotor.BGeom.r_R[i]+rotor.BGeom.r_R[i+1])/2;
    ratint(rotor.BGeom.r_R,rotor.BGeom.c_R,rotor.BGeom.NumStations,r_R,&chord);
    ratint(rotor.BGeom.r_R,rotor.BGeom.twist,rotor.BGeom.NumStations,r_R,&twist);
    if(DEBUG==0) printf(" r/R = %lf chod = %lf twist = %lf\n",r_R, chord, twist);

    sTube=solveStreamTube(rotor,\
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

    fprintf(fm,"%8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f %8.4f %8.4f %8.4f\n",sTube.r_R\
             ,sTube.a,sTube.a_p,sTube.F_n,sTube.F_t,sTube.gamma, sTube.AoA, sTube.Cl, sTube.Cd);
    double dr=(rotor.BGeom.r_R[i+1]-rotor.BGeom.r_R[i])*rotor.R;
    rotor.CT +=dr*sTube.r_R*rotor.NB/(0.5*rotor.Vinf*rotor.Vinf*PI*rotor.R*rotor.R);
    rotor.CP +=dr*sTube.F_n*sTube.a_p*rotor.NB*rotor.R*rotor.omega/(0.5*rotor.Vinf*rotor.Vinf*rotor.Vinf*PI*rotor.R*rotor.R);

  }
  fclose(fm);
  printf("CT = %lf\t CP = %lf\n",rotor.CT, rotor.CP);






/* dealloc */
//  free(polar);
   free_dvector(Polar.AoA,0,numPolarP);
   free_dvector(Polar.Cl,0,numPolarP);
   free_dvector(Polar.Cd,0,numPolarP);

  end=clock();
  cpu_time=get_cpu_time(start,end);

  return 0;

}
