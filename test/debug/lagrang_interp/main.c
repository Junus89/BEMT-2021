# include <math.h>
# include <stdbool.h>
# include <stdlib.h>
# include <stdio.h>
# include <time.h>
# include <string.h>

# include "lagrange_interp_1d.h"
# include "r8lib.h"



int main(){

double *cl=NULL;
double ac[5]={-3, -2.3,-1,0 ,1.2};
double ym[5]={-0.5,-0.2,-0.04, 0.5, 1.5};
double xi[5]={-3.1, -2.1,-1.1,0.1,1.5};
double *yi;

yi=lagrange_basis_1d ( 5,ac,5,xi );


return 0;
}
