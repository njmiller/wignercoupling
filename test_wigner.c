#include <stdio.h>
#include "WignerCoupling.h"

int main(int argc, char *argv[]) {
    double WigVal, NormVal, *WigVal2, ytmp, norm;
    int i, tj1, tj2, tj3;
	int tm1, tm2, tm3;

	//WigVal = Wigner6j(2*3,2*4,2*4,2*2,2*3,2*2);
    //printf("WigVal = %f\n",WigVal);
	
	//These should test all 4 cases (normal and the 3 special cases)
    WigVal = wigner3j(2*1,2*1,0,0,0,0);
    printf("WigVal = %f -0.577350\n",WigVal);
    WigVal = wigner3j(2*4,2*5,2*3,2*2,-2*2,0);
    printf("WigVal = %f 0.0215917\n",WigVal);
    WigVal = wigner3j(2*3,2*3,2*3,1*2,-1*2,0);
    printf("WigVal = %f -0.1543033\n",WigVal);
    WigVal = wigner3j(3,9,6,1,3,-4);
    printf("WigVal = %f 0.146385\n",WigVal);

    return 0;
}
