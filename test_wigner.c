#include <stdio.h>
#include "WignerCoupling.h"

int main(int argc, char *argv[]) {
    double WigVal, NormVal, *WigVal2;
    int i, j1, j2, jmin;

    WigVal = Wigner6j(2*3,2*4,2*4,2*2,2*3,2*2);

    printf("WigVal = %f\n",WigVal);

    return 0;
}
