/* Calculation of the Wigner 9j symbols. Calculation is using definition of 
 * 9j symbols that writes them as a function of the 6j symbols*/

#include <math.h>
#include <stdlib.h>
#include "wignercoupling.h"

double Wigner9j(int tj1, int tj2, int tJ12, int tj3, int tj4, int tJ34, int tJ13, int tJ24, int tJ) {
    //calculate in terms of 6j symbols or 3j symbols
    //looks to be easier with 6j symbols since only one sum and not 6 as for 3j symbols

    //use Wigner6jVect function. Three calls and have all the families
    //that we need to sum over. Then do the sum of the required gs
    double j1, j2, J12, j3, j4, J34, J13, J24, J;
    double g, tempa, tempb, *Fam1, *Fam2, *Fam3, WigVal;
    int i, tmaxg, tming, a[3], NumG;
	double ming1, ming2, ming3, maxg1, maxg2, maxg3, ming, maxg;
	double NormVal1, NormVal2, NormVal3;

	j1 = (double) tj1 / 2.;
	j2 = (double) tj2 / 2.;
	J12 = (double) tJ12 / 2.;
	j3 = (double) tj3 / 2.;
	j4 = (double) tj4 / 2.;
	J34 = (double) tJ34 / 2.;
	J13 = (double) tJ13 / 2.;
	J24 = (double) tJ24 / 2.;
	J = (double) tJ / 2.;

	//9j is a sum over a range of product of 3 6j symbols
	//Find the acceptable range for g in each of the 3 symbols
	//Sum is over the intersection of this range
    tempa = (j1 > j2) ? j1-j2 : j2-j1;
    tempb = (J34 > J) ? J34-J : J-J34;
    ming1 = (tempa > tempb) ? tempa : tempb;
    maxg1 = ((j1+j2) < (J34+J)) ? j1+j2 : J34+J;

    tempa = (j2 > J24) ? j2-J24 : J24-j2;
    tempb = (j3 > J34) ? j3-J34 : J34-j3;
    ming2 = (tempa > tempb) ? tempa : tempb;
    maxg2 = ((j3+J34) < (j2+J24)) ? j3+J34 : j2+J24;

    tempa = (J > J24) ? J-J24 : J24-J;
    tempb = (j3 > j1) ? j3-j1 : j1-j3;
    ming3 = (tempa > tempb) ? tempa : tempb;
    maxg3 = ((J+J24) < (j3+j1)) ? J+J24 : j3+j1;

    ming = (ming1 > ming2) ? ming1 : ming2;
    ming = (ming > ming3) ? ming : ming3;
    maxg = (maxg1 < maxg2) ? maxg1 : maxg2;
    maxg = (maxg < maxg3) ? maxg : maxg3;

    NumG = (int) (maxg - ming);

	//Three different families
    Fam1 = Wigner6jVect(tJ34,tJ,tj1,tj2,tJ12,&NormVal1);
    Fam2 = Wigner6jVect(tj2,tJ24,tj3,tJ34,tj4,&NormVal2);
    Fam3 = Wigner6jVect(tj1,tj3,tJ13,tJ24,tJ,&NormVal3);

    tmaxg = (int) 2*maxg;
    tming = (int) 2*ming;
    for (i=tming;i<=tmaxg;i+=2) {
        g = (double) i / 2.;
		//the a values are the offsets in the Fam* arrays for this g value
        a[0] = (int) g-ming1;
        a[1] = (int) g-ming2;
        a[2] = (int) g-ming3;
        WigVal += pow(-1.,i) * (i+1.) * Fam1[a[0]]*Fam2[a[1]]*Fam3[a[2]];
    }

	WigVal *= NormVal1*NormVal2*NormVal3; //normalize all the Wigner 6j symbols
	
	free(Fam1);
	free(Fam2);
	free(Fam3);

    return WigVal;
}
