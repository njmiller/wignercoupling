/* Calculation of the Wigner 6j symbols. Same routine as Wigner 3j with different 
 * X, Y, and Z functions. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double Zj6j(double j1, double j2, double j3, double k1, double k2, double k3) {
    double tempa, tempb, tempc, tempd, tempe;

    tempa = (j3+1);
    tempb = pow(j3,2) - pow(j2-j1,2);
    tempc = pow(j1+j2+1,2) - pow(j3,2);
    tempd = pow(j3,2) - pow(k1-k2,2);
    tempe = pow(k1+k2+1,2) - pow(j3,2);

    return tempa * sqrt(tempb*tempc*tempd*tempe);
}

double Yj6j(double j1, double j2, double j3, double k1, double k2, double k3) {
    double tempa, tempb, tempc, tempd;

    tempa = 2*j3+1.;
    tempb = j3*(j3+1)*(-j3*(j3+1)+j1*(j1+1)+j2*(j2+1)-2*k3*(k3+1));
    tempc = k1*(k1+1)*(j3*(j3+1)+j1*(j1+1)-j2*(j2+1));
    tempd = k2*(k2+1)*(j3*(j3+1)-j1*(j1+1)+j2*(j2+1));

    return tempa*(tempb+tempc+tempd);
}

double Xj6j(double j1, double j2, double j3, double k1, double k2, double k3) {
    double tempa, tempb, tempc, tempd, tempe;

    tempa = j3;
    tempb = pow(j3+1,2) - pow(j2-j1,2);
    tempc = pow(j1+j2+1,2) - pow(j3+1,2);
    tempd = pow(j3+1,2) - pow(k1-k2,2);
    tempe = pow(k1+k2+1,2) - pow(j3+1,2);

    return tempa * sqrt(tempb*tempc*tempd*tempe);
}

double* SpecCase6j1(double j1, double j2, double k1, double k2, double k3, double *NormVal) {
    double tempa, tempb, tempval2, Xjval, Zjval, j1pj2, jcent, absm;
	int ngo, i, Nj;
    double *wigvals, jmax, jmin, DVal, SVal;
    
	jmax = ((j1+j2) < (k1+k2)) ? j1+j2 : k1+k2; //maximum valid j3
    tempa = ((j1-j2) > 0) ? j1-j2 : j2-j1;
    tempb = ((k1-k2) > 0) ? k1-k2 : k2-k1;
    jmin = (tempa > tempb) ? tempa : tempb; //minimum valid j3

	Nj = (int) (jmax-jmin + 1);

	DVal = 0;

	wigvals = (double *) malloc(sizeof(double)*Nj);

	//if ((j1 == 0 && m1 == 0) || (j2 == 0 && m2 == 0))
	
	wigvals[Nj-1] = 1; //starting at top
	if (Nj > 1) wigvals[Nj-2] = 0;
	DVal = (2*jmax+1);

	for (i=Nj-3;i>=0;i -= 2) {
		Xjval = Xj6j(j1,j2,i+1+jmin,k1,k2,k3);
		Zjval = Zj6j(j1,j2,i+1+jmin,k1,k2,k3);
		wigvals[i] = - Xjval * wigvals[i+2] / Zjval;
		if (i >= 1) wigvals[i-1] = 0;
		DVal += (2*(i+jmin)+1.)*pow(wigvals[i],2);
	}

	DVal = sqrt((2*k3+1.)*DVal);

	SVal = pow(-1.,j1+j2+k1+k2);

	*NormVal = SVal / DVal;

    return wigvals;
}

double* SpecCase6j2(double j1, double j2, double k1, double k2, double k3, double *NormVal) {
    double *s, *um, DVal, SVal, Yjmin, tempa, tempb, jmax, jmin, jm;
    double jval, tempa2;
    int i, k, Nj, a, Nj2;

    jmax = ((j1+j2) < (k1+k2)) ? j1+j2 : k1+k2; //maximum valid j3
    tempa = ((j1-j2) > 0) ? j1-j2 : j2-j1;
    tempb = ((k1-k2) > 0) ? k1-k2 : k2-k1;
    jmin = (tempa > tempb) ? tempa : tempb; //minimum valid j3

    Nj = (int) (jmax-jmin + 1);

    s = (double *) malloc(sizeof(double)*Nj);
    um = (double *) malloc(sizeof(double)*Nj);

    Yjmin = Yj6j(j1,j2,jmin,k1,k2,k3);

    s[0] = -Xj6j(j1,j2,jmin,k1,k2,k3) / Yjmin;
    tempa = s[0];
    i = 0;
	if (Nj == 1) tempa2 = 2; else tempa2 = 0;
    while (tempa < 1) {
        i++;
        s[i] = -Xj6j(j1,j2,jmin+i,k1,k2,k3) / (Yj6j(j1,j2,jmin+i,k1,k2,k3)+Zj6j(j1,j2,jmin+i,k1,k2,k3)*s[i-1]);
        tempa = s[i];
        tempa2 = (tempa > 0) ? tempa : -tempa;
        if (i == Nj-1) tempa2 = 5; //to make sure we don't go out of bounds of our array and hence crashing the program when deleting
    }

    jm = jmin+i;

    um[(int) (jm-jmin)] = 1;

    DVal = (2.*jm+1.);

    for (k=1;k<=(int)(jm-jmin);k++) {
        a = (int) (jm-k-jmin);
        jval = jm-k;
        um[a] = um[a+1]*s[a];
        DVal += (2.*jval+1.)*pow(um[a],2);
    }

    Nj2 = (int) (jmax-jm);
    for (i=1;i<=Nj2;i++) { //redo for non-integer jm's
        a = (int) (jm+i-jmin);
        jval = jm+i;
        tempa = - Yj6j(j1,j2,jval-1,k1,k2,k3)*um[a-1] - Zj6j(j1,j2,jval-1,k1,k2,k3)*um[a-2];
        tempb = Xj6j(j1,j2,jval-1,k1,k2,k3);
        um[a] = tempa / tempb;
        DVal += (2.*jval+1.)*pow(um[a],2);
    }

    SVal = (um[Nj-1] > 0) ? 1 : -1;
    SVal *= tempa * pow(-1.,j1+j2+k1+k2);

    DVal = sqrt((2.*k3+1.)*DVal);
    free(s);

    *NormVal = SVal / DVal;
    return um;
}

double* SpecCase6j3(double j1, double j2, double k1, double k2, double k3, double *NormVal) {
    double *r, *um, DVal, SVal, Yjmax, tempa, tempb, jmax, jmin, jp;
	double jval, tempa2;
	int i, k, Nj, a, Nj2;

	jmax = ((j1+j2) < (k1+k2)) ? j1+j2 : k1+k2; //maximum valid j3
    tempa = ((j1-j2) > 0) ? j1-j2 : j2-j1;
    tempb = ((k1-k2) > 0) ? k1-k2 : k2-k1;
    jmin = (tempa > tempb) ? tempa : tempb; //minimum valid j3

	Nj = (int) (jmax-jmin + 1);

	r = (double *) malloc(sizeof(double)*Nj);
	um = (double *) malloc(sizeof(double)*Nj);

	Yjmax = Yj6j(j1,j2,jmax,k1,k2,k3);
	r[Nj-1] = -Zj6j(j1,j2,jmax,k1,k2,k3) / Yjmax;
	tempa = r[Nj-1];
	i = 0;

	if (Nj == 1) tempa2 = 2; else tempa2 = 0;
	while (tempa2 < 1) {
		i++;
		r[Nj-i-1] = -Zj6j(j1,j2,jmax-i,k1,k2,k3) / (Yj6j(j1,j2,jmax-i,k1,k2,k3)+Xj6j(j1,j2,jmax-i,k1,k2,k3)*r[Nj-i]);
		tempa = r[Nj-i-1];
		tempa2 = (tempa > 0) ? tempa : -tempa;
		if (i == Nj-1) tempa2 = 5;
	}

	jp = jmax-i;

	um[(int) (jp-jmin)] = 1;

	DVal = 2.*jp+1.;

	for (k=1;k<=(int)(jmax-jp);k++) {
		a = (int) (jp+k-jmin);
		jval = jp+k;
		um[a] = um[a-1]*r[a];
		DVal += (2.*jval+1.)*pow(um[a],2);
	}

	Nj2 = (int) (jp-jmin);
	for (i=1;i<=Nj2;i++) {
		a = (int) (jp-i-jmin);
		jval = jp-i;
		tempa = - Yj6j(j1,j2,jval+1,k1,k2,k3)*um[a+1] - Xj6j(j1,j2,jval+1,k1,k2,k3)*um[a+2];
		tempb = Zj6j(j1,j2,jval+1,k1,k2,k3);
		um[a] = tempa / tempb;
		DVal += (2.*jval+1.)*pow(um[a],2);
	}

	SVal = (um[Nj-1] > 0) ? 1 : -1;
	SVal *= pow(-1.,j1+j2+k1+k2);

	DVal = sqrt((2.*k3+1)*DVal);

	free(r);
	*NormVal = SVal / DVal;
	return um;
}

double* wigner6jvect(int tj1, int tj2, int tk1, int tk2, int tk3, double *NormVal) {
    //can do something similar to 3j symbols
    //follow Luscombe and Luban doc for Xj, Yj, and Zj
    //vectorized version which will return all values in the family
    //  (j1 j2 j )
    //  (k1 k2 k3)
	//Uses equation 9 in Mathworld article for normalization

    double j1, j2, j, k1, k2, k3, jmax, jmin, jp, jm, tempa2;
    double Triad[2][3], TriadM, TriadP, tempa, tempb, *r, *s, *um;
    double Yjmin, Yjmax, *WigVal, jval, Sval, D;
    int TriadSum, Nj, i, Nj2, k, a;

    j1 = (double) tj1 / 2.;
    j2 = (double) tj2 / 2.;
    //j3 = (double) tj3 / 2.;
    k1 = (double) tk1 / 2.;
    k2 = (double) tk2 / 2.;
    k3 = (double) tk3 / 2.;

    Triad[0][0] = j1; Triad[0][1] = k2; Triad[0][2] = k3;
    Triad[1][0] = k1; Triad[1][1] = j2; Triad[1][2] = k3;

    //Each triad must satisfy triangle inequalities
    //Each triad must have the sum of its elements be an integer (or 2*sum is even)
    //Not testing two triads with j3, since we will determine and calculate all valid j3s
    for (i=0;i<=1;i++) {
        TriadSum = (int) (2*(Triad[i][0]+Triad[1][1]+Triad[1][2]));
        TriadM = (Triad[i][0] > Triad[i][1]) ? Triad[i][0]-Triad[i][1] : Triad[i][1]-Triad[i][0];
        TriadP = Triad[i][0]+Triad[i][1];
        if ((Triad[i][2] < TriadM) || (Triad[i][2] > TriadP) || ((TriadSum % 2) == 1)) {
            WigVal = (double *) malloc(sizeof(double));
            WigVal[0] = 0;
            return WigVal;
        }
    }
    //end triad testing, if past here should be a valid 6j symbol family

    jmax = ((j1+j2) < (k1+k2)) ? j1+j2 : k1+k2; //maximum valid j3
    tempa = ((j1-j2) > 0) ? j1-j2 : j2-j1;
    tempb = ((k1-k2) > 0) ? k1-k2 : k2-k1;
    jmin = (tempa > tempb) ? tempa : tempb; //minimum valid j3

    Nj = (int) (jmax - jmin + 1);

    Yjmin = Yj6j(j1,j2,jmin,k1,k2,k3);
    Yjmax = Yj6j(j1,j2,jmax,k1,k2,k3);

    //For special cases when either Yjmin or Yjmax is 0
    //maybe I could change around j's and k's to make non-zero Yjs
    //and avoid special cases
	//or just do some as for 3j with different normalization
    if ((Yjmin==0) && (Yjmax==0)) {
        WigVal = SpecCase6j1(j1,j2,k1,k2,k3,NormVal);
        return WigVal;
    } else if (Yjmax==0) {
        WigVal = SpecCase6j2(j1,j2,k1,k2,k3,NormVal);
        return WigVal;
    } else if (Yjmin==0) {
        WigVal = SpecCase6j3(j1,j2,k1,k2,k3,NormVal);
        return WigVal;
    }

	r = (double *) malloc(sizeof(double)*Nj);
    s = (double *) malloc(sizeof(double)*Nj);
    um = (double *) malloc(sizeof(double)*Nj);

    r[Nj-1] = -Zj6j(j1,j2,jmax,k1,k2,k3) / Yjmax;
    tempa = r[Nj-1];
    i = 0;
	if (Nj == 1) tempa2 = 5; else tempa2 = 0;
    while (tempa2 < 1) {
        i++;
        r[Nj-i-1] = -Zj6j(j1,j2,jmax-i,k1,k2,k3) / (Yj6j(j1,j2,jmax-i,k1,k2,k3)+Xj6j(j1,j2,jmax-i,k1,k2,k3)*r[Nj-i]);
        tempa = r[Nj-i-1];
        tempa2 = (tempa > 0) ? tempa : -tempa;
        if (i == Nj-1) tempa2 = 1.5;
    }

    jp = jmax-i;
    jp = (jp > jmin) ? jp : jmin; //is this already taken care of in while statement??? (if (i >=Nj-1) thingy)

    s[0] = -Xj6j(j1,j2,jmin,k1,k2,k3) / Yjmin;
    tempa = s[0];
    i = 0;
	if (Nj == 1) tempa2 = 5; else tempa2 = 0;
    while (tempa2 < 1) {
        i++;
        s[i] = -Xj6j(j1,j2,jmin+i,k1,k2,k3) / (Yj6j(j1,j2,jmin+i,k1,k2,k3)+Zj6j(j1,j2,jmin+i,k1,k2,k3)*s[i-1]);
        tempa = s[i];
        tempa2 = (tempa > 0) ? tempa : -tempa;
        if (i == Nj-1) tempa2 = 1.5;
    }

    jm = jmin+i;
    jm = (jm < jmax) ? jm : jmax;

    um[(int) (jm-jmin)] = 1;

    for (k=1;k<=(int)(jm-jmin);k++) {
        a = (int) (jm-k-jmin);
        um[a] = um[a+1]*s[a];
    }

    Nj2 = (int) (jp-jm);
    for (i=1;i<=Nj2;i++) {
        a = (int) (jm+i-jmin);
        jval = jm+i;
        tempa = - Yj6j(j1,j2,jval-1,k1,k2,k3)*um[a-1] - Zj6j(j1,j2,jval-1,k1,k2,k3)*um[a-2];
        tempb = Xj6j(j1,j2,jval-1,k1,k2,k3);
        um[a] = tempa / tempb;
    }

    for (k=1;k<=(int)(jmax-jp);k++) {
        a = (int) (jp+k-jmin);
        um[a] = um[a-1]*r[a];
    }

    tempa = (um[Nj-1] > 0) ? 1 : -1;
    Sval = tempa * pow(-1.,j1+j2+k1+k2);

    tempa = 0;
    for (i=0;i<Nj;i++) {
        jval = i+jmin;
        tempa += (2.*jval + 1.) * pow(um[i],2);
    }
    D = sqrt((2.*k3+1.)*tempa);


    *NormVal = Sval/D;

    free(r);
    free(s);

    return um;
}

double wigner6j(int tj1, int tj2, int tj3, int tk1, int tk2, int tk3) {

    double Triad[4][3], TriadM, TriadP, tempa, tempb;
    double jmin, *WigVal, WigRet, NormVal, j1, j2, j3, k1, k2, k3;
    int i, TriadSum;

    j1 = (double) tj1 / 2.;
    j2 = (double) tj2 / 2.;
    j3 = (double) tj3 / 2.;
    k1 = (double) tk1 / 2.;
    k2 = (double) tk2 / 2.;
    k3 = (double) tk3 / 2.;

    Triad[0][0] = j1; Triad[0][1] = j2; Triad[0][2] = j3;
    Triad[1][0] = j1; Triad[1][1] = k2; Triad[1][2] = k3;
    Triad[2][0] = k1; Triad[2][1] = j2; Triad[2][2] = k3;
    Triad[3][0] = k1; Triad[3][1] = k2; Triad[3][2] = j3;

    //Each triad must satisfy triangle inequalities
    //Each triad must have the sum of its elements be an integer (or 2*sum is even)
    for (i=0;i<=3;i++) {
        TriadSum = (int) (2*(Triad[i][0]+Triad[1][1]+Triad[1][2]));
        TriadM = (Triad[i][0] > Triad[i][1]) ? Triad[i][0]-Triad[i][1] : Triad[i][1]-Triad[i][0];
        TriadP = Triad[i][0]+Triad[i][1];
        if ((Triad[i][2] < TriadM) || (Triad[i][2] > TriadP) || ((TriadSum % 2) == 1)) {
            return 0;
        }
    }
    //should be good Wigner 6j symbol

    tempa = ((j1-j2) > 0) ? j1-j2 : j2-j1;
    tempb = ((k1-k2) > 0) ? k1-k2 : k2-k1;
    jmin = (tempa > tempb) ? tempa : tempb; //minimum valid j3

    WigVal = wigner6jvect(tj1, tj2, tk1, tk2, tk3, &NormVal);
    WigRet = NormVal*WigVal[(int) (j3-jmin)];
	free(WigVal);

	return WigRet;
}

double racahw(int ta, int tb, int tc, int td, int te, int tf) {
    double WigVal, powval;

    powval = (double) (ta+tb+tc+td) / 2.;

    WigVal = wigner6j(ta,tb,tc,td,te,tf);

    return pow(-1.,powval)*WigVal;
}
