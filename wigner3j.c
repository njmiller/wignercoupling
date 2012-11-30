//Trying to follow Brock's algorithm in calculation of Wigner 3j symbols
//Should reduce numerical overflow
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double Zj(double j, double j1, double j2, double m1, double m2) {
	double tempa, tempb, tempc, tempd;

	tempa = (j+1);
	tempb = pow(j,2) - pow(j1-j2,2);
	tempc = pow(j1+j2+1,2) - pow(j,2);
	tempd = pow(j,2) - pow(m1+m2,2);

	return tempa*sqrt(tempb*tempc*tempd);
}

double Yj(double j, double j1, double j2, double m1, double m2) {
	double tempa, tempb, tempc;

	tempa = (m1+m2)*(j1*(j1+1)-j2*(j2+1));
	tempb = (m1-m2)*j*(j+1);
	tempc = (2*j+1)*(tempa-tempb);
    
	return tempc;
}

double Xj(double j, double j1, double j2, double m1, double m2) {
	double tempa, tempb, tempc;

	tempa = pow(j+1,2.) - pow(j1-j2,2);
	tempb = pow(j1+j2+1,2.) - pow(j+1,2);
	tempc = pow(j+1,2.) - pow(m1+m2,2);

	return j*sqrt(tempa*tempb*tempc);
}


double* SpecCase1(double j1, double j2, double m1, double m2) {
    //uses three term recursion and normalization to do all
    //this is the vectorized version

    double tempa, tempb, tempval2, Xjval, Zjval, j1pj2, jcent, absm;
    int ngo, i, Nj;
    double *wigvals, jmax, jmin, DVal, SVal, norm;

    jmax = j1+j2;
    tempa = (j1-j2 > 0) ? j1-j2 : j2-j1;

    tempb = (m1+m2 > 0) ? m1+m2 : -m1-m2;
    jmin = (tempa > tempb) ? tempa : tempb;

    Nj = (int) (jmax-jmin + 1);

    DVal = 0;

    wigvals = (double *) malloc(sizeof(double)*Nj);

    if ((j1 == 0 && m1 == 0) || (j2 == 0 && m2 == 0)) {
        //since m1 = 0 and/or m2 = 0 then m1+m2 is value of non-zero term
        wigvals[0] = pow(-1.,jmin+m1+m2) / sqrt(2.*jmin+1.);
        return wigvals;
    }

    wigvals[Nj-1] = 1; //starting at top since that always has j1+j2+j3 as even (2*j1+2*j2)
    if (Nj > 1) wigvals[Nj-2] = 0; //only do if exists in sum
    DVal = (2*jmax+1);

	//3-term recursion calculation. In this special case every other term
	//is 0 so it is a 2-term recursion
    for (i=Nj-3;i>=0;i -= 2) {
        Xjval = Xj(i+1+jmin,j1,j2,m1,m2);
        Zjval = Zj(i+1+jmin,j1,j2,m1,m2);
        wigvals[i] = - Xjval * wigvals[i+2] / Zjval;
        if (i >= 1) wigvals[i-1] = 0;
        DVal += (2*(i+jmin)+1.)*pow(wigvals[i],2);
    }

    DVal = sqrt(DVal);

    SVal = pow(-1.,j1-j2+m1+m2);

    norm = SVal / DVal;
	for (i=0;i<Nj;i++) {
		wigvals[i] *= norm;
	}

    return wigvals;

}

double* SpecCase2(double j1, double j2, double m1, double m2) {
    //This is case where Yjmax = 0. Upper end recursion can't get started
    //use two term lower end recursion, then 3-term recursion for all other terms

    double *s, *um, DVal, SVal, Yjmin, tempa, tempb, jmax, jmin, jm;
    double jval, tempa2, norm;
    int i, k, Nj, a, Nj2;

    jmax = j1+j2;
    tempa = (j1-j2 > 0) ? j1-j2 : j2-j1;

    tempb = (m1+m2 > 0) ? m1+m2 : -m1-m2;
    jmin = (tempa > tempb) ? tempa : tempb;

    Nj = (int) (jmax-jmin + 1);

    s = (double *) malloc(sizeof(double)*Nj);
    um = (double *) malloc(sizeof(double)*Nj);

    Yjmin = Yj(jmin,j1,j2,m1,m2);

	//loop is 2-term recursion from lower-end
    s[0] = -Xj(jmin,j1,j2,m1,m2) / Yjmin;
    tempa = s[0];
    i = 0;
    if (Nj == 1) tempa2 = 2; else tempa2 = 0;
    while (tempa2 < 1) {
        i++;
        s[i] = -Xj(jmin+i,j1,j2,m1,m2) / (Yj(jmin+i,j1,j2,m1,m2)+Zj(jmin+i,j1,j2,m1,m2)*s[i-1]);
        tempa = s[i];
        tempa2 = (tempa > 0) ? tempa : -tempa;
        if (i == Nj-1) tempa2 = 5;
    }

    jm = jmin+i;

    um[(int) (jm-jmin)] = 1;

    DVal = (2.*jm+1.);
	
	//calculate values based on recursion ratios and unnormalized choice
    for (k=1;k<=(int)(jm-jmin);k++) {
        a = (int) (jm-k-jmin);
        jval = jm-k;
        um[a] = um[a+1]*s[a];
        DVal += (2.*jval+1.)*pow(um[a],2); //part of normalization
    }

	//3-term recursion, goes from end of 2-term recursion calculations
    Nj2 = (int) (jmax-jm);
    for (i=1;i<=Nj2;i++) { //redo for non-integer jm's
        a = (int) (jm+i-jmin);
        jval = jm+i;
        tempa = - Yj(jval-1,j1,j2,m1,m2)*um[a-1] - Zj(jval-1,j1,j2,m1,m2)*um[a-2];
        tempb = Xj(jval-1,j1,j2,m1,m2);
        um[a] = tempa / tempb;
        DVal += (2.*jval+1.)*pow(um[a],2); //part of normalization
    }

	//normalization
    SVal = (um[Nj-1] > 0) ? 1 : -1;
    SVal *= pow(-1.,j1-j2+m1+m2);

    DVal = sqrt(DVal);

	//free memory
    free(s);

    norm = SVal / DVal;
	for (i=0;i<Nj;i++) {
		um[i] *= norm;
	}

    return um;
}

double* SpecCase3(double j1, double j2, double m1, double m2) {
    //This is the one where Yjmin = 0. Do Same thing as special case 2 except
    //use the higher end instead of lower end

    double *r, *um, DVal, SVal, Yjmax, tempa, tempb, jmax, jmin, jp;
    double jval, tempa2, norm;
    int i, k, Nj, a, Nj2;

    jmax = j1+j2;
    tempa = (j1-j2 > 0) ? j1-j2 : j2-j1;

    tempb = (m1+m2 > 0) ? m1+m2 : -m1-m2;
    jmin = (tempa > tempb) ? tempa : tempb;

    Nj = (int) (jmax-jmin + 1);

    r = (double *) malloc(sizeof(double)*Nj);
    um = (double *) malloc(sizeof(double)*Nj);

    Yjmax = Yj(jmax,j1,j2,m1,m2);

    r[Nj-1] = -Zj(jmax,j1,j2,m1,m2) / Yjmax;
    tempa = r[Nj-1];
    i = 0;

    if (Nj == 1) tempa2 = 2; else tempa2 = 0;
    while (tempa2 < 1) {
        i++;
        r[Nj-i-1] = -Zj(jmax-i,j1,j2,m1,m2) / (Yj(jmax-i,j1,j2,m1,m2)+Xj(jmax-i,j1,j2,m1,m2)*r[Nj-i]);
        tempa = r[Nj-i-1];
        tempa2 = (tempa > 0) ? tempa : -tempa;
        if (Nj-i-1 == 0) tempa2 = 5;
    }

    jp = jmax-i;

    um[(int) (jp-jmin)] = 1;

    DVal = (2.*jp+1.);

    for (k=1;k<=(int)(jmax-jp);k++) {
        a = (int) (jp+k-jmin);
        jval = jp+k;
        um[a] = um[a-1]*r[a];
        DVal += (2.*jval+1.)*pow(um[a],2);
    }

    Nj2 = (int) (jp-jmin);
    for (i=1;i<=Nj2;i++) { //redo for non-integer jm's
        a = (int) (jp-i-jmin);
        jval = jp-i;
        tempa = - Yj(jval+1,j1,j2,m1,m2)*um[a+1] - Xj(jval+1,j1,j2,m1,m2)*um[a+2];
        tempb = Zj(jval+1,j1,j2,m1,m2);
        um[a] = tempa / tempb;
        DVal += (2.*jval+1.)*pow(um[a],2);
    }

    SVal = (um[Nj-1] > 0) ? 1 : -1;
    SVal *= pow(-1.,j1-j2+m1+m2);

    DVal = sqrt(DVal);

    free(r);

    norm = SVal / DVal;
	for (i=0;i<Nj;i++) {
		um[i] *= norm;
	}
    return um;

}

double* wigner3jvect(int tj1, int tj2, int tm1, int tm2) {
    //returns a vector with Wigner 3j values for the family (j1 j2    j  )
    //                                                      (m1 m2 -m1-m2)

    double j1, j2, j3, m1, m2, m3, jmin, jmax, tempa, tempb, Yjmax, Yjmin;
    double *r, *s, *um, jp, jm, jval, *WigVal, D, Sval, tempa2, norm;
    int Nj, i, j, k, Nj2, a, t1, t2;

    j1 = (double)tj1 / 2.;
    j2 = (double)tj2 / 2.;
    m1 = (double)tm1 / 2.;
    m2 = (double)tm2 / 2.;
    m3 = -(m1+m2);

	t1 = (tm1 > 0) ? tm1 : -tm1;
	t2 = (tm2 > 0) ? tm2 : -tm2;
    if (m1 > j1 || m1 < -j1 || m2 > j2 || m2 < -j2 || tj1%2 != t1%2 || tj2%2 != t2%2) {
        WigVal = (double *) malloc(sizeof(double));
        WigVal[0] = 0;
        return WigVal;
    }

    jmax = j1+j2;
    tempa = (j1-j2 > 0) ? j1-j2 : j2-j1;

    tempb = (m1+m2 > 0) ? m1+m2 : -m1-m2;
    jmin = (tempa > tempb) ? tempa : tempb;

    Nj = (int) (jmax-jmin + 1);

    Yjmax = Yj(jmax,j1,j2,m1,m2);
    Yjmin = Yj(jmin,j1,j2,m1,m2);

    //check to see if code already compensates for Cases 2 AND 3. It might since Infinity > 1
    //and would set jp or jm = jmax or jmin as needed. (no since I need good values
	//Special Cases require something different as we can't get the recursion started at one 
	//(or both) ends
    if (Yjmax == 0 && Yjmin == 0) {
        WigVal = SpecCase1(j1,j2,m1,m2);
        return WigVal;
    } else if (Yjmax == 0) {
        WigVal = SpecCase2(j1,j2,m1,m2);
        return WigVal;
    } else if (Yjmin == 0) {
        WigVal = SpecCase3(j1,j2,m1,m2);
        return WigVal;
    }

	r = (double *) malloc(sizeof(double)*Nj);
    s = (double *) malloc(sizeof(double)*Nj);
    um = (double *) malloc(sizeof(double)*Nj);

	//Start at the upper end and do two-term recursion
	//End when the abs(ratios) start to become > 1
	r[Nj-1] = -Zj(jmax,j1,j2,m1,m2) / Yjmax;
    tempa = r[Nj-1];
    i = 0;
	if (Nj == 1) tempa2 = 5; else tempa2 = 0;
    while (tempa2 < 1) {
        i++;
        r[Nj-i-1] = -Zj(jmax-i,j1,j2,m1,m2) / (Yj(jmax-i,j1,j2,m1,m2)+Xj(jmax-i,j1,j2,m1,m2)*r[Nj-i]);
        tempa = r[Nj-i-1];
        tempa2 = (tempa > 0) ? tempa : -tempa;
        if (i == Nj-1) tempa2 = 1.5; //in case of cases where never get > 1
    }
	
    jp = jmax-i;
    jp = (jp > jmin) ? jp : jmin;

	//Start at the lower end and do two-term recursion (recursion of Wigner symbol ratios)
	//Ends when the abs(ratios) start to become > 1. 
    s[0] = -Xj(jmin,j1,j2,m1,m2) / Yjmin;
    tempa = s[0];
    i = 0;
    if (Nj == 1) tempa2 = 5; else tempa2 = 0;
    while (tempa2 < 1) {
        i++;
        s[i] = -Xj(jmin+i,j1,j2,m1,m2) / (Yj(jmin+i,j1,j2,m1,m2)+Zj(jmin+i,j1,j2,m1,m2)*s[i-1]);
        tempa = s[i];
        tempa2 = (tempa > 0) ? tempa : -tempa;
        if (i == Nj-1) tempa2 = 1.5; //in case of cases where never get > 1
    }

    jm = jmin+i;
    jm = (jm < jmax) ? jm : jmax;

	//starting to calculate unnormalized wigner symbols, using the ratios found
	//in the earlier loops and the normal 3-term recursion. 
	//We just pick a certain position to be 1 and spread out from there
    um[(int) (jm-jmin)] = 1;

	//lower end with ratios calculated earlier
    for (k=1;k<=(int)(jm-jmin);k++) {
        a = (int) (jm-k-jmin);
        um[a] = um[a+1]*s[a];
    }

	//normal 3-term recursion
    Nj2 = (int) (jp-jm);
    for (i=1;i<=Nj2;i++) {
        a = (int) (jm+i-jmin);
        jval = jm+i;
        tempa = - Yj(jval-1,j1,j2,m1,m2)*um[a-1] - Zj(jval-1,j1,j2,m1,m2)*um[a-2];
        tempb = Xj(jval-1,j1,j2,m1,m2);
        um[a] = tempa / tempb;
    }

	//upper end with ratios calculated earlier
    for (k=1;k<=(int)(jmax-jp);k++) {
        a = (int) (jp+k-jmin);
        um[a] = um[a-1]*r[a];
        
    }

	//calculation of the normalization constant
    tempa = (um[Nj-1] > 0) ? 1 : -1;
    Sval = tempa * pow(-1.,j1-j2+m1+m2);

    tempa = 0;
    for (i=0;i<Nj;i++) {
        jval = i+jmin;
        tempa += (2.*jval + 1.) * pow(um[i],2);
    }
    D = sqrt(tempa);

    free(r);
    free(s);

    norm = Sval / D;
	
	for (i=0;i<Nj;i++) {
		um[i] *= norm;
	}

    return um;
}

double wigner3j(int tj1, int tj2, int tj3, int tm1, int tm2, int tm3) {

    double *WigVal, WigRet, jmin, tempa, tempb, j1, j2, j3, m1, m2, m3, jmax;
    int temp3;

    j1 = (double) tj1 / 2.;
    j2 = (double) tj2 / 2.;
    j3 = (double) tj3 / 2.;
    m1 = (double) tm1 / 2.;
    m2 = (double) tm2 / 2.;
    m3 = (double) tm3 / 2.;

    //Do all the testing for the Wigner 3j symbols here
	temp3 = (tm1 > 0) ? tm1 : -tm1;
	if ((tj1%2) != (temp3%2)) return 0;
	temp3 = (tm2 > 0) ? tm2 : -tm2;
	if ((tj2%2) != (temp3%2)) return 0;
	temp3 = (tm3 > 0) ? tm3 : -tm3;
	if ((tj3%2) != (temp3%2)) return 0;
    if (m1+m2 != -m3) return 0;
    if (m1 > j1 || m1 < -j1) return 0;
    if (m2 > j2 || m2 < -j2) return 0;
    if (m3 > j3 || m3 < -j3) return 0;

    jmax = j1+j2;

    tempa = (j1-j2 > 0) ? j1-j2 : j2-j1;
    if (j3 < tempa || j3 > jmax) return 0;

    tempb = (m1+m2 > 0) ? m1+m2 : -m1-m2;
    jmin = (tempa > tempb) ? tempa : tempb;

    WigVal = wigner3jvect(tj1, tj2, tm1, tm2);
    temp3 = (int) (j3-jmin);

	WigRet = WigVal[temp3];
	free(WigVal);

    return WigRet;
}

double clebschgordon(int tj1, int tj2, int tj3, int tm1, int tm2, int tm3) {
    double WigVal, j1, j2, m3;

    j1 = (double) tj1 / 2.;
    j2 = (double) tj2 / 2.;
    m3 = (double) tm3 / 2.;

    WigVal = wigner3j(tj1,tj2,tj3,tm1,tm2,tm3);

    return pow(-1.,m3+j1-j2) * sqrt((double)tj3+1.) * WigVal;
}


double racahv(int tj1, int tj2, int tj3, int tm1, int tm2, int tm3) {
    double WigVal, j3, m3;

    j3 = (double) tj3 / 2.;
    m3 = (double) tm3 / 2.;

    WigVal = wigner3j(tj2,tj2,tj3,tm1,tm2,tm3);

    return pow(-1.,j3+m3)*WigVal;
}


