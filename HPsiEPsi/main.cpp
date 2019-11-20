#include <stdio.h>
#include <iostream>
//#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_odeiv2.h>
#include <math.h>
#include <vector>
using namespace std;

class
generalParametersClass{
public:
    generalParametersClass(int N_, double RMIN_, double RMAX_, int LMAX_, double EPSILON_);
    int     N;
    double  RMAX;
    double  RMIN;
    double  H;
    int     LMAX;
    double  EPSILON;
    vector<double> R;
};

generalParametersClass::generalParametersClass(int N_, double RMIN_, double RMAX_, int LMAX_, double EPSILON_){
    N         = N_;
    RMIN      = RMIN_;
    RMAX      = RMAX_;
    LMAX      = LMAX_;
    H         =(RMAX -RMIN) /(N-1);
    EPSILON   = EPSILON_;
    for (int i =0; i <N; i++) {
        R.push_back(H*i+0.001);
    }
}



class
partitionClass {
public:
    partitionClass (bool energyLab, double energy, double m1, double m2, int z1, int z2, int N);
    double  hBarSquared =41.801651165221026;
    bool    energyLab;
    double  energy;
    double  m1, m2;
    int     z1, z2;
    double  m, twoMuDevidedByhBarSquared;
    double  k, kSquared;
    double  n;
    vector<double>  waveFunctionRe;
    vector<double>  waveFunctionIm;
 //   double  interactionKernel (double R);
//    double  getAssimpFunctionAtV0 (double V0);
    bool    boundStateQ();
    void    setWoodsSaxonParameters(double V0_r_, double R0_r, double a0_r_,
      double V0_i_ , double R0_i_, double a0_i_ );
    double  kernelFunction(double R);
    double  V0_r, R0_r, a0_r;
    double  V0_i, R0_i, a0_i;
    int l=0;
private:

};

partitionClass::partitionClass
 (bool energyLabBool,  double energyDouble,
  double m1Double,     double m2Double,
  int z1Int,           int z2Int, int N){
    energyLab   =energyLabBool;
    energy      =energyDouble;
    m1          =m1Double;
    m2          =m2Double;
    z1          =z1Int;
    z2          =z2Int;
    if(energyLab)   energy      = energy *m2 /(m1 +m2);
    m                           = m1 *m2 /(m1 +m2);
    //2*amu/hbar^2 = 0.0478450
    twoMuDevidedByhBarSquared   = 2 *m /hBarSquared;
    kSquared                    = 2 *m *energy /hBarSquared;
    k                           = sqrt(abs(kSquared));
    n                           = z1 *z2 *1.43997 /hBarSquared /k;
    waveFunctionRe.resize(N, 0.0);
    waveFunctionIm.resize(N, 0.0);

}



/*
double
partitionClass::interactionPotential(char potentialChar, double R_){
    double var;
    switch(potentialChar){
         case "WS":

     }
    return var;
}
*/

double
partitionClass::kernelFunction(double R_){
   return  kSquared +twoMuDevidedByhBarSquared *V0_r /( 1 +exp( (R_ - R0_r) /a0_r))
    - l *(l +1) /R_ /R_;
}

void
partitionClass::setWoodsSaxonParameters(double V0_r_, double R0_r_, double a0_r_,
double V0_i_ =0.0, double R0_i_ =1.250, double a0_i_ =0.650){
    V0_r  =V0_r_;
    R0_r  =R0_r_;
    a0_r  =a0_r_;
    V0_i  =V0_i_;
    R0_i  =R0_i_;
    a0_i  =a0_i_;
}

bool
partitionClass::boundStateQ(){
    if (energy<0) return true;
    else    return false;
}

void
applyFiniteDiffMethodFor (partitionClass &partition, generalParametersClass &parameter){
    partition.waveFunctionRe[0]=0.0;
    partition.waveFunctionRe[1]=0.1;
    partition.waveFunctionIm[0]=0.0;
    partition.waveFunctionIm[1]=0.1;
    for (int i=2; i < parameter.N; i++) {
        partition.waveFunctionRe[i] =
        (2 -parameter.H *parameter.H
        *partition.kernelFunction(parameter.R[i-1]))
        *partition.waveFunctionRe[i-1] -partition.waveFunctionRe[i-2]
        -parameter.H *parameter.H *partition.twoMuDevidedByhBarSquared
        *partition.V0_i /(1 +exp(parameter.R[i-1] - partition.R0_i) /partition.a0_i)
        *partition.waveFunctionIm[i-1];

        partition.waveFunctionIm[i] =
        (2 -parameter.H *parameter.H
        *partition.kernelFunction(parameter.R[i-1]))
        *partition.waveFunctionIm[i-1] -partition.waveFunctionIm[i-2]
        -parameter.H *parameter.H *partition.twoMuDevidedByhBarSquared
        *partition.V0_i /(1 +exp(parameter.R[i-1] - partition.R0_i) /partition.a0_i)
        *partition.waveFunctionRe[i-1];
    }
}

void
applyNumerovMethodFor (partitionClass &partition, generalParametersClass &parameter){
    double gamma = parameter.H *parameter.H /12.0;
    partition.waveFunctionRe[0]=0.0;
    partition.waveFunctionRe[1]=0.1;
    for (int i=2; i < parameter.N; i++) {
        partition.waveFunctionRe[i] =
        1. /(  1. +gamma *partition.kernelFunction( parameter.R[i] )  )
        *(
          (2. -10.0 *gamma *partition.kernelFunction(parameter.R[i-1]) )
          *partition.waveFunctionRe[i-1]
          - (1. +gamma *partition.kernelFunction(parameter.R[i-2]))
          *partition.waveFunctionRe[i-2]
          );
    }
}

void
findRoot (double (*f)(double), double argument){
    double rightTest= argument +argument *0.5;
    double leftTest= argument +argument *0.5;
    double middleTest=leftTest;
    int rootCounter=0;
    while ( true ) {
        middleTest =(rightTest +leftTest) /2;
 /*       partition->V0 =middleTest;
        applyNumerovMethodFor(&partition, y, &partition);
        if ( y[N-1] < parameter.EPSILON && y[N-1] > 0.0) break;
        if( y[N-1] > 0) leftTest =middleTest;
        else            rightTest= middleTest;*/
        rootCounter++;
        printf("%d \t %.9f \n",rootCounter, middleTest);
    }

}


int
main (void) {
    generalParametersClass generalParameters(201, 0.001, 40.001, 40, 1.E-6);
    partitionClass firstPartition(false, -2.225, 1.0, 1.0, 0, 0, 201);
    firstPartition.setWoodsSaxonParameters(63.368593, 1.25, 0.65, 30.0, 1.25, 0.65);
    applyNumerovMethodFor(firstPartition, generalParameters);

    partitionClass secondPartition(false, -2.225, 1.0, 1.0, 0, 0, 201);
    secondPartition.setWoodsSaxonParameters(63.368593, 1.25, 0.65);
    applyFiniteDiffMethodFor(secondPartition, generalParameters);

    for (int i=0; i<generalParameters.N; i++) {
        printf("%.3f\t%.5f\t%.5f\n", generalParameters.R[i], firstPartition.waveFunctionIm[i], secondPartition.waveFunctionRe[i]);
    }
   // for (auto i = generalParameters.R.begin(); i != generalParameters.R.end(); ++i)
    //cout << *i << "\n";
    return 0;
}
