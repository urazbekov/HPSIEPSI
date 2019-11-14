#include <stdio.h>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <math.h>
#include <vector>
using namespace std;

struct
generalParameters{
    int     N;
    double  RMAX;
    double  RMIN;
    double  H;
    int     LMAX;
    double  EPSILON;
};


class
boundPartiton {
public:
    boundPartiton (bool energyLab, double energy,
                             double m1, double m2,
                             int z1, int z2);
    double  hBarSquared =41.801651165221026;
    bool    energyLab;
    double  energy;
    double  m1, m2;
    int     z1, z2;
    double  m, twoMuDevidedByhBarSquared;
    double  k, kSquared;
    double  n;
    void    initializeData();
    double  interactionKernel (double R);
    double  V0;
    double  getAssimpFunctionAtV0 (double V0);
    double  waveFunction[];
};

boundPartiton::boundPartiton
 (bool energyLabBool,  double energyDouble,
  double m1Double,     double m2Double,
  int z1Int,           int z2Int){
    energyLab   =energyLabBool;
    energy      =energyDouble;
    m1          =m1Double;
    m2          =m2Double;
    z1          =z1Int;
    z2          =z2Int;
}

void boundPartiton::initializeData(){
    if(energyLab)   energy      = energy *m2 /(m1 +m2);
    m                           = m1 *m2 /(m1 +m2);
    //2*amu/hbar^2 = 0.0478450
    twoMuDevidedByhBarSquared   = 2 *m /hBarSquared;
    kSquared                    = 2 *m *energy /hBarSquared;
    k                           = sqrt(abs(kSquared));
    n                           = z1 *z2 *1.43997 /hBarSquared /k;
}

double boundPartiton::interactionKernel(double R){
    return kSquared +twoMuDevidedByhBarSquared *V0
    /(1.0 +exp((R -1.25)/0.65 ));
    
}

void
applyFiniteDiffMethodFor (boundPartiton partition,
                           double *y,
                           generalParameters parameter){
    y[0]=0.0;
    y[1]=0.1;
    for (int i=2; i < parameter.N; i++) {
        y[i] =
        (2 -parameter.H *parameter.H
        *partition.interactionKernel(parameter.H *(i-1)))
        *y[i-1] -y[i-2];
    }
}


void
applyNumerovMethodFor (boundPartiton *partition,
                       double *y,
                       generalParameters *parameter){
    double gamma = parameter->H *parameter->H /12.0;
    y[0]=0.0;
    y[1]=0.1;
    for (int i=2; i < parameter->N; i++) {
        y[i] =
        1. /(1. +gamma *partition->interactionKernel(parameter->H *(i)))
        *((2. -10. *gamma *partition->interactionKernel(parameter->H *(i-1)) ) *y[i-1]
          - (1. +gamma *partition->interactionKernel(parameter->H *(i-1.))) *y[i-2]);
    }
}

double boundPartiton::getAssimpFunctionAtV0(double testV0){
    applyNumerovMethodFor(<#boundPartiton *partition#>, <#double *y#>, <#generalParameters *parameter#>)
}

void
findRoot (double (*f)(double), double argument){
    double rightTest= argument +argument *0.5;
    double leftTest= argument +argument *0.5;
    double middleTest=leftTest;
    int rootCounter=0;
    while ( true ) {
        middleTest =(rightTest +leftTest) /2;
        partition->V0 =middleTest;
        applyNumerovMethodFor(&partition, y, &partition);
        if ( y[N-1] < parameter.EPSILON && y[N-1] > 0.0) break;
        if( y[N-1] > 0) leftTest =middleTest;
        else            rightTest= middleTest;
        rootCounter++;
        printf("%d \t %.9f \n",rootCounter, middleTest);
    }
    
}


int
main (void) {
    const int N =201;
    generalParameters parameter;
    boundPartiton firstPartition(false, -2.225, 1.0, 1.0, 0, 0);
    firstPartition.initializeData();
    double y[N];
    parameter.N         = N;
    parameter.RMIN      = 0.000;
    parameter.RMAX      =40.000;
    parameter.LMAX      =40;
    parameter.H         =(parameter.RMAX -parameter.RMIN)
                        /(parameter.N-1);
    parameter.EPSILON   =1.E-6;
    firstPartition.V0=63.19828;
    
   
//    applyNumerovMethodFor(&firstPartition, y, &setOfParametersODE);
    for (int i=0; i<N; i++) {
        printf("%f \t %.5e \n", parameter.H *i, y[i]);
    }
    printf("%d\t%.7f \n", rootCounter, middleTest);
    return 0;
}
