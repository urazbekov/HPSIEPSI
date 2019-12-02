#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <string>
#include <gsl/gsl_sf_coulomb.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
using namespace std;

class
        generalParametersClass {
public:
generalParametersClass(int N_, double RMIN_, double RMAX_,
                       int LMIN_, int LMAX_, double EPSILON_);
int N;
double RMAX;
double RMIN;
double H;
int LMIN;
int LMAX;
double EPSILON;
vector<double> R;
};

generalParametersClass::generalParametersClass(int N_, double RMIN_, double RMAX_,
                                               int LMIN_, int LMAX_, double EPSILON_){
        N         = N_;
        RMIN      = RMIN_;
        RMAX      = RMAX_;
        LMIN      = LMIN_;
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
partitionClass ( const char*  name, bool energyLab, double energy,
                 double m1, double m2, int z1, int z2, generalParametersClass &parameters);
string name;
double hBarSquared;
bool energyLab;
double energy;
double m1, m2;
int z1, z2;
double m, twoMuDevidedByhBarSquared;
double k, kSquared;
double n;
vector<double>  coulombF;
vector<double>  coulombG;
vector<double>  coulombFPrime;
vector<double>  coulombGPrime;
vector<double> sMatrixRe;
vector<double> sMatrixIm;
vector<double> rowWaveFunctionRe;
vector<double> rowWaveFunctionIm;
vector<double> rowWaveFunctionPrimeRe;
vector<double> rowWaveFunctionPrimeIm;
vector< vector<double> > waveFunctionRe;
vector< vector<double> > waveFunctionIm;
vector< vector<double> > waveFunctionPrimeRe;
vector< vector<double> > waveFunctionPrimeIm;
bool    boundStateQ();
void    setWoodsSaxonParameters(double V0_r_, double R0_r, double a0_r_,
                                double V0_i_, double R0_i_, double a0_i_, double Rc_ );
double getWoodSaxonVolumeIm(double R);
double getWoodSaxonVolumeRe(double R);
double getCoulombPotential(double R);
double  kernelFunction(double R);
double V0_r, R0_r, a0_r;
double V0_i, R0_i, a0_i;
double Rc;
int l;
int setCoulombWaveFunctions(generalParametersClass &parameter);
void printWaveFunction(generalParametersClass &parameter);
void getSmatrix(generalParametersClass &parameter);

private:

};



partitionClass::partitionClass
        (const char*   nameChar,
        bool energyLabBool,
        double energyDouble,
        double m1Double,
        double m2Double,
        int z1Int,
        int z2Int,
        generalParametersClass &parameters)
// using Initializer List
        : rowWaveFunctionRe(parameters.N,0.0),
        rowWaveFunctionIm(parameters.N,0.0),
        rowWaveFunctionPrimeRe(parameters.N, 0.0),
        rowWaveFunctionPrimeIm(parameters.N, 0.0),
        waveFunctionRe(parameters.LMAX+1,rowWaveFunctionRe),
        waveFunctionIm(parameters.LMAX+1,rowWaveFunctionIm),
        waveFunctionPrimeRe(parameters.LMAX +1, rowWaveFunctionPrimeRe),
        waveFunctionPrimeIm(parameters.LMAX +1, rowWaveFunctionPrimeIm)
{
        hBarSquared =41.801651165221026;
        name =nameChar;
        energyLab   =energyLabBool;
        energy      =energyDouble;
        m1          =m1Double;
        m2          =m2Double;
        z1          =z1Int;
        z2          =z2Int;
        if(energyLab) energy      = energy *m2 /(m1 +m2);
        m                           = m1 *m2 /(m1 +m2);
        //2*amu/hbar^2 = 0.0478450
        twoMuDevidedByhBarSquared   = 2 *m /hBarSquared;
        kSquared                    = 2 *m *energy /hBarSquared;
        k                           = sqrt(abs(kSquared));
        n                           = z1 *z2 *1.43997 /hBarSquared /k;
        coulombF.resize(parameters.LMAX+1, 0.0);
        coulombG.resize(parameters.LMAX+1, 0.0);
        coulombFPrime.resize(parameters.LMAX+1, 0.0);
        coulombGPrime.resize(parameters.LMAX+1, 0.0);
        sMatrixRe.resize(parameters.LMAX+1, 0.0);
        sMatrixIm.resize(parameters.LMAX+1, 0.0);
}

double
partitionClass::getWoodSaxonVolumeRe(double R){
        return V0_r /exp(1+ (R -R0_r)/a0_r);
}


double
partitionClass::getWoodSaxonVolumeIm(double R){
        return V0_i /exp(1+ (R -R0_i)/a0_i);
}

double
partitionClass::getCoulombPotential(double R){
        if (R > Rc)
                return z1 *z2*1.43997 /R;
        else
                return z1 *z2*1.43997 *(3 -R *R /Rc /Rc) /Rc *0.5;
}





void
partitionClass::getSmatrix(generalParametersClass &parameter) {
        setCoulombWaveFunctions(parameter);
        double a, b, c, d;
        double A, B, C, D;
        double AsquaredPlusBsquared;
        for (int l_ = parameter.LMIN; l_ < parameter.LMAX; l_++) {
                a = waveFunctionPrimeRe[l_][parameter.N-1] *coulombF[l_]
                    - waveFunctionRe[l_][parameter.N-1] *coulombFPrime[l_];
                b= waveFunctionIm[l_][parameter.N-1] *coulombGPrime[l_]
                   - waveFunctionPrimeIm[l_][parameter.N-1] *coulombG[l_];
                c= waveFunctionPrimeIm[l_][parameter.N-1] *coulombF[l_]
                   - waveFunctionIm[l_][parameter.N-1] *coulombFPrime[l_];
                d= waveFunctionRe[l_][parameter.N-1] *coulombGPrime[l_]
                   - waveFunctionPrimeRe[l_][parameter.N-1] *coulombG[l_];
                A= b- a;
                B= -c -d;
                C= a +b;
                D= c -d;
                AsquaredPlusBsquared = A *A + B *B;
                sMatrixRe[l_]= (A *C + B *D) /AsquaredPlusBsquared;
                sMatrixIm[l_]= (A *D - B *C) /AsquaredPlusBsquared;
        }


}


void
partitionClass::printWaveFunction(generalParametersClass &parameter){
        FILE * pFile;
        string filename =name+".wfn";
        pFile = fopen (filename.data(),"w");
        fprintf(pFile, "# Wave Function For the ");
        fprintf(pFile, "%s",name.data());
        fprintf(pFile, " partition \n");

        if (boundStateQ()) {
                fprintf(pFile, "# R [fm]     Psi \n");
                for (int i=0; i< parameter.N; i++) {
                        //                      fprintf (pFile, "% .3e\t% .5e\n",parameter.R[i],waveFunctionRe[i]);
                }
        }
        else{
                fprintf(pFile, "# R [fm]       ");
                for (int l_ = parameter.LMIN; l_ <= parameter.LMAX; l_++) {
                        fprintf(pFile, "PsiRe[%d]      PsiIm[%d]        ",l_,l_);
                }
                fprintf(pFile, "\n");
                for (int Ri =0; Ri < parameter.N; Ri++) {
                        fprintf(pFile,"% .3e\t\t",parameter.R[Ri]);

                        for(int l_=parameter.LMIN; l_<=parameter.LMAX; l_++) {
                                fprintf (pFile, "% .5e\t% .5e\t\t",
                                         waveFunctionRe[l_][Ri], waveFunctionIm[l_][Ri]);
                        }
                        fprintf(pFile, "\n");
                }


        }
        fclose (pFile);
}

double
partitionClass::kernelFunction(double R){
        return kSquared - l *(l +1) /R /R - twoMuDevidedByhBarSquared
               *(  getWoodSaxonVolumeRe(R) +getCoulombPotential(R));
}

void
partitionClass::setWoodsSaxonParameters(double V0_r_, double R0_r_, double a0_r_,
                                        double V0_i_ =0.0, double R0_i_ =1.250, double a0_i_ =0.650,
                                        double Rc_ =1.25){
        V0_r  =V0_r_;
        R0_r  =R0_r_;
        a0_r  =a0_r_;
        V0_i  =V0_i_;
        R0_i  =R0_i_;
        a0_i  =a0_i_;
        Rc    =Rc_;
}

bool
partitionClass::boundStateQ(){
        if (energy<0) return true;
        else return false;
}

void
applyFiniteDiffMethodFor (partitionClass &partition, generalParametersClass &parameter){
        double tempKernelFunction;
        double tempKappaSquared = parameter.H *parameter.H;
        for (int l_ = parameter.LMIN; l_ <= parameter.LMAX; l_++) {
                partition.l=l_;
                partition.waveFunctionRe[l_][0]=0.0;
                partition.waveFunctionRe[l_][1]=0.4;
                partition.waveFunctionIm[l_][0]=0.0;
                partition.waveFunctionIm[l_][1]=0.4;
                partition.waveFunctionPrimeRe[l_][0]=0.4;
                partition.waveFunctionPrimeIm[l_][0]=0.4;

                for (int i=1; i < parameter.N; i++) {
                        tempKernelFunction = 2.0 -tempKappaSquared
                                             *partition.kernelFunction(parameter.R[i]);
                        partition.waveFunctionRe[l_][i+1] =
                                tempKernelFunction *partition.waveFunctionRe[l_][i]
                                -partition.waveFunctionRe[l_][i-1]
                                -tempKappaSquared *partition.twoMuDevidedByhBarSquared
                                *partition.getWoodSaxonVolumeIm(parameter.R[i])
                                *partition.waveFunctionIm[l_][i];
                        partition.waveFunctionIm[l_][i+1] =
                                tempKernelFunction *partition.waveFunctionIm[l_][i]
                                -partition.waveFunctionIm[l_][i-1]
                                +tempKappaSquared *partition.twoMuDevidedByhBarSquared
                                *partition.getWoodSaxonVolumeIm(parameter.R[i])
                                *partition.waveFunctionRe[l_][i];
                        partition.waveFunctionPrimeRe[l_][i]=
                                (partition.waveFunctionRe[l_][i+1] -partition.waveFunctionRe[l_][i-1] )
                                /parameter.H;
                        partition.waveFunctionPrimeIm[l_][i]=
                                (partition.waveFunctionIm[l_][i+1] -partition.waveFunctionIm[l_][i-1] )
                                /parameter.H;
                }
        }
}

void
applyNumerovMethodFor (partitionClass &partition, generalParametersClass &parameter){
/*
        double gamma = parameter.H *parameter.H /12.0;
        partition.waveFunctionRe[0]=0.0;
        partition.waveFunctionRe[1]=0.01;
        for (int i=2; i < parameter.N; i++) {
                partition.waveFunctionRe[i] =
                        1. /(  1. +gamma *partition.kernelFunction( parameter.R[i] )  )
 *(
                                (2. -10.0 *gamma *partition.kernelFunction(parameter.R[i-1]) )
 * partition.waveFunctionRe[i-1]
                                - (1. +gamma *partition.kernelFunction(parameter.R[i-2]))
 * partition.waveFunctionRe[i-2]
                                );
        }
 */
}

void
fitDepthOfPotential (partitionClass &partition, generalParametersClass &parameters, double deviation=0.2){
/*
        double leftTest=  partition.V0_r -partition.V0_r *deviation;
        double rightTest= partition.V0_r +partition.V0_r *deviation;

        partition.V0_r= leftTest;
        applyFiniteDiffMethodFor(partition, parameters);
        double leftFunction =partition.waveFunctionRe[parameters.N-1];
        partition.V0_r=rightTest;
        applyFiniteDiffMethodFor(partition, parameters);
        double rightFunction =partition.waveFunctionRe[parameters.N-1];

        if(leftFunction<rightFunction) {
                double temp;
                temp =rightTest;
                rightTest =leftTest;
                leftTest =temp;
        }
        cout<<"FOR THE PARTITION " <<partition.name<<": "<<endl;
        printf("THE DEPTH IN THE RANGE FROM %.3f TO %.3f IS BEING FITTED \n", leftTest, rightTest);
        double waveFunctionAssimp;
        double middleTest=leftTest;
        bool test =false;
        int rootCounter=0;
        while ( !test ) {
                middleTest =(rightTest +leftTest) /2;
                partition.V0_r =middleTest;
                applyFiniteDiffMethodFor(partition, parameters);
                waveFunctionAssimp = partition.waveFunctionRe[parameters.N-1];
                test =  waveFunctionAssimp < parameters.EPSILON
                                             && waveFunctionAssimp > 0.0;
                if( waveFunctionAssimp > 0) leftTest =middleTest;
                else rightTest= middleTest;
                rootCounter++;
                if (rootCounter > 100) break;
                //        printf("%d \t %.9f \n",rootCounter, middleTest );
        }

        if (test) printf("THE DEPTH V0 FITTED INTO %.5f WITH %d TRIAL\n\n", middleTest, rootCounter);

        else printf("WARNING: COULDN'T FIT THE DEPTH OF THE POTENTIAL!\n\n" );
 */
}

int
partitionClass::setCoulombWaveFunctions(generalParametersClass &parameters) {
        double F_exponent, G_exponent;
        return
                gsl_sf_coulomb_wave_FGp_array(parameters.LMIN, parameters.LMAX,
                                              n, k *parameters.R[parameters.N-1],
                                              //1.0, 5.0,
                                              &coulombF[0], &coulombFPrime[0],
                                              &coulombG[0], &coulombGPrime[0],
                                              &F_exponent, &G_exponent);
}

int arraySize =201;
generalParametersClass generalParameters(arraySize, 0.001, 40.001, 0, 1, 1.E-6);
partitionClass firstPartition("Zn67_d", true, 6.0, 2., 67., 1, 30, generalParameters);


int
functionForChi (double t, const double x[], double f[],
                void *params)
{
        double remnantForChi = *(double *)params;
        f[0] = x[1];
        f[1] =-firstPartition.kernelFunction(t) *x[0] +remnantForChi;
        return GSL_SUCCESS;
}


int
functionForUpsilon (double t, const double y[], double f[],
                    void *params)
{
        double remnantForUpsilon = *(double *)params;
        f[0] = y[1];
        f[1] =-firstPartition.kernelFunction(t) *y[0] +remnantForUpsilon;
        return GSL_SUCCESS;
}

void
gsl_solveODE(partitionClass &partition, generalParametersClass &parameter){
        double tChi = parameter.RMIN;
        double tUpsilon = parameter.RMIN;

        double remnantForChi;
        double remnantForUpsilon;

        double y[2] = { 0.0, 0.4 };
        double x[2] = { 0.0, 0.4 };

        gsl_odeiv2_system systemForChi =       {functionForChi, NULL, 2, &remnantForChi};
        gsl_odeiv2_system systemForUpsilon =       {functionForUpsilon, NULL, 2, &remnantForUpsilon};

        gsl_odeiv2_driver * driverForChi =
                gsl_odeiv2_driver_alloc_y_new (&systemForChi, gsl_odeiv2_step_rk4, 1e-6, 1e-6, 0.0);
        gsl_odeiv2_driver * driverForUpsilon =
                gsl_odeiv2_driver_alloc_y_new (&systemForUpsilon, gsl_odeiv2_step_rk4, 1e-6, 1e-6, 0.0);

        int statusForChi;
        int statusForUpsilon;

        for (int i = 1; i < parameter.N; i++)
        {
                remnantForChi= -partition.twoMuDevidedByhBarSquared
                               *partition.V0_i/ exp( (parameter.R[i-1]-partition.R0_i) /partition.a0_i)
                               * partition.waveFunctionIm[0][i-1];
                remnantForUpsilon= partition.twoMuDevidedByhBarSquared
                                   *partition.V0_i/ exp( (parameter.R[i-1]-partition.R0_i) /partition.a0_i)
                                   *partition.waveFunctionRe[0][i-1];
                statusForChi =
                        gsl_odeiv2_driver_apply (driverForChi, &tChi, parameter.R[i], x);
                statusForUpsilon
                        = gsl_odeiv2_driver_apply (driverForUpsilon, &tUpsilon, parameter.R[i], y);
                if (statusForChi && statusForUpsilon != GSL_SUCCESS)  {
                        printf ("error, return value=%d\n", statusForChi);
                        break;
                }
                printf ("% .3f\t % .5e\t% .5e\t \t % .5e\t% .5e\n", parameter.R[i], y[0], y[1], x[0], x[1]);
                partition.waveFunctionRe[0][i]=x[0];
                partition.waveFunctionIm[0][i]=y[0];
                partition.waveFunctionPrimeRe[0][i]=x[1];
                partition.waveFunctionPrimeIm[0][i]=y[1];

        }
        gsl_odeiv2_driver_free (driverForUpsilon);
        gsl_odeiv2_driver_free (driverForChi);
}

int
main (void) {

        firstPartition.setWoodsSaxonParameters(-107.0, 4.32555, 0.86, 0.0, 1.0, 1.0, 4.87386);

        gsl_solveODE(firstPartition, generalParameters);

//        partitionClass secondPartition("Zn67_d", true, 6.0, 2., 67., 1, 30, generalParameters);
//        secondPartition.setWoodsSaxonParameters(-107.0, 4.32555, 0.86, 0.0, 1.0, 1.0, 4.87386);
//        applyFiniteDiffMethodFor(secondPartition, generalParameters);
//    fitDepthOfPotential(firstPartition, generalParameters, 0.5);
//        secondPartition.printWaveFunction(generalParameters);
//        secondPartition.getSmatrix(generalParameters);


        for(int Ri=0; Ri<arraySize; Ri++) {
                //          printf("%.3f\t% .5e\t% .5e\n", generalParameters.R[Ri], firstPartition.waveFunctionRe[0][Ri],
                //               firstPartition.waveFunctionIm[0][Ri]);
        }

        // for (auto i = generalParameters.R.begin(); i != generalParameters.R.end(); ++i)
        //cout << *i << "\n";
        return 0;
}
