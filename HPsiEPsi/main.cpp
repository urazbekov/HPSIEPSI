#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <string>
#include <gsl/gsl_sf_coulomb.h>
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
vector<double>  waveFunctionRe;
vector<double>  waveFunctionIm;
vector<double> rowWaveFunctionRe_;
vector<double> rowWaveFunctionIm_;
vector< vector<double> > waveFunctionRe_;
vector< vector<double> > waveFunctionIm_;
//    vector< vector <double> >  myVector;
//    vector< vector<int> >
//   double  interactionKernel (double R);
//    double  getAssimpFunctionAtV0 (double V0);
bool    boundStateQ();
void    setWoodsSaxonParameters(double V0_r_, double R0_r, double a0_r_,
                                double V0_i_, double R0_i_, double a0_i_, double Rc_ );
double  kernelFunction(double R);
double V0_r, R0_r, a0_r;
double V0_i, R0_i, a0_i;
double Rc;
int l;
void printWaveFunction(generalParametersClass &parameter);

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
        : rowWaveFunctionRe_(parameters.N,-1.0),
        rowWaveFunctionIm_(parameters.N,-1.0),
        waveFunctionRe_(parameters.LMAX+1,rowWaveFunctionRe_),
        waveFunctionIm_(parameters.LMAX+1,rowWaveFunctionIm_)
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
        waveFunctionRe.resize(parameters.LMAX+1, 0.0);
        waveFunctionIm.resize(parameters.LMAX+1, 0.0);
        coulombF.resize(parameters.LMAX+1, 0.0);
        coulombG.resize(parameters.LMAX+1, 0.0);
        coulombFPrime.resize(parameters.LMAX+1, 0.0);
        coulombGPrime.resize(parameters.LMAX+1, 0.0);

        l =0;

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
                        fprintf (pFile, "% .3e\t% .5e\n",parameter.R[i],waveFunctionRe[i]);
                }
        }
        else{
                fprintf(pFile, "# R [fm]       ");
                for (int l_ = 0; l_ <= parameter.LMAX; l_++) {
                        fprintf(pFile, "PsiRe[%d]      PsiIm[%d]        ",l_,l_);
                }
                fprintf(pFile, "\n");
                for (int Ri =0; Ri < parameter.N; Ri++) {
                        fprintf(pFile,"% .3e\t\t",parameter.R[Ri]);

                        for(int l_=0; l_<=parameter.LMAX; l_++) {
                                fprintf (pFile, "% .5e\t% .5e\t\t",
                                         waveFunctionRe_[l_][Ri], waveFunctionIm_[l_][Ri]);
                        }
                        fprintf(pFile, "\n");
                }


        }
        fclose (pFile);
}

double
partitionClass::kernelFunction(double R_){
        double coulumbPart;
        double volumePart;
        double centriFugalPart;
        volumePart = twoMuDevidedByhBarSquared *V0_r /( 1 +exp( (R_ - R0_r) /a0_r));
        centriFugalPart = l *(l +1) /R_ /R_;
        if (R_ > Rc) coulumbPart =1.  /R_;
        else coulumbPart =1.  *(3 -R_ *R_ /Rc /Rc) /Rc *0.5;
        return kSquared -volumePart -twoMuDevidedByhBarSquared *z1 *z2*1.43997 *coulumbPart
               -centriFugalPart;
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

        for (int l_ = 0; l_ <= parameter.LMAX; l_++) {
                partition.l=l_;
                partition.waveFunctionRe_[l_][0]=0.0;
                partition.waveFunctionRe_[l_][1]=0.01;
                partition.waveFunctionIm_[l_][0]=0.0;
                partition.waveFunctionIm_[l_][1]=0.01;
                for (int i=2; i < parameter.N; i++) {
                        partition.waveFunctionRe_[l_][i] =
                                (2 -parameter.H *parameter.H
                                 *partition.kernelFunction(parameter.R[i-1]))
                                *partition.waveFunctionRe_[l_][i-1] -partition.waveFunctionRe_[l_][i-2]
                                -parameter.H *parameter.H *partition.twoMuDevidedByhBarSquared
                                *partition.V0_i /(1 +exp(parameter.R[i-1] - partition.R0_i) /partition.a0_i)
                                *partition.waveFunctionIm_[l_][i-1];

                        partition.waveFunctionIm_[l_][i] =
                                (2 -parameter.H *parameter.H
                                 *partition.kernelFunction(parameter.R[i-1]))
                                *partition.waveFunctionIm_[l_][i-1] -partition.waveFunctionIm_[l_][i-2]
                                +parameter.H *parameter.H *partition.twoMuDevidedByhBarSquared
                                *partition.V0_i /(1 +exp(parameter.R[i-1] - partition.R0_i) /partition.a0_i)
                                *partition.waveFunctionRe_[l_][i-1];
                }
        }
}

void
applyNumerovMethodFor (partitionClass &partition, generalParametersClass &parameter){
        double gamma = parameter.H *parameter.H /12.0;
        partition.waveFunctionRe[0]=0.0;
        partition.waveFunctionRe[1]=0.01;
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
fitDepthOfPotential (partitionClass &partition, generalParametersClass &parameters, double deviation=0.2){
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
}

int
setCoulombWaveFunctionsFor(partitionClass &partition, generalParametersClass &parameters) {
        double F_exponent, G_exponent;
        return
                gsl_sf_coulomb_wave_FGp_array(0.0, parameters.LMAX,
                                              partition.n, partition.k *parameters.R.back(),
                                              //1.0, 5.0,
                                              &partition.coulombF[0], &partition.coulombFPrime[0],
                                              &partition.coulombG[0], &partition.coulombGPrime[0],
                                              &F_exponent, &G_exponent);
}


int
main (void) {
        int arraySize =201;
        generalParametersClass generalParameters(arraySize, 0.001, 40.001, 0, 30, 1.E-6);
        partitionClass firstPartition("16O_d", false, 17.526, 16, 2, 8, 1, generalParameters);
        firstPartition.setWoodsSaxonParameters(-27.0, 3.15, 0.65, 0.0, 1.25, 0.65, 3.15);
        applyFiniteDiffMethodFor(firstPartition, generalParameters);
//    fitDepthOfPotential(firstPartition, generalParameters, 0.5);
        firstPartition.printWaveFunction(generalParameters);

        //  cout<<firstPartition.myVector[0][0]<<endl;

        printf("%f\n", firstPartition.waveFunctionIm_[4][201]);


        //      if (setCoulombWaveFunctionsFor(firstPartition, generalParameters)!=0)
//                cout<< "error!"<<endl;

        for (int l = 0; l <= generalParameters.LMAX; l++) {
                printf("%d\t% .5e\t% .5e\n", l, firstPartition.waveFunctionRe_[l][4], firstPartition.waveFunctionRe_[l][10]);
        }


        // for (auto i = generalParameters.R.begin(); i != generalParameters.R.end(); ++i)
        //cout << *i << "\n";
        return 0;
}
