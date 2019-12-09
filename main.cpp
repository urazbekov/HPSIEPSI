#include <fstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_coulomb.h>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string>
#include <vector>
using namespace std;

class generalParametersClass {
public:
generalParametersClass(int N_, double RMIN_, double RMAX_, int LMIN_,
                       int LMAX_, double EPSILON_);
int N;
double RMAX;
double RMIN;
double H;
int LMIN;
int LMAX;
double EPSILON;
vector<double> R;
};

generalParametersClass::generalParametersClass(int N_, double RMIN_,
                                               double RMAX_, int LMIN_,
                                               int LMAX_, double EPSILON_){
        N = N_;
        RMIN = RMIN_;
        RMAX = RMAX_;
        LMIN = LMIN_;
        LMAX = LMAX_;
        H = (RMAX - RMIN) / (N - 1.0 );
        EPSILON = EPSILON_;
        for (double iteratorR = RMIN; iteratorR <= RMAX + H; iteratorR = iteratorR + H) {
                R.push_back(iteratorR);
        }
}

class partitionClass {
public:
partitionClass(const char * name, bool energyLab, double energy, double m1,
               double m2, int z1, int z2, generalParametersClass &parameters);
string name;
double hBarSquared;
bool energyLab;
double energy;
double m1, m2;
int z1, z2;
double m, twoMuDevidedByhBarSquared;
double k, kSquared;
double n;
vector<double> coulombF;
vector<double> coulombG;
vector<double> coulombFPrime;
vector<double> coulombGPrime;
vector<double> sMatrixRe;
vector<double> sMatrixIm;
vector<double> rowWaveFunctionRe;
vector<double> rowWaveFunctionIm;
vector<double> rowWaveFunctionPrimeRe;
vector<double> rowWaveFunctionPrimeIm;
vector<vector<double> > waveFunctionRe;
vector<vector<double> > waveFunctionIm;
vector<vector<double> > waveFunctionPrimeRe;
vector<vector<double> > waveFunctionPrimeIm;
bool
boundStateQ();

void
setWoodsSaxonRe(double V0_, double R0_, double a0_);
void
setWoodsSaxonIm(double W0_, double R0_, double a0_);
void
setCoulombR(double Rc_);

double
getWoodSaxonVolumeIm(double R);
double
getWoodSaxonVolumeRe(double R);
double
getCoulombPotential(double R);
double
kernelFunction(double R);
double V0_r, R0_r, a0_r;
double V0_i, R0_i, a0_i;
double Rc;
int l;
int
setCoulombWaveFunctions(generalParametersClass &parameter);
void
printWaveFunction(generalParametersClass &parameter);
void
getSmatrix(generalParametersClass &parameter);

private:
};




partitionClass::partitionClass(const char * nameChar, bool energyLabBool,
                               double energyDouble, double m1Double,
                               double m2Double, int z1Int, int z2Int,
                               generalParametersClass &parameters)
// using Initializer List
        : rowWaveFunctionRe(parameters.N, 0.0),
        rowWaveFunctionIm(parameters.N, 0.0),
        rowWaveFunctionPrimeRe(parameters.N, 0.0),
        rowWaveFunctionPrimeIm(parameters.N, 0.0),
        waveFunctionRe(parameters.LMAX + 1, rowWaveFunctionRe),
        waveFunctionIm(parameters.LMAX + 1, rowWaveFunctionIm),
        waveFunctionPrimeRe(parameters.LMAX + 1, rowWaveFunctionPrimeRe),
        waveFunctionPrimeIm(parameters.LMAX + 1, rowWaveFunctionPrimeIm){
        hBarSquared =  41.801651165221026;
        name = nameChar;
        energyLab = energyLabBool;
        energy = energyDouble;
        m1 = m1Double;
        m2 = m2Double;
        z1 = z1Int;
        z2 = z2Int;
        if (energyLab)
                energy = energy * m2 / (m1 + m2);
        m = m1 * m2 / (m1 + m2);
        // 2*amu/hbar^2 = 0.0478450
        twoMuDevidedByhBarSquared = 2 * m / hBarSquared;
        kSquared = 2 * m * energy / hBarSquared;
        k = sqrt(abs(kSquared));
        n = z1 * z2 * 1.43997 *m / hBarSquared / k;
        coulombF.resize(parameters.LMAX + 1, 0.0);
        coulombG.resize(parameters.LMAX + 1, 0.0);
        coulombFPrime.resize(parameters.LMAX + 1, 0.0);
        coulombGPrime.resize(parameters.LMAX + 1, 0.0);
        sMatrixRe.resize(parameters.LMAX + 1, 0.0);
        sMatrixIm.resize(parameters.LMAX + 1, 0.0);
}

void
partitionClass::setWoodsSaxonRe(double V0_, double R0_, double a0_){
        V0_r =V0_;
        R0_r=R0_;
        a0_r=a0_;
}

void
partitionClass::setWoodsSaxonIm(double W0_, double R0_, double a0_){
        V0_i =W0_;
        R0_i=R0_;
        a0_i=a0_;
}

void
partitionClass::setCoulombR(double Rc_){
        Rc =Rc_;
}

double
partitionClass::getWoodSaxonVolumeRe(double R){
        return -twoMuDevidedByhBarSquared * V0_r / (1 + exp((R - R0_r) / a0_r));
}

double
partitionClass::getWoodSaxonVolumeIm(double R){
        return -twoMuDevidedByhBarSquared * V0_i / (1 + exp((R - R0_i) / a0_i));
}

double
partitionClass::getCoulombPotential(double R){
        if (R > Rc)
                return twoMuDevidedByhBarSquared *
                       z1 * z2 * 1.43997 / R;
        else
                return twoMuDevidedByhBarSquared *
                       z1 * z2 * 1.43997 * (3 - R * R / Rc / Rc) / Rc * 0.5;
}

void
partitionClass::getSmatrix(generalParametersClass &parameter){
        setCoulombWaveFunctions(parameter);
        double a, b, c, d;
        double A, B, C, D;
        double AsquaredPlusBsquared;
        for (int l_ = parameter.LMIN; l_ <= parameter.LMAX; l_++) {
                a = waveFunctionPrimeRe[l_][parameter.N - 2] * coulombF[l_]
                    - waveFunctionRe[l_][parameter.N - 2] * coulombFPrime[l_];
                b = waveFunctionIm[l_][parameter.N - 2] * coulombGPrime[l_]
                    - waveFunctionPrimeIm[l_][parameter.N - 2] * coulombG[l_];
                c = waveFunctionPrimeIm[l_][parameter.N - 2] * coulombF[l_]
                    - waveFunctionIm[l_][parameter.N - 2] * coulombFPrime[l_];
                d = waveFunctionRe[l_][parameter.N - 2] * coulombGPrime[l_]
                    - waveFunctionPrimeRe[l_][parameter.N - 2] * coulombG[l_];
                A = b - a;
                B = -c - d;
                C = a + b;
                D = c - d;
                AsquaredPlusBsquared = A * A + B * B;
                sMatrixRe[l_] = (A * C + B * D) / AsquaredPlusBsquared;
                sMatrixIm[l_] = (A * D - B * C) / AsquaredPlusBsquared;
        }
} // partitionClass::getSmatrix

void
partitionClass::printWaveFunction(generalParametersClass &parameter){
        FILE * pFile;
        string filename = name + ".wfn";

        pFile = fopen(filename.data(), "w");
        //      fprintf(pFile, "# Wave Function For the ");
        //      fprintf(pFile, "%s", name.data());
        //      fprintf(pFile, " partition \n");

        if (boundStateQ()) {
                fprintf(pFile, "# R [fm]     Psi \n");
                for (int i = 0; i < parameter.N; i++) {
                        //                      fprintf (pFile, "% .3e\t%
                        //                      .5e\n",parameter.R[i],waveFunctionRe[i]);
                }
        } else {
                //          fprintf(pFile, "# R [fm]       ");
                for (int l_ = parameter.LMIN; l_ <= parameter.LMAX; l_++) {
                        //                  fprintf(pFile, "PsiRe[%d]      PsiIm[%d]        ", l_, l_);
                }
                //      fprintf(pFile, "\n");
                for (int Ri = 0; Ri < parameter.N; Ri++) {
                        fprintf(pFile, "% .3e\t", parameter.R[Ri]);

                        for (int l_ = parameter.LMIN; l_ <= parameter.LMAX; l_++) {
                                fprintf(pFile, "% .5e\t% .5e\t", waveFunctionRe[l_][Ri],
                                        waveFunctionIm[l_][Ri]);
                        }
                        fprintf(pFile, "\n");
                }
        }
        fclose(pFile);
} // partitionClass::printWaveFunction

double
partitionClass::kernelFunction(double R){
        return kSquared - l * (l + 1) / R / R
               - getWoodSaxonVolumeRe(R) -  getCoulombPotential(R);
}



bool
partitionClass::boundStateQ(){
        if (energy < 0)
                return true;
        else
                return false;
}

void
applyFiniteDiffMethodFor(partitionClass &partition,
                         generalParametersClass                &parameter){
        double subKernelFunction;
        double kappaSquared = parameter.H * parameter.H;

        for (int l_ = parameter.LMIN; l_ <= parameter.LMAX; l_++) {
                partition.l = l_;
                partition.waveFunctionRe[l_][0] = 0.0;
                partition.waveFunctionRe[l_][1] = 0.01;
                partition.waveFunctionIm[l_][0] = 0.0;
                partition.waveFunctionIm[l_][1] = 0.01;
                partition.waveFunctionPrimeRe[l_][0] = 0.01;
                partition.waveFunctionPrimeIm[l_][0] = 0.01;
                for (int i = 1; i < parameter.N; i++) {
                        subKernelFunction =
                                2.0 - kappaSquared * partition.kernelFunction(parameter.R[i]);
                        partition.waveFunctionRe[l_][i + 1] =
                                subKernelFunction * partition.waveFunctionRe[l_][i]
                                - partition.waveFunctionRe[l_][i - 1]
                                - kappaSquared
                                * partition.getWoodSaxonVolumeIm(parameter.R[i])
                                * partition.waveFunctionIm[l_][i];
                        partition.waveFunctionIm[l_][i + 1] =
                                subKernelFunction * partition.waveFunctionIm[l_][i]
                                - partition.waveFunctionIm[l_][i - 1]
                                + kappaSquared
                                * partition.getWoodSaxonVolumeIm(parameter.R[i])
                                * partition.waveFunctionRe[l_][i];
                        partition.waveFunctionPrimeRe[l_][i] =
                                (partition.waveFunctionRe[l_][i + 1] - partition.waveFunctionRe[l_][i - 1])
                                / parameter.H;
                        partition.waveFunctionPrimeIm[l_][i] =
                                (partition.waveFunctionIm[l_][i + 1] - partition.waveFunctionIm[l_][i - 1])
                                / parameter.H;
                }
        }
} // applyFiniteDiffMethodFor

void
applyNumerovMethodFor(partitionClass &partition,
                      generalParametersClass &parameter){

        double gamma = parameter.H * parameter.H / 12.0;
        double gammaSquared = gamma * gamma;
        double Wip1, Wi, Wim1;
        double Vip1, Vi, Vim1;
        double xim1, xi;
        double yim1, yi;

        double A1, B1, C1, D1, E;
        double A2, B2, C2, D2;

        for (int l_ = 0; l_ <= parameter.LMAX; l_++) {
                partition.l = l_;
                partition.waveFunctionRe[l_][0] = 0.0;
                partition.waveFunctionRe[l_][1] = 0.01;
                partition.waveFunctionIm[l_][0] = 0.0;
                partition.waveFunctionIm[l_][1] = 0.01;
                partition.waveFunctionPrimeRe[l_][0]=0.01;
                partition.waveFunctionPrimeIm[l_][0]=0.01;
                for (int i = 1; i < parameter.N; i++) {

                        Wip1 = partition.kernelFunction(parameter.R[i + 1]);
                        Wi = partition.kernelFunction(parameter.R[i]);
                        Wim1 = partition.kernelFunction(parameter.R[i - 1]);

                        Vip1 = partition.getWoodSaxonVolumeIm(parameter.R[i + 1]);
                        Vi = partition.getWoodSaxonVolumeIm(parameter.R[i]);
                        Vim1 = partition.getWoodSaxonVolumeIm(parameter.R[i - 1]);

                        xim1 = partition.waveFunctionRe[l_][i - 1];
                        xi = partition.waveFunctionRe[l_][i];

                        yim1 = partition.waveFunctionIm[l_][i - 1];
                        yi = partition.waveFunctionIm[l_][i];

                        A1 = xi * (2 - 10 * gammaSquared * Vi * Vip1 - 10 * gamma * Wi + 2 * gamma * Wip1
                                   - 10 * gammaSquared * Wi * Wip1);
                        B1 = xim1 * (-1 - gammaSquared * Vim1 * Vip1 - gamma * Wim1 - gamma * Wip1
                                     -  gammaSquared * Wim1 * Wip1);
                        C1 = yi * (-10 * gamma * Vi - 2 * gamma * Vip1 + 10 * gammaSquared * Vip1 * Wi
                                   - 10 * gammaSquared * Vi * Wip1);
                        D1 = yim1 * (-gamma * Vim1 + gamma * Vip1 + gammaSquared * Vip1 * Wim1
                                     - gammaSquared * Vim1 * Wip1);
                        E = 1 + gammaSquared * Vip1 * Vip1 + 2 * gamma * Wip1 + gammaSquared * Wip1 * Wip1;


                        partition.waveFunctionRe[l_][i + 1] = (A1 + B1 + C1 + D1) / E;

                        A2 = xi * (10 * gamma * Vi + 2 * gamma * Vip1 - 10 * gammaSquared * Vip1 * Wi
                                   + 10 * gammaSquared * Vi * Wip1);
                        B2 = xim1 * (gamma * Vim1 - gamma * Vip1 - gammaSquared * Vip1 * Wim1
                                     + gammaSquared * Vim1 * Wip1);
                        C2 = yi * (2 - 10 * gammaSquared * Vi * Vip1 - 10 * gamma * Wi + 2 * gamma * Wip1
                                   - 10 * gammaSquared * Wi * Wip1);
                        D2 = yim1 * (-1 - gammaSquared * Vim1 * Vip1 - gamma * Wim1 - gamma * Wip1
                                     - gammaSquared * Wim1 * Wip1);

                        partition.waveFunctionIm[l_][i + 1] = (A2 + B2 + C2 + D2) / E;

                        partition.waveFunctionPrimeRe[l_][i]= (partition.waveFunctionRe[l_][i + 1] - xim1) /parameter.H;
                        partition.waveFunctionPrimeIm[l_][i]= (partition.waveFunctionIm[l_][i + 1] - yim1) /parameter.H;

                }
        }

} // applyNumerovMethodFor

void
fitDepthOfPotential(partitionClass &partition,
                    generalParametersClass           &parameters,
                    double deviation = 0.2){
        /*
         *      double leftTest=  partition.V0_r -partition.V0_r *deviation;
         *      double rightTest= partition.V0_r +partition.V0_r *deviation;
         *
         *      partition.V0_r= leftTest;
         *      applyFiniteDiffMethodFor(partition, parameters);
         *      double leftFunction =partition.waveFunctionRe[parameters.N-1];
         *      partition.V0_r=rightTest;
         *      applyFiniteDiffMethodFor(partition, parameters);
         *      double rightFunction =partition.waveFunctionRe[parameters.N-1];
         *
         *      if(leftFunction<rightFunction) {
         *              double temp;
         *              temp =rightTest;
         *              rightTest =leftTest;
         *              leftTest =temp;
         *      }
         *      cout<<"FOR THE PARTITION " <<partition.name<<": "<<endl;
         *      printf("THE DEPTH IN THE RANGE FROM %.3f TO %.3f IS BEING FITTED \n",
         * leftTest, rightTest); double waveFunctionAssimp; double
         * middleTest=leftTest; bool test =false; int rootCounter=0; while ( !test ) {
         *              middleTest =(rightTest +leftTest) /2;
         *              partition.V0_r =middleTest;
         *              applyFiniteDiffMethodFor(partition, parameters);
         *              waveFunctionAssimp = partition.waveFunctionRe[parameters.N-1];
         *              test =  waveFunctionAssimp < parameters.EPSILON
         *                                           && waveFunctionAssimp > 0.0;
         *              if( waveFunctionAssimp > 0) leftTest =middleTest;
         *              else rightTest= middleTest;
         *              rootCounter++;
         *              if (rootCounter > 100) break;
         *              //        printf("%d \t %.9f \n",rootCounter, middleTest );
         *      }
         *
         *      if (test) printf("THE DEPTH V0 FITTED INTO %.5f WITH %d TRIAL\n\n",
         * middleTest, rootCounter);
         *
         *      else printf("WARNING: COULDN'T FIT THE DEPTH OF THE POTENTIAL!\n\n" );
         */
}

int
partitionClass::setCoulombWaveFunctions(
        generalParametersClass &parameters){
        double F_exponent, G_exponent;

        return gsl_sf_coulomb_wave_FGp_array(
                parameters.LMIN, parameters.LMAX, n, k * parameters.R[parameters.N - 2],
                // 1.0, 5.0,
                &coulombF[0], &coulombFPrime[0], &coulombG[0], &coulombGPrime[0],
                &F_exponent, &G_exponent);
}

int
functionForChi(double t, const double x[], double f[], void * params){
        double temp1 = ((double *)params)[0];
        double temp2 = ((double *)params)[1];

        //  printf("%f\t%f\n", ((double *)params)[0], ((double *)params)[1]);
        f[0] = x[1];
        f[1] =  temp1 * x[0] + temp2;
        return GSL_SUCCESS;
}


int
functionForUpsilon(double t, const double y[], double f[], void * params){
        double temp1 = ((double *)params)[0];
        double temp2 = ((double *)params)[1];

        //  printf("%f\t%f\n", ((double *)params)[0], ((double *)params)[1]);
        f[0] = y[1];
        f[1] =  temp1 * y[0] + temp2;
        return GSL_SUCCESS;
}




void
gsl_solveODE(partitionClass &partition,
             generalParametersClass    &parameter){
        int l_ = partition.l;

        double tChi = parameter.RMIN;
        double tUpsilon = parameter.RMIN;

        double remnantForChi[2] = { 0.0, 0.0 };
        double remnantForUpsilon[2] = { 0.0, 0.0 };

        void * pointerRemnantForChi = &remnantForChi;
        void * pointerRemnantForUpsilon = &remnantForUpsilon;

        double y[2];
        double x[2];

        gsl_odeiv2_system systemForChi = { functionForChi, NULL, 2,
                                           pointerRemnantForChi };
        gsl_odeiv2_system systemForUpsilon = { functionForUpsilon, NULL, 2,
                                               pointerRemnantForUpsilon };

        gsl_odeiv2_driver * driverForChi = gsl_odeiv2_driver_alloc_y_new(
                &systemForChi, gsl_odeiv2_step_rk2, 1e-6, 1e-6, 0.0);
        gsl_odeiv2_driver * driverForUpsilon = gsl_odeiv2_driver_alloc_y_new(
                &systemForUpsilon, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);

        int statusForChi;
        int statusForUpsilon;

        //    for (int l_ = 0; l_ < parameter.LMAX; l_++) {
        x[0] = 0.0;
        x[1] = 0.001;
        y[0] = 0.0;
        y[1] = 0.001;
        partition.waveFunctionRe[l_][0] = 0.0;
        partition.waveFunctionIm[l_][0] = 0.0;
        partition.waveFunctionPrimeRe[l_][0] = 0.001;
        partition.waveFunctionPrimeIm[l_][0] = 0.001;
        for (int i = 1; i < parameter.N; i++) {
                remnantForChi[0] = -partition.kernelFunction(parameter.R[i]);
                remnantForChi[1] =  -partition.getWoodSaxonVolumeIm(parameter.R[i])
                                   * partition.waveFunctionIm[l_][i - 1];
                remnantForUpsilon[0] = -partition.kernelFunction(parameter.R[i]);
                remnantForUpsilon[1] = partition.getWoodSaxonVolumeIm(parameter.R[i])
                                       * partition.waveFunctionRe[l_][i - 1];
                statusForChi =  gsl_odeiv2_driver_apply(driverForChi, &tChi, parameter.R[i], x);
                statusForUpsilon =  gsl_odeiv2_driver_apply(driverForUpsilon, &tUpsilon, parameter.R[i], y);

                //  printf("%f\t%f\n", remnantForChi[1], remnantForUpsilon[1]);
                if (statusForChi && statusForUpsilon != GSL_SUCCESS) {
                        printf("error, return value=%d\n", statusForChi);
                }
                partition.waveFunctionRe[l_][i] = x[0];
                partition.waveFunctionIm[l_][i] = y[0];
                partition.waveFunctionPrimeRe[l_][i] = x[1];
                partition.waveFunctionPrimeIm[l_][i] = y[1];
                //        printf("%f\t% .5e\t% .5e\n", parameter.R[i], x[0], y[0]);
        }

//        }
        gsl_odeiv2_driver_free(driverForUpsilon);
        gsl_odeiv2_driver_free(driverForChi);
} // gsl_solveODE


int
main(void){

        int arraySize = 801;
        generalParametersClass generalParameters(arraySize, 0.001, 40.0, 0, 50, 1.E-6);
        partitionClass firstPartition("Zn67_d", true, 6.0, 2., 67., 1, 30,
                                      generalParameters);

        firstPartition.setWoodsSaxonRe(107.,4.26463,0.86);
        firstPartition.setWoodsSaxonIm(40.24,6.22635,0.884);
        firstPartition.setCoulombR(5.28001);



        //  gsl_solveODE(firstPartition, generalParameters);

        //        partitionClass secondPartition("Zn67_d", true, 6.0, 2., 67., 1, 30,
        //        generalParameters);
        //        secondPartition.setWoodsSaxonParameters(-107.0, 4.32555, 0.86,
        //        0.0, 1.0, 1.0, 4.87386);
        // applyFiniteDiffMethodFor(firstPartition, generalParameters);
        //        fitDepthOfPotential(firstPartition, generalParameters, 0.5);
        // firstPartition.printWaveFunction(generalParameters);
        applyNumerovMethodFor(firstPartition, generalParameters);
        firstPartition.getSmatrix(generalParameters);

        for (int l_ = 0; l_ <= generalParameters.LMAX; l_++) {
                printf("%d\t% .5f\t% .5f\n", l_, firstPartition.sMatrixRe[l_], firstPartition.sMatrixIm[l_]);
        }

        printf("%f\n",firstPartition.n);
        return 0;
        for (int i = 0; i < arraySize; i++) {
                printf("%.3f\t% .5e\t% .5e\n", generalParameters.R[i],
                       firstPartition.waveFunctionPrimeRe[firstPartition.l][i],
                       firstPartition.waveFunctionPrimeIm[firstPartition.l][i]);
        }
        //    printf("%f\n", generalParameters.R[0]);
        // for (auto i = generalParameters.R.begin(); i != generalParameters.R.end();
        // ++i)
        // cout << *i << "\n";
        return 0;
} // main
