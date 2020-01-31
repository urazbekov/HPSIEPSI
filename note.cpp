
struct structParams{
  double *a;
  double *b;
  double *c;
  double *d;
  double *e;
};

double
f(double x, void * structParams)
{
  double *eta = ((struct structParams *) structParams) -> a;
  double f = gsl_pow_int(*eta, 2) + 1.0;
  return f;
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
                                / parameter.H /2.;
                        partition.waveFunctionPrimeIm[l_][i] =
                                (partition.waveFunctionIm[l_][i + 1] - partition.waveFunctionIm[l_][i - 1])
                                / parameter.H /2.;
                }
        }
} // applyFiniteDiffMethodFor




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
