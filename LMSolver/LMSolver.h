//
//  LMSolver.h
//  LMSolver
//
//  Created by Claude Rogers on 7/9/11.
//  Copyright 2011 California Institute of Technology. All rights reserved.
//

#ifndef __LMSOLVER__
#define __LMSOLVER__

#define NORM_G 1e-10
#define NORM_H 1e-15
#define MAX_ITER 500

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multifit_nlin.h>

enum Models {
    boltzmann,
    expdecay,
    gaussian,
    hill,
    ic50, 
    mm,
    modsin
};

class LMSolver {
  public:
    LMSolver (char *filename, int myModel, 
              double var1, double var2, double var3, double var4);
    ~LMSolver ();
    void   equation ();
    void   levenberg_marquardt ();
  private:
    char   *filename;
    int    model;
    int    varlen;
    int    datalen;
    double param[4];
    double var[4];
    double *x;
    double *y;
    double *m;
    double *d1;
    double *d2;
    double *d3;
    double *d4;
    void   derivatives ();
    void   get_f (gsl_vector *fvect);
    void   get_jac (double *jac);
    void   solve_h (double *a, gsl_matrix *muImat, double *g, double *h);
    double get_rho (gsl_vector *fvect, gsl_vector *newfvect, 
                    double mu, double *h, double *g);
    double get_mu (double *a);
    
};

#endif // __LMSOLVER__ 
