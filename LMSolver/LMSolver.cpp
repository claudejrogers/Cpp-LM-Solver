//
//  LMSolver.cpp
//  LMSolver
//
//  Created by Claude Rogers on 7/9/11.
//  Copyright 2011 California Institute of Technology. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
using namespace std;

#include "LMSolver.h"

LMSolver::LMSolver (char *filename, int myModel, 
                    double var1, double var2, double var3, double var4)
{
    vector<double> xvals;
    vector<double> yvals;
    int i;
    double tempX, tempY;
    ifstream infile;
    string line, xstr, ystr;
    infile.open(filename);
    while (infile.good()) {
        getline(infile, line);
        istringstream iss(line);
        
        getline(iss, xstr, '\t');
        getline(iss, ystr, '\t');
        
        tempX = strtod(xstr.c_str(), NULL);
        tempY = strtod(ystr.c_str(), NULL);
        
        xvals.push_back(tempX);
        yvals.push_back(tempY);
    }
    datalen = (int) xvals.size();
    param[0] = var1;
    param[1] = var2;
    param[2] = var3;
    param[3] = var4;
    var[0] = var1;
    var[1] = var2;
    var[2] = var3;
    var[3] = var4;
    x  = new double[datalen];
    y  = new double[datalen];
    m  = new double[datalen];
    d1 = new double[datalen];
    d2 = new double[datalen];
    d3 = new double[datalen];
    d4 = new double[datalen];
    model = myModel;
    switch (model) {
        case boltzmann:
        case gaussian:
            varlen = 4;
            break;
        case expdecay:
        case hill:
        case ic50:
        case modsin:
            varlen = 3;
            break;
        case mm:
            varlen = 2;
            break;
        default:
            exit(0);
            break;
    }
    vector<double>::iterator it;
    for (it = xvals.begin(), i = 0; it < xvals.end(); it++, i++)
        x[i] = *it;
    for (it = yvals.begin(), i = 0; it < yvals.end(); it++, i++) 
        y[i] = *it;
}

void LMSolver::equation ()
{
    int i;
    for (i = 0; i < datalen; i++) {
        double xi = x[i];
        switch (model) {
            case boltzmann:
                m[i] = var[0] + ((var[1] - var[0]) / 
                                 (1 + exp((var[2] - xi) / var[3])));
                break;
            case expdecay:
                m[i] = var[0] + var[1] + exp(-var[2] * xi);
                break;
            case gaussian:
                m[i] = var[0] + var[1] * exp(-(xi - var[2]) * (xi - var[2]) / 
                                             (var[3] * var[3]));
                break;
            case hill:
                m[i] = (var[0] / (1 + pow((var[1] / xi), var[2])));
                break;
            case ic50:
                m[i] = (1 - (var[0]/(1 + pow((var[1]/xi), var[2]))));
                break;
            case mm:
                m[i] = ((var[0] * xi) / (var[1] + xi));
                break;
            case modsin:
                m[i] = var[0] * sin(M_PI * (xi - var[1]) / var[2]);
                break;
            default:
                exit(1);
                break;
        }
    }
}

void LMSolver::derivatives ()
{
    int i;
    for (i = 0; i < datalen; i++) {
        double xi = x[i];
        switch (model) {
            case boltzmann:
                d1[i] = 1 - (1 / (1 + exp((var[2] - xi) / var[3])));
                d2[i] = 1 / (1 + exp((var[2] - xi) / var[3]));
                d3[i] = (((var[0] - var[1]) * exp((var[2] + xi) / var[3])) / 
                         (var[3] * pow((exp(var[2] / var[3]) + 
                                        exp(xi / var[3])), 2)));
                d4[i] = (((var[1] - var[0]) * (var[2] - xi) * 
                          exp((var[2] - xi) / var[3])) / (var[3] * var[3] * 
                          pow((exp((var[2] - xi) / var[3]) + 1), 2)));
                break;
            case expdecay:
                d1[i] = 1.0;
                d2[i] = exp(-var[2] * xi);
                d3[i] = -xi*var[1]*exp(-var[2]*xi);
                break;
            case gaussian:
                d1[i] = 1.0;
                d2[i] = exp(-(xi - var[2]) * (xi - var[2])/
                            (var[3] * var[3]));
                d3[i] = 2 * (((xi - var[2]) / (var[3] * var[3])) * var[1] *
                             exp(-(xi - var[2]) * (xi - var[2]) / 
                                 (var[3]*var[3])));
                d4[i] = 2 * ((xi - var[2])*(xi - var[2])/
                             (var[3]*var[3]*var[3]))*var[1]*\
                exp(-(xi - var[2])*(xi - var[2]) / (var[3] * var[3]));
                break;
            case hill:
                d1[i] = (1/(1 + pow((var[1]/xi), var[2])));
                d2[i] = ((-var[0] * var[2] * pow(var[1], (var[2] - 1)) *
                          pow(xi, var[2])) / 
                         pow((pow(var[1], var[2]) + pow(xi, var[2])), 2));
                d3[i] = ((var[0]*pow((var[1] * xi), var[2]) * 
                          log(xi/var[1]))/pow((pow(var[1], var[2]) + 
                                               pow(xi, var[2])), 2));
                break;
            case ic50:
                d1[i] = (-(1/(1 + pow((var[1]/xi), var[2]))));
                d2[i] = ((var[0] * var[2] * pow((var[1]/xi), (var[2] - 1))) / 
                         (xi * pow((1 + (pow((var[1]/xi), var[2]))), 2)));
                d3[i] = ((var[0] * pow((var[1]/xi), var[2]) * 
                          log((var[1]/xi))) / 
                         (pow((1 + pow((var[1]/xi), var[2])), 2)));
                break;
            case mm:
                d1[i] = (xi / (var[1] + xi));
                d2[i] = (-(var[0] * xi) / pow((var[1] + xi), 2.0));
                break;
            case modsin:
                d1[i] = sin(M_PI * (xi - var[1]) / var[2]);
                d2[i] = (-var[0] * M_PI * 
                         cos(M_PI * (xi - var[1]) / var[2])) / var[2];
                d3[i] = ((var[0] * M_PI * (var[1] - xi) 
                          * cos(M_PI * (xi - var[1]) / var[2])) / 
                         pow(var[2], 2));
                break;
            default:
                exit(2);
                break;
        }
    }
}

void LMSolver::get_f (gsl_vector *fvect)
{
    equation();
    double fv;
    int i;
    for (i = 0; i < datalen; i++) {
        fv = y[i] - m[i];
        gsl_vector_set(fvect, i, fv);
    }
}

void LMSolver::get_jac (double *jac)
{
    derivatives();
    int i, j, k, l;
    i = 0;
    j = datalen;
    k = 2 * datalen;
    l = 3 * datalen;
    for (i = 0; i < datalen; i++, j++, k++, l++) {
        switch (model) {
            case boltzmann:
            case gaussian:
                jac[i] = -d1[i];
                jac[j] = -d2[i];
                jac[k] = -d3[i];
                jac[l] = -d4[i];
                break;
            case expdecay:
            case hill:
            case ic50:
            case modsin:
                jac[i] = -d1[i];
                jac[j] = -d2[i];
                jac[k] = -d3[i];
                break;
            case mm:
                jac[i] = -d1[i];
                jac[j] = -d2[i];
                break;
            default:
                exit(3);
                break;
        }
    }
}

void LMSolver::solve_h (double *a, gsl_matrix *muImat, double *g, double *h)
{
    int i, s, l, dl;
    l = varlen;
    dl = l * l;
    double *na, *ng;
    
    na = new double[dl];
    ng = new double[l];
    
    for (i = 0; i < dl; i++)
        na[i] = a[i];
    for (i = 0; i < l; i++)
        ng[i] = g[i];
    
    gsl_matrix_view nA = gsl_matrix_view_array(na, l, l);
    gsl_vector_view nG = gsl_vector_view_array(ng, l);
    gsl_vector_view nH = gsl_vector_view_array(h, l);
    gsl_permutation *p = gsl_permutation_alloc(l);
    
    gsl_matrix_add(&nA.matrix, muImat);
    gsl_vector_scale(&nG.vector, -1);
    gsl_linalg_LU_decomp(&nA.matrix, p, &s);
    gsl_linalg_LU_solve(&nA.matrix, p, &nG.vector, &nH.vector);
    
    delete [] na;
    delete [] ng;
    gsl_permutation_free(p);
}

double LMSolver::get_rho(gsl_vector *fvect, gsl_vector *newfvect, 
                         double mu, double *h, double *g)
{
    int i;
    double dF, dF1, dF2, dL;
    double *c;
    
    c = new double[varlen];
    
    gsl_vector_view new_h = gsl_vector_view_array(h, varlen);
    
    gsl_blas_ddot(fvect, fvect, &dF1);
    gsl_blas_ddot(newfvect, newfvect, &dF2);
    
    dF = 0.5 * dF1 - 0.5 * dF2;
    
    for (i = 0; i < varlen; i++)
        c[i] = mu * h[i] - g[i];
    
    gsl_vector_view c_v = gsl_vector_view_array(c, varlen);
    gsl_blas_ddot(&new_h.vector, &c_v.vector, &dL);
    dL *= 0.5;
    
    delete [] c;
    
    return dF/dL;
}

double LMSolver::get_mu (double *a)
{
    int i;
    int alen = varlen * varlen;
    double max = a[0];
    
    for (i = 0; i < alen; i += (varlen + 1))
        max = (max < a[i]) ? a[i] : max;
    
    return 1e-3*max;
}

void LMSolver::levenberg_marquardt ()
{
    int i, k, v;
    double mu, rho, normH;
    double *j;
    double *a;
    double *g;
    double *h;
    
    j = new double[datalen * varlen];
    a = new double[varlen * varlen];
    g = new double[varlen];
    h = new double[varlen];
    
    k = 0;
    v = 2;
    
    gsl_vector *f      = gsl_vector_alloc(datalen);
    gsl_vector *new_f  = gsl_vector_alloc(datalen);
    gsl_matrix *muI    = gsl_matrix_alloc(varlen, varlen);
    
    gsl_matrix_view JT = gsl_matrix_view_array(j, varlen, datalen);
    gsl_matrix_view A  = gsl_matrix_view_array(a, varlen, varlen);
    
    gsl_vector_view G  = gsl_vector_view_array(g, varlen);
    gsl_vector_view H  = gsl_vector_view_array(h, varlen);
    
    get_f(f);
    get_jac(j);
    
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, 
                   &JT.matrix, &JT.matrix, 0.0, &A.matrix);
    gsl_blas_dgemv(CblasNoTrans, 1.0, &JT.matrix, f, 0.0, &G.vector);
    
    mu = get_mu(a);
    
    while ((fabs(g[gsl_blas_idamax(&G.vector)]) >= NORM_G)) {
        k++;
        gsl_matrix_set_identity(muI);
        gsl_matrix_scale(muI, mu);
        
        solve_h(a, muI, g, h);
        
        for (i = 0; i < varlen; i++) {
            param[i] = var[i];
            var[i] += h[i];
        }
        
        if ((normH = gsl_blas_dnrm2(&H.vector)) <= NORM_H) {
            cout << "var converged" << endl;
            break;
        }
        
        get_f(new_f);
        
        rho = get_rho(f, new_f, mu, h, g);
        
        if (rho > 0) {
            gsl_vector_swap(f, new_f);
            get_jac(j);
            gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, 
                           &JT.matrix, &JT.matrix, 0.0, &A.matrix);
            gsl_blas_dgemv(CblasNoTrans, 1.0, &JT.matrix, f, 0.0, &G.vector);
            mu *= GSL_MAX_DBL((1.0/3.0), (1 - pow((2*rho - 1), 3.0)));
            v = 2;
        } else {
            for (i = 0; i < varlen; i++) 
                var[i] = param[i];
            mu *= v;
            v *= 2;
        }
        cout << "iter " << k << ": var = ";
        for (i = 0; i < varlen; i++) 
            cout << " " << var[i];
        cout << ", |f(x)| = " << gsl_blas_dnrm2(f) << endl;
    }
    
    delete [] j;
    delete [] a;
    delete [] g;
    delete [] h;
    
    gsl_vector_free(f);
    gsl_vector_free(new_f);
    gsl_matrix_free(muI);
}

LMSolver::~LMSolver() 
{
    delete [] x;
    delete [] y;
    delete [] m;
    delete [] d1;
    delete [] d2;
    delete [] d3;
    delete [] d4;
}