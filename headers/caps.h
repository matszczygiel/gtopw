#ifndef CAPS_H
#define CAPS_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <complex>

#include "../headers/gtopw.h"

void SphericalBesselJ( const int nmax, const double z, double* jn );
void BuildAlphaIJ( const usint ij_max, double** v );

inline double 
    CartSolidAng( const usint i, const usint j, 
                  const usint k, double** a ) {
    /* x^i y^j z^k integrated over the solid angle */
    if( ( i + j ) % 2 == 1 ) return 0.0;
    return 2.0 * a[i+j+1][k] * a[j][i];
};

void CartSolidAng( const usint i, const usint j, const usint k, 
                   const usint l, double** a, double* proj );

void CompSolidHarm( const double x, 
                    const double y,
                    const double z,
                    const usint l, double* ylm );

void CompAngFac( const double kx, 
                 const double ky,
                 const double kz,
                 const usint lP, const usint lG, 
                 const usint L , double** aij,
                 double* ang_fac );

void CAPAngFac( const double kx, const double ky, const double kz,
                const usint lG, const usint Lmax , 
                double** aij  , double** ang_fac );

void Faddeeva( const double xi, const double yi,
               double & u, double & v, bool flag );
                
void Faddeeva( const quad xi  , const quad yi  ,
               quad   & u, quad   & v, bool flag );
               
void CAPRadFac( const int lmax, const int nmax,
                const double k, const double p,
                double** rad_fac );

#endif //CAPS_H
