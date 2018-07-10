#include <cmath>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <ctime>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <complex>
#include <limits>
#include "omp.h"
#include "../headers/gtopw.h"
#include "../headers/caps.h"
#include "../headers/auxfun1.h"

using namespace std;

extern "C" { 
    #include <quadmath.h>
}

void SphericalBesselJ( const int nmax, const double z, double* jn ) {
    /* spherical Bessel functions, j_n(z), z>0 */
    if( nmax > 25 ) {
        cout << " The code has never been tested for n>25!" << endl;
        cout << " Emergency halt." << endl;
        exit( EXIT_FAILURE );
    };
    
    if( z < numeric_limits<double>::epsilon() ) {
        jn[0] = 1.0;
        return ;
    };
    
    if( z < 25.0 ) {
        usint n_strt,n;
        n_strt = 51;
        double* tn = new double[n_strt];
        n_strt--;
        
        tn[n_strt] = 0.0;
        for(n=n_strt; n>=1; n--) tn[n-1] = z / ( 2*n + 1 - z * tn[n] );
    
        jn[0] = sin(z) / z;
        for(n=1; n<=nmax; n++) jn[n] = tn[n-1] * jn[n-1];
        
        delete []tn;
    } else {
        usint n;
        jn[0] = sin(z) / z;
        jn[1] = ( jn[0] - cos(z) ) / z;
        for(n=1; n<nmax; n++) jn[n+1] = ( 2*n + 1 ) * jn[n] / z - jn[n-1];
    };
    ///for(n=0; n<=nmax; n++) cout << n << " " << jn[n] << endl;
    /* */
};

void BuildAlphaIJ( const usint ij_max, double** v ) {
    /* \alpha_{ij} = \int_0^\pi d\phi \sin^i \phi \cos^j \phi */
    usint i,j;
    
    v[0][0] = M_PI;
    v[1][0] = 2.0;
    for(i=2; i<=ij_max; i++) v[i][0] = v[i-2][0] * ( i - 1 ) / i;
    
    for(j=2; j<=ij_max; j+=2)
        for(i=0; i<=ij_max-j; i++)
            v[i][j] = v[i][j-2] - v[i+2][j-2];
    /* */
};

void CartSolidAng( const usint i, const usint j, const usint k, 
                   const usint l, double** a, double* proj ) {
    /* \int d\Omega x^i y^j z^k Y_{lm}(r) / r^(i+j+k) */
    double sum;
    int m;
    usint x,y,z;
    
    for(m=-l; m<=l; m++) {
        sum = 0.0;
        for(x=0; x<=l; x++)
            for(y=0; y<=l-x; y++) 
                for(z=0; z<=l-x-y; z++)
            sum += NoNormCalcClmR( l, m, x, y, z )
                 * CartSolidAng( i + x, j + y, k + z, a );
        proj[m+l] = sum;
    };
    
    /* */
};

void CompSolidHarm( const double x, 
                    const double y,
                    const double z,
                    const usint l, double* ylm ) {
    /* real solid harmonics (no Racah) at [x,y,z] */
    double sum;
    double xp,yp,zp,rr;
    double xr,yr,zr;
    int m;
    usint i,j,k;
    
    rr = x * x + y * y + z * z;
    rr = sqrt( rr );
    
    if( rr < numeric_limits<double>::epsilon() ) {
        memset( ylm, 0, (l+l+1) * sizeof(double) );
        if( l == 0 ) ylm[0] = 1.0 / sqrt_4pi; // by convention
        return ;
    };
    
    xr = x / rr;
    yr = y / rr;
    zr = z / rr;
    
    for(m=-l; m<=l; m++) {
        sum = 0.0;
        xp  = 1.0;
        for(i=0; i<=l; i++) {
            yp = 1.0;
            for(j=0; j<=l-i; j++) {
                zp = 1.0;
                for(k=0; k<=l-i-j; k++) {
                    sum += NoNormCalcClmR( l, m, i, j, k )
                         * xp * yp * zp;
                    zp = zp * zr;
                };
                yp = yp * yr;
            };
            xp = xp * xr;
        };
        ylm[m+l] = sum;
    };
    
    /* */
};

void CompAngFac( const double kx, 
                 const double ky,
                 const double kz,
                 const usint lP, const usint lG, 
                 const usint L , double** aij,
                 double* ang_fac ) {
    /* angular factors for the projector operators */
    int M;
    usint a,b,c;
    usint ia,ja,ka,mi;
    double s1,s2;
    arr1d <double> ylm( 2*L + 1 );
    arr1d <double> ang( 2*L + 1 );
    
    CompSolidHarm( kx, ky, kz, L, ylm.v );
    
    ia = 0;
    ja = 0;
    for(mi=0; mi<crt_siz[lG]; mi++) {
        ka = lG - ia - ja;
        
        s1 = 0.0;
        for(a=0; a<=lP; a++) for(b=0; b<=lP-a; b++) for(c=0; c<=lP-a-b; c++) {
            s2 = 0.0;
            CartSolidAng( ia + a, ja + b, ka + c, L, aij, ang.v );
            for(M=-L; M<=L; M++) s2 += ang[M+L] * ylm[M+L];
            s1 += s2 * NoNormCalcClmR( lP, 0, a, b, c );
        };
        
        ang_fac[mi] = s1;
            
        ja++;
        if( ja > lG - ia ) { ja = 0; ia++; };
    };

    /* */
};

void CAPAngFac( const double kx, const double ky, const double kz,
                const usint lG, const usint Lmax , 
                double** aij  , double** ang_fac ) {
    /* angular factors for the atomic CAP integrals */
    int M;
    usint ia,ja,ka,mi;
    usint L,twoL;
    twoL = 2*Lmax + 1;
    
    arr1d <double> ylm( twoL );
    arr1d <double> ang( twoL );
    
    for(L=0; L<=Lmax; L++) {
        ylm.zero(); ang.zero();
        CompSolidHarm( kx, ky, kz, L, ylm.v );
        
        ia = 0;
        ja = 0;
        for(mi=0; mi<crt_siz[lG]; mi++) {
            ka = lG - ia - ja;
            
            CartSolidAng( ia, ja, ka, L, aij, ang.v );
            for(M=0; M<twoL; M++) 
                ang_fac[L][mi] += ang[M] * ylm[M];
            
            ja++;
            if( ja > lG - ia ) { ja = 0; ia++; };
        };
    };
    /* */
};

void Faddeeva( const double xi, const double yi,
               double & u, double & v, bool flag ) {
    /* Faddeeva function for complex arguments;
     * based on: G. M. Poppe and C. M. J. Wijers,
     * ACM Trans. Math. Soft. 16, 38 (1990).
     */
    double factor,rmaxreal,rmaxexp;
    double xabs,yabs,x,y,qrho,xabsq,xquad,yquad;
    double xsum,ysum,xaux,daux,u1,u2,v1,v2;
    double h,h2,qlambda,rx,sx,ry,sy,tx,ty,c,w1;
    bool a,b;
    int n,j,kapn,nu,np1;
    
    factor   = 2.0 / sqrt( M_PI );
    rmaxreal = sqrt( numeric_limits<double>::max() );
    rmaxexp  = 2.0 * log( rmaxreal ) - M_LN2;
    
    flag = false;
    xabs = abs( xi );
    yabs = abs( yi );
    x    = xabs / 6.3;
    y    = yabs / 4.4;
    
    if( xabs > rmaxreal || yabs > rmaxreal ) {
        flag = true; return ;
    };
    
    qrho  = x * x + y * y;
    xabsq = xabs * xabs;
    xquad = xabsq - yabs * yabs;
    yquad = 2.0 * xabs * yabs;
    a     = qrho < 0.085264;
    
    if( a ) {
        qrho = ( 1.0 - 0.85 * y ) * sqrt( qrho );
        n    = (int) nearbyint( 6.0 + 72.0 * qrho );
        j    = 2*n + 1;
        xsum = 1.0 / j;
        ysum = 0.0;
        
        for(int i=n; i>=1; i--) {
            j    = j - 2;
            xaux = ( xsum * xquad - ysum * yquad ) / i;
            ysum = ( xsum * yquad + ysum * xquad ) / i;
            xsum = xaux + 1.0 / j;
        };
        
        u1   = -factor * ( xsum * yabs + ysum * xabs ) + 1.0;
        v1   = +factor * ( xsum * xabs - ysum * yabs );
        daux = +exp( -xquad );
        u2   = +daux * cos( yquad );
        v2   = -daux * sin( yquad );
        
        u    = u1 * u2 - v1 * v2;
        v    = u1 * v2 + v1 * u2;
    } else {
        h2      = 0.0;
        qlambda = 0.0;
        if( qrho > 1.0 ) {
            h    = 0.0;
            kapn = 0;
            qrho = sqrt( qrho );
            nu   = (int) nearbyint( 3.0 + ( 1442.0 / ( 26.0 * qrho + 77.0 ) ) );
        } else {
            qrho = ( 1.0 - y ) * sqrt( 1.0 - qrho );
            h    = 1.88 * qrho;
            h2   = 2.0 * h;
            kapn = (int) nearbyint(  7.0 + 34.0 * qrho );
            nu   = (int) nearbyint( 16.0 + 26.0 * qrho );
        };
        
        b = h > 0.0;
        
        if( b ) qlambda = pow( h2, kapn );
        
        rx = 0.0;
        ry = 0.0;
        sx = 0.0;
        sy = 0.0;
        
        for(n=nu; n>=0; n--) {
            np1 = n + 1;
            tx  = yabs + h + np1 * rx;
            ty  = xabs - np1 * ry;
            c   = 0.5 / ( tx * tx + ty * ty );
            rx  = c * tx;
            ry  = c * ty;
            if( b && ( n <= kapn ) ) {
                tx = qlambda + sx;
                sx = rx * tx - ry * sy;
                sy = ry * tx + rx * sy;
                qlambda = qlambda / h2;
            };
        };
        
        if( abs( h ) < numeric_limits<double>::epsilon() ) {
            u = factor * rx;
            v = factor * ry;
        } else {
            u = factor * sx;
            v = factor * sy;
        };
        
        if( abs( yabs ) < numeric_limits<double>::epsilon() ) 
            u = exp( -xabs * xabs );
    };
    
    if( yi < 0.0 ) {
        if( a ) {
            u2 = 2.0 * u2;
            v2 = 2.0 * v2;
        } else {
            xquad = -xquad;
            if( xquad > rmaxexp ) {
                flag = true; return ;
            };
            w1 = +2.0 * exp( xquad );
            u2 = +w1 * cos( yquad );
            v2 = -w1 * sin( yquad );
        };
        
        u = u2 - u;
        v = v2 - v;
        if( xi > 0.0 ) v = -v;
    } else {
        if( xi < 0.0 ) v = -v;
    };
    /* */
};

void Faddeeva( const quad xi, const quad yi,
               quad & u, quad & v, bool flag ) {
    /* Faddeeva function for complex arguments,
     * quadruple precision subroutine version;
     * based on: G. M. Poppe and C. M. J. Wijers,
     * ACM Trans. Math. Soft. 16, 38 (1990).
     */
    quad factor,rmaxreal,rmaxexp;
    quad xabs,yabs,x,y,qrho,xabsq,xquad,yquad;
    quad xsum,ysum,xaux,daux,u1,u2,v1,v2;
    quad h,h2,qlambda,rx,sx,ry,sy,tx,ty,c,w1;
    bool a,b;
    int n,j,kapn,nu,np1;
    
    factor   = 2.0q / sqrtq( M_PIq );
    rmaxreal = sqrt( numeric_limits<double>::max() );
    rmaxexp  = 2.0q * logq( rmaxreal ) - M_LN2q;
    
    flag = false;
    xabs = fabsq( xi );
    yabs = fabsq( yi );
    x    = xabs / 6.3q;
    y    = yabs / 4.4q;
    
    if( xabs > rmaxreal || yabs > rmaxreal ) {
        flag = true; return ;
    };
    
    qrho  = x * x + y * y;
    xabsq = xabs * xabs;
    xquad = xabsq - yabs * yabs;
    yquad = 2.0q * xabs * yabs;
    a     = qrho < 0.085264q;
    
    PrintQuadMath( qrho );
    
    if( a ) {
        qrho = ( 1.0q - 0.85q * y ) * sqrtq( qrho );
        n    = (int) nearbyint( 16.0 + 132.0 * (double) qrho );
        j    = 2*n + 1;
        xsum = 1.0q / j;
        ysum = 0.0q;
        
        for(int i=n; i>=1; i--) {
            j    = j - 2;
            xaux = ( xsum * xquad - ysum * yquad ) / i;
            ysum = ( xsum * yquad + ysum * xquad ) / i;
            xsum = xaux + 1.0q / j;
        };
        
        u1   = -factor * ( xsum * yabs + ysum * xabs ) + 1.0q;
        v1   = +factor * ( xsum * xabs - ysum * yabs );
        daux = +expq( -xquad );
        u2   = +daux * cosq( yquad );
        v2   = -daux * sinq( yquad );
        
        u    = u1 * u2 - v1 * v2;
        v    = u1 * v2 + v1 * u2;
    } else {
        h2      = 0.0q;
        qlambda = 0.0q;
        if( qrho > 1.0q ) {
            h    = 0.0q;
            kapn = 0;
            qrho = sqrtq( qrho );
            nu   = (int) nearbyint( 12.0 + ( 2800.0 / ( 26.0 * (double) qrho + 77.0 ) ) );
        } else {
            qrho = ( 1.0q - y ) * sqrtq( 1.0q - qrho );
            h    = 1.88q * qrho;
            h2   = 2.0q * h;
            kapn = (int) nearbyint(  7.0 + 34.0 * (double) qrho );
            nu   = (int) nearbyint( 16.0 + 26.0 * (double) qrho );
        };
        
        b = h > 0.0;
        
        if( b ) qlambda = powq( h2, kapn );
        
        rx = 0.0;
        ry = 0.0;
        sx = 0.0;
        sy = 0.0;
        
        for(n=nu; n>=0; n--) {
            np1 = n + 1;
            tx  = yabs + h + np1 * rx;
            ty  = xabs - np1 * ry;
            c   = 0.5q / ( tx * tx + ty * ty );
            rx  = c * tx;
            ry  = c * ty;
            if( b && ( n <= kapn ) ) {
                tx = qlambda + sx;
                sx = rx * tx - ry * sy;
                sy = ry * tx + rx * sy;
                qlambda = qlambda / h2;
            };
        };
        
        if( fabsq( h ) < numeric_limits<double>::epsilon()
                       * numeric_limits<double>::epsilon() ) {
            u = factor * rx;
            v = factor * ry;
        } else {
            u = factor * sx;
            v = factor * sy;
        };
        
        if( fabsq( yabs ) < numeric_limits<double>::epsilon()
                          * numeric_limits<double>::epsilon() ) 
            u = expq( -xabs * xabs );
    };
    
    if( yi < 0.0q ) {
        if( a ) {
            u2 = 2.0q * u2;
            v2 = 2.0q * v2;
        } else {
            xquad = -xquad;
            if( xquad > rmaxexp ) {
                flag = true; return ;
            };
            w1 = +2.0 * expq( xquad );
            u2 = +w1  * cosq( yquad );
            v2 = -w1  * sinq( yquad );
        };
        
        u = u2 - u;
        v = v2 - v;
        if( xi > 0.0q ) v = -v;
    } else {
        if( xi < 0.0q ) v = -v;
    };
    /* */
};

void CAPRadFac( const int lmax, const int nmax,
                const double k, const double p,
                double** rad_fac ) {
    /* radial factors for the atomic CAP integrals */
    if( nmax == 0 ) return ;
    bool flag;
    int l,n,ln;
    double sq_p,j01,j02,j03,fac,last;
    double fed_re,fed_im,arg_re,arg_im;
    cdouble fed,dfed,tmp,face;
    
    /* special case for small k < 1 */
    if( k * k / p < 1.0 ) {
        sq_p   = sqrt( p );
        fed_re = exp( -p );
        fed_im = 2.0 * p;
        arg_re = -k * k / 2.0;
        
        arr1d <double> gn( lmax + nmax + 1 );
        gn[0] = 0.25 * sqrt_4pi * erfc( sq_p ) / sq_p;
        gn[1] = fed_re / fed_im;
        
        for(n=2; n<=nmax+lmax; n++)
            gn[n] = ( ( n - 1 ) * gn[n-2] + fed_re ) / fed_im;
        
        for(l=0; l<=lmax; l++) for(n=l+1; n<=nmax; n++) {
            ln     = l + n;
            last   = 0.0;
            arg_im = 1.0;
            fac    = gn[ln];
            
            for(int m=0; m<=50; m++) {
                j01    = arg_im * fac / dfact[l+m+1] / fact[m];
                last  += j01;
                if( abs( j01 / last ) < numeric_limits<double>::epsilon() ) break;
                fac    = ( ( ln + m + m + 1 ) * fac + fed_re ) / fed_im;
                arg_im = arg_im * arg_re;
            };
            
            rad_fac[l][n] = last * pow( k, l );
        };
        
        return ;
    };
    
    /* part 1: J_{01} and J_{02} from scratch */
    flag   = false;
    sq_p   = sqrt( p );
    arg_re = k / ( 2.0 * sq_p );
    arg_im = sq_p;
    
    Faddeeva( arg_re, arg_im, fed_re, fed_im, flag );
    fed  = fed_re + I * fed_im;
    dfed = 2.0 * ( 2.0 * I / sqrt_4pi - fed * ( arg_re + I * arg_im ) );
    
    face = exp( I * k - p );
    tmp  = fed * face;
    fac  = 0.25 * sqrt_4pi / sq_p / k;
    j01  = fac * tmp.imag();
    
    tmp  = face * ( I * fed + 0.5 / sq_p * dfed );
    j02  = -fac * tmp.real();
    
    /* part 2: recursion to increase n */
    rad_fac[0][0] = 0.0;
    rad_fac[0][1] = j01;
    if( nmax == 1 ) return ;
    
    rad_fac[0][2] = j02;
    if( nmax == 2 ) return ;
    
    fed_re = exp( -p ) / ( 4.0 * p * p );
    fed_im = 1.0 / ( 2.0 * p );
    arg_im = sin( k );
    arg_re = ( cos( k ) + arg_im / k / fed_im ) * fed_re;
    arg_im = -arg_im * fed_re / k;
    fed_re = k * fed_im;
    fed_re = fed_re  * fed_re;
    
    j03 = ( 0.5 / p - fed_re ) * j01 + arg_re;
    rad_fac[0][3] = j03;
    
    for(n=1; n<=nmax-3; n++) 
        rad_fac[0][n+3] = ( fed_im * ( n + n + 1 ) - fed_re ) * rad_fac[0][n+1]
                        - n * ( n - 1 ) * rad_fac[0][n-1] * fed_im * fed_im
                        + arg_re + arg_im * n;
    n = nmax - 2;
    last = ( fed_im * ( n + n + 1 ) - fed_re ) * rad_fac[0][n+1]
         - n * ( n - 1 ) * rad_fac[0][n-1] * fed_im * fed_im
         + arg_re + arg_im * n;
    
    /* part 3: build l at the cost of n */
    fed_re = 2.0 * p;
    fed_im = exp( -p ) * sin( k ) / k;
    
    if( lmax > 0 ) {
        for(n=2; n<nmax; n++)
            rad_fac[1][n] = ( fed_im + n * rad_fac[0][n-1] 
                          -   fed_re * rad_fac[0][n+1] ) / k;
    
        rad_fac[1][nmax]  = ( nmax * rad_fac[0][nmax-1] 
                          -   fed_re * last + fed_im ) / k;
    };
    
    for(l=1; l<lmax; l++) for(n=l+2; n<=nmax; n++)
        rad_fac[l+1][n] = ( l + l + 1 ) * rad_fac[l][n-1] / k - rad_fac[l-1][n];
    
    /* */
};











