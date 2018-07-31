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
#include "../headers/auxfun1.h"

using namespace std;

void WriteDown( cdouble* data, const int len, std::ofstream & ofs ) {
    /* write down a matrix to the disk */
    double val;
    
    for(int i=0; i<len; i++) {
        val = data[i].real();
        ofs.write(reinterpret_cast<char*>(&val),sizeof(double));
    };
    
    for(int i=0; i<len; i++) {
        val = data[i].imag();
        ofs.write(reinterpret_cast<char*>(&val),sizeof(double));
    };
    /* */
};

void TransToSpher1D( cdouble* shl_crt, 
                     cdouble* shl_sph, 
                     const int li ) {
    /* transform to the spherical representation */
    int shi,shgi;
    int mi,mig;
    int ia,ja,ka;
    cdouble tmp;
    
    shi  = sph_siz[li];
    shgi = crt_siz[li];
    
    for(mi=0; mi<shi; mi++) {
        tmp = 0;
        
        ia  = 0;
        ja  = 0;
        for(mig=0; mig<shgi; mig++) {
            ka = li - ia - ja;
            tmp += shl_crt[ mig ]
                 *  clmr[li][mi][mig];
        };
        
        shl_sph[ mi ] = tmp;
        
        ja++;
        if( ja > li - ia ) { ja = 0; ia++; };
    };
    /* */
};

void TransToSpher( cdouble* shl_crt, cdouble* shl_sph, 
                   const int li, const int lj ) {
    /* transform to the spherical representation */
    int shi,shj,shgi,shgj;
    int mi,mj,mig,mjg;
    int ia,ib,ja,jb,ka,kb;
    cdouble tmp;
    
    shi  = sph_siz[li];
    shgi = crt_siz[li];
    
    shj  = sph_siz[lj];
    shgj = crt_siz[lj];
    
    for(mi=0; mi<shi; mi++) for(mj=0; mj<shj; mj++) {
        ia  = 0;
        ja  = 0;
        tmp = 0;
        for(mig=0; mig<shgi; mig++) {
            ka = li - ia - ja;
            ib = 0;
            jb = 0;
            for(mjg=0; mjg<shgj; mjg++) {
                kb = lj - ib - jb;
                tmp += shl_crt[ mig * shgj + mjg ]
                     *  clmr[li][mi][mig]
                     *  clmr[lj][mj][mjg];
                jb++;
                if( jb > lj - ib ) { jb = 0; ib++; };
            };
            ja++;
            if( ja > li - ia ) { ja = 0; ia++; };
        };
        shl_sph[ mi * shj + mj ] = tmp;
    };
    /* */
};

void GamessHolyOrder( cdouble* shl_crt, cdouble* shl_crt_dum, 
                      const int li, const int lj ) {
    /* transform to the Gamess shell indexing */
    int shgi,shgj;
    int mi,mj,mi_gam,mj_gam;
    
    shgi = crt_siz[li];
    shgj = crt_siz[lj];
    
    for(mi=0; mi<shgi; mi++) {
        mi_gam = pos_change_gamess( li, mi );
        for(mj=0; mj<shgj; mj++) {
            mj_gam = pos_change_gamess( lj, mj );
            shl_crt_dum[ mi_gam * shgj + mj_gam ] = shl_crt[ mi * shgj + mj ];
        };
    };
    /* */
};

void GamessHolyOrder2E( FullC <cdouble> & full_trans, 
                        FullC <cdouble> & full_trans_gam,
                        const int li, const int lj,
                        const int lk, const int ll,
                        const bool overwrite ) {
    /* transform to the Gamess shell indexing */
    int shgi,shgj,shgk,shgl;
    int mi,mj,mk,ml;
    int mi_gam,mj_gam,mk_gam,ml_gam;
    
    shgi = crt_siz[li];
    shgj = crt_siz[lj];
    shgk = crt_siz[lk];
    shgl = crt_siz[ll];
    
    for(mi=0; mi<shgi; mi++) {
        mi_gam = pos_change_gamess( li, mi );
        for(mj=0; mj<shgj; mj++) {
            mj_gam = pos_change_gamess( lj, mj );
            for(mk=0; mk<shgk; mk++) {
                mk_gam = pos_change_gamess( lk, mk );
                for(ml=0; ml<shgl; ml++) {
                    ml_gam = pos_change_gamess( ll, ml );
                    full_trans_gam.v[ mi_gam ][ mj_gam ][ mk_gam ][ ml_gam ] = 
                        full_trans.v[ mi ][ mj ][ mk ][ ml ];
                };
            };
        };
    };
    
    if( !overwrite ) return ;
    
    for(mi=0; mi<shgi; mi++) for(mj=0; mj<shgj; mj++) 
        for(mk=0; mk<shgk; mk++) for(ml=0; ml<shgl; ml++)
            full_trans.v[ mi ][ mj ][ mk ][ ml ] = 
                full_trans_gam.v[ mi ][ mj ][ mk ][ ml ];
    /* */
};

usint FindPos( usint ia_, usint ja_, usint ka_ ) {
    /* find the position in the internal ordering */
    usint ia,ja,ka,la,shgi,mi;
    la   = ia_ + ja_ + ka_;
    shgi = ( la + 1 ) * ( la + 2 ) / 2;
    
    ia = 0;
    ja = 0;
    for(mi=0; mi<shgi; mi++) {
        ka = la - ia - ja;
        if( ia == ia_ &&
            ja == ja_ &&
            ka == ka_  ) return mi;
        
        ja++;
        if( ja > la - ia ) { ja = 0; ia++; };
    };
    /* */
};

void CalcRNucA( RInts1E <cdouble> & R ) {
    /* R_tuv for nuclear attraction integrals */
    cdouble ex1,ex2,fac0,bigT,twop;
    cdouble fax,fay,faz,fix,fiy,fiz;
    
    fax  = R.Px - R.Cx + I * R.kPx / ( 2.0 * R.aP );
    fay  = R.Py - R.Cy + I * R.kPy / ( 2.0 * R.aP );
    faz  = R.Pz - R.Cz + I * R.kPz / ( 2.0 * R.aP );
    
    bigT = fax * fax + fay * fay + faz * faz;
    bigT = bigT * R.aP;
    twop = -2.0 * R.aP;
    
    if( bigT.real() > 0.0 ) {
        ex1  = R.kPx * R.kPx + R.kPy * R.kPy + R.kPz * R.kPz;
        ex1  = exp( -ex1 / ( 4.0 * R.aP ) );
        ex2  = exp( I * ( R.kPx * R.Px + R.kPy * R.Py + R.kPz * R.Pz ) );
        fac0 = +2.0 * M_PI * ex1 * ex2 / R.aP;
    } else {
        ex1  = ( R.Px - R.Cx ) * ( R.Px - R.Cx )
             + ( R.Py - R.Cy ) * ( R.Py - R.Cy )
             + ( R.Pz - R.Cz ) * ( R.Pz - R.Cz );
        ex1  = exp( -R.aP * ex1 );
        ex2  = exp( I * ( R.kPx * R.Cx + R.kPy * R.Cy + R.kPz * R.Cz ) );
        fac0 = +2.0 * M_PI * ex1 * ex2 / R.aP;
    };
    
    fix  = I * R.kPx;
    fiy  = I * R.kPy;
    fiz  = I * R.kPz;

    BoysFunction <cdouble> boys( R.lsum, bigT );
    BoysFn( boys );
    
    /* start with R_{000j} - scaled Boys function */
    for(int j=0; j<=R.lsum; j++) {
        R.rtuvj[0][0][0][j] = fac0 * boys.fn[j];
        fac0 = fac0 * twop;
    };
    
    /* R_{00vj}: grow v at cost of j */
    if( R.lsum > 0 ) {
        for(int j=0; j<R.lsum; j++) 
            R.rtuvj[0][0][1][j] = faz * R.rtuvj[0][0][0][j+1]
                                + fiz * R.rtuvj[0][0][0][j  ];
    };
    
    for(int v=1; v<R.lsum; v++) for(int j=0; j<R.lsum-v; j++) 
        R.rtuvj[0][0][v+1][j] = faz * R.rtuvj[0][0][v][j+1]
                              + fiz * R.rtuvj[0][0][v][j  ]
                              + (double)v * R.rtuvj[0][0][v-1][j+1];
    
    /* R_{0uvj}: grow u at cost of j, constant v */
    if( R.lsum > 0 ) {
        for(int v=0; v<=R.lsum; v++) for(int j=0; j<R.lsum-v; j++) 
            R.rtuvj[0][1][v][j] = fay * R.rtuvj[0][0][v][j+1]
                                + fiy * R.rtuvj[0][0][v][j  ];
    };
        
    for(int u=1; u<R.lsum; u++) for(int v=0; v<=R.lsum; v++) 
        for(int j=0; j<R.lsum-v-u; j++) 
            R.rtuvj[0][u+1][v][j] = fay * R.rtuvj[0][u][v][j+1]
                                  + fiy * R.rtuvj[0][u][v][j  ]
                                  + (double)u * R.rtuvj[0][u-1][v][j+1];
        
    /* R_{tuvj}: grow t at cost of j, constant u and v */
    if( R.lsum > 0 ) {
        for(int u=0; u<=R.lsum; u++) for(int v=0; v<=R.lsum; v++) 
            for(int j=0; j<R.lsum-v-u; j++) 
                R.rtuvj[1][u][v][j] = fax * R.rtuvj[0][u][v][j+1]
                                    + fix * R.rtuvj[0][u][v][j  ];
    };
        
    for(int t=1; t<R.lsum; t++) for(int u=0; u<=R.lsum; u++) 
        for(int v=0; v<=R.lsum; v++) for(int j=0; j<R.lsum-v-u-t; j++) 
            R.rtuvj[t+1][u][v][j] = fax * R.rtuvj[t][u][v][j+1]
                                  + fix * R.rtuvj[t][u][v][j  ]
                                  + (double)t * R.rtuvj[t-1][u][v][j+1];
    R.strip();
    /* */
};

void CalcROvrl( RInts1E <cdouble> & R, const cdouble fac ) {
    /* R_tuv for overlap integrals */
    if( R.lsum == 0 ) {
        R.rtuv[0][0][0] = fac;
        return ;
    };
    
    cdouble* kpxt = new cdouble[R.lsum+1];
    cdouble* kpyu = new cdouble[R.lsum+1];
    cdouble* kpzv = new cdouble[R.lsum+1];
    
    memset( kpxt, 0, ( R.lsum + 1 ) * sizeof( cdouble ) );
    memset( kpyu, 0, ( R.lsum + 1 ) * sizeof( cdouble ) );
    memset( kpzv, 0, ( R.lsum + 1 ) * sizeof( cdouble ) );
    
    kpxt[0] = 1.0;
    kpyu[0] = 1.0;
    kpzv[0] = 1.0;
    
    kpxt[1] = I * R.kPx;
    kpyu[1] = I * R.kPy;
    kpzv[1] = I * R.kPz;
    
    for(int l=2; l<=R.lsum; l++) {
        kpxt[l] = kpxt[l-1] * kpxt[1];
        kpyu[l] = kpyu[l-1] * kpyu[1];
        kpzv[l] = kpzv[l-1] * kpzv[1];
    };
    
    for(int t=0; t<=R.lsum; t++) for(int u=0; u<=R.lsum; u++) {
        for(int v=0; v<=R.lsum-t-u; v++)
            R.rtuv[t][u][v] = fac * kpxt[t] * kpyu[u] * kpzv[v];
    };
    
    /* */
    delete []kpxt;
    delete []kpyu;
    delete []kpzv;
};

void CalcEijt( ECoefs <double> & E ) {
    /* calculate the momentum transfer coefficients */
    int i1,j1;
    double abp,ap,bp;
    abp = E.a * E.b / E.p;
    ap  = E.a / E.p;
    bp  = E.b / E.p;
    
    E.v[0][0][0] = 1.0;
    
    if( E.jmax > 0 ) E.v[0][1][0] = ap * E.AB * E.v[0][0][0];
    for(int j=2; j<=E.jmax; j++) {
        j1 = j - 1;
        E.v[0][j][0] = ap * E.AB * E.v[0][j1][0] + j1 * E.v[0][j-2][0] / E.twop;
    };
    
    for(int j=1; j<=E.jmax; j++) for(int k=1; k<=j; k++)
        E.v[0][j][k] = j * E.v[0][j-1][k-1] / ( E.twop * k );
    
    if( E.imax > 0 ) {
        E.v[1][0][0] = - bp * E.AB * E.v[0][0][0]; 
        for(int j=1; j<=E.jmax; j++)
            E.v[1][j][0] = -bp * E.AB * E.v[0][j][0] + j * E.v[0][j-1][0] / E.twop;
    };
    
    for(int i=2; i<=E.imax; i++) {
        i1 = i - 1;
        E.v[i][0][0] = - bp * E.AB * E.v[i1][0][0] + i1 * E.v[i-2][0][0] / E.twop;
        for(int j=1; j<=E.jmax; j++) 
            E.v[i][j][0] = -bp * E.AB * E.v[i1][j][0] + ( i1 * E.v[i-2][j][0] + j * E.v[i1][j-1][0] ) / E.twop;
    };
    
    for(int i=1; i<=E.imax; i++) {
        i1 = i - 1;
        for(int k=1; k<=i; k++)
            E.v[i][0][k] = i * E.v[i1][0][k-1] / ( E.twop * k );
    };
    
    for(int i=1; i<=E.imax; i++) for(int j=1; j<=E.jmax; j++) for(int k=1; k<=i+j; k++)
          E.v[i][j][k] = ( i * E.v[i-1][j][k-1] + j * E.v[i][j-1][k-1] ) / ( E.twop * k );
    
    /* */
};

void RenormContr( GaussC& gau_tmp, const string norm_name ) {
    /* renormalise the Gaussian contraction */
    double num,den,al,fra,norm,tot_norm;
    double d1_re,d2_re,d1_im,d2_im;
    num = 0;
    den = 0;
    fra = 0;
    
    for(int i=0; i<gau_tmp.clen; i++) {
        for( int j=0; j<gau_tmp.clen; j++) {
            d1_re = gau_tmp.dA_re[i] * pow( 2.0 * gau_tmp.alphaA[i] / M_PI, 0.75 );
            d1_re = d1_re * pow( 2.0 * sqrt( gau_tmp.alphaA[i] ), gau_tmp.lA ) / sqrt( dfact[gau_tmp.lA] );
            
            d2_re = gau_tmp.dA_re[j] * pow( 2.0 * gau_tmp.alphaA[j] / M_PI, 0.75 );
            d2_re = d2_re * pow( 2.0 * sqrt( gau_tmp.alphaA[j] ), gau_tmp.lA ) / sqrt( dfact[gau_tmp.lA] );
            
            d1_im = gau_tmp.dA_im[i] * pow( 2.0 * gau_tmp.alphaA[i] / M_PI, 0.75 );
            d1_im = d1_im * pow( 2.0 * sqrt( gau_tmp.alphaA[i] ), gau_tmp.lA ) / sqrt( dfact[gau_tmp.lA] );
            
            d2_im = gau_tmp.dA_im[j] * pow( 2.0 * gau_tmp.alphaA[j] / M_PI, 0.75 );
            d2_im = d2_im * pow( 2.0 * sqrt( gau_tmp.alphaA[j] ), gau_tmp.lA ) / sqrt( dfact[gau_tmp.lA] );
            
            num = d1_re * d2_re + d1_im * d2_im;
            den = pow( gau_tmp.alphaA[i] + gau_tmp.alphaA[j], 1.5 + gau_tmp.lA );
            fra += num / den;
        };
    };
    
    norm = pow( M_PI, -0.75 ) * pow( 2.0, 0.5 * gau_tmp.lA ) / sqrt( fra );
    for(int i=0; i<gau_tmp.clen; i++) {
        gau_tmp.dA_re[i] = gau_tmp.dA_re[i] * norm * pow( 2.0 * gau_tmp.alphaA[i] / M_PI, 0.75 );
        gau_tmp.dA_re[i] = gau_tmp.dA_re[i] * pow( 2.0 * sqrt( gau_tmp.alphaA[i] ), gau_tmp.lA ) / sqrt( dfact[gau_tmp.lA] );
        
        gau_tmp.dA_im[i] = gau_tmp.dA_im[i] * norm * pow( 2.0 * gau_tmp.alphaA[i] / M_PI, 0.75 );
        gau_tmp.dA_im[i] = gau_tmp.dA_im[i] * pow( 2.0 * sqrt( gau_tmp.alphaA[i] ), gau_tmp.lA ) / sqrt( dfact[gau_tmp.lA] );
    };
    
    ofstream ofs_norm( norm_name, std::ios::out|std::ios::binary|std::ios_base::app );
    if(!ofs_norm.is_open())
    {
        cout << "Cannot open norm1E!\n";
        exit(EXIT_FAILURE);
    }

    for(int mi=0; mi<crt_siz[gau_tmp.lA]; mi++) {
        fra = norm / xyz_norm[gau_tmp.lA][mi];
        ofs_norm.write(reinterpret_cast<char*>(&fra),sizeof(double));
    };
    ofs_norm.close();
    
    /* */
};

double CalcClmR( const int l, const int m, const int lx, const int ly, const int lz ) {
    /* cartesian to spherical transformation coefficients */
    if( lx + ly + lz != l ) return 0.0;
    if( (lx+ly-m) % 2 == 0 && (lx+ly) >= abs(m) ) {
        double phi, theta, coef;
        int ll1, ll2, am;
        phi   = fact[2*lx] * fact[2*ly] / fact[lx] / fact[ly] / fact[lx+ly] / pow ( 4.0, lx + ly );
        theta = 0.0;
        for(int n=0; n<=lx+ly; n++) theta += belt[n] * binom[lx+ly][n] / (lz + n + 0.5);
        phi  *= theta;
        phi   = sqrt(phi) * omega[l][abs(m)] / ( pow( 2.0, l ) * fact[l] );
        am    = abs(m);
        ll1   = (lx+ly-am)/2;
        coef  = 0.0;
        
        if( m < 0  ) {
            if( (lx-m) % 2 == 0 ) return 0.0;
            for(int a = ll1; a <= (l-am)/2; a++) {
                for(int c = max( 0, (lx-am+1)/2 ); c <= min(lx/2,ll1); c++) {
                    coef += belt[abs(c+a+(am-lx-1)/2)] *
                            binom[l][a] * binom[a][ll1] * binom[ll1][c] * binom[am][2*c+am-lx] *
                            fact[2*l-2*a] / fact[l-am-2*a];
                };
            };
            return coef * phi * M_SQRT2;
        };
        
        if( m == 0 ) {
            if(lx % 2 == 0) ll2 = lx;
            else            ll2 = lx+1;
            for(int a = ll1; a <= l/2; a++) {
                for(int c = max( 0, ll2/2 ); c <= min(lx/2,ll1); c++) {
                    coef += belt[abs(c+a+(ll2)/2)] *
                            binom[l][a] * binom[a][ll1] * binom[ll1][c] * binom[am][2*c-lx] *
                            fact[2*l-2*a] / fact[l-2*a];
                };
            };
            return coef * phi;
        };
        
        if( m > 0  ) {
            if( (lx-m+1) % 2 == 0 ) return 0.0;
            for(int a = ll1; a <= (l-am)/2; a++) {
                for(int c = max( 0, (lx-am)/2 ); c <= min(lx/2,ll1); c++) {
                    coef += belt[abs(c+a+(am-lx)/2)] *
                            binom[l][a] * binom[a][ll1] * binom[ll1][c] * binom[am][2*c+am-lx] *
                            fact[2*l-2*a] / fact[l-am-2*a];
                };
            };
            return coef * phi * M_SQRT2;
        };
        
    };
    return 0.0;
    /* */
};

double NoNormCalcClmR( const int l, const int m, const int lx, const int ly, const int lz ) {
    /* cartesian to spherical transformation coefficients */
    if( lx + ly + lz != l ) return 0.0;
    if( m > 0 ) {
        int bs,am,lima,lim1b,lim2b,twol,twoa,twob;
        double val,fac,faci;
        am  = abs( m );
        bs  = l - lz - am;
        twol= 2 * l;
        if( bs % 2 == 1 || bs < 0 ) return 0.0;
        bs  = bs / 2;
        fac = M_SQRT2 * omega[l][m] * belt[am] / ( pow( 2.0, l ) * fact[l] );

        lima  = l - am;
        lima  = ( lima  % 2 == 0 ) ? lima / 2 : ( lima - 1 ) / 2;
        lim1b = lx - am;
        lim1b = ( lim1b % 2 == 0 ) ? lim1b / 2 : ( lim1b + 1 ) / 2;
        lim1b = max( lim1b, 0  );
        lim2b = ( lx    % 2 == 0 ) ? lx / 2 : ( lx - 1 ) / 2;
        lim2b = min( lim2b, bs );

        val = 0.0;
        for(int a=max(bs,0); a<=lima; a++) {
            twoa = 2 * a;
            faci = binom[l][a] * belt[a] * fact[twol-twoa]
                 * binom[a][bs] / fact[l-am-twoa];
            for(int b=lim1b; b<=lim2b; b++) {
                twob = 2 * b + am - lx;
            if( twob % 2 == 1 ) continue;
                val += binom[bs][b] * binom[am][twob]
                     * belt[twob/2] * faci;
            };
        };
        return belt[am] * val * fac / sqrt( 2.0 * M_PI );
    };
    if( m == 0 ) {
        int lxx,lxy,lim1a,lim2a,twoa,twol;
        lxx = lx;
        lxy = lx + ly;
        if( lxx % 2 == 1 ) return 0.0;
        if( lxy % 2 == 1 ) return 0.0;
        twol = 2 * l;
        lxx /= 2;
        lxy /= 2;
        double val,fac;
        fac = sqrt((2*l+1)/2.)* binom[lxy][lxx] / ( pow( 2.0, l ) * fact[l] );
        lim1a = lxy;
        lim2a = ( l % 2 == 0 ) ? l / 2 : ( l - 1 ) / 2;

        val = 0.0;
        for(int a=lim1a; a<=lim2a; a++) {
            twoa = 2 * a;
            val += binom[l][a] * belt[a] * fact[twol-twoa]* binom[a][lxy] / fact[l-twoa];
        };
        return  val * fac / sqrt( 2.0 * M_PI );
    };
    if( m < 0 ) {
        int bs,am,lima,lim1b,lim2b,twol,twoa,twob;
        double val,fac,faci;
        am  = abs( m );
        bs  = l - lz - am;
        twol= 2 * l;
        if( bs % 2 == 1 || bs < 0 ) return 0.0;
        bs  = bs / 2;
        fac = - M_SQRT2 * omega[l][am] * belt[am] / ( pow( 2.0, l ) * fact[l] );
   
        lima  = l - am;
        lima  = ( lima  % 2 == 0 ) ? lima / 2 : ( lima - 1 ) / 2;
        lim1b = lx - am;
        lim1b = ( lim1b % 2 == 0 ) ? lim1b / 2 : ( lim1b + 1 ) / 2;
        lim1b = max( lim1b, 0  );
        lim2b = ( lx    % 2 == 0 ) ? lx / 2 : ( lx - 1 ) / 2;
        lim2b = min( lim2b, bs );
   
        val = 0.0;
        for(int a=max(bs,0); a<=lima; a++) {
            twoa = 2 * a;
            faci = binom[l][a] * belt[a] * fact[twol-twoa] * binom[a][bs] / fact[l-am-twoa];
            for(int b=lim1b; b<=lim2b; b++) {
                twob = 2 * b + am - lx + 1;
                if( twob % 2 == 1 ) continue;
                val += binom[bs][b] * binom[am][twob-1]* belt[twob/2] * faci;
            };
        };
        return belt[am] * val * fac / sqrt( 2.0 * M_PI );
    };
    return 0;
    /* */
};

void BoysFn( BoysFunction <cdouble> & v ) {
    /* complex Boys function, based mostly on: 
     * K. Ishida, J. Comp. Chem. 25, 739 (2004) */
    double rez,imz,abz,ab2,del;
    rez = v.z.real();
    imz = v.z.imag();
    abz = abs( v.z );
    ab2 = abz * abz;
    del = numeric_limits<double>::epsilon();
    
    if( rez >= 0.0 ) {
        if( imz >= 0.0 ) { // quadrant x >= 0. && y >= 0.
            if( rez + imz <= 36.0 && rez <= 20.0 ) { // Luke's rational appr.
                int nfmx,nfmx1;
                nfmx  = 47;
                nfmx1 = nfmx - 1;
                
                double* AnR = new double[nfmx];
                double* AnI = new double[nfmx];
                
                double* BnR = new double[nfmx];
                double* BnI = new double[nfmx];
                
                double x2,y2,xy,den,x2y2,x3y,y3x;
                double t1,t2,t3,t4,twom;
                double F1,F2,F3;
                
                x2 = rez * rez;
                y2 = imz * imz;
                xy = 2.0 * rez * imz;
                
                x3y  = x2 - 3.0 * y2;
                y3x  = 3.0 * x2 - y2;
                x2y2 = x2 - y2;
                
                BnR[0] = 1.0;
                BnI[0] = 0.0;

                BnR[1] = 1.0 + 0.5 * rez;
                BnI[1] = 0.5 * imz;
                
                BnR[2] = BnR[1] + ( x2 - y2 ) / 12.0;
                BnI[2] = BnI[1] + xy / 12.0;
                
                for(int n=3; n<nfmx; n++) {
                    den = 4.0 * ( n + n - 1 ) * ( n + n - 3 );
                    BnR[n] = BnR[n-1] + ( x2y2 * BnR[n-2] - xy * BnI[n-2] ) / den;
                    BnI[n] = BnI[n-1] + ( x2y2 * BnI[n-2] + xy * BnR[n-2] ) / den;
                };
                
                for(int m=0; m<=v.nmax; m++) {
                    AnR[0] = 1.0;
                    AnI[0] = 0.0;
                    twom   = m + m + 1;
                    
                    t1 = twom / ( twom + 2 );
                    t2 = imz * t1;
                    t1 = rez * t1;
                    
                    AnR[1] = BnR[1] - t1;
                    AnI[1] = BnI[1] - t2;
                    
                    t3 = twom / ( twom + 2 ) / ( twom + 4 );
                    t4 = xy * t3;
                    t3 = x2y2 * t3;
                    
                    AnR[2] = BnR[2] - t1 - t3;
                    AnI[2] = BnI[2] - t2 - t4;
                    
                    for(int n=3; n<nfmx; n++) {
                        F1 = 0.5 * ( 2*n - 2*m - 5 ) / ( ( 2*n - 3 ) * ( 2*n + 2*m + 1 ) );
                        F2 = 0.25 / ( 2*n - 1 ) / ( 2*n - 3 );
                        F3 = -0.25 * F1 / ( 2*n - 3 ) / ( 2*n - 5 );
                        
                        AnR[n] = ( 1.0 + rez * F1 ) * AnR[n-1] - F1 * imz * AnI[n-1]
                               + ( x2y2 * F2 - F1 * rez ) * AnR[n-2] - ( xy * F2 - F1 * imz ) * AnI[n-2]
                               + F3 * rez * x3y * AnR[n-3] - F3 * imz * y3x * AnI[n-3];
                        
                        AnI[n] = ( 1.0 + rez * F1 ) * AnI[n-1] + F1 * imz * AnR[n-1]
                               + ( x2y2 * F2 - F1 * rez ) * AnI[n-2] + ( xy * F2 - F1 * imz ) * AnR[n-2]
                               + F3 * rez * x3y * AnI[n-3] + F3 * imz * y3x * AnR[n-3];
                    };
                    
                    den = twom * ( BnR[nfmx1] * BnR[nfmx1] + BnI[nfmx1] * BnI[nfmx1] );
                    v.fn[m] = ( AnR[nfmx1] * BnR[nfmx1] + AnI[nfmx1] * BnI[nfmx1] ) / den
                        + I * ( AnI[nfmx1] * BnR[nfmx1] - AnR[nfmx1] * BnI[nfmx1] ) / den;
                };
                
                delete []BnI; delete []BnR;
                delete []AnI; delete []AnR;
                
                return ;
            } else { // large-|z| asymptotic formula
                int n1,nf;
                double c,s,d,phir,phii,ccut,scut;
                double ex,cosy,siny,twon;
                
                ex   = exp( -rez );
                cosy = cos( +imz );
                siny = sin( +imz );
                
                phir = 0.5 * ( imz * ex * siny - rez * ex * cosy ) / ab2;
                phii = 0.5 * ( imz * ex * cosy + rez * ex * siny ) / ab2;
                
                d = 8.0 * ab2 / M_PI;
                c = ( abz + rez ) / d;
                s = ( abz - rez ) / d;
                c = sqrt( c ); s = sqrt( s );
                
                ccut = del * c;
                scut = del * s;
                
                int nfmx = 50;
                double* anR = new double[nfmx];
                double* anI = new double[nfmx];
                
                anR[0] = phir;
                anI[0] = phii;
                
                nf = nfmx - 1;
                for(int n=1; n<nfmx; n++) {
                    n1 = n - 1;
                    twon = 0.5 * ( n + n - 1 );
                    anR[n] = -twon * ( rez * anR[n1] + imz * anI[n1] ) / ab2;
                    anI[n] = +twon * ( imz * anR[n1] - rez * anI[n1] ) / ab2;
                    if( abs( anR[n] ) < ccut && abs( anI[n] ) < scut ) {
                        nf = n;
                        break ;
                    };
                };
                
                double* fnR = new double[v.nmax+1];
                double* fnI = new double[v.nmax+1];
                
                fnR[0] = +c;
                fnI[0] = -s;
                for(int n=0; n<=nf; n++) {
                    fnR[0] += anR[n];
                    fnI[0] += anI[n];
                };
                
                for(int n=1; n<=v.nmax; n++) {
                    n1 = n - 1;
                    twon = 0.5 * ( n + n - 1 );
                    fnR[n] = +twon * ( rez * fnR[n1] + imz * fnI[n1] ) / ab2;
                    fnI[n] = +twon * ( rez * fnI[n1] - imz * fnR[n1] ) / ab2;
                };
                
                double* bnR = new double[v.nmax+1];
                double* bnI = new double[v.nmax+1];
                
                bnR[0] = 0.0;
                bnI[0] = 0.0;
                
                if( v.nmax > 0 ) {
                    bnR[1] = phir;
                    bnI[1] = phii;
                };
                
                for(int n=2; n<=v.nmax; n++) {
                    n1 = n - 1;
                    twon = 0.5 * ( n + n - 1 );
                    bnR[n] = phir + twon * ( rez * bnR[n1] + imz * bnI[n1] ) / ab2;
                    bnI[n] = phii + twon * ( rez * bnI[n1] - imz * bnR[n1] ) / ab2;
                };
                
                for(int n=0; n<=v.nmax; n++)
                    v.fn[n] = fnR[n] + bnR[n] + I * ( fnI[n] + bnI[n] );
                
                delete []bnI; delete []bnR;
                delete []fnI; delete []fnR;
                delete []anI; delete []anR;
                return ;
            };
        } else { // quadrant x >= 0. && y >= 0., reflection
            v.z = conj( v.z ); BoysFn( v );
            for(int n=0; n<=v.nmax; n++) v.fn[n] = conj( v.fn[n] );
            v.z = conj( v.z );
            return ;
        };
    } else {
        if( imz <= 0.0 ) { // quadrant x <= 0. && y <= 0.
            if( abs( rez ) + abs( imz ) <= 51.0 && 
                abs( rez ) <= 36.0 && abs( imz ) <= 36.0 ) { // Luke's rational appr.
                rez = -rez;
                imz = -imz;
                
                int nfmx,nfmx1;
                nfmx  = 39;
                nfmx1 = nfmx - 1;
                
                double* AnR = new double[nfmx];
                double* AnI = new double[nfmx];
                
                double* BnR = new double[nfmx];
                double* BnI = new double[nfmx];
                
                double x2,y2,xy,x2y2;
                double t1,twom,F1,F2,den;
                int n4,tmp;
                cdouble et;
                
                x2 = rez * rez;
                y2 = imz * imz;
                xy = 2.0 * rez * imz;
                
                x2y2 = x2 - y2;
                
                for(int m=0; m<=v.nmax; m++) {
                    twom   = m + m + 1;
                    
                    BnR[0] = 1.0;
                    BnI[0] = 0.0;
                    
                    t1 = 2.0 / ( twom + 4 );
                
                    BnR[1] = 1.0 + t1 * rez;
                    BnI[1] =       t1 * imz;
                    
                    AnR[0] = 1.0;
                    AnI[0] = 0.0;
                    
                    t1 = 2.0 / ( twom + 2 );
                
                    AnR[1] = BnR[1] - t1 * rez;
                    AnI[1] = BnI[1] - t1 * imz;
                
                    for(int n=2; n<nfmx; n++) {
                        n4 = n  + n ;
                        n4 = n4 + n4;
                        tmp= n4 + twom - 4;
                        
                        F1 = 2.0 * twom / ( n4 + twom ) / tmp;
                        F2 = 8.0 * ( n - 1 ) * ( n + n + twom - 2 ) 
                           / ( tmp + 2 ) / tmp / tmp / ( tmp - 2 );
                        
                        AnR[n] = ( 1.0 + F1 * rez ) * AnR[n-1] - F1 * imz * AnI[n-1]
                               + F2 * ( x2y2 * AnR[n-2] - xy * AnI[n-2] );
                        AnI[n] = ( 1.0 + F1 * rez ) * AnI[n-1] + F1 * imz * AnR[n-1]
                               + F2 * ( x2y2 * AnI[n-2] + xy * AnR[n-2] );
                        
                        BnR[n] = ( 1.0 + F1 * rez ) * BnR[n-1] - F1 * imz * BnI[n-1]
                               + F2 * ( x2y2 * BnR[n-2] - xy * BnI[n-2] );
                        BnI[n] = ( 1.0 + F1 * rez ) * BnI[n-1] + F1 * imz * BnR[n-1]
                               + F2 * ( x2y2 * BnI[n-2] + xy * BnR[n-2] );
                    };
                    
                    den = twom * ( BnR[nfmx1] * BnR[nfmx1] + BnI[nfmx1] * BnI[nfmx1] );
                    v.fn[m] = ( AnR[nfmx1] * BnR[nfmx1] + AnI[nfmx1] * BnI[nfmx1] ) / den
                        + I * ( AnI[nfmx1] * BnR[nfmx1] - AnR[nfmx1] * BnI[nfmx1] ) / den;
                };
                
                delete []BnI; delete []BnR;
                delete []AnI; delete []AnR;
                
                /*rez = -rez;
                imz = -imz;
                
                et = exp(-v.z);
                for(int m=0; m<=v.nmax; m++) v.fn[m] = v.fn[m] * et;*/
                
                return ;
            } else { // large-|z| asymptotic formula
                rez = -rez;
                imz = -imz;
                
                int n1,nf;
                double c,s,d,SnR,SnI;
                double ex,cosy,siny,twon;
                cdouble asym,et;
                
                ex   = exp( -rez );
                cosy = cos( +imz );
                siny = sin( +imz );
                
                d = 8.0 * ab2 / M_PI;
                c = ( abz + rez ) / d;
                s = ( abz - rez ) / d;
                c = sqrt( c ); s = sqrt( s );
                
                SnR  = c * ex * siny + s * ex * cosy;
                SnI  = c * ex * cosy - s * ex * siny;
                
                asym = SnR + I * SnI;
                
                int nfmx = 50;
                double* anR = new double[nfmx];
                double* anI = new double[nfmx];
                
                anR[0] = +0.5 * rez / ab2;
                anI[0] = -0.5 * imz / ab2;
                
                SnR += anR[0];
                SnI += anI[0];
                
                nf = nfmx - 1;
                for(int n=1; n<nfmx; n++) {
                    n1 = n - 1;
                    twon = 0.5 * ( n + n - 1 );
                    
                    anR[n] = +twon * ( rez * anR[n1] + imz * anI[n1] ) / ab2;
                    anI[n] = +twon * ( rez * anI[n1] - imz * anR[n1] ) / ab2;
                    
                    SnR += anR[n];
                    SnI += anI[n];
                    
                    if( abs( anR[n] ) < del * abs( SnR ) && 
                        abs( anI[n] ) < del * abs( SnI ) ) {
                        nf = n;
                        break ;
                    };
                };
                
                double* gnR = new double[v.nmax+1];
                double* gnI = new double[v.nmax+1];
                
                gnR[0] = +SnR;
                gnI[0] = +SnI;
                
                for(int n=1; n<=v.nmax; n++) {
                    n1 = n - 1;
                    twon = 0.5 * ( n + n - 1 );
                    gnR[n] = -twon * ( rez * gnR[n1] + imz * gnI[n1] ) / ab2;
                    gnI[n] = +twon * ( imz * gnR[n1] - rez * gnI[n1] ) / ab2;
                };
                
                double* bnR = new double[v.nmax+1];
                double* bnI = new double[v.nmax+1];
                
                bnR[0] = 0.0;
                bnI[0] = 0.0;
                
                if( v.nmax > 0 ) {
                    bnR[1] = anR[0];
                    bnI[1] = anI[0];
                };
                
                for(int n=2; n<=v.nmax; n++) {
                    n1 = n - 1;
                    twon = 0.5 * ( n + n - 1 );
                    bnR[n] = anR[0] - twon * ( rez * bnR[n1] + imz * bnI[n1] ) / ab2;
                    bnI[n] = anI[0] + twon * ( imz * bnR[n1] - rez * bnI[n1] ) / ab2;
                };
                
                for(int n=0; n<=v.nmax; n++)
                    v.fn[n] = gnR[n] + bnR[n] + I * ( gnI[n] + bnI[n] );
                
                delete []bnI; delete []bnR;
                delete []gnI; delete []gnR;
                delete []anI; delete []anR;
                
                /*rez = -rez;
                imz = -imz;
                
                et = exp(-v.z);
                for(int m=0; m<=v.nmax; m++) v.fn[m] = v.fn[m] * et;*/
                
                return ;
            };
        } else { // quadrant x <= 0. && y >= 0., reflection
            v.z = conj( v.z ); BoysFn( v );
            for(int n=0; n<=v.nmax; n++) v.fn[n] = conj( v.fn[n] );
            v.z = conj( v.z );
            return ;
        };
    };
    /* */
};
