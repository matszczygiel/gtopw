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
#include "../headers/auxfun2.h"

using namespace std;

void CalcRERI( RInts2E <cdouble> & R ) {
    /* R_{tuv,t'u'v'} for two-electron ERI */
    cdouble kP2,kQ2,eP1,eP2,eQ1,eQ2;
    cdouble Ppx,Ppy,Ppz,Qpx,Qpy,Qpz,fac0;
    cdouble PQx,PQy,PQz,bigT,xi,twoxi;
    cdouble ikPx,ikPy,ikPz,ikQx,ikQy,ikQz;
    
    Ppx = R.Px + I * R.kPx / ( 2.0 * R.aP );
    Ppy = R.Py + I * R.kPy / ( 2.0 * R.aP );
    Ppz = R.Pz + I * R.kPz / ( 2.0 * R.aP );
    
    Qpx = R.Qx + I * R.kQx / ( 2.0 * R.aQ );
    Qpy = R.Qy + I * R.kQy / ( 2.0 * R.aQ );
    Qpz = R.Qz + I * R.kQz / ( 2.0 * R.aQ );
    
    PQx = Ppx - Qpx;
    PQy = Ppy - Qpy;
    PQz = Ppz - Qpz;
    
    xi   = R.aP * R.aQ / ( R.aP + R.aQ );
    bigT = PQx * PQx + PQy * PQy + PQz * PQz;
    bigT = bigT * xi;
    twoxi= -2.0 * xi;
    
    if( bigT.real() > 0.0 ) {
        kP2 = R.kPx * R.kPx + R.kPy * R.kPy + R.kPz * R.kPz;
        eP1 = exp( -kP2 / ( 4.0 * R.aP ) );
        eP2 = exp( I * ( R.kPx * R.Px + R.kPy * R.Py + R.kPz * R.Pz ) );
    
        kQ2 = R.kQx * R.kQx + R.kQy * R.kQy + R.kQz * R.kQz;
        eQ1 = exp( -kQ2 / ( 4.0 * R.aQ ) );
        eQ2 = exp( I * ( R.kQx * R.Qx + R.kQy * R.Qy + R.kQz * R.Qz ) );
        
        fac0= 2.0 * pi52 * eP1 * eP2 * eQ1 * eQ2 / ( R.aP * R.aQ * sqrt( R.aP + R.aQ ) );
    } else {
        kP2 = R.kPx * R.kPx + R.kPy * R.kPy + R.kPz * R.kPz;
        eP1 = exp( -kP2 / ( 4.0 * R.aP ) );
        eP2 = exp( I * ( R.kPx * R.Px + R.kPy * R.Py + R.kPz * R.Pz ) );
    
        kQ2 = R.kQx * R.kQx + R.kQy * R.kQy + R.kQz * R.kQz;
        eQ1 = exp( -kQ2 / ( 4.0 * R.aQ ) );
        eQ2 = exp( I * ( R.kQx * R.Qx + R.kQy * R.Qy + R.kQz * R.Qz ) );
        
        /*eP1 = ( R.kPx - R.kQx ) * ( R.kPx - R.kQx )
            + ( R.kPy - R.kQy ) * ( R.kPy - R.kQy )
            + ( R.kPz - R.kQz ) * ( R.kPz - R.kQz );
        eP1 = eP1 / 4.0 / ( R.aP + R.aQ );
        eP1 = exp( -eP1 );
        
        eP2 = ( R.Px - R.Qx ) * ( R.Px - R.Qx )
            + ( R.Py - R.Qy ) * ( R.Py - R.Qy )
            + ( R.Pz - R.Qz ) * ( R.Pz - R.Qz );
        eP2 = xi * eP2;
        eP2 = exp( -eP2 );
        
        eQ1 = R.aQ * ( R.kPx + R.kQx ) * R.Qx
            + R.aQ * ( R.kPy + R.kQy ) * R.Qy
            + R.aQ * ( R.kPz + R.kQz ) * R.Qz;
        eQ1 = eQ1 / ( R.aP + R.aQ );
        eQ1 = exp( I * eQ1 );
        
        eQ2 = R.aP * ( R.kPx + R.kQx ) * R.Px
            + R.aP * ( R.kPy + R.kQy ) * R.Py
            + R.aP * ( R.kPz + R.kQz ) * R.Pz;
        eQ2 = eQ2 / ( R.aP + R.aQ );
        eQ2 = exp( I * eQ2 );
        
        cout << "eP1 " << eP1 << endl;
        cout << "eP2 " << eP2 << endl;
        cout << "eQ1 " << eQ1 << endl;
        cout << "eQ2 " << eQ2 << endl;*/
        
        fac0= 2.0 * pi52 * eP1 * eP2 * eQ1 * eQ2 / ( R.aP * R.aQ * sqrt( R.aP + R.aQ ) ) * exp( -bigT );
    };
    
    ikPx = I * R.kPx; ikPy = I * R.kPy; ikPz = I * R.kPz;
    ikQx = I * R.kQx; ikQy = I * R.kQy; ikQz = I * R.kQz;
    
    BoysFunction <cdouble> boys( R.lsum_PQ, bigT );
    BoysFn( boys );
    
    /* R_{000,000,n} */
    register int t1,u1,v1,t2,u2,v2,j;
    
    for(j=0; j<=R.lsum_PQ; j++) {
        R.R2En[0][0][0][0][0][0][j] = fac0 * boys.fn[j];
        fac0 = fac0 * twoxi;
    };
    
    if( R.lsum_PQ == 0 ) {
        R.strip();
        return ;
    };
    
    /* R_{000,00v',n} */
    for(j=0; j<R.lsum_PQ; j++) 
        R.R2En[0][0][0][0][0][1][j] = ikQz * R.R2En[0][0][0][0][0][0][j  ]
                                    -  PQz * R.R2En[0][0][0][0][0][0][j+1];
    
    for(v2=1; v2<R.lsum_PQ; v2++) for(j=0; j<R.lsum_PQ-v2; j++) 
         R.R2En[0][0][0][0][0][v2+1][j] = ikQz * R.R2En[0][0][0][0][0][v2][j  ]
                                        -  PQz * R.R2En[0][0][0][0][0][v2][j+1]
                                        + (double)v2 * R.R2En[0][0][0][0][0][v2-1][j+1];
    
    /* R_{000,0u'v',n} */
    for(v2=0; v2<R.lsum_PQ; v2++) for(j=0; j<R.lsum_PQ-v2; j++) 
        R.R2En[0][0][0][0][1][v2][j] = ikQy * R.R2En[0][0][0][0][0][v2][j  ]
                                     -  PQy * R.R2En[0][0][0][0][0][v2][j+1];
    
    for(u2=1; u2<R.lsum_PQ; u2++) for(v2=0; v2<R.lsum_PQ-u2; v2++)
        for(j=0; j<R.lsum_PQ-v2-u2; j++) 
        R.R2En[0][0][0][0][u2+1][v2][j] = ikQy * R.R2En[0][0][0][0][u2][v2][j  ]
                                        -  PQy * R.R2En[0][0][0][0][u2][v2][j+1]
                                        + (double)u2 * R.R2En[0][0][0][0][u2-1][v2][j+1];
    /* R_{000,t'u'v',n} */
    for(u2=0; u2<R.lsum_PQ; u2++) for(v2=0; v2<R.lsum_PQ-u2; v2++) 
        for(j=0; j<R.lsum_PQ-v2-u2; j++) 
        R.R2En[0][0][0][1][u2][v2][j] = ikQx * R.R2En[0][0][0][0][u2][v2][j  ]
                                      -  PQx * R.R2En[0][0][0][0][u2][v2][j+1];
    
    for(t2=1; t2<R.lsum_PQ; t2++) for(u2=0; u2<R.lsum_PQ-t2; u2++)
        for(v2=0; v2<R.lsum_PQ-t2-u2; v2++) for(j=0; j<R.lsum_PQ-t2-v2-u2; j++) 
        R.R2En[0][0][0][t2+1][u2][v2][j] = ikQx * R.R2En[0][0][0][t2][u2][v2][j  ]
                                         -  PQx * R.R2En[0][0][0][t2][u2][v2][j+1]
                                         + (double)t2 * R.R2En[0][0][0][t2-1][u2][v2][j+1];
    /* R_{00v,t'u'v',n} */
    for(t2=0; t2<R.lsum_PQ; t2++) for(u2=0; u2<R.lsum_PQ-t2; u2++) 
        for(j=0; j<R.lsum_PQ-u2-t2; j++) 
        R.R2En[0][0][1][t2][u2][0][j] = ikPz * R.R2En[0][0][0][t2][u2][0][j  ]
                                      +  PQz * R.R2En[0][0][0][t2][u2][0][j+1];
    
    for(t2=0; t2<R.lsum_PQ; t2++) for(u2=0; u2<R.lsum_PQ-t2; u2++) 
        for(v2=1; v2<R.lsum_PQ-u2-t2; v2++) for(j=0; j<R.lsum_PQ-v2-u2-t2; j++) 
        R.R2En[0][0][1][t2][u2][v2][j] = ikPz * R.R2En[0][0][0][t2][u2][v2][j  ]
                                       +  PQz * R.R2En[0][0][0][t2][u2][v2][j+1]
                                       - (double)v2 * R.R2En[0][0][0][t2][u2][v2-1][j+1];
    
    for(v1=1; v1<R.lsum_PQ; v1++) for(t2=0; t2<R.lsum_PQ-v1; t2++)
        for(u2=0; u2<R.lsum_PQ-t2-v1; u2++) for(j=0; j<R.lsum_PQ-u2-t2-v1; j++) 
            R.R2En[0][0][v1+1][t2][u2][0][j] = ikPz * R.R2En[0][0][v1][t2][u2][0][j  ]
                                       +  PQz * R.R2En[0][0][v1][t2][u2][0][j+1]
                                       + (double)v1 * R.R2En[0][0][v1-1][t2][u2][0][j+1];
    
    for(v1=1; v1<R.lsum_PQ; v1++) for(t2=0; t2<R.lsum_PQ-v1; t2++)
        for(u2=0; u2<R.lsum_PQ-t2-v1; u2++) for(v2=1; v2<R.lsum_PQ-u2-t2-v1; v2++)
            for(j=0; j<R.lsum_PQ-v2-u2-t2-v1; j++) 
            R.R2En[0][0][v1+1][t2][u2][v2][j] = ikPz * R.R2En[0][0][v1][t2][u2][v2][j  ]
                                       +  PQz * R.R2En[0][0][v1][t2][u2][v2][j+1]
                                       - (double)v2 * R.R2En[0][0][v1][t2][u2][v2-1][j+1]
                                       + (double)v1 * R.R2En[0][0][v1-1][t2][u2][v2][j+1];
    /* R_{0uv,t'u'v',n} */
    for(v1=0; v1<R.lsum_PQ; v1++) for(t2=0; t2<R.lsum_PQ-v1; t2++)
        for(v2=0; v2<R.lsum_PQ-t2-v1; v2++) for(j=0; j<R.lsum_PQ-v1-t2-v2; j++) 
        R.R2En[0][1][v1][t2][0][v2][j] = ikPy * R.R2En[0][0][v1][t2][0][v2][j  ]
                                       +  PQy * R.R2En[0][0][v1][t2][0][v2][j+1];
    
    for(v1=0; v1<R.lsum_PQ; v1++) for(t2=0; t2<R.lsum_PQ-v1; t2++)
        for(u2=1; u2<R.lsum_PQ-t2-v1; u2++) for(v2=0; v2<R.lsum_PQ-t2-v1-u2; v2++) 
            for(j=0; j<R.lsum_PQ-v1-t2-v2-u2; j++) 
            R.R2En[0][1][v1][t2][u2][v2][j] = ikPy * R.R2En[0][0][v1][t2][u2][v2][j  ]
                                            +  PQy * R.R2En[0][0][v1][t2][u2][v2][j+1]
                                            - (double)u2 * R.R2En[0][0][v1][t2][u2-1][v2][j+1];
    
    for(u1=1; u1<R.lsum_PQ; u1++) for(v1=0; v1<R.lsum_PQ-u1; v1++) 
        for(t2=0; t2<R.lsum_PQ-v1-u1; t2++) for(v2=0; v2<R.lsum_PQ-t2-v1-u1; v2++)
            for(j=0; j<R.lsum_PQ-v1-t2-v2-u1; j++) 
        R.R2En[0][u1+1][v1][t2][0][v2][j] = ikPy * R.R2En[0][u1][v1][t2][0][v2][j  ]
                                          +  PQy * R.R2En[0][u1][v1][t2][0][v2][j+1]
                                          + (double)u1 * R.R2En[0][u1-1][v1][t2][0][v2][j+1];
    
    for(u1=1; u1<R.lsum_PQ; u1++) for(v1=0; v1<R.lsum_PQ-u1; v1++) 
        for(t2=0; t2<R.lsum_PQ-v1-u1; t2++) for(u2=1; u2<R.lsum_PQ-v1-u1; u2++) 
            for(v2=0; v2<R.lsum_PQ-t2-v1-u1-u2; v2++) for(j=0; j<R.lsum_PQ-v1-t2-v2-u1-u2; j++) 
        R.R2En[0][u1+1][v1][t2][u2][v2][j] = ikPy * R.R2En[0][u1][v1][t2][u2][v2][j  ]
                                           +  PQy * R.R2En[0][u1][v1][t2][u2][v2][j+1]
                                           - (double)u2 * R.R2En[0][u1][v1][t2][u2-1][v2][j+1]
                                           + (double)u1 * R.R2En[0][u1-1][v1][t2][u2][v2][j+1];
    /* R_{tuv,t'u'v',n} */
    for(u1=0; u1<R.lsum_PQ; u1++) for(v1=0; v1<R.lsum_PQ-u1; v1++)
        for(u2=0; u2<R.lsum_PQ-u1-v1; u2++) for(v2=0; v2<R.lsum_PQ-u1-v1-u2; v2++) 
            for(j=0; j<R.lsum_PQ-u1-v1-u2-v2; j++) 
            R.R2En[1][u1][v1][0][u2][v2][j] = ikPx * R.R2En[0][u1][v1][0][u2][v2][j  ]
                                            +  PQx * R.R2En[0][u1][v1][0][u2][v2][j+1];
    
    for(u1=0; u1<R.lsum_PQ; u1++) for(v1=0; v1<R.lsum_PQ-u1; v1++)
        for(t2=1; t2<R.lsum_PQ-v1-u1; t2++) for(u2=0; u2<R.lsum_PQ-u1-v1-t2; u2++) 
            for(v2=0; v2<R.lsum_PQ-u1-v1-t2-u2; v2++) for(j=0; j<R.lsum_PQ-u1-v1-t2-u2-v2; j++) 
            R.R2En[1][u1][v1][t2][u2][v2][j] = ikPx * R.R2En[0][u1][v1][t2][u2][v2][j  ]
                                             +  PQx * R.R2En[0][u1][v1][t2][u2][v2][j+1]
                                             - (double)t2 * R.R2En[0][u1][v1][t2-1][u2][v2][j+1];
    
    for(t1=1; t1<R.lsum_PQ; t1++) for(u1=0; u1<R.lsum_PQ-t1; u1++) 
        for(v1=0; v1<R.lsum_PQ-u1-t1; v1++) for(u2=0; u2<R.lsum_PQ-u1-v1-t1; u2++)
            for(v2=0; v2<R.lsum_PQ-u1-v1-u2-t1; v2++) for(j=0; j<R.lsum_PQ-u1-v1-u2-v2-t1; j++) 
                R.R2En[t1+1][u1][v1][0][u2][v2][j] = ikPx * R.R2En[t1][u1][v1][0][u2][v2][j  ]
                                                   +  PQx * R.R2En[t1][u1][v1][0][u2][v2][j+1]
                                                   + (double)t1 * R.R2En[t1-1][u1][v1][0][u2][v2][j+1];
    
    for(t1=1; t1<R.lsum_PQ; t1++) for(u1=0; u1<R.lsum_PQ-t1; u1++) 
        for(v1=0; v1<R.lsum_PQ-u1-t1; v1++) for(t2=1; t2<R.lsum_PQ-t1-u1-v1; t2++)
         for(u2=0; u2<R.lsum_PQ-u1-v1-t1-t2; u2++) for(v2=0; v2<R.lsum_PQ-u1-v1-t1-t2-u2; v2++) 
            for(j=0; j<R.lsum_PQ-u1-v1-t1-t2-u2-v2; j++) 
                R.R2En[t1+1][u1][v1][t2][u2][v2][j] = ikPx * R.R2En[t1][u1][v1][t2][u2][v2][j  ]
                                                    +  PQx * R.R2En[t1][u1][v1][t2][u2][v2][j+1]
                                                    + (double)t1 * R.R2En[t1-1][u1][v1][t2][u2][v2][j+1]
                                                    - (double)t2 * R.R2En[t1][u1][v1][t2-1][u2][v2][j+1];
    R.strip();
    /* */
};
