#ifndef AUXFUN1_H
#define AUXFUN1_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <complex>

#include "../headers/gtopw.h"
#include "../headers/auxfun2.h"

void WriteDown( cdouble* data, const int len, std::ofstream & ofs );

void TransToSpher1D( cdouble* shl_crt, 
                     cdouble* shl_sph, 
                     const int li );

void TransToSpher   ( cdouble* shl_crt, cdouble* shl_sph    , const int li, const int lj );
void GamessHolyOrder( cdouble* shl_crt, cdouble* shl_crt_dum, const int li, const int lj );

void GamessHolyOrder2E( FullC <cdouble> & full_trans, 
                        FullC <cdouble> & full_trans_gam,
                        const int li, const int lj,
                        const int lk, const int ll,
                        const bool overwrite );

usint FindPos( usint ia_, usint ja_, usint ka_ );

inline int pos_change_gamess( const int L, const int pos ) {
    /* change ordering within the shell to the Gamess convention */
    switch( L ) {
        case 0: // s-type
        return 0;
        case 1: // p-type
        switch( pos ) {
            case 0: return 2;
            case 1: return 1;
            case 2: return 0;
            default: 
            std::cout << " Error in pos_change!" << std::endl;
            return EXIT_FAILURE;
        }
        case 2: // d-type
        switch( pos ) {
            case 0: return 2;
            case 1: return 5;
            case 2: return 1;
            case 3: return 4;
            case 4: return 3;
            case 5: return 0;
            default:
            std::cout << " Error in pos_change!" << std::endl;
            return EXIT_FAILURE;
        }
        case 3: // f-type
        switch( pos ) {
            case 0: return 2;
            case 1: return 8;
            case 2: return 6;
            case 3: return 1;
            case 4: return 7;
            case 5: return 9;
            case 6: return 5;
            case 7: return 4;
            case 8: return 3;
            case 9: return 0;
            default:
            std::cout << " Error in pos_change!" << std::endl;
            return EXIT_FAILURE;
        }
        case 4: // g-type
        switch( pos ) {
            case 0 : return 2 ;
            case 1 : return 8 ;
            case 2 : return 11;
            case 3 : return 6 ;
            case 4 : return 1 ;
            case 5 : return 7 ;
            case 6 : return 14;
            case 7 : return 13;
            case 8 : return 5 ;
            case 9 : return 10;
            case 10: return 12;
            case 11: return 9 ;
            case 12: return 4 ;
            case 13: return 3 ;
            case 14: return 0 ;
            default:
            std::cout << " Error in pos_change!" << std::endl;
            return EXIT_FAILURE;
        }
        case 5: // h-type
        switch( pos ) {
            case 0 : return 2 ;
            case 1 : return 8 ;
            case 2 : return 14;
            case 3 : return 12;
            case 4 : return 6 ;
            case 5 : return 1 ;
            case 6 : return 7 ;
            case 7 : return 17;
            case 8 : return 20;
            case 9 : return 16;
            case 10: return 5 ;
            case 11: return 13;
            case 12: return 19;
            case 13: return 18;
            case 14: return 11;
            case 15: return 10;
            case 16: return 15;
            case 17: return 9 ;
            case 18: return 4 ;
            case 19: return 3 ;
            case 20: return 0 ;
            default:
            std::cout << " Error in pos_change!" << std::endl;
            return EXIT_FAILURE;
        }
        case 6: // i-type
        switch( pos ) {
            case 0 : return 2 ;
            case 1 : return 8 ;
            case 2 : return 14;
            case 3 : return 20;
            case 4 : return 12;
            case 5 : return 6 ;
            case 6 : return 1 ;
            case 7 : return 7 ;
            case 8 : return 17;
            case 9 : return 26;
            case 10: return 24;
            case 11: return 16;
            case 12: return 5;
            case 13: return 13;
            case 14: return 25;
            case 15: return 27;
            case 16: return 23;
            case 17: return 11;
            case 18: return 19;
            case 19: return 22;
            case 20: return 21;
            case 21: return 18;
            case 22: return 10;
            case 23: return 15;
            case 24: return 9 ;
            case 25: return 4 ;
            case 26: return 3 ;
            case 27: return 0 ;
            default:
            std::cout << " Error in pos_change!" << std::endl;
            return EXIT_FAILURE;
        }
        default: // higher L not supported in Gamess
        return pos;
    };
};

template <class type> 
    class BoysFunction {
    public:
    
    static const int ncap = 32;
    int nmax;
    type  z;
    type* fn;
    
    BoysFunction( const int nmax_, const type z_ ) : nmax(nmax_), z(z_) { 
        if( nmax > ncap ) {
            std::cout<< " " << "Boys function is implemented only up to n = 32!"<<std::endl;
            std::cout<< " " << "Emergency halt."<<std::endl;
            std::exit( EXIT_FAILURE );
        };
        fn = new type[ncap+1];
    };
    
    void print() {
        for(int n=0; n<=nmax; n++) 
            std::cout<< n << " " << fn[n] << std::endl;
    };
    
   ~BoysFunction() { delete []fn; };
};

void BoysFn( BoysFunction <cdouble> & v );

template <class type> class RInts1E {
    public:
    int lsum;
    
    type kPx;
    type kPy;
    type kPz;
    
    type Px;
    type Py;
    type Pz;
    
    type Cx;
    type Cy;
    type Cz;
    
    type aP;
    
    type***  rtuv;
    type**** rtuvj;
    
    void strip() {
        for(int t=0; t<=lsum; t++) for(int u=0; u<=lsum; u++) for(int v=0; v<=lsum; v++) 
            rtuv[t][u][v] = rtuvj[t][u][v][0];
    };
    
    
    void load( const type kPx_, const type kPy_, const type kPz_, 
               const type Px_ , const type Py_ , const type Pz_ ,
               const type aP_ ) {
        kPx = kPx_; kPy = kPy_; kPz = kPz_;
        Px  = Px_ ; Py  = Py_ ; Pz  = Pz_ ;
        aP  = aP_ ;
    };
    
    void load( const type kPx_, const type kPy_, const type kPz_, 
               const type Px_ , const type Py_ , const type Pz_ ,
               const type Cx_ , const type Cy_ , const type Cz_ ,
               const type aP_ ) {
        kPx = kPx_; kPy = kPy_; kPz = kPz_;
        Px  = Px_ ; Py  = Py_ ; Pz  = Pz_ ;
        Cx  = Cx_ ; Cy  = Cy_ ; Cz  = Cz_ ;
        aP  = aP_ ;
    };
    
    void zero() {
        for(int t=0; t<=lsum; t++) for(int u=0; u<=lsum; u++) 
            for(int v=0; v<=lsum; v++) 
                memset(rtuvj[t][u][v],0,(lsum+1)*sizeof(type));
        for(int t=0; t<=lsum; t++) for(int u=0; u<=lsum; u++) 
                memset(rtuv[t][u],0,(lsum+1)*sizeof(type));
    };
    
    void print() {
        for(int t=0; t<=lsum; t++) for(int u=0; u<=lsum; u++) for(int v=0; v<=lsum-t-u; v++)
            std::cout<< t << " " << u << " " << v << " " << rtuv[t][u][v] << std::endl;
    };
    
    RInts1E( const int lsum_ ) : lsum(lsum_) {
        rtuv = new type**[lsum+1];
        for(int t=0; t<=lsum; t++) {
            rtuv[t] = new type*[lsum+1];
            for(int u=0; u<=lsum; u++) {
                rtuv[t][u] = new type[lsum+1];
                memset(rtuv[t][u],0,(lsum+1)*sizeof(type));
            };
        };
    
        rtuvj = new type***[lsum+1];
        for(int t=0; t<=lsum; t++) {
            rtuvj[t] = new type**[lsum+1];
            for(int u=0; u<=lsum; u++) {
                rtuvj[t][u] = new type*[lsum+1];
                for(int v=0; v<=lsum; v++) {
                    rtuvj[t][u][v] = new type[lsum+1];
                    memset(rtuvj[t][u][v],0,(lsum+1)*sizeof(type));
                };
            };
        };  
    };
    
   ~RInts1E() {
        for(int t=0; t<=lsum; t++) {
            for(int u=0; u<=lsum; u++) {
                for(int v=0; v<=lsum; v++) delete []rtuvj[t][u][v];
                delete []rtuvj[t][u];
            };
            delete []rtuvj[t];
        };
        delete []rtuvj;
    
        for(int t=0; t<=lsum; t++) {
            for(int u=0; u<=lsum; u++) delete []rtuv[t][u];
            delete []rtuv[t];
        };
        delete []rtuv;
    };
};

template <class type> 
    class ECoefs {
    public:
    
    int imax;
    int jmax;
    int umax;
    
    type a;
    type b;
    type p;
    type twop;
    
    type AB;
    
    type*** v;
    
    void load( const type a_, const type b_ ) { 
        a = a_;
        b = b_;
        p = a + b;
        twop = p + p;
    };
    
    void zero() { 
        for(int i=0; i<=imax; i++) for(int j=0; j<=jmax; j++) 
            memset(v[i][j],0,(umax+1)*sizeof(type));
    };
    
    ECoefs( const int imax_, const int jmax_, const type a_, const type b_, 
        const type AB_ ) : imax(imax_), jmax(jmax_), a(a_), b(b_), AB(AB_) {
    
        umax = imax + jmax;
        p    = a + b;
        twop = p + p;
    
        v = new type**[imax+1];
        for(int i=0; i<=imax; i++) {
            v[i] = new type*[jmax+1];
            for(int j=0; j<=jmax; j++) {
                v[i][j] = new type[umax+1];
                memset(v[i][j],0,(umax+1)*sizeof(type));
            };
        };
    };
    
   ~ECoefs() {
        for(int i=0; i<=imax; i++) {
            for(int j=0; j<=jmax; j++) delete []v[i][j];
            delete []v[i];
        };
        delete []v;
    };
    
    void print() {
        for(int i=0; i<=imax; i++) for(int j=0; j<=jmax; j++) for(int u=0; u<=umax; u++) 
            std::cout << i << " " << j << " " << u << " " << v[i][j][u] << std::endl;
    };
};

void CalcROvrl( RInts1E <cdouble> & R, const cdouble fac );
void CalcRNucA( RInts1E <cdouble> & R );

void CalcEijt ( ECoefs  <double>  & E );

//void RenormContr( GaussC& gau_tmp );
void RenormContr( GaussC& gau_tmp, const std::string norm_name );

double CalcClmR( const int l, const int m, const int lx, const int ly, const int lz );
double NoNormCalcClmR( const int l, const int m, const int lx, const int ly, const int lz );

#endif //AUXFUN1_H
