#ifndef AUXFUN2_H
#define AUXFUN2_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <complex>
#include <fstream>

#include "../headers/gtopw.h"


inline bool is_unique(const int& i, const int& j, const int& k, const int& l) {
    if(i+j-k-l >= 0 && i-j+k-l >= 0 && i-j-k+l >= 0)
        return true;
    if(i+j-k-l < 0 && i-j+k-l < 0 && i-j-k+l < 0)
        return true;

    return false;
}


template <class type> class FullC {
    public:
    int lA;
    int lB;
    int lC;
    int lD;
    
    int shgA;
    int shgB;
    int shgC;
    int shgD;
    
    int num;
    
    int posA_strt;
    int posB_strt;
    int posC_strt;
    int posD_strt;
    
    int Ash;
    int Bsh;
    int Csh;
    int Dsh;
    
    type**** v;
    
    FullC( const int lA_, const int lB_, const int lC_, const int lD_, const int num_ ) :
        lA( lA_ ), lB( lB_ ), lC( lC_ ), lD( lD_ ), num(num_) {
        
        shgA = crt_siz[lA];
        shgB = crt_siz[lB];
        shgC = crt_siz[lC];
        shgD = crt_siz[lD];
        
        v = new type***[shgA];
        for(int mi=0; mi<shgA; mi++) {
            v[mi] = new type**[shgB];
            for(int mj=0; mj<shgB; mj++) {
                v[mi][mj] = new type*[shgC];
                for(int mk=0; mk<shgC; mk++) {
                    v[mi][mj][mk] = new type[shgD];
                    memset( v[mi][mj][mk], 0, shgD * sizeof( type ) );
                };
            };
        };
        
    };
    
    ~FullC( ) {
        for(int mi=0; mi<shgA; mi++) {
            for(int mj=0; mj<shgB; mj++) {
                for(int mk=0; mk<shgC; mk++) delete []v[mi][mj][mk];
                delete []v[mi][mj];
            };
            delete []v[mi];
        };
        delete []v;
    };
    
    void zero( ) {
        for(int mi=0; mi<shgA; mi++)
            for(int mj=0; mj<shgB; mj++)
                for(int mk=0; mk<shgC; mk++)
                    memset( v[mi][mj][mk], 0, shgD * sizeof( type ) );
    };
    
    void load_pos_mx( const int pA , const int pB , const int pC , const int pD, 
                      const int Ash_, const int Bsh_, const int Csh_, const int Dsh_ ) {
        posA_strt = pA;
        posB_strt = pB;
        posC_strt = pC;
        posD_strt = pD;
        
        Ash = Ash_;
        Bsh = Bsh_;
        Csh = Csh_;
        Dsh = Dsh_;
    };
    
    void to_disk( llint & tot_2E, const double thrsh, std::ofstream & ofs ) {
        usint indi,indj,indk,indl;
        double re_data,im_data,znorm;
        usint mA,mB,mC,mD;
        
        for(mA=0; mA<shgA; mA++) {
            indi = posA_strt + mA + 1;
            for(mB=0; mB<shgB; mB++) {
                indj = posB_strt + mB + 1;
                for(mC=0; mC<shgC; mC++) {
                    indk = posC_strt + mC + 1;
                    for(mD=0; mD<shgD; mD++) {
                        indl = posD_strt + mD + 1;

                        if(!is_unique(indi, indj, indk, indl)) continue;

                        re_data = ( v[ mA ][ mB ][ mC ][ mD ] ).real();
                        im_data = ( v[ mA ][ mB ][ mC ][ mD ] ).imag();
                        znorm   = std::abs( v[ mA ][ mB ][ mC ][ mD ] );
                        
                        if( znorm < thrsh ) continue;
                        
                        if( indk > indi ) continue ;
                        if( indj > indi ) continue ;
                        
                        ofs.write(reinterpret_cast<char*>(&indi),sizeof(usint));
                        ofs.write(reinterpret_cast<char*>(&indj),sizeof(usint));
                        ofs.write(reinterpret_cast<char*>(&indk),sizeof(usint));
                        ofs.write(reinterpret_cast<char*>(&indl),sizeof(usint));
                        
                        ofs.write(reinterpret_cast<char*>(&re_data),sizeof(double));
                        ofs.write(reinterpret_cast<char*>(&im_data),sizeof(double));
                        
                        tot_2E++;
                    };
                };
            };
        };
    };
    
    void print( const double thrsh ) {
        usint indi,indj,indk,indl;
        double re_data,im_data,znorm;
        usint mA,mB,mC,mD;
        
        for(mA=0; mA<shgA; mA++) {
            indi = posA_strt + mA + 1;
            for(mB=0; mB<shgB; mB++) {
                indj = posB_strt + mB + 1;
                for(mC=0; mC<shgC; mC++) {
                    indk = posC_strt + mC + 1;
                    for(mD=0; mD<shgD; mD++) {
                        indl = posD_strt + mD + 1;
                        
                        re_data = ( v[ mA ][ mB ][ mC ][ mD ] ).real();
                        im_data = ( v[ mA ][ mB ][ mC ][ mD ] ).imag();
                        znorm   = std::abs( v[ mA ][ mB ][ mC ][ mD ] );
                        
                        if( znorm < thrsh ) continue;
                        
                        if( indk > indi ) continue ;
                        if( indj > indi ) continue ;
                        
                        std::cout << indi << " " << indj << " " << indk << " " << indl << " " << re_data << std::endl;
                    };
                };
            };
        };
    };
    
};

template <class type> class HalfC {
    public:
    int lsum_AB;
    int lC;
    int lD;
    
    int shgC;
    int shgD;
    
    type***** v;
    
    HalfC( const int lsum_AB_, const int lC_, const int lD_ ) :
        lsum_AB( lsum_AB_ ), lC( lC_ ), lD( lD_ ) {
        
        shgC = crt_siz[lC];
        shgD = crt_siz[lD];
        
        v = new type****[lsum_AB+1];
        for(int t1=0; t1<=lsum_AB; t1++) {
            v[t1] = new type***[lsum_AB+1];
            for(int u1=0; u1<=lsum_AB; u1++) {
                v[t1][u1] = new type**[lsum_AB+1];
                for(int v1=0; v1<=lsum_AB; v1++) {
                    v[t1][u1][v1] = new type*[shgC];
                    for(int mk=0; mk<shgC; mk++) {
                        v[t1][u1][v1][mk] = new type[shgD];
                        memset( v[t1][u1][v1][mk], 0, shgD * sizeof( type ) );
                    };
                };
            };
        };
        
    };
    
    ~HalfC( ) {
        
        for(int t1=0; t1<=lsum_AB; t1++) {
            for(int u1=0; u1<=lsum_AB; u1++) {
                for(int v1=0; v1<=lsum_AB; v1++) {
                    for(int mk=0; mk<shgC; mk++) delete []v[t1][u1][v1][mk];
                    delete []v[t1][u1][v1];
                };
                delete []v[t1][u1];
            };
            delete []v[t1];
        };
        delete []v;
        
    };
    
    void zero( ) {
        
        for(int t1=0; t1<=lsum_AB; t1++) {
            for(int u1=0; u1<=lsum_AB; u1++)
                for(int v1=0; v1<=lsum_AB; v1++)
                    for(int mk=0; mk<shgC; mk++)
                        memset( v[t1][u1][v1][mk], 0, shgD * sizeof( type ) );
        };
        
    };
    
};

template <class type> class RInts2E {
    public:
    int lsum_AB;
    int lsum_CD;
    int lsum_PQ;
    
    type kPx;
    type kPy;
    type kPz;
    
    type kQx;
    type kQy;
    type kQz;
    
    type Px;
    type Py;
    type Pz;
    
    type Qx;
    type Qy;
    type Qz;
    
    type aP;
    type aQ;
    
    type******  R2E;
    type******* R2En;
    
    void strip() {
        register int t1,u1,v1,t2,u2,v2;
        for(t1=0; t1<=lsum_AB; t1++) for(u1=0; u1<=lsum_AB-t1; u1++) for(v1=0; v1<=lsum_AB-t1-u1; v1++) 
            for(t2=0; t2<=lsum_CD; t2++) for(u2=0; u2<=lsum_CD-t2; u2++) for(v2=0; v2<=lsum_CD-t2-u2; v2++) 
                R2E[t1][u1][v1][t2][u2][v2] = R2En[t1][u1][v1][t2][u2][v2][0];
    };
    
    
    void load( const type kPx_, const type kPy_, const type kPz_, 
               const type kQx_, const type kQy_, const type kQz_, 
               const type Px_ , const type Py_ , const type Pz_ ,
               const type Qx_ , const type Qy_ , const type Qz_ ,
               const type aP_ , const type aQ_ ) {
        kPx = kPx_; kPy = kPy_; kPz = kPz_;
        kQx = kQx_; kQy = kQy_; kQz = kQz_;
        Px  = Px_ ; Py  = Py_ ; Pz  = Pz_ ;
        Qx  = Qx_ ; Qy  = Qy_ ; Qz  = Qz_ ;
        aP  = aP_ ; aQ  = aQ_ ;
    };
    
    void zero() {
        register int t1,u1,v1,t2,u2,v2;
        for(t1=0; t1<=lsum_AB; t1++) for(u1=0; u1<=lsum_AB; u1++) for(v1=0; v1<=lsum_AB; v1++) 
            for(t2=0; t2<=lsum_CD; t2++) for(u2=0; u2<=lsum_CD; u2++) {
                memset(R2E[t1][u1][v1][t2][u2],0,(lsum_CD+1)*sizeof(type));
                for(v2=0; v2<=lsum_CD; v2++) 
                    memset(R2En[t1][u1][v1][t2][u2][v2],0,(lsum_PQ+1)*sizeof(type));
            };
    };
    
    void print() {
        register int t1,u1,v1,t2,u2,v2;
        for(t1=0; t1<=lsum_AB; t1++) for(u1=0; u1<=lsum_AB-t1; u1++) for(v1=0; v1<=lsum_AB-t1-u1; v1++) 
            for(t2=0; t2<=lsum_CD; t2++) for(u2=0; u2<=lsum_CD-t2; u2++) for(v2=0; v2<=lsum_CD-t2-u2; v2++) 
            std::cout<< t1 << " " << u1 << " " << v1 << " " << t2 << " " << u2 << " " << v2 << "   " << R2E[t1][u1][v1][t2][u2][v2] << std::endl;
    };
    
    RInts2E( const int lsum_AB_, const int lsum_CD_ ) : 
            lsum_AB( lsum_AB_ ), lsum_CD( lsum_CD_ ) {
        lsum_PQ = lsum_AB + lsum_CD;
        register int t1,u1,v1,t2,u2,v2;
        
        R2E  = new type***** [lsum_PQ+1];
        R2En = new type******[lsum_PQ+1];
        for(t1=0; t1<=lsum_PQ; t1++) {
            R2E [t1] = new type**** [lsum_PQ+1];
            R2En[t1] = new type*****[lsum_PQ+1];
            for(u1=0; u1<=lsum_PQ; u1++) {
                R2E [t1][u1] = new type*** [lsum_PQ+1];
                R2En[t1][u1] = new type****[lsum_PQ+1];
                for(v1=0; v1<=lsum_PQ; v1++) {
                    R2E [t1][u1][v1] = new type** [lsum_PQ+1];
                    R2En[t1][u1][v1] = new type***[lsum_PQ+1];
                    for(t2=0; t2<=lsum_PQ; t2++) {
                        R2E [t1][u1][v1][t2] = new type*[lsum_PQ+1];
                        R2En[t1][u1][v1][t2] = new type**[lsum_PQ+1];
                        for(u2=0; u2<=lsum_PQ; u2++) {
                            R2E [t1][u1][v1][t2][u2] = new type [lsum_PQ+1];
                            R2En[t1][u1][v1][t2][u2] = new type*[lsum_PQ+1];
                            memset( R2E[t1][u1][v1][t2][u2], 0, (lsum_PQ+1) * sizeof(type) );
                            for(v2=0; v2<=lsum_PQ; v2++) {
                                R2En[t1][u1][v1][t2][u2][v2] = new type[lsum_PQ+1];
                                memset( R2En[t1][u1][v1][t2][u2][v2], 0, (lsum_PQ+1) * sizeof(type) );
                            };
                        };
                    };
                };
            };
        };
    };
    
   ~RInts2E() {
        register int t1,u1,v1,t2,u2,v2;
        for(t1=0; t1<=lsum_PQ; t1++) {
            for(u1=0; u1<=lsum_PQ; u1++) {
                for(v1=0; v1<=lsum_PQ; v1++) {
                    for(t2=0; t2<=lsum_PQ; t2++) {
                        for(u2=0; u2<=lsum_PQ; u2++) {
                            for(v2=0; v2<=lsum_PQ; v2++) delete []R2En[t1][u1][v1][t2][u2][v2];
                            delete []R2E [t1][u1][v1][t2][u2];
                            delete []R2En[t1][u1][v1][t2][u2];
                        };
                        delete []R2E [t1][u1][v1][t2];
                        delete []R2En[t1][u1][v1][t2];
                    };
                    delete []R2E [t1][u1][v1];
                    delete []R2En[t1][u1][v1];
                };
                delete []R2E [t1][u1];
                delete []R2En[t1][u1];
            };
            delete []R2E [t1];
            delete []R2En[t1];
        };
        delete []R2E;
        delete []R2En;
    };
};

void CalcRERI( RInts2E <cdouble> & R );




#endif //AUXFUN2_H
