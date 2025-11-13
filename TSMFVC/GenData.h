#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/matrix.h>
#include <NTL/vec_vec_ZZ.h>


using namespace std;
using namespace NTL;
struct Data {
    vec_ZZ X;
    vec_ZZ F;
    Vec<vec_ZZ_pX> C_X;
    Vec<vec_ZZ_pX> C_F;
    Vec<vec_ZZ_pX> PRF;
};

void GenData(Data &data, PKE_Para pkePara, vec_ZZ_pX pkePk) {
    data.X.SetLength(pkePara.num_data);
    vec_ZZ_pX C_x, prf;
    C_x.SetLength(4);
    prf.SetLength(2);
    ZZ_pXModulus modulus(pkePara.xN);
    for (int i = 0; i < pkePara.num_data; i++) {
        RandomBits(data.X[i], pkePara.msg_bit);
        Probgen(C_x, pkePara, modulus, pkePk, data.X[i]);
        data.C_X.append(C_x);
    }
    for (int i = 0; i < 10; i++) {
        Random_ZZ_pX(prf[0], pkePara.N, pkePara.q_bit);
        Random_ZZ_pX(prf[1], pkePara.N, pkePara.q_bit);
        data.PRF.append(prf);
    }
}

void GenDataF(Data &data, PKE_Para pkePara, vec_ZZ_pX pkePk) {
    ZZ_pXModulus modulus(pkePara.xN);
    int a=1,b=1;
    for(int i=1;i<=pkePara.num_data;i++)
    {
        a=a*i;
        b=b*(i+pkePara.d);
    }
    a=b/a;
    //cout<<"items:"<<a<<endl;
    data.F.SetLength(10);
    vec_ZZ_pX C_f;
    C_f.SetLength(4);
    for (int i = 0; i < 10; i++) {
        RandomBits(data.F[i], pkePara.msg_bit);
        Keygen(C_f, pkePara, modulus, pkePk, data.F[i]);
        data.C_F.append(C_f);
    }
    // for (int i = 0; i < a; i++) {
    //     data.F[i]=data.F[i%10];
    // }
}
