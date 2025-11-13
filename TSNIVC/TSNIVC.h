#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include "PKE.h"

using namespace std;
using namespace NTL;


void Setup(vec_ZZ_pX &hssEk_1, vec_ZZ_pX &hssEk_2,
    PKE_Para &pkePara, vec_ZZ_pX &pkePk, vec_ZZ_pX &pkeSk,int eval_poly) {

    PKE_Gen(pkePara, pkePk, pkeSk,eval_poly);
    ZZ_pXModulus modulus(pkePara.xN);
    ZZ_p::init(pkePara.q);
    ZZ_pX temp;
    MulMod(temp, pkeSk[1], pkeSk[1], modulus);
    
    Random_ZZ_pX(hssEk_1[0], pkePara.N, pkePara.q_bit);
    Random_ZZ_pX(hssEk_1[1], pkePara.N, pkePara.q_bit);
    Random_ZZ_pX(hssEk_1[2], pkePara.N, pkePara.q_bit);
    Random_ZZ_pX(hssEk_1[3], pkePara.N, pkePara.q_bit);

    hssEk_2[0] = pkeSk[0] - hssEk_1[0];
    hssEk_2[1] = pkeSk[1] - hssEk_1[1];
    hssEk_2[2] = pkeSk[1] - hssEk_1[2];
    hssEk_2[3] = temp - hssEk_1[3];
    //cout<<"sk0:"<<pkeSK[0]<<endl;
}

void Keygen(vec_ZZ_pX &C, const PKE_Para &pkePara, ZZ_pXModulus &modulus, vec_ZZ_pX &pkePk, const ZZ &x) {
    PKE_Enc(C, pkePara, modulus, pkePk, x);
}
//Load
void Load(vec_ZZ_pX &db, const PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX pkeSk, vec_ZZ_pX C) {
    PKE_DDec(db, pkePara, modulus, pkeSk, C);
}

//Output
void Output(vec_ZZ_pX &db, const PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX pkeSk, vec_ZZ_pX C) {
    PKE_DDec2(db, pkePara, modulus, pkeSk, C);
}

int prfkey = 1;
ZZ static_temp;
vec_ZZ_pX tb_static, F_static;

void f(Vec<vec_ZZ_pX> &tb, int b, int d, int num_data, int loop, int beg_ind, int *ind_var,
       PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX ek, Vec<vec_ZZ> X, Vec<vec_ZZ_pX> C_F, Vec<vec_ZZ_pX> PRF, int &count, ZZ r) {
    if (loop == d) {
        Load(F_static, pkePara, modulus, ek, C_F[count%10]);
        tb_static=F_static;
        if (d == 1) {
            for(int j=0;j<pkePara.d+1;j++)
            {
                static_temp= X[j][ind_var[0]];
                for (long i = 0; i <= deg(F_static[0]); i++) {
                    tb_static[0].rep[i] = F_static[0].rep[i]*conv<ZZ_p>(static_temp);
                    tb_static[1].rep[i] = F_static[1].rep[i]*conv<ZZ_p>(static_temp);
                }
                tb[j][0] = tb[j][0] + tb_static[0];
                tb[j][1] = tb[j][1] + tb_static[1];
            } 
            count++;
        } else {
            for(int j=0;j<pkePara.d+1;j++)
            {
                static_temp= X[j][ind_var[0]];
                for (int i = 1; i < d; i++) {
                    static_temp = (static_temp * X[j][ind_var[i]])%r;
                }
                for (long i = 0; i <= deg(F_static[0]); i++) {
                    tb_static[0].rep[i] = F_static[0].rep[i]*conv<ZZ_p>(static_temp);
                    tb_static[1].rep[i] = F_static[1].rep[i]*conv<ZZ_p>(static_temp);
                }
                tb[j][0] = tb[j][0] + tb_static[0];
                tb[j][1] = tb[j][1] + tb_static[1];
            }
            count++;
        }
    } else {
        loop = loop + 1;
        for (int i = beg_ind; i < num_data; i++) {
            ind_var[loop - 1] = i;
            f(tb, b, d, num_data, loop, i, ind_var,
              pkePara, modulus, ek, X, C_F, PRF, count, r);
        }
    }
}


void Compute(vec_ZZ &tb_y, ZZ &tb_a, vec_ZZ &tb_pi, int b, PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX ek, Vec<vec_ZZ> X, Vec<vec_ZZ_pX> C_F,
              Vec<vec_ZZ_pX> PRF, vec_ZZ_pX C1, vec_ZZ_pX Ca, ZZ r, vec_ZZ v) {
    tb_static.SetLength(2);
    F_static.SetLength(2);
    int count = 0;
    int loop = 0;
    int beg_ind = 0;
    int *ind_var = (int *) malloc(sizeof(int) * pkePara.d);
    Vec<vec_ZZ_pX> tb, tb_temp;
    tb.SetLength(pkePara.d+1);
    tb_temp.SetLength(pkePara.d+1);
    for(int i=0;i<pkePara.d+1;i++)
    {
        tb[i].SetLength(2);
        tb_temp[i].SetLength(2);
    }
    ZZX temp;
    vec_ZZ Inv_v;
    Inv_v.SetLength(pkePara.d+1);
    for (int i = 1; i < pkePara.d + 1; i++) {
        f(tb_temp, b, i, pkePara.num_data, loop, beg_ind, ind_var, pkePara, modulus, ek, X, C_F, PRF, count, r);
        prfkey = 1;
        for(int j=0;j<pkePara.d+1;j++)
        {
            tb[j][0] = tb[j][0] + tb_temp[j][0];
            tb[j][1] = tb[j][1] + tb_temp[j][1];
            tb_temp[j][0] = 0;
            tb_temp[j][1] = 0;
        }
        cout<<"here"<<endl;
    }
    cout<<"items:"<<count<<endl;
    //y
    
    
    ZZ_pX temp1;
    ZZ_pX e1;
    GaussRand(e1, pkePara.N);
    MulMod(temp1, ek[1], Ca[1], modulus);
    temp1=e1-temp1;
    conv(temp,temp1);
    GetCoeff(tb_a, temp, 0);

    for(int j=0;j<pkePara.d+1;j++)
    {
        v[j]%=r;
        InvMod(Inv_v[j], v[j], r);
        GaussRand(e1, pkePara.N);
        MulMod(temp1, tb[j][1], Ca[1], modulus);
        temp1=e1+temp1;
        conv(temp,temp1);
        GetCoeff(tb_pi[j], temp, 0);
        tb_pi[j]*=Inv_v[j];
    }

    for(int j=0;j<pkePara.d+1;j++)
    {
        Output(tb[j], pkePara, modulus, tb[j], C1);
        conv(temp,tb[j][0]);
        GetCoeff(tb_y[j], temp, 0);
        tb_y[j]*=Inv_v[j];
    }
}