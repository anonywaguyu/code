#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include "PKE.h"
#include <vector>

using namespace std;
using namespace NTL;


void Setup(vec_ZZ_pX &hssEk_1, vec_ZZ_pX &hssEk_2,
    PKE_Para &pkePara, vec_ZZ_pX &pkePk, vec_ZZ_pX &pkeSk,int eval_poly) {
    
    PKE_Gen(pkePara, pkePk, pkeSk,eval_poly);
    ZZ_pXModulus modulus(pkePara.xN);
    Random_ZZ_pX(hssEk_1[0], pkePara.N, pkePara.q_bit);
    Random_ZZ_pX(hssEk_1[1], pkePara.N, pkePara.q_bit);

    hssEk_2[0] = pkeSk[0] - hssEk_1[0];
    hssEk_2[1] = pkeSk[1] - hssEk_1[1];
    //cout<<"sk0:"<<pkeSK[0]<<endl;

}

void Keygen(vec_ZZ_pX &C, const PKE_Para &pkePara, ZZ_pXModulus &modulus, vec_ZZ_pX &pkePk, const ZZ &x) {
    PKE_OKDM(C, pkePara, modulus, pkePk, x);
}

void Probgen(vec_ZZ_pX &C, const PKE_Para &pkePara, ZZ_pXModulus &modulus, vec_ZZ_pX &pkePk, const ZZ &x) {
    PKE_OKDM(C, pkePara, modulus, pkePk, x);
}
//Mult&Load
void Mult(vec_ZZ_pX &db, const PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX pkeSk, vec_ZZ_pX C) {
    PKE_DDec(db, pkePara, modulus, pkeSk, C);
}

//Output
void Output(vec_ZZ_pX &db, const PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX pkeSk, vec_ZZ_pX C) {
    PKE_DDec2(db, pkePara, modulus, pkeSk, C);
}

int prfkey = 1;
ZZ static_temp;
vec_ZZ_pX tb_static, F_static;

// void f(vec_ZZ_pX &tb, int b, int d, int num_data, int loop, int beg_ind, int *ind_var,
//        PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX ek, Vec<vec_ZZ_pX> C_X, Vec<vec_ZZ_pX> C_F, Vec<vec_ZZ_pX> PRF, int &count) {
//     if (loop == d) {
//         if (d == 1) {
//             Mult(tb_static, pkePara, modulus, ek, C_X[ind_var[0]]);
//             Mult(tb_static, pkePara, modulus, tb_static, C_F[count%10]);
//             tb[0] = tb[0] + tb_static[0];
//             tb[1] = tb[1] + tb_static[1];
//             count++;
//         } 
//         else {
//             Mult(tb_static, pkePara, modulus, ek, C_X[ind_var[0]]);
//             for (int i = 1; i < d; i++) {
//                 Mult(tb_static, pkePara, modulus, tb_static, C_X[ind_var[i]]);
//             }
//             Mult(tb_static, pkePara, modulus, tb_static, C_F[count%10]);
//             tb[0] = tb[0] + tb_static[0];
//             tb[1] = tb[1] + tb_static[1];
//             count++;
//         }
//     } else {
//         loop = loop + 1;
//         for (int i = beg_ind; i < num_data; i++) {
//             ind_var[loop - 1] = i;
//             f(tb, b, d, num_data, loop, i, ind_var,
//               pkePara, modulus, ek, C_X, C_F, PRF, count);
//         }
//     }
// }


// void Compute(ZZ &tb_y, ZZ &tb_a, ZZ &tb_pi, int b, PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX ek, Vec<vec_ZZ_pX> C_X, Vec<vec_ZZ_pX> C_F,
//               Vec<vec_ZZ_pX> PRF, vec_ZZ_pX C1, vec_ZZ_pX Ca) {
//     int count = 0;
//     int loop = 0;
//     int beg_ind = 0;
//     int *ind_var = (int *) malloc(sizeof(int) * pkePara.d);
//     vec_ZZ_pX tb, tb_temp;
//     tb_temp.SetLength(2);
//     ZZX temp;
//     tb_static.SetLength(2);
//     tb.SetLength(2);
//     cout<<"daole"<<endl;
//     for (int i = 1; i < pkePara.d + 1; i++) {
//         f(tb_temp, b, i, pkePara.num_data, loop, beg_ind, ind_var, pkePara, modulus, ek, C_X, C_F, PRF, count);
//         prfkey = 1;
//         tb[0] = tb[0] + tb_temp[0];
//         tb[1] = tb[1] + tb_temp[1];
//         tb_temp[0] = 0;
//         tb_temp[1] = 0;
//         cout<<"here"<<endl;
//     }
//     cout<<"items:"<<count<<endl;
//     //y
//     Mult(tb, pkePara, modulus, tb, C1);
//     conv(temp,tb[0]);
//     GetCoeff(tb_y, temp, 0);
    
//     ZZ_pX temp1;
//     ZZ_pX e1;
//     GaussRand(e1, pkePara.N);
//     MulMod(temp1, ek[1], Ca[1], modulus);
//     temp1=e1-temp1;
//     conv(temp,temp1);
//     GetCoeff(tb_a, temp, 0);

//     GaussRand(e1, pkePara.N);
//     MulMod(temp1, tb[1], Ca[1], modulus);
//     temp1=e1+temp1;
//     conv(temp,temp1);
//     GetCoeff(tb_pi, temp, 0);
// }

void Compute(ZZ &tb_y, ZZ &tb_a, ZZ &tb_pi, int b, PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX ek, Vec<vec_ZZ_pX> C_X, Vec<vec_ZZ_pX> C_F,
    Vec<vec_ZZ_pX> PRF, vec_ZZ_pX C1, vec_ZZ_pX Ca, vector<vector<int>> &F, int item) {
    int count = 0;
    vec_ZZ_pX tb, tb_temp;
    tb_temp.SetLength(2);
    tb.SetLength(2);

    for (int i = 0; i < item; i++)
    {
        Mult(tb_temp, pkePara, modulus, ek, C_F[i%10]);
        if(i%10==0)cout<<i<<endl;
        for (int j = 0; j < 4; j++)
        {
            while(F[i][j]!=0)
            {
                Mult(tb_temp, pkePara, modulus, tb_temp, C_X[j]);
                F[i][j]--;
            }
        }
        tb[0] = tb[0] + tb_temp[0];
        tb[1] = tb[1] + tb_temp[1];
    }
    cout<<"items:"<<item<<endl;
    //y
    ZZX temp;
    ZZ_pX temp1;
    ZZ_pX e1;
    GaussRand(e1, pkePara.N);
    MulMod(temp1, ek[1], Ca[1], modulus);
    temp1=e1-temp1;
    conv(temp,temp1);
    GetCoeff(tb_a, temp, 0);

    GaussRand(e1, pkePara.N);
    MulMod(temp1, tb[1], Ca[1], modulus);
    temp1=e1+temp1;
    conv(temp,temp1);
    GetCoeff(tb_pi, temp, 0);

    Output(tb, pkePara, modulus, tb, C1);
    conv(temp,tb[0]);
    GetCoeff(tb_y, temp, 0);
}