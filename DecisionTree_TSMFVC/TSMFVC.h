#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include "PKE.h"

using namespace std;
using namespace NTL;



void Setup(vec_ZZ_pX &hssEk_1, vec_ZZ_pX &hssEk_2,
             PKE_Para&pkePara, vec_ZZ_pX &pkePk, vec_ZZ_pX &pkeSk) {

    PKE_Gen(pkePara, pkePk, pkeSk);
    ZZ_pXModulus modulus(pkePara.xN);
    
    Random_ZZ_pX(hssEk_1[0], pkePara.N, pkePara.q_bit);
    Random_ZZ_pX(hssEk_1[1], pkePara.N, pkePara.q_bit);

    hssEk_2[0] = pkeSk[0] - hssEk_1[0];
    hssEk_2[1] = pkeSk[1] - hssEk_1[1];
}

void Keygen(vec_ZZ_pX &C, const PKE_Para &pkePara, ZZ_pXModulus &modulus, vec_ZZ_pX &pkePk, const ZZ &x) {
    PKE_OKDM(C, pkePara, modulus, pkePk, x);
}

void Probgen(vec_ZZ_pX &C, const PKE_Para &pkePara, ZZ_pXModulus &modulus, vec_ZZ_pX &pkePk, const ZZ &x) {
    PKE_OKDM(C, pkePara, modulus, pkePk, x);
}

//Load&Mult
void Mult(vec_ZZ_pX &db, const PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX pkeSk, vec_ZZ_pX C) {
    PKE_DDec(db, pkePara, modulus, pkeSk, C);
}

void Add(vec_ZZ_pX &out, vec_ZZ_pX in1,vec_ZZ_pX in2) {
    out[0]=in1[0]+in2[0];
    out[1]=in1[1]+in2[1];
}

void Sub(vec_ZZ_pX &out, vec_ZZ_pX in1,vec_ZZ_pX in2) {
    out[0]=in1[0]-in2[0];
    out[1]=in1[1]-in2[1];
}

int prfkey = 1;

void f(vec_ZZ_pX &tb, int b, int d, int num_data, int loop, int beg_ind, int *ind_var,
       PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX ek, Vec<vec_ZZ_pX> C_X, Vec<vec_ZZ_pX> C_F, Vec<vec_ZZ_pX> PRF, int &count) {
    if (loop == d) {
        vec_ZZ_pX tb_temp;
        tb_temp.SetLength(2);
        if (d == 1) {
            Mult(tb_temp, pkePara, modulus, ek, C_X[ind_var[0]]);
            Mult(tb_temp, pkePara, modulus, tb_temp, C_F[count%10]);

            tb[0] = tb[0] + tb_temp[0];
            tb[1] = tb[1] + tb_temp[1];
            count++;

        } else {
            Mult(tb_temp, pkePara, modulus, ek, C_X[ind_var[0]]);
            for (int i = 1; i < d; i++) {
                Mult(tb_temp, pkePara, modulus, tb_temp, C_X[ind_var[i]]);
            }
            Mult(tb_temp, pkePara, modulus, tb_temp, C_F[count%10]);

            tb[0] = tb[0] + tb_temp[0];
            tb[1] = tb[1] + tb_temp[1];
            count++;
        }
    } else {
        loop = loop + 1;
        for (int i = beg_ind; i < num_data; i++) {
            ind_var[loop - 1] = i;
            f(tb, b, d, num_data, loop, i, ind_var,
              pkePara, modulus, ek, C_X, C_F, PRF, count);
        }
    }
}

void Compute(ZZ &tb_y, ZZ &tb_a, ZZ &tb_pi, int b, PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX ek, Vec<vec_ZZ_pX> C_X, Vec<vec_ZZ_pX> C_F,
              Vec<vec_ZZ_pX> PRF, vec_ZZ_pX C1, vec_ZZ_pX Ca) {
    int count = 0;
    int loop = 0;
    int beg_ind = 0;
    int *ind_var = (int *) malloc(sizeof(int) * pkePara.d);
    vec_ZZ_pX tb, tb_temp;
    ZZX temp;
    tb.SetLength(2);
    tb_temp.SetLength(2);
    
    for (int i = 1; i < pkePara.d + 1; i++) {
        f(tb_temp, b, i, pkePara.num_data, loop, beg_ind, ind_var, pkePara, modulus, ek, C_X, C_F, PRF, count);
        prfkey = 1;
        tb[0] = tb[0] + tb_temp[0];
        tb[1] = tb[1] + tb_temp[1];
        tb_temp[0] = 0;
        tb_temp[1] = 0;
        cout<<"here"<<endl;
    }
    cout<<"items:"<<count<<endl;

    //y
    Mult(tb, pkePara, modulus, tb, C1);
    conv(temp,tb[0]);
    GetCoeff(tb_y, temp, 0);
    
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
}

void Verify(ZZ y1, ZZ y2, ZZ tba1, ZZ tba2, ZZ tbpi1, ZZ tbpi2, PKE_Para pkePara, int &flag, ZZ &y, ZZ &a, ZZ& c_1, ZZ&pi, ZZ& ay, ZZ&pi_)
{
    y = (y1 + y2)%pkePara.p;
    a=RandomBits_ZZ(128);
    c_1=(tba1+tba2+(a<<(pkePara.lift_bit)));
    pi = (((y*c_1+tbpi1+tbpi2)*pkePara.twice_p+pkePara.q)>>(1+pkePara.q_bit))%pkePara.p;
    ay=a*y%pkePara.p;
    pi_ = ((((y*c_1+tbpi1)*pkePara.twice_p+pkePara.q)>>(1+pkePara.q_bit))+
        ((tbpi2*pkePara.twice_p+pkePara.q)>>(1+pkePara.q_bit)))%pkePara.p;
    if (pi==ay && pi==pi_) ;
    else {
        printf("******************** ERROR ********************\n");
        flag=1;
        cout << "The value of y   is  :" << y << endl;
        printf("\n\n");
        cout << "THE value of pi  is  :" << pi <<endl;
    }
}