#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include "tool.h"
using namespace std;
using namespace NTL;
struct PKE_Para {
    int N,msg_bit,p_bit,q_bit,lift_bit;
    ZZ msg,p,q,coeff,twice_p,twice_q,half_p,conversion_modulus,mod_p;
    ZZ_pX xN;
    int hsk = 64;
    int d;
    int num_data;
};

//eval_poly=1, eval a poly
void SetPara(PKE_Para & pkePara);

void PKE_Gen(PKE_Para &pkePara, vec_ZZ_pX &pkePk, vec_ZZ_pX &pkeSk) {
    //initialize the parameters
    SetPara(pkePara);

    NTL::power(pkePara.p, 2, pkePara.p_bit);
    NTL::power(pkePara.q, 2, pkePara.q_bit);
    ZZ_p::init(pkePara.q);

    SetCoeff(pkePara.xN, 0, 1);
    SetCoeff(pkePara.xN, pkePara.N, 1);
    ZZ_pXModulus modulus(pkePara.xN);

    pkePara.twice_p=2*pkePara.p;
    pkePara.twice_q=2*pkePara.q;
    pkePara.half_p=pkePara.p/2;
    ZZ constant1=conv<ZZ>("1");
    pkePara.mod_p=pkePara.p-constant1;

    //gen sk
    ZZ_pX hat_s, e;
    Random_ZZ_pX(pkePk[0], pkePara.N, pkePara.q_bit);
    SecretKey(hat_s, pkePara.N, pkePara.hsk);
    GaussRand(e, pkePara.N);

    MulMod(pkePk[1], pkePk[0], hat_s, modulus);
    pkePk[1] = pkePk[1] + e;
    //gen sk
    SetCoeff(pkeSk[0], 0, 1);
    pkeSk[1] = hat_s;
    
}

void PKE_Enc(vec_ZZ_pX &c, const PKE_Para pkePara, ZZ_pXModulus modulus,vec_ZZ_pX pkePk, const ZZ &x) {
    ZZ_pX v, e1, e2, x_ZZ_pX;
    ZZ q_div_p = pkePara.q / pkePara.p;
    ZZ_p coeff;
    SecretKey(v, pkePara.N, pkePara.hsk);
    GaussRand(e1, pkePara.N);
    GaussRand(e2, pkePara.N);
    MulMod(c[0], pkePk[1], v, modulus);
    c[0] = c[0] + e1;
    conv(coeff,q_div_p*x);
    SetCoeff(x_ZZ_pX, 0, coeff);
    c[0] = c[0] + x_ZZ_pX;
    MulMod(c[1], pkePk[0], v, modulus);
    c[1] = e2 - c[1];
}

void PKE_OKDM(vec_ZZ_pX &C, const PKE_Para &pkePara, ZZ_pXModulus& modulus,vec_ZZ_pX &pkePk, const ZZ &x) {
    ZZ zero;
    zero = 0;
    ZZ q_div_p = pkePara.q / pkePara.p;
    vec_ZZ_pX c_xs1, c_xs2;
    ZZ_p coeff;
    ZZ_pX x_ZZ_pX;
    c_xs1.SetLength(2);
    c_xs2.SetLength(2);

    PKE_Enc(c_xs1, pkePara, modulus,pkePk, x);
    PKE_Enc(c_xs2, pkePara, modulus,pkePk, zero);
    conv(coeff,q_div_p*x);
    SetCoeff(x_ZZ_pX, 0, coeff);
    c_xs2[1] = c_xs2[1] + x_ZZ_pX;
    C[0] = c_xs1[0];
    C[1] = c_xs1[1];
    C[2] = c_xs2[0];
    C[3] = c_xs2[1];
}

void PKE_DDec(vec_ZZ_pX &db, const PKE_Para pkePara, ZZ_pXModulus modulus,vec_ZZ_pX pkeSk, vec_ZZ_pX C) {
    ZZ coeff;
    ZZX temp;
    ZZ_pX temp1, temp2;
    MulMod(temp1, pkeSk[0], C[0], modulus);
    MulMod(temp2, pkeSk[1], C[1], modulus);
    temp1 = temp1 + temp2;
    conv(temp,temp1);
    for (int i = 0; i < pkePara.N; i++) {
        GetCoeff(coeff, temp, i);
        coeff = (coeff*pkePara.twice_p+pkePara.q)/(pkePara.twice_q);
        coeff=coeff%pkePara.p;
        if(coeff>pkePara.half_p)
        {
            coeff-=pkePara.p;
        }
        SetCoeff(temp, i, coeff);
    }
    conv(db[0],temp);
    MulMod(temp1, pkeSk[0], C[2], modulus);
    MulMod(temp2, pkeSk[1], C[3], modulus);
    temp1 = temp1 + temp2;
    conv(temp,temp1);
    for (int i = 0; i < pkePara.N; i++) {
        GetCoeff(coeff, temp, i);
        coeff = (coeff*pkePara.twice_p+pkePara.q)/(pkePara.twice_q);
        coeff=coeff%pkePara.p;
        if(coeff>pkePara.half_p)
        {
            coeff-=pkePara.p;
        }
        SetCoeff(temp, i, coeff);
    }
    conv(db[1],temp);
}

void SetPara(PKE_Para & pkePara) {
    pkePara.conversion_modulus=conv<ZZ>("340282366920938463463374607431768211297");
    pkePara.N=8192;
    pkePara.p_bit=62;
    pkePara.q_bit=146;
    pkePara.lift_bit=pkePara.q_bit-pkePara.p_bit;
}