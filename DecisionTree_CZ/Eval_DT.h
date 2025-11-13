#include <gmp.h>
extern "C" {
#include <relic/relic.h>
}

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>

#include <pbc/pbc.h>

void Depth2Tree(vec_ZZ_pX &res, vec_ZZ_pX p, vec_ZZ_pX lc, vec_ZZ_pX rc,
                vec_ZZ_pX tb_one, PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX ek);

void Depth2Tree_Type1(vec_ZZ_pX &res, vec_ZZ_pX p, vec_ZZ_pX lc, vec_ZZ_pX rc,
                      vec_ZZ_pX tb_one, PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX ek);

void Depth2Tree_Type2(vec_ZZ_pX &res, vec_ZZ_pX p, vec_ZZ_pX lc, vec_ZZ_pX rc,
                      vec_ZZ_pX tb_one, PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX ek);

void Depth2Tree_Type3(vec_ZZ_pX &res, vec_ZZ_pX p, vec_ZZ_pX lc, vec_ZZ_pX rc,
                      vec_ZZ_pX tb_one, PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX ek);

void Eval_DT(ZZ_pX &tby, int B, ep_t &g1T1, ep2_t &g2T2, int model, int d, CZ_Para pvhssPara, PKE_Para pkePara,
             ZZ_pXModulus modulus,
             vec_ZZ_pX ek,
             vec_ZZ_pX C_alpha, Vec<vec_ZZ_pX> b, Vec<vec_ZZ_pX> c, vec_ZZ_pX C_one) {
    vec_ZZ_pX temp1, temp2, temp3, tb_one, tb;
    temp1.SetLength(2);
    temp2.SetLength(2);
    temp3.SetLength(2);
    tb_one.SetLength(2);
    tb.SetLength(2);
    ZZ_pX tb0_ZZ_pX;
    if (model == 2) {
        switch (d) {
            case 4:
                Mult(tb_one, pkePara, modulus, ek, C_one);//load 1
                Depth2Tree(temp1, b[3], c[2], c[3], tb_one, pkePara, modulus, ek);
                Depth2Tree(temp2, b[4], c[4], c[5], tb_one, pkePara, modulus, ek);
                Depth2Tree_Type3(temp1, b[1], temp1, temp2, tb_one, pkePara, modulus, ek);
                Depth2Tree(temp2, b[2], c[0], c[1], tb_one, pkePara, modulus, ek);
                Depth2Tree_Type3(temp1, b[0], temp1, temp2, tb_one, pkePara, modulus, ek);
                break;
            case 6:
                Mult(tb_one, pkePara, modulus, ek, C_one);//load 1
                Depth2Tree(temp1, b[5], c[5], c[6], tb_one, pkePara, modulus, ek);
                Depth2Tree_Type2(temp1, b[4], c[4], temp1, tb_one, pkePara, modulus, ek);
                Depth2Tree_Type2(temp1, b[3], c[3], temp1, tb_one, pkePara, modulus, ek);
                Depth2Tree(temp2, b[2], c[1], c[2], tb_one, pkePara, modulus, ek);
                Depth2Tree_Type3(temp3, b[1], temp2, temp1, tb_one, pkePara, modulus, ek);
                Depth2Tree_Type1(temp1, b[0], temp3, c[0], tb_one, pkePara, modulus, ek);
                break;
            case 8:
                Mult(tb_one, pkePara, modulus, ek, C_one);//load 1
                Depth2Tree(temp1, b[8], c[8], c[9], tb_one, pkePara, modulus, ek);
                Depth2Tree_Type2(temp1, b[7], c[7], temp1, tb_one, pkePara, modulus, ek);
                Depth2Tree_Type2(temp1, b[6], c[6], temp1, tb_one, pkePara, modulus, ek);
                Depth2Tree_Type2(temp1, b[5], c[5], temp1, tb_one, pkePara, modulus, ek);
                Depth2Tree_Type1(temp1, b[3], temp1, c[2], tb_one, pkePara, modulus, ek);
                Depth2Tree(temp2, b[4], c[3], c[4], tb_one, pkePara, modulus, ek);
                Depth2Tree_Type1(temp2, b[2], temp2, c[1], tb_one, pkePara, modulus, ek);
                Depth2Tree_Type3(temp1, b[1], temp2, temp1, tb_one, pkePara, modulus, ek);
                Depth2Tree_Type1(temp1, b[0], temp1, c[0], tb_one, pkePara, modulus, ek);
                break;
            default:
                printf("error\n");
                break;
        }
    } else {
        switch (d) {
            case 4:
                Mult(tb_one, pkePara, modulus, ek, C_one);//load 1
                Depth2Tree(temp1, b[3], c[3], c[4], tb_one, pkePara, modulus, ek);
                Depth2Tree_Type1(temp1, b[2], temp1, c[2], tb_one, pkePara, modulus, ek);
                Depth2Tree(temp2, b[1], c[0], c[1], tb_one, pkePara, modulus, ek);
                Depth2Tree_Type3(temp1, b[0], temp2, temp1, tb_one, pkePara, modulus, ek);
                break;
            case 6:
                Mult(tb_one, pkePara, modulus, ek, C_one);//load 1
                Depth2Tree(temp1, b[7], c[6], c[7], tb_one, pkePara, modulus, ek);
                Depth2Tree_Type2(temp1, b[5], c[4], temp1, tb_one, pkePara, modulus, ek);

                Depth2Tree(temp2, b[8], c[8], c[9], tb_one, pkePara, modulus, ek);
                Depth2Tree_Type1(temp2, b[6], temp2, c[5], tb_one, pkePara, modulus, ek);

                Depth2Tree_Type3(temp1, b[3], temp1, temp2, tb_one, pkePara, modulus, ek);
                Depth2Tree(temp2, b[4], c[2], c[3], tb_one, pkePara, modulus, ek);
                Depth2Tree_Type3(temp1, b[2], temp1, temp2, tb_one, pkePara, modulus, ek);
                Depth2Tree(temp2, b[1], c[0], c[1], tb_one, pkePara, modulus, ek);
                Depth2Tree_Type3(temp1, b[0], temp2, temp1, tb_one, pkePara, modulus, ek);

                break;
            case 8:
                Mult(tb_one, pkePara, modulus, ek, C_one);//load 1
                Depth2Tree(temp1, b[12], c[11], c[12], tb_one, pkePara, modulus, ek);
                Depth2Tree_Type1(temp1, b[10], temp1, c[9], tb_one, pkePara, modulus, ek);
                Depth2Tree_Type1(temp1, b[8], temp1, c[7], tb_one, pkePara, modulus, ek);

                Depth2Tree(temp2, b[13], c[13], c[14], tb_one, pkePara, modulus, ek);
                Depth2Tree_Type1(temp2, b[11], temp2, c[10], tb_one, pkePara, modulus, ek);
                Depth2Tree_Type1(temp2, b[9], temp2, c[8], tb_one, pkePara, modulus, ek);

                Depth2Tree_Type3(temp1, b[6], temp1, temp2, tb_one, pkePara, modulus, ek);
                Depth2Tree(temp2, b[7], c[5], c[6], tb_one, pkePara, modulus, ek);
                Depth2Tree_Type3(temp1, b[5], temp1, temp2, tb_one, pkePara, modulus, ek);
                Depth2Tree(temp2, b[4], c[3], c[4], tb_one, pkePara, modulus, ek);
                Depth2Tree_Type3(temp1, b[2], temp2, temp1, tb_one, pkePara, modulus, ek);/////

                Depth2Tree(temp2, b[3], c[1], c[2], tb_one, pkePara, modulus, ek);
                Depth2Tree_Type1(temp2, b[1], temp1, c[0], tb_one, pkePara, modulus, ek);
                Depth2Tree_Type3(temp1, b[0], temp2, temp1, tb_one, pkePara, modulus, ek);

                break;
            default:
                printf("error\n");
                break;
        }
    }


    tby = temp1[0];


    //output
    ZZX tb0;
    bn_t Tb_bn;
    bn_new(Tb_bn);

    ZZ Tb;
    Mult(tb, pkePara, modulus, temp1, C_alpha);
    conv(tb0, tb[0]);
    eval_ZZX(Tb, tb0);
    if (B == 1) {
        Tb = Tb % pvhssPara.g1_order_ZZ;
        ZZ2bn(Tb_bn, Tb);
        ep_mul_gen(g1T1, Tb_bn);
    } else {
        Tb = Tb % pvhssPara.g2_order_ZZ;
        ZZ2bn(Tb_bn, Tb);
        ep2_mul_gen(g2T2, Tb_bn);
    }
}

void Depth2Tree(vec_ZZ_pX &res, vec_ZZ_pX p, vec_ZZ_pX lc, vec_ZZ_pX rc,
                vec_ZZ_pX tb_one, PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX ek) {
    vec_ZZ_pX temp, temp1;
    temp.SetLength(2);
    temp1.SetLength(2);
    Mult(temp, pkePara, modulus, ek, p);//load p
    Sub(temp1, tb_one, temp);//1-p
    Mult(temp1, pkePara, modulus, temp1, lc);//(1-p)*lc
    Mult(res, pkePara, modulus, temp, rc);//p*rc
    Add(res, res, temp1);
}

//lc 2 rc 4
void Depth2Tree_Type1(vec_ZZ_pX &res, vec_ZZ_pX p, vec_ZZ_pX lc, vec_ZZ_pX rc,
                      vec_ZZ_pX tb_one, PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX ek) {
    vec_ZZ_pX temp, temp1;
    temp.SetLength(2);
    temp1.SetLength(2);
    Mult(temp, pkePara, modulus, lc, p);//temp=p*lc
    Sub(temp, lc, temp);//temp=lc-p*lc

    Mult(temp1, pkePara, modulus, ek, p);//load p
    Mult(temp1, pkePara, modulus, temp1, rc);//p*rc
    Add(res, temp, temp1);
}

//lc 4 rc 2
void Depth2Tree_Type2(vec_ZZ_pX &res, vec_ZZ_pX p, vec_ZZ_pX lc, vec_ZZ_pX rc,
                      vec_ZZ_pX tb_one, PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX ek) {
    vec_ZZ_pX temp, temp1;
    temp.SetLength(2);
    temp1.SetLength(2);
    Mult(temp, pkePara, modulus, ek, p);//load p
    Sub(temp, tb_one, temp);//1-p
    Mult(temp, pkePara, modulus, temp, lc);//(1-p)*lc

    Mult(temp1, pkePara, modulus, rc, p);//p*rc
    Add(res, temp, temp1);
}

//lc 2 rc 2
void Depth2Tree_Type3(vec_ZZ_pX &res, vec_ZZ_pX p, vec_ZZ_pX lc, vec_ZZ_pX rc,
                      vec_ZZ_pX tb_one, PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX ek) {
    vec_ZZ_pX temp, temp1;
    temp.SetLength(2);
    temp1.SetLength(2);
    Mult(temp, pkePara, modulus, lc, p);//temp=p*lc
    Sub(temp, lc, temp);//temp=lc-p*lc

    Mult(temp1, pkePara, modulus, rc, p);//temp1=p*rc
    Add(res, temp, temp1);
}

void DataProcess(double &mean, double &stdev, double *Time, int cyctimes) {
    double temp;
    double sum = 0;
    for (int i = 0; i < cyctimes; i++) {
        sum = sum + Time[i];
    }
    mean = sum / cyctimes;

    double temp_sum = 0;
    for (int i = 0; i < cyctimes; i++) {
        temp = mean - Time[i];
        temp = temp * temp;
        temp_sum = temp_sum + temp;
    }

    stdev = sqrt(temp_sum / cyctimes);
}

void EvalDT(int model, int d, int cyctimes)
{
    PKE_Para pkePara;
    CZ_Para pvhssPara;
    vec_ZZ_pX pkePk, pkeSk, hssEk_1, hssEk_2, C_alpha;
    pkePk.SetLength(2);
    pkeSk.SetLength(2);
    hssEk_1.SetLength(2);
    hssEk_2.SetLength(2);
    C_alpha.SetLength(4);
    
    Setup(hssEk_1, hssEk_2, C_alpha, pvhssPara, pkePara, pkePk, pkeSk);
    ZZ_pXModulus modulus(pkePara.xN);

    int num_node = 20;
    int num_class = 20;
    vec_ZZ_pX p, c;
    p.SetLength(num_node);
    c.SetLength(num_class);
    for(int i=0;i<num_node;i++)
    {
        p[i]= random()%2;
        c[i]= 1;
    }
    Vec<vec_ZZ_pX> C_p, C_c;
    C_p.SetLength(num_node);
    C_c.SetLength(num_class);
    vec_ZZ_pX cipher;
    cipher.SetLength(4);
    for (int i = 0; i < num_node; i++) {
        Keygen(cipher, pkePara, modulus, pkePk, p[i]);
        C_p[i] = cipher;
    }
    for (int i = 0; i < num_class; i++) {
        Keygen(cipher, pkePara, modulus, pkePk, c[i]);
        C_c[i] = cipher;
    }
    ZZ_pX one;
    vec_ZZ_pX C_one;
    one=1;
    C_one.SetLength(4);
    Keygen(C_one, pkePara, modulus, pkePk, one);

    ZZ_pX y1, y2;
    ep_t g1T1;
    ep2_t g2T2;
    ep_new(g1T1);
    ep2_new(g2T2);
    auto Time=new double[cyctimes];
    double time,mean,stdev;
    cout << "******************** Server 1 Evaluating ********************" << "\n";
    for(int i=0;i<cyctimes;i++)
    {
        y1=0;
        time= GetTime();
        Eval_DT(y1,1,g1T1,g2T2,
                model,d,
                pvhssPara,
                pkePara,modulus,hssEk_1,C_alpha,C_p,C_c,C_one);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean,stdev,Time,cyctimes);
    cout << "running time Server 1: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";

    cout << "******************** Server 2 Evaluating ********************" << "\n";
    for(int i=0;i<cyctimes;i++)
    {
        y2=0;
        time = GetTime();
        Eval_DT(y2,2,g1T1,g2T2,
                model,d,
                pvhssPara,
                pkePara,modulus,hssEk_2,C_alpha,C_p,C_c,C_one);
        Time[i] = GetTime() - time;
    }
    DataProcess(mean,stdev,Time,cyctimes);
    cout << "running time Server 2: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";

    int flag=0;
    cout << "******************** Verifying ********************" << "\n";
    for(int i=0;i<cyctimes;i++)
    {   
        time = GetTime();
        Verify(y1, y2, g1T1, g2T2, pvhssPara, pkePara, modulus, flag);
        Time[i] = GetTime() - time;
    }
    if (flag==0) cout << "Verification passed!\n";
    DataProcess(mean,stdev,Time,cyctimes);
    cout << "running time Client: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";
}