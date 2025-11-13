#include <gmp.h>
extern "C" {
#include <relic/relic.h>
}

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include "TSMFVC.h"
#include "time.h"
#include "sys/time.h"


void Depth2Tree(vec_ZZ_pX &res, vec_ZZ_pX p, vec_ZZ_pX lc, vec_ZZ_pX rc,
                vec_ZZ_pX tb_one, PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX ek);

void Depth2Tree_Type1(vec_ZZ_pX &res, vec_ZZ_pX p, vec_ZZ_pX lc, vec_ZZ_pX rc,
                      vec_ZZ_pX tb_one, PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX ek);

void Depth2Tree_Type2(vec_ZZ_pX &res, vec_ZZ_pX p, vec_ZZ_pX lc, vec_ZZ_pX rc,
                      vec_ZZ_pX tb_one, PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX ek);

void Depth2Tree_Type3(vec_ZZ_pX &res, vec_ZZ_pX p, vec_ZZ_pX lc, vec_ZZ_pX rc,
                      vec_ZZ_pX tb_one, PKE_Para pkePara, ZZ_pXModulus modulus, vec_ZZ_pX ek);

void Eval_DT(ZZ &tb_y, ZZ &tb_a, ZZ &tb_pi, int B, int model, int d, PKE_Para pkePara,
             ZZ_pXModulus modulus,
             vec_ZZ_pX ek, Vec<vec_ZZ_pX> b, Vec<vec_ZZ_pX> c, vec_ZZ_pX C_one, vec_ZZ_pX Ca) {
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
    //y
    ZZX temp;
    Mult(tb, pkePara, modulus, temp1, C_one);
    conv(temp,tb[0]);
    GetCoeff(tb_y, temp, 0);
    
    ZZ_pX temp12;
    ZZ_pX e1;
    GaussRand(e1, pkePara.N);
    MulMod(temp12, ek[1], Ca[1], modulus);
    temp12=e1-temp12;
    conv(temp,temp12);
    GetCoeff(tb_a, temp, 0);

    GaussRand(e1, pkePara.N);
    MulMod(temp12, tb[1], Ca[1], modulus);
    temp12=e1+temp12;
    conv(temp,temp12);
    GetCoeff(tb_pi, temp, 0);
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

void DataProcess(double &mean, double &stdev,  double *Time,int cyctimes)
{
    double temp;
    double sum=0;
    for(int i=0;i<cyctimes;i++)
    {
        sum=sum+Time[i];
    }
    mean=sum/cyctimes;

    double temp_sum=0;
    for(int i=0;i<cyctimes;i++)
    {
        temp=mean-Time[i];
        temp=temp*temp;
        temp_sum=temp_sum+temp;
    }

    stdev=sqrt(temp_sum/cyctimes);
    stdev=stdev/mean;
}

void EvalDT(int model, int d, int cyctimes)
{
    PKE_Para pkePara;
    vec_ZZ_pX pkePk, pkeSk, hssEk_1, hssEk_2;
    pkePk.SetLength(2);
    pkeSk.SetLength(2);
    hssEk_1.SetLength(2);
    hssEk_2.SetLength(2);
    
    Setup(hssEk_1, hssEk_2, pkePara, pkePk, pkeSk);
    ZZ_pXModulus modulus(pkePara.xN);

    int num_node = 20;
    int num_class = 20;
    vec_ZZ p, c;
    p.SetLength(num_node);
    c.SetLength(num_class);

    for(int i=0;i<num_node;i++)
    {
        if(random()%2==1)p[i]= conv<ZZ>("1");
        else p[i]= conv<ZZ>("0");
        c[i]= conv<ZZ>("1");
    }
    Vec<vec_ZZ_pX> C_p, C_c;
    C_p.SetLength(num_node);
    C_c.SetLength(num_class);
    vec_ZZ_pX cipher;
    cipher.SetLength(4);
    for (int i = 0; i < num_class; i++) {
        Keygen(cipher, pkePara, modulus, pkePk, c[i]);
        C_c[i] = cipher;
    }
    for (int i = 0; i < num_node; i++) {
        Probgen(cipher, pkePara, modulus, pkePk, p[i]);
        C_p[i] = cipher;
    }
    //C^1
    ZZ one=conv<ZZ>("1");
    vec_ZZ_pX C_one;
    C_one.SetLength(4);
    Keygen(C_one, pkePara, modulus, pkePk, one);
    ZZ y, y1, y2, a;
    //C^a
    vec_ZZ_pX Ca;
    Ca.SetLength(2);
    Random_ZZ_pX(Ca[1], pkePara.N, 128);

    ZZ tba1, tba2, tbpi1, tbpi2, coeff, c_1, pi, ay, pi_;
    ZZX ytemp;
    auto Time1=new double[cyctimes];
    auto Time2=new double[cyctimes];
    auto Time3=new double[cyctimes];
    double time,mean,stdev;
    int flag=0;
    struct timeval start;
    struct timeval end; 
    cout << "******************** Server 1 Evaluating ********************" << "\n";

    for(int i=0;i<cyctimes;i++)
    {
        gettimeofday(&start,NULL);
        Eval_DT(y1, tba1, tbpi1, 1, model, d, pkePara, modulus, hssEk_1, C_p,C_c,C_one,Ca);
        gettimeofday(&end,NULL);
        Time1[i]=((end.tv_sec-start.tv_sec)*1000000+(end.tv_usec-start.tv_usec));//us
    }
    DataProcess(mean,stdev,Time1,cyctimes);
    cout << "running time Server 1: " << mean/1000.0  << " ms  RSD: "<<stdev*100<<"%\n";


    cout << "******************** Server 2 Evaluating ********************" << "\n";
    for(int i=0;i<cyctimes;i++)
    {
        gettimeofday(&start,NULL);
        Eval_DT(y2, tba2, tbpi2, 2, model, d, pkePara, modulus, hssEk_2, C_p,C_c,C_one,Ca);
        gettimeofday(&end,NULL);
        Time2[i]=((end.tv_sec-start.tv_sec)*1000000+(end.tv_usec-start.tv_usec));//us
    }
    DataProcess(mean,stdev,Time2,cyctimes);
    cout << "running time Server 2: " << mean/1000.0 << " ms  RSD: "<<stdev*100<<"%\n";
    
    cout << "************************* Verifying *************************" << "\n";
    
    for(int i=0;i<cyctimes;i++)
    {
        // Eval_DT(y1, tba1, tbpi1, 1, model, d, pkePara, modulus, hssEk_1, C_p,C_c,C_one,Ca);
        // Eval_DT(y2, tba2, tbpi2, 2, model, d, pkePara, modulus, hssEk_2, C_p,C_c,C_one,Ca);
        gettimeofday(&start,NULL);
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
        //Verify(y1, y2, tba1, tba2, tbpi1, tbpi2, pkePara, flag, y, a, c_1, pi, ay, pi_);
        gettimeofday(&end,NULL);
        Time3[i]=((end.tv_sec-start.tv_sec)*1000000+(end.tv_usec-start.tv_usec));//us
    }
    if (flag==0) cout << "Verification passed!\n";
    DataProcess(mean,stdev,Time3,cyctimes);
    cout << "running time Client: " << mean/1000.0 << " ms  RSD: "<<stdev*100<<"%\n";
}