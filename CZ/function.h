#include "CZ.h"
#include "GenData.h"
#include <algorithm>
#include <stdlib.h>
using namespace std;

void Eval_Poly(int deg, int num_data) {
    //para
    PKE_Para pkePara;
    CZ_Para pvhssPara;
    vec_ZZ_pX pkePk, pkeSk, hssEk_1, hssEk_2, C_alpha;
    pkePk.SetLength(2);
    pkeSk.SetLength(2);
    hssEk_1.SetLength(2);
    hssEk_2.SetLength(2);
    C_alpha.SetLength(4);
    pkePara.msg_bit = 32;
    pkePara.d = deg;
    pkePara.num_data = num_data;
    PKE_Gen(pkePara, pkePk, pkeSk,1);
    ZZ_pXModulus modulus(pkePara.xN);
    Setup(hssEk_1, hssEk_2, C_alpha, pvhssPara, pkePara, modulus, pkePk, pkeSk);

    //data
    Data data;
    GenData(data, pkePara, modulus, pkePk);
    GenDataF(data, pkePara, pkePk);
    ZZ T1, T2;
    ep_t g1T1;
    ep2_t g2T2;
    ep_new(g1T1);
    ep2_new(g2T2);
    ZZ_pX t1y, t2y;

    double time = GetTime();
    ZZ y_ZZ;
    NativeEval(y_ZZ, pkePara.d, pkePara.num_data, data.X_decimal,data.F);
    time = GetTime() - time;
    cout << "running time native: " << time * 1000 << "ms\n";
    //eval
    cout << "******************** Server 1 Evaluating ********************" << "\n";
    time = GetTime();
    Compute(t1y,  g1T1, g2T2, 1, pvhssPara, pkePara,
               modulus, hssEk_1, C_alpha, data.C_X, data.F, data.PRF);
    time = GetTime() - time;
    cout << "running time Server 1: " << time * 1000 << "ms\n";

    cout << "******************** Server 2 Evaluating ********************" << "\n";
    time = GetTime();
    Compute(t2y, g1T1, g2T2, 2, pvhssPara, pkePara,
               modulus, hssEk_2, C_alpha, data.C_X, data.F, data.PRF);
    time = GetTime() - time;
    cout << "running time Server 2: " << time * 1000 << "ms\n";

    int flag=1;
    cout << "******************** Client Verifying ********************" << "\n";
    time = GetTime();
    Verify(t1y, t2y, g1T1, g2T2, pvhssPara, pkePara, flag);
    time = GetTime() - time;
    if ((flag==1))
    {
        cout << "Verification passed!\n";
    }
    cout << "running time Client: " << time * 1000 << "ms\n";
}

void Run_Setup(int msg_bit) {
    //para
    PKE_Para pkePara;
    CZ_Para pvhssPara;
    vec_ZZ_pX pkePk, pkeSk, hssEk_1, hssEk_2, C_alpha;
    pkePk.SetLength(2);
    pkeSk.SetLength(2);
    hssEk_1.SetLength(2);
    hssEk_2.SetLength(2);
    C_alpha.SetLength(4);
    pkePara.msg_bit = msg_bit;
    pkePara.num_data = 2;
    PKE_Gen(pkePara, pkePk, pkeSk,0);
    ZZ_pXModulus modulus(pkePara.xN);
    Setup(hssEk_1, hssEk_2, C_alpha, pvhssPara, pkePara, modulus, pkePk, pkeSk);
}

void Time_Setup(int msg_bit,int cyctimes) {
    auto *Time = new double[cyctimes];
    double time,mean,stdev;
    stdev=1;
    while(100*stdev>10)
    {
        for (int i = 0; i < cyctimes; i++) {
        time = GetTime();
        Run_Setup(msg_bit);
        Time[i] = GetTime() - time;
        }
        DataProcess(mean,stdev,Time,cyctimes);
    }
    cout << "Setup algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";
}

void Time_Probgen(int msg_bit,  int cyctimes) {
    //para
    PKE_Para pkePara;
    CZ_Para pvhssPara;
    vec_ZZ_pX pkePk, pkeSk, hssEk_1, hssEk_2, C_alpha;
    pkePk.SetLength(2);
    pkeSk.SetLength(2);
    hssEk_1.SetLength(2);
    hssEk_2.SetLength(2);
    C_alpha.SetLength(4);
    pkePara.msg_bit = msg_bit;
    pkePara.num_data = 2;
    PKE_Gen(pkePara, pkePk, pkeSk,0);
    ZZ_pXModulus modulus(pkePara.xN);
    Setup(hssEk_1, hssEk_2, C_alpha, pvhssPara, pkePara, modulus, pkePk, pkeSk);

    vec_ZZ_pX C;
    ZZ_pX x;
    C.SetLength(4);
    auto *Time = new double[cyctimes];
    double time,mean,stdev;
    stdev=1;
    while(100*stdev>10)
    {
        for (int i = 0; i < cyctimes; i++) {
        Random_ZZ_pX(x, pkePara.msg_bit, 1);
        time = GetTime();
        Enc(C, pkePara, modulus, pkePk, x);
        Time[i] = GetTime() - time;
        }
        DataProcess(mean,stdev,Time,cyctimes);
    }
    cout << "Probgen algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";
}

void Time_Ver(int msg_bit,  int cyctimes) {
    //para
    PKE_Para pkePara;
    CZ_Para pvhssPara;
    vec_ZZ_pX pkePk, pkeSk, hssEk_1, hssEk_2, C_alpha;
    pkePk.SetLength(2);
    pkeSk.SetLength(2);
    hssEk_1.SetLength(2);
    hssEk_2.SetLength(2);
    C_alpha.SetLength(4);
    pkePara.msg_bit = msg_bit;
    pkePara.d = 2;
    pkePara.num_data = 2;
    PKE_Gen(pkePara, pkePk, pkeSk,0);
    ZZ_pXModulus modulus(pkePara.xN);
    Setup(hssEk_1, hssEk_2, C_alpha, pvhssPara, pkePara, modulus, pkePk, pkeSk);
    
    //data
    Data data;
    GenData(data, pkePara, modulus, pkePk);
    GenDataF(data, pkePara, pkePk);
    ZZ T1, T2;
    ep_t g1T1;
    ep2_t g2T2;
    ep_new(g1T1);
    ep2_new(g2T2);
    ZZ_pX t1y, t2y;
    ZZ y_ZZ;
    Compute(t1y,  g1T1, g2T2, 1, pvhssPara, pkePara,
               modulus, hssEk_1, C_alpha, data.C_X, data.F, data.PRF);
    Compute(t2y, g1T1, g2T2, 2, pvhssPara, pkePara,
               modulus, hssEk_2, C_alpha, data.C_X, data.F, data.PRF);
    auto *Time = new double[cyctimes];
    double time,mean,stdev;
    stdev=1;
    int flag=1;
    while(100*stdev>10)
    {
        for (int i = 0; i < cyctimes; i++) {
            time = GetTime();
            Verify(t1y, t2y, g1T1, g2T2, pvhssPara, pkePara, flag);
            Time[i]  = GetTime() - time;
        }
        DataProcess(mean,stdev,Time,cyctimes);
    }
    if ((flag==1)) cout << "Ver algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";
}

void Time_Eval_Subalgo(int msg_bit, int cyctimes) {
    //para
    PKE_Para pkePara;
    CZ_Para pvhssPara;
    vec_ZZ_pX pkePk, pkeSk, hssEk_1, hssEk_2, C_alpha;
    pkePk.SetLength(2);
    pkeSk.SetLength(2);
    hssEk_1.SetLength(2);
    hssEk_2.SetLength(2);
    C_alpha.SetLength(4);
    pkePara.msg_bit = msg_bit;
    pkePara.num_data = 2;
    PKE_Gen(pkePara, pkePk, pkeSk,0);
    ZZ_pXModulus modulus(pkePara.xN);
    Setup(hssEk_1, hssEk_2, C_alpha, pvhssPara, pkePara, modulus, pkePk, pkeSk);
    //data
    Data data;
    GenData(data, pkePara, modulus, pkePk);

    //load
    vec_ZZ_pX db1, db2;
    db1.SetLength(2);
    db2.SetLength(2);
    auto *Time = new double[cyctimes];
    double time,mean,stdev;
    stdev=1;
    while(100*stdev>10)
    {
        for (int i = 0; i < cyctimes; i++) {
        int index = rand() % 10;
        time = GetTime();
        Mult(db1, pkePara, modulus, hssEk_1, data.C_X[rand() % pkePara.num_data]);
        db1[0] = db1[0] + data.PRF[index][0];
        db1[1] = db1[1] + data.PRF[index][1];
        Time[i] = GetTime() - time;
        }
        DataProcess(mean,stdev,Time,cyctimes);
    }
    cout << "Load algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";

    //add1
    stdev=1;
    while(100*stdev>10)
    {
        for (int i = 0; i < cyctimes; i++) {
        Random_ZZ_pX(db1[0], pkePara.N, pkePara.q_bit);
        Random_ZZ_pX(db1[1], pkePara.N, pkePara.q_bit);
        Random_ZZ_pX(db2[0], pkePara.N, pkePara.q_bit);
        Random_ZZ_pX(db2[1], pkePara.N, pkePara.q_bit);
        int index = rand() % 10;
        time = GetTime();
        db1[0] = db1[0] + db2[0] + data.PRF[index][0];
        db1[1] = db1[1] + db2[1] + data.PRF[index][1];
        Time[i] = GetTime() - time;
        }
        DataProcess(mean,stdev,Time,cyctimes);
    }
    cout << "Add1 algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";

    //add2
    vec_ZZ_pX C;
    C.SetLength(4);
    stdev=1;
    while(100*stdev>10)
    {
        for (int i = 0; i < cyctimes; i++) {
        time = GetTime();
        C[0] = data.C_X[0][0] + data.C_X[1][0];
        C[1] = data.C_X[0][1] + data.C_X[1][1];
        C[2] = data.C_X[0][2] + data.C_X[1][2];
        C[3] = data.C_X[0][3] + data.C_X[1][3];
        Time[i] = GetTime() - time;
        }
        DataProcess(mean,stdev,Time,cyctimes);
    }
    cout << "Add2 algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";

    //cmult
    ZZ constant= RandomBits_ZZ(msg_bit);
    stdev=1;
    while(100*stdev>10)
    {
        for (int i = 0; i < cyctimes; i++) {
        Mult(db1, pkePara, modulus, hssEk_1, data.C_X[rand() % pkePara.num_data]);
        int index = rand() % 10;
        time = GetTime();
        ZZ_pX_ScaleMul_ZZ(db1[0],db1[0],constant);
        ZZ_pX_ScaleMul_ZZ(db1[1],db1[1],constant);
        db1[0] = db1[0] + data.PRF[index][0];
        db1[1] = db1[1] + data.PRF[index][1];
        Time[i] = GetTime() - time;
        }
        DataProcess(mean,stdev,Time,cyctimes);
    }
    cout << "cMult algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";

    //mult
    stdev=1;
    while(100*stdev>10)
    {
        for (int i = 0; i < cyctimes; i++) {
        Mult(db1, pkePara, modulus, hssEk_1, data.C_X[rand() % pkePara.num_data]);
        int index = rand() % 10;
        time = GetTime();
        Mult(db1, pkePara, modulus, db1, data.C_X[rand() % pkePara.num_data]);
        db1[0] = db1[0] + data.PRF[index][0];
        db1[1] = db1[1] + data.PRF[index][1];
        Time[i] = GetTime() - time;
        }
        DataProcess(mean,stdev,Time,cyctimes);
    }
    cout << "Mult algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";

    //output
    int b = 1;
    ep_t g1T1;
    ep2_t g2T2;
    ep_new(g1T1);
    ep2_new(g2T2);
    ZZ_p yb;
    vec_ZZ_pX tb;
    tb.SetLength(2);
    Mult(tb, pkePara, modulus, hssEk_1, data.C_X[rand() % pkePara.num_data]);
    stdev=1;
    while(100*stdev>10)
    {
        for (int i = 0; i < cyctimes; i++) {
        time = GetTime();
        ZZX tb0;
        ZZ two_ZZ, Tb;
        ZZ_p two_ZZ_p;
        bn_t Tb_bn;
        two_ZZ_p = 2;
        bn_new(Tb_bn);
        eval(yb, tb[0], two_ZZ_p);
        Mult(tb, pkePara, modulus, tb, C_alpha);
        conv(tb0, tb[0]);
        eval_ZZX(Tb, tb0);
        if (b == 1) {
            Tb = Tb % pvhssPara.g1_order_ZZ;
            ZZ2bn(Tb_bn, Tb);
            ep_mul_gen(g1T1, Tb_bn);
        } else {
            Tb = Tb % pvhssPara.g2_order_ZZ;
            ZZ2bn(Tb_bn, Tb);
            ep2_mul_gen(g2T2, Tb_bn);
        }
        Time[i] = GetTime() - time;
        }
        DataProcess(mean,stdev,Time,cyctimes);
    }
    cout << "Output algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";
}