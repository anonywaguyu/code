#include "TSMFVC.h"
#include "GenData.h"
#include <algorithm>
#include <cstdlib>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include <vector>

using namespace std;

int total_solutions; 

void generate_leq_d(vector<vector<int>> &F,int current_var, int *current_solution, int remaining_d, int *index) {
    if (current_var == 3) {
        // 遍历所有可能的 x4 值，使 x1 + x2 + x3 + x4 <= d
        for (int i = 0; i <= remaining_d; ++i) {
            current_solution[current_var] = i;
            // 存入全局数组 F
            for (int j = 0; j < 4; ++j) {
                F[*index][j] = current_solution[j];
            }
            (*index)++;
        }
        return;
    }
    for (int i = 0; i <= remaining_d; ++i) {
        current_solution[current_var] = i;
        generate_leq_d(F, current_var + 1, current_solution, remaining_d - i, index);
    }
}

void generateF(vector<vector<int>> &F, int d) {
    total_solutions = 0;
    int current_solution[4] = {0}; // 临时存储解
    generate_leq_d(F, 0, current_solution, d, &total_solutions);
}

void Set_pX(ZZ_pX &a, int N, int q_bit) {
    ZZ_p coeff;
    for (int i = 0; i < N; i++) {
        conv(coeff, RandomBits_ZZ(q_bit)<<213);
        SetCoeff(a, i, coeff);
    }
}

void Eval_Poly(int deg, int num_data) {
    cout << "*************************** Setup -***************************" << "\n";
    int item=1,b=1;
    for(int i=1;i<=num_data;i++)
    {
        item=item*i;
        b=b*(i+deg);
    }
    item=b/item;
    vector<vector<int>> F(item, vector<int>(4, 0)); // 初始化为 0
    vector<vector<int>> F_temp(item, vector<int>(4, 0)); // 初始化为 0
    //cout<<item<<endl;
    generateF(F,deg);
    // for (int i = 0; i < item; i++)
    // {
    //     for (int j = 0; j < 4; j++)
    //     {
    //         cout<<F[i][j]<<" ";
    //     }
    //     cout<<"\n";
        
    // }
    
    //para
    PKE_Para pkePara;
    vec_ZZ_pX pkePk, pkeSk, hssEk_1, hssEk_2;
    vec_ZZ_pX C1;
    C1.SetLength(4);
    pkePk.SetLength(2);
    pkeSk.SetLength(2);
    hssEk_1.SetLength(2);
    hssEk_2.SetLength(2);
    pkePara.msg_bit = 32;
    pkePara.d = deg;
    pkePara.num_data = num_data;
    Setup(hssEk_1, hssEk_2, pkePara, pkePk, pkeSk, 1);
    ZZ_pXModulus modulus(pkePara.xN);
    //C^1
    ZZ constant1=conv<ZZ>("1");
    PKE_OKDM(C1, pkePara, modulus, pkePk, constant1);
    //data
    Data data;
    cout << "*************************** Keygen ***************************" << "\n";
    GenDataF(data, pkePara, pkePk);
    cout << "************************** Probgen ***************************" << "\n";
    GenData(data, pkePara, pkePk);

    ZZ y, y1, y2;
    //C^a
    vec_ZZ_pX Ca;
    Ca.SetLength(2);
    Random_ZZ_pX(Ca[1], pkePara.N, 128);
    ZZ ay, tba1, tba2, pi, tbpi1, tbpi2, coeff, y_native, c_1, pi_;
    ZZX ytemp;
    ZZ a=RandomBits_ZZ(128);
    for (int i = 0; i < item; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            F_temp[i][j]=F[i][j];
        }
    }
    cout << "******************** Server 1 Evaluating *********************" << "\n";
    double time = GetTime();
    Compute(y1, tba1, tbpi1, 1, pkePara, modulus, hssEk_1, data.C_X, data.C_F, data.PRF, C1, Ca, F_temp, item);
    time = GetTime() - time;
    cout << "running time Server 1: " << time * 1000 << "ms\n";

    for (int i = 0; i < item; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            F_temp[i][j]=F[i][j];
        }
    }
    cout << "******************** Server 2 Evaluating ********************" << "\n";
    time = GetTime();
    Compute(y2, tba2, tbpi2, 2, pkePara, modulus, hssEk_2, data.C_X, data.C_F, data.PRF, C1, Ca, F_temp, item);
    time = GetTime() - time; 
    cout << "running time Server 2: " << time * 1000 << "ms\n";

    cout << "********************** Client Verifying **********************" << "\n";
    time = GetTime();
    y = (y1 + y2)%pkePara.p;
    c_1=(tba1+tba2+(a<<(pkePara.lift_bit)));
    pi = (((y*c_1+tbpi1+tbpi2)*pkePara.twice_p+pkePara.q)>>(1+pkePara.q_bit))%pkePara.p;
    ay=a*y%pkePara.p;
    pi_ = ((((y*c_1+tbpi1)*pkePara.twice_p+pkePara.q)>>(1+pkePara.q_bit))+
            ((tbpi2*pkePara.twice_p+pkePara.q)>>(1+pkePara.q_bit)))%pkePara.p;
    if (pi==ay && pi==pi_) {
        cout << "Verification passed!";
    } else {
        printf("******************** ERROR ********************\n");
        cout << "The value of y   is  :" << y << endl;
        printf("\n\n");
        cout << "THE value of pi  is  :" << pi <<endl;
    }
    time = GetTime() - time;
    cout << "running time Client: " << time * 1000 << "ms\n";
    
    cout << "********************** native computing **********************" << "\n";
    time = GetTime();
    NativeEval(y_native, pkePara.d, pkePara.num_data, data.X, data.F, F, item);
    time = GetTime()-time;
    cout << "native computing time: " << time * 1000 << "ms\n";
    cout << "The value of native y:" << y_native << endl;
    cout << "The value of y   is  :" << y << endl;
    cout << "THE value of pi  is  :" << pi <<endl;
    cout << "THE value of pi_ is  :" << pi_ <<endl;
    cout << "THE value of ay  is  :" << ay <<endl;
}

void Run_Setup(int msg_bit, int num_data) {
    PKE_Para pkePara;
    vec_ZZ_pX pkePk, pkeSk, hssEk_1, hssEk_2;
    pkePk.SetLength(2);
    pkeSk.SetLength(2);
    hssEk_1.SetLength(2);
    hssEk_2.SetLength(2);
    pkePara.msg_bit = msg_bit;
    pkePara.num_data = 2;
    Setup(hssEk_1, hssEk_2, pkePara, pkePk, pkeSk, 0);
}

void Time_Setup(int msg_bit, int cyctimes) {
    auto *Time = new double[cyctimes];
    double time,mean,stdev=1;
    while(stdev*100>10)
    {
        for (int i = 0; i < cyctimes; i++) {
        time = GetTime();
        Run_Setup(msg_bit, 2);
        Time[i] = GetTime() - time;
        }
        DataProcess(mean,stdev,Time,cyctimes);
    }
    cout << "Setup algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";
}

void Time_Keygen(int msg_bit,  int cyctimes) {
    //para
    PKE_Para pkePara;
    vec_ZZ_pX pkePk, pkeSk, hssEk_1, hssEk_2;
    pkePk.SetLength(2);
    pkeSk.SetLength(2);
    hssEk_1.SetLength(2);
    hssEk_2.SetLength(2);
    pkePara.msg_bit = msg_bit;
    pkePara.num_data = 2;
    Setup(hssEk_1, hssEk_2, pkePara, pkePk, pkeSk, 0);
    ZZ_pXModulus modulus(pkePara.xN);

    vec_ZZ_pX C;
    ZZ x;
    C.SetLength(4);
    auto *Time = new double[cyctimes];
    double time,mean,stdev=1;
    while(stdev*100>10)
    {
        for (int i = 0; i < cyctimes; i++) {
        x = RandomBits_ZZ(pkePara.msg_bit);
        time = GetTime();
        Keygen(C, pkePara, modulus, pkePk, x);
        Time[i] = GetTime() - time;
        }
        DataProcess(mean,stdev,Time,cyctimes);
    }
    cout << "Keygen algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";
}

void Time_Probgen(int msg_bit,  int cyctimes) {
    //para
    PKE_Para pkePara;
    vec_ZZ_pX pkePk, pkeSk, hssEk_1, hssEk_2;
    pkePk.SetLength(2);
    pkeSk.SetLength(2);
    hssEk_1.SetLength(2);
    hssEk_2.SetLength(2);
    pkePara.msg_bit = msg_bit;
    pkePara.num_data = 2;
    Setup(hssEk_1, hssEk_2, pkePara, pkePk, pkeSk, 0);
    ZZ_pXModulus modulus(pkePara.xN);

    vec_ZZ_pX C;
    ZZ x;
    C.SetLength(4);
    auto *Time = new double[cyctimes];
    double time,mean,stdev=1;
    while(stdev*100>10)
    {
        for (int i = 0; i < cyctimes; i++) {
        x = RandomBits_ZZ(pkePara.msg_bit);
        time = GetTime();
        Probgen(C, pkePara, modulus, pkePk, x);
        Time[i] = GetTime() - time;
        }
        DataProcess(mean,stdev,Time,cyctimes);
    }
    cout << "Probgen algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";
}

void Time_Verify(int msg_bit, int cyctimes){
    //para
    PKE_Para pkePara;
    vec_ZZ_pX pkePk, pkeSk, hssEk_1, hssEk_2,tb1,tb2;
    pkePk.SetLength(2);
    pkeSk.SetLength(2);
    hssEk_1.SetLength(2);
    hssEk_2.SetLength(2);
    tb1.SetLength(2);
    tb2.SetLength(2);
    pkePara.msg_bit = msg_bit;
    pkePara.num_data = 2;
    Setup(hssEk_1, hssEk_2, pkePara, pkePk, pkeSk, 0);
    ZZ_pXModulus modulus(pkePara.xN);
    ZZX temp;
    ZZ_pX temp1;
    ZZ_pX e1;
    ZZ x,y,y1,y2,c_1,ay,pi,pi_;
    ZZ tba1,tba2,tbpi1,tbpi2,a;
    vec_ZZ_pX Ca,C;
    Ca.SetLength(2);
    C.SetLength(4);
    Random_ZZ_pX(Ca[1], pkePara.N, pkePara.q_bit);
    auto *Time = new double[cyctimes];
    double time,mean,stdev=1;
    while(stdev*100>30)
    {
        for (int i = 0; i < cyctimes; i++) {
            x = RandomBits_ZZ(pkePara.msg_bit);
            Probgen(C, pkePara, modulus, pkePk, x);
            Mult(tb1, pkePara, modulus, hssEk_1, C);
            conv(temp,tb1[0]);
            GetCoeff(y1, temp, 0);
            
            GaussRand(e1, pkePara.N);
            MulMod(temp1, hssEk_1[1], Ca[1], modulus);
            temp1=e1-temp1;
            conv(temp,temp1);
            GetCoeff(tba1, temp, 0);

            GaussRand(e1, pkePara.N);
            MulMod(temp1, tb1[1], Ca[1], modulus);
            temp1=e1+temp1;
            conv(temp,temp1);
            GetCoeff(tbpi1, temp, 0);

            Mult(tb2, pkePara, modulus, hssEk_2, C);
            conv(temp,tb2[0]);
            GetCoeff(y2, temp, 0);

            GaussRand(e1, pkePara.N);
            MulMod(temp1, hssEk_2[1], Ca[1], modulus);
            temp1=e1-temp1;
            conv(temp,temp1);
            GetCoeff(tba2, temp, 0);

            GaussRand(e1, pkePara.N);
            MulMod(temp1, tb2[1], Ca[1], modulus);
            temp1=e1+temp1;
            conv(temp,temp1);
            GetCoeff(tbpi2, temp, 0);

            time = GetTime();
            y = (y1 + y2)%pkePara.p;
            a=RandomBits_ZZ(128);
            c_1=(tba1+tba2+(a<<(pkePara.lift_bit)));
            pi = (((y*c_1+tbpi1+tbpi2)*pkePara.twice_p+pkePara.q)>>(1+pkePara.q_bit))%pkePara.p;
            ay=a*y%pkePara.p;
            pi_ = ((((y*c_1+tbpi1)*pkePara.twice_p+pkePara.q)>>(1+pkePara.q_bit))+
                ((tbpi2*pkePara.twice_p+pkePara.q)>>(1+pkePara.q_bit)))%pkePara.p;
            if (pi==ay && pi==pi_) {
            } else {
                printf("******************** ERROR ********************\n");
                cout << "The value of y   is  :" << y << endl;
                printf("\n\n");
                cout << "THE value of pi  is  :" << pi <<endl;
            }
            Time[i] = GetTime() - time;
        }
        DataProcess(mean,stdev,Time,cyctimes);
    }
    cout << "Verify algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";
}

void Time_Eval_Subalgo(int msg_bit,  int cyctimes) {
    //para
    PKE_Para pkePara;
    vec_ZZ_pX pkePk, pkeSk, hssEk_1, hssEk_2,tb1;
    pkePk.SetLength(2);
    pkeSk.SetLength(2);
    hssEk_1.SetLength(2);
    hssEk_2.SetLength(2);
    pkePara.msg_bit = msg_bit;
    pkePara.num_data = 2;
    Setup(hssEk_1, hssEk_2, pkePara, pkePk, pkeSk, 0);
    ZZ_pXModulus modulus(pkePara.xN);
    //data
    Data data;
    GenData(data, pkePara, pkePk);

    //load
    vec_ZZ_pX db1,db2;
    db1.SetLength(2);
    db2.SetLength(2);
    auto *Time = new double[cyctimes];
    double time,mean,stdev=1;
    while(stdev*100>10)
    {
        for (int i = 0; i < cyctimes; i++) {
        int index=rand()%10;
        time = GetTime();
        Mult(db1,pkePara,modulus,hssEk_1,data.C_X[rand()%pkePara.num_data]);
        db1[0]=db1[0]+data.PRF[index][0];
        db1[1]=db1[1]+data.PRF[index][1];
        Time[i] = GetTime() - time;
        }
        DataProcess(mean,stdev,Time,cyctimes);
    }
    cout << "Load algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";

    //add1
    stdev=1;
    while(stdev*100>10)
    {
        for (int i = 0; i < cyctimes; i++) {
        Random_ZZ_pX(db1[0],pkePara.N,pkePara.q_bit);
        Random_ZZ_pX(db1[1],pkePara.N,pkePara.q_bit);
        Random_ZZ_pX(db2[0],pkePara.N,pkePara.q_bit);
        Random_ZZ_pX(db2[1],pkePara.N,pkePara.q_bit);
        int index=rand()%10;
        time = GetTime();
        db1[0]=db1[0]+db2[0]+data.PRF[index][0];
        db1[1]=db1[1]+db2[1]+data.PRF[index][1];
        Time[i] = GetTime() - time;
        }
        DataProcess(mean,stdev,Time,cyctimes);
    }
    cout << "Add1 algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";

    //mult
    stdev=1;
    while(stdev*100>10)
    {
        for (int i = 0; i < cyctimes; i++) {
        Mult(db1,pkePara,modulus,hssEk_1,data.C_X[rand()%pkePara.num_data]);
        int index=rand()%10;
        time = GetTime();
        Mult(db1,pkePara,modulus,db1,data.C_X[rand()%pkePara.num_data]);
        db1[0]=db1[0]+data.PRF[index][0];
        db1[1]=db1[1]+data.PRF[index][1];
        Time[i] = GetTime() - time;
        }
        DataProcess(mean,stdev,Time,cyctimes);
    }
    cout << "Mult algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";

    //out
    ZZX temp;
    ZZ_pX temp1;
    ZZ_pX e1;
    ZZ x,y,y1,tba1,tbpi1;
    vec_ZZ_pX Ca,C;
    Ca.SetLength(2);
    C.SetLength(4);
    tb1.SetLength(2);
    Random_ZZ_pX(Ca[1], pkePara.N, pkePara.q_bit);
    stdev=1;
    while(stdev*100>10)
    {
        for (int i = 0; i < cyctimes; i++) {
            x = RandomBits_ZZ(pkePara.msg_bit);
            Probgen(C, pkePara, modulus, pkePk, x);
            
            time = GetTime();
            Output(tb1, pkePara, modulus, hssEk_1, C);
            conv(temp,tb1[0]);
            GetCoeff(y1, temp, 0);
            
            GaussRand(e1, pkePara.N);
            MulMod(temp1, tb1[1], Ca[1], modulus);
            temp1=e1+temp1;
            conv(temp,temp1);
            GetCoeff(tbpi1, temp, 0);
            Time[i] = GetTime() - time;
        }
        DataProcess(mean,stdev,Time,cyctimes);
    }
    cout << "Output algo time: " << mean * 1000 << " ms  RSD: "<<stdev*100<<"%\n";
}

