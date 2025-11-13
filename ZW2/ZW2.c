#include <stdio.h>
#include <stdlib.h>
#include <sodium.h>
#include <gmp.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "time.h"
#include "sys/time.h"

#define M 4
#define degree 4 //degree of polynomial F(x)
#define N  70//composition number C(m+d d)
#define zero 0 //
#define experiment 10
#define ttt 1
mpz_t ZERO;//0
mpz_t PRIME_Q;//q
mpz_t F_coefficient[N];
mpz_t x[M],TEMP_512,TEMP_256;//coefficients of the polynomial F(X)
int F[N][M]={0}, F_variable_degree[N][M]={0}, count=0;//degree of every variable
int total_terms = 0;
int total_solutions; 

void generate_leq_d(int current_var, int *current_solution, int remaining_d, int *index) {
    if (current_var == M - 1) {
        // 遍历所有可能的 x4 值，使 x1 + x2 + x3 + x4 <= d
        for (int i = 0; i <= remaining_d; ++i) {
            current_solution[current_var] = i;
            // 存入全局数组 F
            for (int j = 0; j < M; ++j) {
                F[*index][j] = current_solution[j];
            }
            (*index)++;
        }
        return;
    }
    for (int i = 0; i <= remaining_d; ++i) {
        current_solution[current_var] = i;
        generate_leq_d(current_var + 1, current_solution, remaining_d - i, index);
    }
}

void generateF(int d) {
    total_solutions = 0;
    int current_solution[M] = {0}; // 临时存储解
    generate_leq_d(0, current_solution, d, &total_solutions);
}

//generate F
void random_generate_F()
{
    //generate random seed
    gmp_randstate_t grt;
    clock_t time=clock();
    gmp_randinit_default(grt);
    gmp_randseed_ui(grt,time);
    //generate F
    for(int i=0;i<N;i++)
    {
        mpz_init2(F_coefficient[i],256);
        mpz_urandomb(F_coefficient[i],grt,256);//random
        mpz_tdiv_r(F_coefficient[i],F_coefficient[i],PRIME_Q);
        while(!mpz_cmp(F_coefficient[i],ZERO))
        {
            mpz_urandomb(F_coefficient[i],grt,256);//random
            mpz_tdiv_r(F_coefficient[i],F_coefficient[i],PRIME_Q);
        }
    }
    generateF(degree);
    for(int i=0;i<N;i++)
        for(int j=0;j<M;j++)
        F_variable_degree[i][j]=F[i][j];
    for(int i=0;i<M;i++)
    {
        mpz_init2(x[i],256);
        mpz_urandomb(x[i],grt,256);//random
        mpz_tdiv_r(x[i],x[i],PRIME_Q);
        while(!mpz_cmp(x[i],ZERO))
        {
            mpz_urandomb(x[i],grt,256);//random
            mpz_tdiv_r(x[i],x[i],PRIME_Q);
        }
    }
}
//mpz_t to libsoidium's string
void mpz_convert_str(mpz_t t, unsigned char*t_c)
{
    for(int i=0;i<32;i++)t_c[i]=0x00;
    int count=0;
    char Temp_C[64]={0x00};
    unsigned char temp_c[32]={0x00};
    mpz_get_str(Temp_C,16,t);
    
    for(int i=0;i<64;i++)
    {
        if(Temp_C[i]>='0'&&Temp_C[i]<='9')count++;
        else if (Temp_C[i]>='a'&&Temp_C[i]<='f')count++;
    }
    
    int flag=(count+1)/2;
    
    if(count%2==0)
    {
        for(int j=0;j<flag;j++)
        {
            int temp_left=0,temp_right=0;
            if(Temp_C[2*j]>='0'&&Temp_C[2*j]<='9')temp_left=Temp_C[2*j]-'0';
            else if(Temp_C[2*j]>='a'&&Temp_C[2*j]<='f')temp_left=Temp_C[2*j]-'a'+10;

            if(Temp_C[2*j+1]>='0'&&Temp_C[2*j+1]<='9')temp_right=Temp_C[2*j+1]-'0';
            else if(Temp_C[2*j+1]>='a'&&Temp_C[2*j+1]<='f')temp_right=Temp_C[2*j+1]-'a'+10;

            temp_c[j]=(temp_left*16+temp_right)&0xFF;
        }
        for(int j=0;j<flag;j++)
        t_c[j]=temp_c[flag-1-j];
    }
    if(count%2==1)
    {
        for(int j=0;j<flag;j++)
        {
            int temp_left=0,temp_right=0;
            if(j==0)
            {
                if(Temp_C[j]>='0'&&Temp_C[j]<='9')temp_left=Temp_C[j]-'0';
                else if(Temp_C[j]>='a'&&Temp_C[j]<='f')temp_left=Temp_C[j]-'a'+10;

                temp_c[j]=temp_left&0xFF;
            }
            else{
                if(Temp_C[2*j-1]>='0'&&Temp_C[2*j-1]<='9')temp_left=Temp_C[2*j-1]-'0';
                else if(Temp_C[2*j-1]>='a'&&Temp_C[2*j-1]<='f')temp_left=Temp_C[2*j-1]-'a'+10;

                if(Temp_C[2*j]>='0'&&Temp_C[2*j]<='9')temp_right=Temp_C[2*j]-'0';
                else if(Temp_C[2*j]>='a'&&Temp_C[2*j]<='f')temp_right=Temp_C[2*j]-'a'+10;

                temp_c[j]=(temp_left*16+temp_right)&0xFF;
            }
        }

        for(int j=0;j<flag;j++)
        {
            t_c[j]=temp_c[flag-1-j];
        }
        for(int j=flag;j<32;j++)
        {
            t_c[j]=0x00;
        }
    }
}

//ProbGen
int ProbGen(mpz_t*x, mpz_t*c,mpz_t*Inter_value, unsigned char*g_alpha)
{
    gmp_randstate_t grt;
    clock_t time=clock();
    gmp_randinit_default(grt);
    gmp_randseed_ui(grt,time);
    mpz_t random[ttt*(M+1)],alpha;   
    for(int i=0;i<ttt*(M+1);i++)
    {
        mpz_init2(random[i],256);
        mpz_urandomb(random[i],grt,256);//random
        mpz_tdiv_r(random[i],random[i],PRIME_Q);
    }
    mpz_init2(alpha,256);
    mpz_urandomb(alpha,grt,256);//random
    mpz_tdiv_r(alpha,alpha,PRIME_Q);

    for (int i = 0; i <(degree+1)*ttt+1; i++)
    {
        mpz_init2(Inter_value[i],256);
        mpz_set_str(Inter_value[i],"0",10);
        for(int j=0;j<ttt;j++)
        {
            mpz_mul_ui(TEMP_512,random[j*(M+1)],pow(i+1,j+1));
            mpz_tdiv_r(TEMP_512,TEMP_512,PRIME_Q);
            mpz_add(Inter_value[i],Inter_value[i],TEMP_512);
            mpz_tdiv_r(Inter_value[i],Inter_value[i],PRIME_Q);
        }
        mpz_add(Inter_value[i],Inter_value[i],alpha);
        mpz_tdiv_r(Inter_value[i],Inter_value[i],PRIME_Q);
    }
    //c(i)
    for (int i = 0; i <(degree+1)*ttt+1; i++)
    {
        for (int j = 0; j <M; j++)
        {
            mpz_init2(c[i*M+j],256);
            mpz_set_str(c[i*M+j],"0",10);
            for (int k = 0; k < ttt; k++)
            {
                mpz_mul_ui(TEMP_512,random[j+1+k*(M+1)],pow(i+1,k+1));
                mpz_tdiv_r(TEMP_512,TEMP_512,PRIME_Q);
                mpz_add(c[i*M+j],c[i*M+j],TEMP_512);
                mpz_tdiv_r(c[i*M+j],c[i*M+j],PRIME_Q);
            }
            mpz_add(c[i*M+j],c[i*M+j],x[j]);
            mpz_tdiv_r(c[i*M+j],c[i*M+j],PRIME_Q);
        }
    }
    //g^alpha
    unsigned char tt[crypto_core_ristretto255_SCALARBYTES];
    mpz_convert_str(alpha,tt);
    crypto_scalarmult_ristretto255_base(g_alpha, tt);
    return 0;
}

void Compute_Original(mpz_t*x)
{
    mpz_t y;
    mpz_init2(y,256);
    mpz_set_str(y,"0",10);
    for(int i=0;i<N;i++)
    {
        mpz_set(TEMP_256,F_coefficient[i]);
        for(int j=0;j<M;j++)
        {
            while (F_variable_degree[i][j]!=0)
            {
                mpz_mul(TEMP_256,TEMP_256,x[j]);
                mpz_tdiv_r(TEMP_256,TEMP_256,PRIME_Q);
                F_variable_degree[i][j]--;
            } 
        }
        mpz_add(y,TEMP_256,y);
        mpz_tdiv_r(y,y,PRIME_Q);
    }
    mpz_tdiv_r(y,y,PRIME_Q);
    for(int i=0;i<N;i++)
        for(int j=0;j<M;j++)F_variable_degree[i][j]=F[i][j];
    unsigned char Y_temp[32];
    mpz_convert_str(y,Y_temp);
    printf("original result:\n");
    for(int i=0;i<32;i++)printf("%02x ",Y_temp[i]);
    printf("\n");
    return ;
}
//F(c(i))
void Compute_partial(mpz_t*x,mpz_t y)
{
    mpz_set_str(y,"0",10);
    for(int i=0;i<N;i++)
    {
        mpz_set(TEMP_256,F_coefficient[i]);
        for(int j=0;j<M;j++)
        {
            while (F_variable_degree[i][j]!=0)
            {
                mpz_mul(TEMP_256,TEMP_256,x[j]);
                mpz_tdiv_r(TEMP_256,TEMP_256,PRIME_Q);
                F_variable_degree[i][j]--;
            } 
        }
        mpz_add(y,TEMP_256,y);
        mpz_tdiv_r(y,y,PRIME_Q);
    }
    mpz_tdiv_r(y,y,PRIME_Q);
    for(int i=0;i<N;i++)
        for(int j=0;j<M;j++)F_variable_degree[i][j]=F[i][j];
    return ;
}
//Compute
void Compute(mpz_t*c,mpz_t*Inter_value,mpz_t*v,mpz_t*w)
{
    mpz_t temp[M];
    for (int i = 0; i <M; i++)mpz_init2(temp[i],256);
    for(int i=0;i<((degree+1)*ttt+1);i++)
    {
        Compute_partial(c+i*M,v[i]);
        mpz_mul(TEMP_512,Inter_value[i],v[i]);
        mpz_tdiv_r(w[i],TEMP_512,PRIME_Q);
    }
}
//verify
int verify(mpz_t*v,mpz_t*w, unsigned char*g_alpha)
{
    mpz_t point[(degree+1)*ttt+1];
    for(int i=0;i<(degree+1)*ttt+1;i++)
    {
        mpz_init2(point[i],256);
        mpz_set_ui(point[i],(unsigned long)(i+1));
    }
    mpz_t temp,y_temp;
    mpz_init2(temp,256);
    mpz_init2(y_temp,256);
    mpz_set_ui(y_temp,0);
    for(int i=0;i<(degree+1)*ttt;i++)
    {
        mpz_set(temp,v[i]);
        for(int j=0;j<(degree+1)*ttt;j++)
        {
            if(j==i)continue;
            else
            {
                mpz_sub(TEMP_256,PRIME_Q,point[j]);
                mpz_tdiv_r(TEMP_256,TEMP_256,PRIME_Q);
                mpz_mul(TEMP_512,TEMP_256,temp);
                mpz_tdiv_r(temp,TEMP_512,PRIME_Q);

                if(j<i)
                {
                    mpz_set_ui(TEMP_256,(unsigned long)(i-j));
                }
                else{
                    mpz_set_ui(TEMP_256,(unsigned long)(j-i));
                    mpz_sub(TEMP_256,PRIME_Q,TEMP_256);
                    mpz_tdiv_r(TEMP_256,TEMP_256,PRIME_Q);
                }
                mpz_invert(TEMP_256,TEMP_256,PRIME_Q);

                mpz_mul(TEMP_512,TEMP_256,temp);
                mpz_tdiv_r(temp,TEMP_512,PRIME_Q);
            }
        }
        mpz_add(TEMP_512,y_temp,temp);
        mpz_tdiv_r(y_temp,TEMP_512,PRIME_Q);
    }
    unsigned char Y_temp[32];
    mpz_convert_str(y_temp,Y_temp);
    printf("reconstruct_value:\n");
    for(int i=0;i<32;i++)printf("%02x ",Y_temp[i]);
    printf("\n");
    unsigned char g_1[crypto_core_ristretto255_BYTES];
    if (crypto_scalarmult_ristretto255(g_1, Y_temp, g_alpha)!=0)
        return -1;
    mpz_set_ui(y_temp,0);
    for(int i=0;i<(degree+1)*ttt+1;i++)
    {
        mpz_set(temp,w[i]);
        for(int j=0;j<(degree+1)*ttt+1;j++)
        {
            if(j==i)continue;
            else
            {
                mpz_sub(TEMP_256,PRIME_Q,point[j]);
                mpz_tdiv_r(TEMP_256,TEMP_256,PRIME_Q);
                mpz_mul(TEMP_512,TEMP_256,temp);
                mpz_tdiv_r(temp,TEMP_512,PRIME_Q);

                if(j<i) mpz_set_ui(TEMP_256,(unsigned long)(i-j));
                else
                {
                    mpz_set_ui(TEMP_256,(unsigned long)(j-i));
                    mpz_sub(TEMP_256,PRIME_Q,TEMP_256);
                    mpz_tdiv_r(TEMP_256,TEMP_256,PRIME_Q);
                }
                mpz_invert(TEMP_256,TEMP_256,PRIME_Q);

                mpz_mul(TEMP_512,TEMP_256,temp);
                mpz_tdiv_r(temp,TEMP_512,PRIME_Q);
            }
        }
        mpz_add(TEMP_512,y_temp,temp);
        mpz_tdiv_r(y_temp,TEMP_512,PRIME_Q);
    }
    mpz_convert_str(y_temp,Y_temp);
    unsigned char g_2[crypto_core_ristretto255_BYTES];
    crypto_scalarmult_ristretto255_base(g_2, Y_temp);
    if (memcmp(g_1, g_2, crypto_core_ristretto255_BYTES) == 0) 
        return 1;
    else 
    {
        return -1;
    }
}

int main()
{
    mpz_init2(ZERO,256);
    mpz_set_str(ZERO,"0",10);
    mpz_init2(PRIME_Q,256);
    mpz_init_set_str(PRIME_Q,"7237005577332262213973186563042994240857116359379907606001950938285454250989",10);

    unsigned char g_alpha[crypto_core_ristretto255_SCALARBYTES];
    double time_original_use=0,time_probgen_use=0,time_compute_use=0,time_verify_use=0;
    struct timeval start;
    struct timeval end; 

    mpz_t v[(degree+1)*ttt+1],w[(degree+1)*ttt+1];
    mpz_t Inter_value[(degree+1)*ttt+1],c[((degree+1)*ttt+1)*M];
    for (int i = 0; i <(degree+1)*ttt+1; i++)
        for (int j = 0; j <M; j++)
            mpz_init2(c[i*M+j],256);
    for (int i = 0; i <(degree+1)*ttt+1; i++)
    {
        mpz_init2(Inter_value[i],256);
        mpz_init2(v[i],256);
        mpz_init2(w[i],256);
    }
    printf("M:%d Degree:%d Items:%d \n",M,degree,N);
    printf("start:\n");
    for(int i=0;i<experiment;i++)
    {
        printf("第%d/%d次:\n", i+1,experiment);
        gmp_randstate_t grt;
        clock_t time=clock();
        gmp_randinit_default(grt);
        gmp_randseed_ui(grt,time);
        random_generate_F();

        gettimeofday(&start,NULL);
        Compute_Original(x);
        gettimeofday(&end,NULL);
        time_original_use+=((end.tv_sec-start.tv_sec)*1000000+(end.tv_usec-start.tv_usec));//us
        //printf("Original Compute is ok!\n");

        gettimeofday(&start,NULL);
        ProbGen(x,c,Inter_value,g_alpha);
        gettimeofday(&end,NULL);
        time_probgen_use+=((end.tv_sec-start.tv_sec)*1000000+(end.tv_usec-start.tv_usec));//us
        printf("ProbGen is ok!\n");

        gettimeofday(&start,NULL);
        Compute(c,Inter_value,v,w);
        gettimeofday(&end,NULL);
        time_compute_use+=((end.tv_sec-start.tv_sec)*1000000+(end.tv_usec-start.tv_usec));//us
        printf("Compute is ok!\n");

        gettimeofday(&start,NULL);
        if(verify(v,w,g_alpha))printf("Verify is passed!\n");
        else printf("Verify not passed!\n");
        gettimeofday(&end,NULL);
        time_verify_use+=((end.tv_sec-start.tv_sec)*1000000+(end.tv_usec-start.tv_usec));//us
    }
    time_original_use/=1000;//ms
    time_probgen_use/=1000;//ms
    time_compute_use/=1000;//ms
    time_verify_use/=1000;//ms

    printf("M:%d Degree:%d Items:%d \n",M,degree,N);
    printf("\ntotal_time_compute_use:%lfms\n",ttt*time_compute_use/experiment);
    printf("\nClien_time_use:%lfms\n",(time_probgen_use+time_verify_use)/experiment);
    return 0;
}



