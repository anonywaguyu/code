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
#define N  70 //composition number C(m+d d)
#define zero 0 //
#define experiment 10
#define ttt 1

mpz_t ZERO;//0
mpz_t PRIME_Q;//q
mpz_t F_coefficient[N],f_shares[(degree+1)];
mpz_t x[M], TEMP_512,TEMP_256,alpha;//coefficients of the polynomial F(X)
int F[N][M]={0},F_variable_degree[N][M]={0}, count=0;//degree of every variable
gmp_randstate_t grt;

int total_solutions;

void generate_leq_d(int current_var, int *current_solution, int remaining_d, int *index) {
    if (current_var == M - 1) {
        for (int i = 0; i <= remaining_d; ++i) {
            current_solution[current_var] = i;
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
    int current_solution[M] = {0};
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
//1-degree multiplication
void multiplication(mpz_t*FFF, mpz_t f_0, mpz_t f_1)
{
    mpz_t temp_multiplication_1[degree+1];
    mpz_t temp_multiplication_2[degree+1];
    for(int i=0;i<degree+1;i++)
    {
        mpz_init2(temp_multiplication_1[i],512);
        mpz_init2(temp_multiplication_2[i],512);
    }
    
    for(int i=0;i<degree;i++)
    {
        if(!mpz_cmp(ZERO,FFF[i]));
        else{
            mpz_mul(temp_multiplication_1[i],FFF[i],f_0);
            mpz_tdiv_r(temp_multiplication_1[i],temp_multiplication_1[i],PRIME_Q);
            mpz_mul(temp_multiplication_2[i+1],FFF[i],f_1);
            mpz_tdiv_r(temp_multiplication_2[i+1],temp_multiplication_2[i+1],PRIME_Q);
        }
    }
    for(int i=0;i<degree+1;i++)
    {
        mpz_add(TEMP_512,temp_multiplication_1[i],temp_multiplication_2[i]);
        mpz_tdiv_r(FFF[i],TEMP_512,PRIME_Q);
    }
}

//KeyGen
void KeyGen(mpz_t*l_0, mpz_t*l_1, mpz_t*f_shares)
{
    //initialize the polynomial l(u)
    for(int i=0;i<M;i++)
    {
        mpz_urandomb(l_0[i],grt,256);//random
        mpz_tdiv_r(l_0[i],l_0[i],PRIME_Q);
        mpz_urandomb(l_1[i],grt,256);//random
        mpz_tdiv_r(l_1[i],l_1[i],PRIME_Q);
    }
    
    for (int j = 0; j <(degree+1); j++)
        mpz_set_str(f_shares[j],"0",10);
    
    mpz_t temp_keyegn[degree+1];
    for(int i=0;i<(degree+1);i++)mpz_init2(temp_keyegn[i],256);

    
    for(int i=0;i<N;i++)
    {
        //set 0 to temp
        for(int j=0;j<degree+1;j++)mpz_set_str(temp_keyegn[j],"0",10);
        //compute every item
        if(!mpz_cmp(F_coefficient[i],ZERO));
        else {
            mpz_set(temp_keyegn[0],F_coefficient[i]);
            for(int j=0;j<M;j++)
            {
                //x_j^
                while(F_variable_degree[i][j]!=0)
                {
                    multiplication(temp_keyegn,l_0[j],l_1[j]);
                        F_variable_degree[i][j]-=1;
                }
            }
        }
        //add every item
        for(int j=0;j<degree+1;j++)
        {
            mpz_add(f_shares[j],temp_keyegn[j],f_shares[j]);
            mpz_tdiv_r(f_shares[j],f_shares[j],PRIME_Q);
        }
            
    }
    
    for(int i=0;i<degree+1;i++)mpz_tdiv_r(f_shares[i],f_shares[i],PRIME_Q);
    for(int i=0;i<N;i++)
        for(int j=0;j<M;j++)
            F_variable_degree[i][j]=F[i][j];
}

//ProbGen
int ProbGen(mpz_t*l_0, mpz_t*l_1,mpz_t*x, mpz_t*c, unsigned char**g_alpha, unsigned char**g_Alpha)
{
    mpz_t random[(ttt+1)*M];
    for(int i=0;i<(ttt+1)*M;i++)
    {
        mpz_init2(random[i],256);
    }
    for(int i=0;i<ttt*M;i++)
    {
        mpz_urandomb(random[i],grt,256);//random
        mpz_tdiv_r(random[i],random[i],PRIME_Q);
    }
    mpz_t b;   
    mpz_init2(b,256);
    mpz_urandomb(b,grt,256);//random
    mpz_tdiv_r(b,b,PRIME_Q);

    mpz_init2(alpha,256);
    mpz_urandomb(alpha,grt,256);//random
    mpz_tdiv_r(alpha,alpha,PRIME_Q);

    mpz_t alpha_vector[degree*(ttt+1)],b_vector[degree];
    for(int i=0;i<degree*(ttt+1);i++)mpz_init2(alpha_vector[i],256);
    for(int i=0;i<degree;i++)mpz_init2(b_vector[i],256);

    mpz_set(alpha_vector[0],alpha);
    for(int i=1;i<degree*(ttt+1);i++)
    {
        mpz_mul(alpha_vector[i],alpha_vector[i-1],alpha);
        mpz_tdiv_r(alpha_vector[i],alpha_vector[i],PRIME_Q);
    }
    mpz_set(b_vector[0],b);
    for(int i=1;i<degree;i++)
    {
        mpz_mul(b_vector[i],b_vector[i-1],b);
        mpz_tdiv_r(b_vector[i],b_vector[i],PRIME_Q);
    }
    
    //generate B=l(b)
    mpz_t B_vector[M];
    for(int i=0;i<M;i++)
    {
        mpz_init2(B_vector[i],256);
        mpz_mul(TEMP_512,b,l_1[i]);
        mpz_tdiv_r(TEMP_512,TEMP_512,PRIME_Q);
        mpz_add(TEMP_512,TEMP_512,l_0[i]);
        mpz_tdiv_r(B_vector[i],TEMP_512,PRIME_Q);
    }
    ////////////
    mpz_t temp[ttt*M];
    for(int i=0;i<ttt*M;i++){
        mpz_init2(temp[i],256);
        mpz_set_str(temp[i],"0",10);
    }
    for (int s = 0; s < ttt; s++)
    {
        for (int j = 0; j < M; j++)
        {
            mpz_mul(temp[s*M+j],random[s*M+j],alpha_vector[s]);
            mpz_tdiv_r(temp[s*M+j],temp[s*M+j],PRIME_Q);
        }
    }
    for (int i = 0; i < M; i++)
    {
        mpz_set_str(TEMP_512,"0",10);
        for (int j = 0; j < ttt; j++)
        {
            mpz_add(TEMP_512,TEMP_512,temp[j*M+i]);\
            mpz_tdiv_r(TEMP_512,TEMP_512,PRIME_Q);
        }
        mpz_tdiv_r(temp[i],TEMP_512,PRIME_Q);
        mpz_add(temp[i],temp[i],x[i]);
        mpz_tdiv_r(temp[i],temp[i],PRIME_Q);
    }
    for(int i=0;i<M;i++)
    {
        if(mpz_cmp(B_vector[i],temp[i]))
        {
            mpz_sub(temp[i],B_vector[i],temp[i]);
            mpz_tdiv_r(temp[i],temp[i],PRIME_Q);
        }
        else 
        {
            mpz_add(B_vector[i],B_vector[i],PRIME_Q);
            mpz_tdiv_r(B_vector[i],B_vector[i],PRIME_Q);
            mpz_sub(temp[i],B_vector[i],temp[i]);
            mpz_tdiv_r(temp[i],temp[i],PRIME_Q);
        }
        mpz_invert(TEMP_256,alpha_vector[ttt],PRIME_Q);
        mpz_mul(temp[i],temp[i],TEMP_256);
        mpz_tdiv_r(random[i+ttt*M],temp[i],PRIME_Q);
    }
    //c(i)
    for (int i = 0; i <(degree)*(ttt+1)+1; i++)
    {
        for (int j = 0; j <M; j++)
        {
            mpz_init2(c[i*M+j],256);
            mpz_set_str(c[i*M+j],"0",10);
            for (int k = 0; k < ttt+1; k++)
            {
                mpz_mul_ui(TEMP_512,random[j+k*(M)],pow(i+1,k+1));
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
    for (int i = 0; i < degree; i++)
    {
        mpz_convert_str(b_vector[i],tt);
        crypto_scalarmult_ristretto255_base(g_Alpha[i], tt);
    }
    for (int i = 0; i < degree*(ttt+1); i++)
    {
        mpz_convert_str(alpha_vector[i],tt);
        crypto_scalarmult_ristretto255_base(g_alpha[i], tt);
    }
    return 0;
}

void Compute_Original(mpz_t*x)
{
    for(int i=0;i<N;i++)
        for(int j=0;j<M;j++)F_variable_degree[i][j]=F[i][j];
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
void Compute(mpz_t*c,mpz_t*v)
{
    mpz_t temp[M];
    for (int i = 0; i <M; i++)mpz_init2(temp[i],256);
    for(int i=0;i<((degree)*(ttt+1)+1);i++)
    {
        Compute_partial(c+i*M,v[i]);
    }
}
//verify
int verify(mpz_t*v, mpz_t*f_shares,unsigned char**g_alpha,unsigned char**g_Alpha)
{
    mpz_t point[(degree)*(ttt+1)+1];
    for(int i=0;i<(degree)*(ttt+1)+1;i++)
    {
        mpz_init2(point[i],256);
        mpz_set_ui(point[i],(unsigned long)(i+1));
    }
    unsigned char tt[crypto_core_ristretto255_SCALARBYTES];
    unsigned char Temp[degree+1][crypto_core_ristretto255_SCALARBYTES];
    mpz_convert_str(f_shares[0],tt);
    crypto_scalarmult_ristretto255_base(Temp[0], tt);
    for(int i=0;i<degree;i++)
    {
        mpz_convert_str(f_shares[i+1],tt);
        if (crypto_scalarmult_ristretto255(Temp[i+1], tt, g_Alpha[i])!=0)return -1;
    }
    unsigned char g_f_b[crypto_core_ristretto255_SCALARBYTES];
    for(int i=0;i<32;i++)g_f_b[i]=0x00;
    for(int i=0;i<degree+1;i++)crypto_core_ristretto255_add(g_f_b,g_f_b,Temp[i]);
    
    mpz_t temp,y_temp;
    mpz_init2(temp,256);
    mpz_init2(y_temp,256);
    mpz_set_ui(y_temp,0);
    unsigned char Y_temp[32];
    for(int i=0;i<(degree)*(ttt+1)+1;i++)
    {
        mpz_set(temp,v[i]);
        for(int j=0;j<(degree)*(ttt+1)+1;j++)
        {
            if(j==i)continue;
            else
            {
                //
                if(mpz_cmp(alpha,point[j]))
                {
                    mpz_sub(TEMP_256,alpha,point[j]);
                    mpz_tdiv_r(TEMP_256,TEMP_256,PRIME_Q);
                }
                else{
                    mpz_add(TEMP_512,alpha,PRIME_Q);
                    mpz_tdiv_r(TEMP_512,TEMP_512,PRIME_Q);
                    mpz_sub(TEMP_256,TEMP_512,point[j]);
                    mpz_tdiv_r(TEMP_256,TEMP_256,PRIME_Q);
                }
                
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
    mpz_convert_str(y_temp,Y_temp);

    crypto_scalarmult_ristretto255_base(Temp[0], Y_temp);
    
    if(memcmp(Temp[0],g_f_b,32)!=0) return -1;

    mpz_set_ui(y_temp,0);
    for(int i=0;i<(degree)*(ttt+1)+1;i++)
    {
        mpz_set(temp,v[i]);
        for(int j=0;j<(degree)*(ttt+1)+1;j++)
        {
            if(j==i)continue;
            else
            {
                //
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
    mpz_convert_str(y_temp,Y_temp);
    printf("reconstruct_value:\n");
    for(int i=0;i<32;i++)printf("%02x ",Y_temp[i]);
    printf("\n");
    return 1;
}

int main()
{
    mpz_init2(ZERO,256);
    mpz_set_str(ZERO,"0",10);
    mpz_init2(PRIME_Q,256);
    mpz_init_set_str(PRIME_Q,"7237005577332262213973186563042994240857116359379907606001950938285454250989",10);

    unsigned char**g_alpha=(unsigned char**)malloc((degree*(ttt+1))*sizeof(unsigned char*));
    for(int i=0;i<degree*(ttt+1);i++)
        g_alpha[i]=(unsigned char*)malloc(crypto_core_ristretto255_BYTES*sizeof(unsigned char));
    unsigned char**g_Alpha=(unsigned char**)malloc((degree)*sizeof(unsigned char*));
    for(int i=0;i<degree;i++)
        g_Alpha[i]=(unsigned char*)malloc(crypto_core_ristretto255_BYTES*sizeof(unsigned char));
    double time_keygen=0,time_original_use=0,time_probgen_use=0,time_compute_use=0,time_verify_use=0;

    struct timeval start;
    struct timeval end; 
    mpz_t l_0[M],l_1[M];
    mpz_t v[(degree)*(ttt+1)+1];
    mpz_t c[(degree*(ttt+1)+1)*M];
    for (int i = 0; i <M; i++)
    {
        mpz_init2(l_0[i],256);
        mpz_init2(l_1[i],256);
    }
    for (int i = 0; i <degree+1; i++)
    {
        mpz_init2(f_shares[i],256);
    }
    for (int i = 0; i <degree*(ttt+1)+1; i++)
        for (int j = 0; j <M; j++)
            mpz_init2(c[i*M+j],256);
    for (int i = 0; i <degree*(ttt+1)+1; i++)
        mpz_init2(v[i],256);
    
    //random seed
    clock_t time=clock();
    gmp_randinit_default(grt);
    gmp_randseed_ui(grt,time);
    printf("M:%d Degree:%d Items:%d \n",M,degree,N);
    
    for(int i=0;i<experiment;i++)
    {
        printf("%d/%d:\n",i+1,experiment);
        random_generate_F();
        gettimeofday(&start,NULL);
        Compute_Original(x);
        gettimeofday(&end,NULL);
        time_original_use+=((end.tv_sec-start.tv_sec)*1000000+(end.tv_usec-start.tv_usec));//us
        //printf("Original Compute is ok!\n");

        //keygen
        gettimeofday(&start,NULL);
        KeyGen(l_0,l_1,f_shares);
        gettimeofday(&end,NULL);
        time_keygen +=((end.tv_sec-start.tv_sec)*1000000+(end. tv_usec-start. tv_usec));//us
        printf("Keygen is ok!\n");

        gettimeofday(&start,NULL);
        ProbGen(l_0,l_1,x,c,g_alpha,g_Alpha);
        gettimeofday(&end,NULL);
        time_probgen_use+=((end.tv_sec-start.tv_sec)*1000000+(end.tv_usec-start.tv_usec));//us
        printf("ProbGen is ok!\n");

        gettimeofday(&start,NULL);
        Compute(c,v);
        gettimeofday(&end,NULL);
        time_compute_use+=((end.tv_sec-start.tv_sec)*1000000+(end.tv_usec-start.tv_usec));//us
        printf("Compute is ok!\n");

        gettimeofday(&start,NULL);
        if(verify(v,f_shares,g_alpha,g_Alpha))printf("Verify is passed!\n");
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



