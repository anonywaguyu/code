#include <iostream>
#include <gmp.h>
#include <gmpxx.h>
extern "C" {
#include <relic/relic.h>
}
#include <cstring>
#include "function.h"
using namespace std;

int main() {
    cout << "**************************Benchmarks**************************" <<"\n";
    // the size of the input (i.e. the value of Log(B) in the paper),
    // and it set to 1,16,32,64,128,256,512,1024 in the experiments in our paper.
    int msg_bit=1;
    // the number of times the program is repeated
    int cyctimes=10;
    // test the running time of generation algorithm
    Time_Setup(msg_bit,  cyctimes);
    Time_Keygen(msg_bit,  cyctimes);
    Time_Probgen(msg_bit,  cyctimes);
    // test the running time of each subroutine in evaluation algorithm
    Time_Eval_Subalgo(msg_bit,  cyctimes);
    Time_Verify(msg_bit,  cyctimes);

    cout << "**************Multivariate polynomial evaluation**************" <<"\n";
    // the degree of the polynomial, and it set to 2,4,...,20 in the experiments in our paper.
    int deg=2;
    // the number of the variables, and it set to 4 in the experiments in our paper.
    int num_data=4;
    // evaluating a polynomial with the input size 32bit
    Eval_Poly(deg, num_data);
}
