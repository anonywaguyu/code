#include <iostream>
#include <gmp.h>
#include <gmpxx.h>
extern "C" {
#include <relic/relic.h>
}
#include <cstring>
#include "Eval_DT.h"
using namespace std;

int main() {

    // model=1: EEG Eye State Data Set
    // model=2: Bank Marketing Data Set
    int model=1;
    // the depth of decision tree
    int d=4;
    // the number of times the program is repeated
    int cyctimes=100;

    cout<<"Mode:";
    if(model==1)cout<<"EEG Eye State Data Set\n";
    else cout<<"Bank Marketing Data Set\n";

    cout<<"Depth of the decision tree: "<<d<<endl;

    EvalDT(model,d,cyctimes);
}
