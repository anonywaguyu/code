#include <iostream>
#include <helib/helib.h>
#include "Eval_DT.h"
using namespace std;
using namespace helib;
using namespace NTL;
int main(int argc, char *argv[]) {
    cout << "***********************************************************" <<"\n";
    // FHE=1: BGV
    // FHE=2: CKKS
    int FHE=1;

    // model=1: EEG Eye State Data Set
    // model=2: Bank Marketing Data Set
    int model=1;

    // the depth of decision tree
    int d=4;

    // the number of times the program is repeated
    int cyctimes=100;

    cout<<"FHE:";
    if(FHE==1)cout<<"BGV\n";
    else cout<<"CKKS\n";
    cout<<"Mode:";
    if(model==1)cout<<"EEG Eye State Data Set\n";
    else cout<<"Bank Marketing Data Set\n";

    cout<<"Depth of the decision tree: "<<d<<endl;

    Eval_DT(FHE,model,d,cyctimes);

    return 0;
}