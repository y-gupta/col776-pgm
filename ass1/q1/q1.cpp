#include <iostream>
#include "bayes.h"
#include <iostream>
#include <cassert>
using namespace std;
int main(int argc, char **argv){
  BayesNet bn;
  if(argc<4)
  {
    cout<<"Usage: ./a.out <g|q> <n|graph_file> <k|query_file>";
    return 1;
  }
  if(argv[1][0]=='g'){
    int n=atoi(argv[2]);
    int k=atoi(argv[3]);
    assert(n>=0);
    assert(k>=0);
    bn.generate(n,k);
    bn.save("./q1a/bn.txt");
  }else if(argv[1][0]=='q'){
    bn.load(argv[2]);
    bn.query(argv[3],"./q1b/out.txt");
  }
  return 0;
}
