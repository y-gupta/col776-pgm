#include <vector>
#include <iostream>

#include "general_net.h"
#define DATA "insurance"
int main(int argc, char **argv){
  GeneralNet gnet;
  gnet.use_markov = true;
  cout<<"Loading...\n";
  gnet.load("data2/" DATA ".bif");
  cout<<"Training...\n";
  gnet.train("data2/" DATA "_small.dat");
  cout<<"Testing...\n";
  gnet.test("data2/" DATA "_test.dat","data2/" DATA "_TrueValues.dat");
  cout<<"Done!\n";
  return 0;
}
