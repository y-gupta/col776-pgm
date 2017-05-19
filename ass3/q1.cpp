#include <vector>
#include <iostream>

#include "word_net.h"

int main(int argc, char **argv){
  WordNet wnet;
  wnet.mode = 3;
  if(argc>=3){
    wnet.mode = argv[1][0]-'0';
  }
  wnet.read_phi("potentials/ocr.dat","potentials/trans.dat");
  // wnet.test_inference();
  wnet.process_word_models("data/data-loops.dat", "data/truth-loops.dat");
  return 0;
}
