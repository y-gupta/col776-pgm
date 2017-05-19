#include <vector>
#include <iostream>

#define MAP_INFERENCE 1

#include "word_net.h"

int main(int argc, char **argv){
  WordNet wnet;
  wnet.mode = 3;
  wnet.loopy = 1;
  if(argc>=3){
    wnet.mode = argv[1][0]-'0';
    wnet.loopy = argv[2][0]-'0';
  }
  wnet.read_phi("potentials/ocr.dat","potentials/trans.dat");
  // wnet.test_inference();
  wnet.process_word_models("data/data-loopsWS.dat", "data/truth-loopsWS.dat");
  return 0;
}
