#include <vector>
#include <iostream>
#include "markov.h"
int main(int argc, char **argv){
  MarkovWord mv;
  int mode = 0; //0 -> ocr, 1 -> transition+ocr, 2 -> trans+ocr+skip
  mv.read_psi("potentials/ocr.dat","potentials/trans.dat");
  FILE *fimg=fopen("data/small/images.dat","r");
  FILE *fword = fopen("data/small/words.dat","r");
  assert(fimg && fword);
  long chars_correct=0, chars_total=0;
  long words_correct=0, words_total=0;
  double log_likelihood=0;
  char line[1025];
  while(1){
    vector<int> img;
    vector<char> word;
    if(!fgets(line, 1024, fimg))
      break;
    char *remaining = &line[0];
    int bytes_read=0,i=0;
    while(1==sscanf(remaining,"%d\t%n",&i,&bytes_read)){
      remaining+=bytes_read;
      img.push_back(i);
    }
    fgets(line, 1024, fword);
    for(i=0;i<img.size();i++){
      assert(line[i] != '0');
      word.push_back(line[i]);
    }
    auto word_chars_logp = mv.infer(img, word, mode);

    // auto predicted_word = get<0>(word_chars_logp);
    if(get<1>(word_chars_logp)==img.size())
    {
      words_correct++;
    }
    chars_correct += get<1>(word_chars_logp);
    chars_total+=img.size();
    words_total++;
    log_likelihood += get<2>(word_chars_logp);

    cout<<".";
  }
  cout<<endl;
  cout<<"Word accuracy: "<<float(words_correct)*100.f/words_total<<"% of "<<words_total<<endl;
  cout<<"Char accuracy: "<<float(chars_correct)*100.f/chars_total<<"% of "<<chars_total<<endl;
  cout<<"Avg. log likelihood: "<<log_likelihood/double(words_total)<<endl;

  return 0;
}