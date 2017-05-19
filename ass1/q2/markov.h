#include <map>
#include <cstdio>
#include <cassert>
#include <iostream>
#include <cmath>
#include <vector>
using namespace std;
class MarkovWord{
public:
  double psi_ocr[1000][10],psi_trans[10][10];
  double psi_skip;
  const char *chars = "etaoinshrd";
  int ctoi(char c){
    int i=0;
    for(const char *p=chars;*p!=0;p++)
    {
      if(c == *p)
        return i;
      i++;
    }
    assert(false);
    return 10;
  }
  void read_psi(const char *ocr,const char *trans){
    FILE *fs=fopen(ocr,"r");
    assert(fs);
    int img;
    double p;
    char c,c1,c2;

    while(3==fscanf(fs,"%d\t%c\t%lf\n",&img,&c,&p)){
      psi_ocr[img][ctoi(c)] = log(p);
    }
    fclose(fs);
    fs = fopen(trans,"r");
    assert(fs);
    while(3==fscanf(fs,"%c\t%c\t%lf\n",&c1,&c2,&p)){
      psi_trans[ctoi(c1)][ctoi(c2)] = log(p);
      // cout<<c1<<" "<<c2<<" "<<p<<" "<<psi_trans[ctoi(c1)][ctoi(c2)]<<endl;
    }
    fclose(fs);
    psi_skip = log(5);
  }
  double calc_score(const vector<int> &img,const vector<int> &word, int mode){
    double S=0;
    for(int i=0;i<word.size();i++){
      S += psi_ocr[img[i]][word[i]];
      if(mode>0 && i>0)
        S += psi_trans[word[i-1]][word[i]];
      if(mode>1)
      for(int j=i+1;j<word.size();j++){
        if(img[j]==img[i] && word[i]==word[j])
          S += psi_skip;
      }
    }
    return S;
  }
  double best_score;
  vector<int> best_word;
  void init(){
    best_word = {};
    best_score = 1;
  }
  double calc_z(const vector<int> &img, int mode, const vector<int> &word = {}, int depth=0){
    if(depth == img.size()){
      double score = calc_score(img, word, mode);
      if(best_word.size()==0 || score>=best_score)
      {
        best_score = score;
        best_word = word;
      }
      return exp(score);
    }
    vector<int> new_word(word);
    new_word.push_back(0);
    double z=0;
    for(int c=0;c<10;c++){
      new_word[depth]=c;
      z += calc_z(img, mode, new_word, depth+1);
    }
    return z;
  }
  tuple<vector<char>, int, double> infer(const vector<int> &img, const vector<char> _correct_word, int mode){
    init();
    double logz = log(calc_z(img, mode));

    vector<int> correct_word;
    for(int i=0;i<img.size();i++)
    {
      correct_word.push_back(ctoi(_correct_word[i]));
    }
    double correct_logp = calc_score(img, correct_word, mode) - logz;
    int chars_correct=0;
    vector<char> predicted_word;
    for(int i=0;i<img.size();i++){
      predicted_word.push_back(chars[best_word[i]]);
      // cerr<<predicted_word[i]<<"-"<<_correct_word[i]<<" ";
      if(correct_word[i]==best_word[i])
        chars_correct++;
    }
    // cerr<<endl;
    return make_tuple(predicted_word, chars_correct, correct_logp);
  }
  // double calc_logp(const vector<int> &img,const vector<int> &word, int mode){
  //   return calc_score(img, word, mode)-log(calc_z(img, mode));
  // }
  // double calc_p(const vector<int> &img,const vector<int> &word, int mode){
  //   return exp(calc_logp(img, word, mode));
  // }
};