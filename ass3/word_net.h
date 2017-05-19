#pragma once

#include <map>
#include <cstdio>
#include <cassert>
#include <iostream>
#include <cmath>
#include <vector>
#include <cstring>

#include "factor_table.h"
#include "markov_net.h"
#include "gibbs_infer.h"
#include "timer.h"

using namespace std;

class WordNet{
public:
  FactorTable phi_ocr, phi_trans, phi_skip, phi_pair_skip;
  // double phi_ocr[1000][10],phi_trans[10][10];
  // double phi_skip;
  // double phi_pair_skip;
  vector<FactorTable> phi_ocr_reduced;
  Timer timer;
  int mode;

  const char *chars = "etaoinshrd";
  map<char, int> char_map;
  WordNet(){
    for(int i=0;i < strlen(chars); i++){ // initialize char-int mapping.
      char_map[chars[i]]=i;
    }

    mode = 3;//0: ocr, 1: +trans, 2: +skip, 3: +pair_skip
  }
  void read_phi(const char *ocr,const char *trans){
    FILE *fs;
    fs=fopen(ocr,"r");
    assert(fs);
    int img;
    double p;
    char c,c1,c2;
    phi_ocr.init({RandomVariable(0,1000), RandomVariable(1,10)});
    while(3==fscanf(fs,"%d\t%c\t%lf\n",&img,&c,&p)){
      phi_ocr.set({img, char_map[c]}, log(p));
      // phi_ocr[img][char_map[c]] = log(p);
    }
    fclose(fs);
    fs=fopen(trans,"r");
    assert(fs);
    phi_trans.init({RandomVariable(0,10), RandomVariable(1,10)});
    while(3==fscanf(fs,"%c\t%c\t%lf\n",&c1,&c2,&p)){
      phi_trans.set({char_map[c1], char_map[c2]}, log(p));
      // phi_trans[char_map[c1]][char_map[c2]] = log(p);
    }
    fclose(fs);
    phi_skip.init({RandomVariable(0,10), RandomVariable(1,10)});
    for(int i=0;i<10;i++){
      for(int j=0;j<10;j++){
        phi_skip.set({i,j}, (i==j)?log(5):log(1));
      }
    }
    phi_pair_skip = phi_skip;
    // phi_skip = log(5);
    // phi_pair_skip = log(5);
    phi_ocr_reduced.clear();
    for(int i=0;i<1000;i++){
      phi_ocr_reduced.push_back(phi_ocr.reduce({{0,i}}));
    }
  }
  vector<int> read_img(char *line){
    vector<int> img;
    char *remaining = &line[0];
    int bytes_read=0,i=0;
    while(1==sscanf(remaining,"%d\t%n",&i,&bytes_read)){
      remaining+=bytes_read;
      img.push_back(i);
    }
    return img;
  }
  vector<char> read_word(char *line){
    vector<char> word;
    while(*line != 0 && *line!='\n'){
      word.push_back(*line);
      line++;
    }
    return word;
  }
  void process_word_models(const char *fname_img, const char *fname_word){
    FILE *fimg=fopen(fname_img,"r");
    FILE *fword=fopen(fname_word,"r");
    assert(fimg && fword);
    long chars_correct=0, chars_total=0;
    long words_correct=0, words_total=0;
    double log_likelihood=0;
    char line[1025];
    double total_time = 0;
    while(1){
      vector<int> img1,img2;
      vector<char> word1,word2;
      if(!fgets(line, 1024, fimg))
        break;
      img1 = read_img(line);
      fgets(line, 1024, fword);
      word1 = read_word(line);
      assert(img1.size() == word1.size());

      fgets(line, 1024, fimg);
      img2 = read_img(line);
      fgets(line, 1024, fword);
      word2 = read_word(line);
      assert(img2.size() == word2.size());

      fgets(line,1024, fimg); //skip empty line b/w pairs
      fgets(line,1024, fword);

      auto net = make_markov(img1,img2);

      auto words = word1;
      words.insert(words.end(),word2.begin(), word2.end()); //concat word1,word2
      timer.reset();
      auto words_logp = infer(net, words); //inferred_word1 and word2 (concat'd), #logp for given words
      total_time += timer.elapsed()/1000.0; //change ms to s.
      auto predicted_words = words_logp.first;

      bool all_correct = true;
      for(int i=0;i<img1.size();i++)
      {
        if(predicted_words[i] == word1[i])
          chars_correct++;
        else
          all_correct = false;
      }
      if(all_correct)
        words_correct++;

      all_correct = true;
      for(int i=0;i<img2.size();i++)
      {
        if(predicted_words[i+img1.size()] == word2[i])
          chars_correct++;
        else
          all_correct = false;
      }
      if(all_correct)
        words_correct++;

      chars_total += img1.size() + img2.size();
      words_total += 2;
      log_likelihood += words_logp.second;

      cerr<<".";
    }

    cout<<endl;
    cout<<"Word accuracy: "<<float(words_correct)*100.f/words_total<<"% of "<<words_total<<endl;
    cout<<"Char accuracy: "<<float(chars_correct)*100.f/chars_total<<"% of "<<chars_total<<endl;
    cout<<"Avg. log likelihood: "<<log_likelihood/double(words_total)<<endl;
    cout<<"Time taken (inference only): "<<total_time<<endl;
  }
  MarkovNet make_markov(const vector<int> &img1, const vector<int> &img2){
    MarkovNet mnet;
    for(int i=0;i<img1.size();i++){
      auto var = RandomVariable(i,10);

      auto ocr_factor = phi_ocr_reduced[img1[i]];
      mnet.add_var(var);
      ocr_factor.remap_vars({{RandomVariable(1), var}});
      mnet.add_factor(ocr_factor);

      if(mode >= 1 && i+1<img1.size()){
        auto trans_factor = phi_trans;
        auto var2 = RandomVariable(i+1, 10);
        trans_factor.remap_vars({{RandomVariable(0), var}, {RandomVariable(1), var2}});
        mnet.add_factor(trans_factor);
      }

      if(mode >= 2)
      for(int j=i+1;j<img1.size();j++){
        if(img1[i]==img1[j]){
          auto skip_factor = phi_skip;
          auto var2= RandomVariable(j, 10);
          skip_factor.remap_vars({{RandomVariable(0), var}, {RandomVariable(1), var2}});
          mnet.add_factor(skip_factor);
        }
      }

      if(mode >= 3)
      for(int j=0;j<img2.size();j++){
        if(img1[i]==img2[j]){
          auto pair_skip_factor = phi_pair_skip;
          auto var2= RandomVariable(img1.size()+j, 10);
          pair_skip_factor.remap_vars({{RandomVariable(0), var}, {RandomVariable(1), var2}});
          mnet.add_factor(pair_skip_factor);
        }
      }
    }
    for(int i=0;i<img2.size();i++){
      auto var = RandomVariable(i+img1.size(), 10);

      auto factor = phi_ocr_reduced[img2[i]];
      mnet.add_var(var);
      factor.remap_vars({{RandomVariable(1), var}});
      mnet.add_factor(factor);

      if(mode >= 1 && i+1<img2.size()){
        auto trans_factor = phi_trans;
        auto var2 = RandomVariable(i+1+img1.size(), 10);
        trans_factor.remap_vars({{RandomVariable(0), var}, {RandomVariable(1), var2}});
        mnet.add_factor(trans_factor);
      }

      if(mode >= 2)
      for(int j=i+1;j<img2.size();j++){
        if(img2[i]==img2[j]){
          auto skip_factor = phi_skip;
          auto var2= RandomVariable(j+img1.size(), 10);
          skip_factor.remap_vars({{RandomVariable(0), var}, {RandomVariable(1), var2}});
          mnet.add_factor(skip_factor);
        }
      }
    }
    return mnet;
  }
  pair<vector<char>, double> infer(MarkovNet &mnet, vector<char> &words){
    pair<vector<char>, double> res;
    map<int, int> predicted;
    map<int, int> correct;
    for(int i=0;i<words.size();i++){
      correct[i] = char_map[words[i]];
    }
    GibbsInferer<MarkovNet> inferer;
    inferer.run(mnet);
    res.second = inferer.infer(predicted);
    res.second = inferer.logLikelihood(correct);
    res.first.resize(words.size());
    for(auto &var_val: predicted)
        res.first[var_val.first] = chars[var_val.second];
    return res;
  }

  void fill_ftable(FactorTable &ft){
    for(size_t i=0;i<ft.table.size();i++){
      ft.table[i] = log(rand()%100+1);
    }
  }
  void test_ftable(){
    printf("---factor_table---\n");
    RandomVariable v1(0,2),v2(1,5),v3(2,4),v4(3,2),v5(4,4);
    FactorTable t1,t2,t3;
    t1.init({v1,v2,v3});
    t2.init({v3,v5,v1,v4});

    t3 = t2.sum({v5,v3});
    printf("table ez: t3(%zu) = sum[t2(%zu)] \n",t3.table.size(), t2.table.size());

    t3 = t1.product(t2);
    printf("table ez: t3(%zu) = t1(%zu) x t2(%zu)\n",t3.table.size(), t1.table.size(), t2.table.size());

    t3 = t2.reduce({{0,3},{2,1}});
  }
  MarkovNet test_mnet(){
    test_ftable();
    MarkovNet mnet;
    RandomVariable v1(0,5),v2(1,2),v3(2,3),v4(3,2),v5(4,10),v6(5,5);
    FactorTable ftable;
    mnet.add_var(v1);
    mnet.add_var(v2);
    mnet.add_var(v3);
    mnet.add_var(v4);
    mnet.add_var(v5);
    mnet.add_var(v6);

    ftable.init({v1,v2});
    fill_ftable(ftable);
    mnet.add_factor(ftable);

    ftable.init({v2,v3,v4});
    fill_ftable(ftable);
    mnet.add_factor(ftable);

    ftable.init({v4,v5});
    fill_ftable(ftable);
    mnet.add_factor(ftable);

    ftable.init({v5,v1});
    fill_ftable(ftable);
    mnet.add_factor(ftable);

    ftable.init({v6});
    fill_ftable(ftable);
    mnet.add_factor(ftable);

    mnet.print();
    return mnet;
  }
  void test_inference(){
    auto mnet = test_mnet();
    printf("---inference---\n");

    GibbsInferer<MarkovNet> inferer;
    inferer.run(mnet);
    map<int, int> predicted;
    double ll = inferer.infer(predicted);
    printf("inferred ll: %f\n",ll);
    double ll2 = inferer.logLikelihood(predicted);
    printf("checked ll: %f\n",ll2);
    for(auto &var_val: predicted){
      printf("  %d: %d\n",var_val.first, var_val.second);
    }
    printf("\n");
  }
};