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
#include "bayes_net.h"
#include "gibbs_infer.h"
#include "timer.h"

using namespace std;

class GeneralNet{
public:
  Timer timer;
  bool use_markov;

  BayesNet bnet;
  MarkovNet mnet;
  GeneralNet(){
    use_markov = false;
  }
  void load(const char* fname){
    bnet = BayesNet::load(fname);
    if(use_markov){
      mnet = bnet.moralize();
    }else{
      bnet.init();
    }
  }
  void train(const char* fname_train){
    FILE *ftrain = fopen(fname_train, "r");
    assert(ftrain);
    char line[10241], tmp[1025];

    fgets(line, 10240, ftrain);
    char *remaining;
    int bytes_read = 0;
    remaining = line;

    vector<int> var_order;

    while(sscanf(remaining, "%s%n", tmp, &bytes_read)==1){
      remaining += bytes_read;
      var_order.push_back(bnet.name_id.at(tmp));
    }
    assert(var_order.size() == bnet.vars.size());

    Sample samp;

    while(fgets(line, 10240, ftrain)){
      remaining = line;
      bytes_read = 0;

      int i=0;
      while(sscanf(remaining, "%s%n", tmp, &bytes_read)==1){
        remaining += bytes_read;
        int var = var_order[i];
        int val = bnet.id_name_value.at(var).at(tmp);
        samp.var_vals[var] = val;
        i++;
        assert(i <= var_order.size());
      }
      if(use_markov){
        mnet.learn(samp);
      }else{
        bnet.learn(samp);
      }
    }
    if(use_markov)
      mnet.finalize_learn();
    else
      bnet.normalize();

    string fname_bif(fname_train);
    if(use_markov == false)
    {  fname_bif += ".bn.bif";
      bnet.output(fname_bif.c_str());
    }else
    { fname_bif += ".mn.bif";
      bnet.output_markov(fname_bif.c_str(), mnet);
    }
  }
  void test(const char* fname_test, const char* fname_true){
    string fname_out(fname_test);
    if(use_markov ==false)
      fname_out += ".bn.out";
    else
      fname_out += ".mn.out";
    FILE *fout = fopen(fname_out.c_str(),"w");
    FILE *ftest = fopen(fname_test, "r"), *ftrue = fopen(fname_true, "r");
    assert(fout && ftest && ftrue);

    char line[10241], tmp[1025];

    fgets(line, 10240, ftrue);
    fgets(line, 10240, ftest);
    char *remaining;
    int bytes_read = 0;
    remaining = line;

    vector<int> var_order;

    while(sscanf(remaining, "%s%n", tmp, &bytes_read)==1){
      remaining += bytes_read;
      var_order.push_back(bnet.name_id.at(tmp));
    }
    assert(var_order.size() == bnet.vars.size());

    Sample samp, true_samp;

    size_t total=0, correct=0, total_queries = 0;
    double total_ll = 0;

    while(fgets(line, 10240, ftest)){
      remaining = line;
      bytes_read = 0;

      int i=0;
      samp.var_vals.clear();
      true_samp.var_vals.clear();
      while(sscanf(remaining, "%s%n", tmp, &bytes_read)==1){
        remaining += bytes_read;
        if(tmp[0]!='?')
        {
          int var = var_order[i];
          int val = bnet.id_name_value.at(var).at(tmp);
          samp.var_vals[var] = val;
        }
        i++;
        assert(i <= var_order.size());
      }
      i=0;
      fgets(line, 10240, ftrue);
      remaining = line;
      bytes_read = 0;
      while(sscanf(remaining, "%s%n", tmp, &bytes_read)==1){
        remaining += bytes_read;
        int var = var_order[i];
        int val = bnet.id_name_value.at(var).at(tmp);
        true_samp.var_vals[var] = val;
        i++;
        assert(i <= var_order.size());
      }
      if(use_markov == false){
        auto reduced_bnet = bnet.reduce(samp.var_vals);
        GibbsInferer<BayesNet> inferer;
        inferer.run(reduced_bnet);
        map<int,int> prediction;
        inferer.infer(prediction);
        map<int, int> unknown_trues = true_samp.var_vals;
        for(auto &id_val: samp.var_vals)
          unknown_trues.erase(id_val.first);
        total_ll += inferer.logLikelihood(unknown_trues);
        total_queries++;
        for(auto id_true: unknown_trues){
          auto id = id_true.first;
          total++;
          if(prediction[id] == true_samp.var_vals[id])
            correct++;
          fprintf(fout,"%s ", bnet.id_name.at(id).c_str());
          for(int i=0;i<bnet.get_var(id).card;i++){
            double p = inferer.beliefs.at(id).logp(id, i);
            p = exp(p);
            fprintf(fout, "%s:%lf ",bnet.id_value_name.at(id).at(i).c_str(),p);
          }
          fprintf(fout, "\n");
        }
      }else{
        auto reduced_mnet = mnet.reduce(samp.var_vals);
        GibbsInferer<MarkovNet> inferer;
        inferer.run(reduced_mnet);
        map<int,int> prediction;
        inferer.infer(prediction);
        map<int, int> unknown_trues = true_samp.var_vals;
        for(auto &id_val: samp.var_vals)
          unknown_trues.erase(id_val.first);
        total_ll += inferer.logLikelihood(unknown_trues);
        total_queries++;
        for(auto id_true: unknown_trues){
          auto id = id_true.first;
          total++;
          if(prediction[id] == true_samp.var_vals[id])
            correct++;
          fprintf(fout,"%s ", bnet.id_name.at(id).c_str());
          for(int i=0;i<bnet.get_var(id).card;i++){
            double p = inferer.beliefs.at(id).logp(id, i);
            p = exp(p);
            fprintf(fout, "%s:%lf ",bnet.id_value_name.at(id).at(i).c_str(),p);
          }
          fprintf(fout, "\n");
        }
      }
      fprintf(fout, "\n");
    }
    fclose(fout);
    cout<<"Average Accuracy: "<<float(correct)/total<<endl;
    cout<<"Average LL: " << total_ll/total_queries<<endl;
  }
};