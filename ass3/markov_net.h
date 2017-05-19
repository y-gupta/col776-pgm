#pragma once

#include <set>
#include <vector>
#include <map>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <cfloat>

#include "rand_var.h"
#include "factor_table.h"
#include "sample.h"

using namespace std;
class MarkovNet{
public:
  vector<RandomVariable> vars;
  vector<FactorTable> factors;
  map<int, vector<int> > factors_map;
  size_t num_samples;
  MarkovNet(){
    num_samples = 0;
  }
  void add_var(const RandomVariable &var){
    vars.push_back(var);
  }

  void add_factor(const FactorTable &factor){
    for(const auto &var: factor.vars){
      factors_map[var.id].push_back(factors.size());
    }
    factors.push_back(factor);
  }

  MarkovNet reduce(const map<int, int> &var_vals){
    MarkovNet res;
    for(auto &var: vars){
      if(var_vals.find(var.id) == var_vals.end())
        res.add_var(var);
    }
    for(auto &factor: factors){
      res.add_factor(factor.reduce(var_vals));
    }
    return res;
  }
  Sample random_sample(){
    Sample samp;
    for(auto &var: vars){
      samp.var_vals[var.id] = rand() % var.card;
    }
    return samp;
  }

  void next_gibbs_sample(Sample &samp){
    // int idx = rand() % vars.size();
    static int idx = 0;
    idx = (idx+1) % vars.size();
    auto &var = vars[idx];
    samp.var_vals.erase(var.id);
    FactorTable dist;
    dist.init({var});
    for(int factor_idx: factors_map[var.id]){
      // auto tmp = factors[factor_idx].reduce(samp.var_vals).table;
      auto tmp = factors[factor_idx].reduce_to_one(var.id, samp.var_vals);
      for(int i=0;i<var.card;i++){
        dist.table[i] += tmp[i];
      }
    }
    assert(dist.vars.size() == 1);
    double sum = 0;
    for(int i=0;i < var.card; i++){
      sum += exp(dist.table[i]);
      dist.table[i] = sum; //dist.table is now cumulative p's
    }
    auto f = sum * double(rand())/RAND_MAX;
    auto it = lower_bound(dist.table.begin(), dist.table.end(), f); //binary search ftw!
    samp.var_vals[var.id] = it - dist.table.begin();
  }

  void finalize_learn(double base_eta=1, double C=1){
    auto data_factors = std::move(factors);
    factors.clear();

    for(auto &factor: data_factors){
      factors.push_back(FactorTable());
      factors.back().init(factor.vars);
      for(int i=0;i<factor.table.size();i++){
        factor.table[i] = (factor.table[i]+1)/double(num_samples);
      }
    }

    vector<FactorTable> initial_factors = factors;
    vector<FactorTable> samp_factors;
    double prev_total_delta = DBL_MAX;
    int iteration = 0;
    double eta = base_eta;
    while(1){
      Sample samp = random_sample();
      size_t count = 0;
      while(count < 8000){ //burn-in
        next_gibbs_sample(samp);
        ++count;
      }
      samp_factors = initial_factors;
      count = 0;
      while(1){
        for(int i=0;i<20;i++)
          next_gibbs_sample(samp);
        for(auto &factor: samp_factors){
          auto idx = factor.get_idx(samp.var_vals);
          ++factor.table[idx];
        }
        ++count;
        if(count==8000)
          break;
      }

      double total_delta = 0;

      for(int f=0;f<factors.size();f++){
        auto &samp_factor = samp_factors[f];
        auto &data_factor = data_factors[f];
        auto &factor = factors[f];
        for(int i=0;i<factor.table.size();i++){
          factor.table[i] = exp(factor.table[i]);
          auto delta = data_factor.table[i] - (samp_factor.table[i])/double(count) + C*factor.table[i]/double(count);
          delta *= eta;
          total_delta += abs(delta/factor.table[i]);
          factor.table[i] += delta;
          if(factor.table[i]<=0)
            factor.table[i] = -DBL_MAX;
          factor.table[i] = log(factor.table[i]);
        }
      }
      cout<< total_delta<<" eta:"<<eta<<endl;
      if(prev_total_delta < total_delta){
        eta = eta/1.5;
      }
      prev_total_delta = total_delta;
      iteration++;
      // eta = base_eta*10/(iteration+10);
      // if(total_delta < 0.1)
      //   break;
      if(iteration > 50)
        break;
    }
  }
  void learn(const Sample &samp){
    ++num_samples;
    for(auto &factor: factors){
      auto idx = factor.get_idx(samp.var_vals);
      ++factor.table[idx];
    }
  }
  void print(){
    int idx=0;
    printf("---markov_net---\nvars:");
    for(auto &var: vars)
        printf("%d ",var.id);
    printf("\n");
    for(auto &f: factors){
      printf("%d: < ",idx);
      for(auto &var: f.vars)
        printf("%d ",var.id);
      printf("> phi sz: %zu\n",f.table.size());
      idx++;
    }
    printf("\n");
  }
};

