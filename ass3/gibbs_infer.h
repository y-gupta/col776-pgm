#pragma once

#include <cstdlib>
#include <cmath>
#include <iostream>

#include "markov_net.h"

using namespace std;

template<typename Net>
class GibbsInferer{
public:
  Net *mnet;
  map<int, FactorTable> beliefs;
  void run(Net &_mnet){
    mnet = &_mnet;
    auto samp = mnet->random_sample();
    size_t vars_sz = mnet->vars.size();
    map<int, vector<int> > counts;
    map<int, double> expectations;
    int total = 0;

    for(auto &var: mnet->vars){
      counts[var.id].resize(var.card, 0);
      expectations[var.id] = 0;
    }
    srand(0);

    while(1){ //burn-in samples
      mnet->next_gibbs_sample(samp);
      // for(auto &var_val: samp.var_vals){
      //   counts[var_val.first][var_val.second]++;
      // }
      total++;
      if(total > 5000)
        break;
    }
    double delta_normalizer = 1;
    for(auto &var: mnet->vars){
      counts[var.id].clear();
      counts[var.id].resize(var.card, 0);
      expectations[var.id] = 0;
      delta_normalizer += var.card;
    }
    total = 0;
    while(1){ //actual samples
      for(int i=0;i<vars_sz;++i)
        mnet->next_gibbs_sample(samp);
      for(auto &var_val: samp.var_vals){
        ++counts[var_val.first][var_val.second];
      }
      ++total;

      if(total % 256 == 0){
        double delta = 0;
        for(auto &var: mnet->vars){
          auto prev_expectation = expectations[var.id];
          expectations[var.id] = 0;
          for(int i=0;i<var.card;i++)
            expectations[var.id] += i*double(counts[var.id][i])/total;
          delta += abs(expectations[var.id] - prev_expectation);
        }
        delta = delta/delta_normalizer;
        if(delta < 0.001)
          break;
      }
    }
    // cout << total <<endl;

    for(auto &var: mnet->vars){
      beliefs[var.id].init({RandomVariable(var.id, var.card)});
      for(int i=0;i<var.card;i++)
        beliefs[var.id].table[i] = log(counts[var.id][i]+1) - log(total+var.card);
    }
  }
  double infer(map<int, int> &res) const { //var.id, value
    double ll = 0;
    res.clear();
    for(const auto &var_belief: beliefs){
      const auto &id = var_belief.first;
      const auto &belief = var_belief.second;
      pair<double, int> logp_val;
      logp_val = belief.maximal_val(id);
      ll += logp_val.first;
      res[id] = logp_val.second;
    }
    return ll;
  }
  double logLikelihood(const map<int, int> &var_vals) const {
    double ll = 0;
    for(const auto &var_val: var_vals){
      const auto &id = var_val.first;
      const auto &belief = beliefs.at(id);
      auto logp = belief.logp(id, var_val.second);
      ll += logp;
    }
    return ll;
  }
};