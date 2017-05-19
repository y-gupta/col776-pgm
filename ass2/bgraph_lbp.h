#pragma once
#include "cluster_graph.h"
#include "factor_table.h"

class BetheGraphLBPInferer:public ClusterGraphInferer{
  ClusterGraph *G;
  vector<FactorTable> psis;
  vector<FactorTable> beliefs;
  map<RandomVariable, int> var_beliefs; //mapping from var to it's cheapest belief
  map<pair<int,int>, FactorTable> msgs;
public:
  void run(ClusterGraph &cgraph){
    G = &cgraph;

    initialize();
    int cluster_idx;
    map<pair<int,int>, FactorTable> next_msgs;
    int iteration = 0;
    double change;
    while(iteration < 50){
      cluster_idx = 0;
      for(const auto &C: G->clusters){
        for(auto nbr: C.nbrs){
          FactorTable del = calcMessage(cluster_idx, nbr);
          next_msgs[make_pair(cluster_idx, nbr)] = del;
        }
        cluster_idx++;
      }
      auto it=next_msgs.begin();
      change = 0;
      for(auto &msg: msgs){
        change += it->second.delta(msg.second);
        ++it;
      }
      msgs = next_msgs;
      iteration++;
      if(change < 0.000001){
        // printf("converge: %d\n",iteration);
        break;
      }
    }

    beliefs.resize(cluster_idx);

    while(cluster_idx-- != 0){
      beliefs[cluster_idx] = calcBelief(cluster_idx);
    }
    // assert(var_beliefs.size() == G->clusters.size());
    psis.clear(); //no longer need psis!!!
  }
  FactorTable calcBelief(int i){
    FactorTable b = psis[i];
    for(auto nbr: G->clusters[i].nbrs){
      b = b.product(msgs.at(make_pair(nbr,i)));
    }

    for(const auto &var: G->clusters[i].vars){ //assign cheapest beliefs to vars
      if(var_beliefs.find(var) != var_beliefs.end()
       && beliefs[var_beliefs[var]].vars.size() <= G->clusters[i].vars.size()){
        continue;
      }
      var_beliefs[var] = i;
    }
    return b;
  }
  FactorTable calcMessage(int i, int j){
    FactorTable psi = psis[i];
    for(auto nbr: G->clusters[i].nbrs){
      if(nbr==j)
        continue;
      psi = psi.product(msgs.at(make_pair(nbr,i)));
    }
    psi = psi.marginalize(RandomVariable::set_difference(G->clusters[i].vars, G->clusters[j].vars));
    return psi;
  }

  void initialize(){
    int cluster_idx=0;
    for(const auto &C: G->clusters){
      FactorTable psi;
      psi.init(C.vars);
      for(const auto &phi: C.alpha){
        psi = psi.product(phi);
      }
      psis.push_back(psi);
      for(auto nbr: C.nbrs){
        msgs[make_pair(cluster_idx,nbr)] = FactorTable();
      }
      cluster_idx++;
    }
  }

  double infer(map<RandomVariable, int> &res) const {
    double ll = 0;
    res.clear();
    for(const auto &var_belief: var_beliefs){
      const auto &var = var_belief.first;
      const auto &belief = beliefs[var_belief.second];
      pair<double, int> logp_val;
      if(MAP_INFERENCE)
        logp_val = belief.reduce(res).maximal_val(var);
      else
        logp_val = belief.maximal_val(var);
      ll += logp_val.first;
      res[var] = logp_val.second;
    }
    return ll;
  }
  double logLikelihood(const map<RandomVariable, int> &var_vals) const {
    double ll = 0;
    for(const auto &var_val: var_vals){
      const auto &var = var_val.first;
      const auto &belief = beliefs[var_beliefs.at(var)];
      auto logp = belief.logp(var, var_val.second);
      ll += logp;
    }
    return ll;
  }
};