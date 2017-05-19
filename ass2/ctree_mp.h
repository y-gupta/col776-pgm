#pragma once
#include "cluster_graph.h"
#include "factor_table.h"

class CliqueTreeMPInferer:public ClusterGraphInferer{
  ClusterGraph *G;
  vector<FactorTable> psis;
  vector<FactorTable> beliefs;
  map<RandomVariable, int> var_beliefs; //mapping from var to it's cheapest belief
  map<pair<int,int>, FactorTable> msgs;
public:
  void run(ClusterGraph &cgraph){
    G = &cgraph;

    initialize();
    int cluster_idx = 0;
    for(const auto &C: G->clusters){ //topological ordering
      if(C.up_nbr != -1) //root nodes have just one var!
      {
        FactorTable del = calcMessage(cluster_idx, C.up_nbr);
        msgs[make_pair(cluster_idx, C.up_nbr)] = del;
      }
      cluster_idx++;
    }

    beliefs.resize(cluster_idx);

    while(cluster_idx-- != 0){ //reverse-topological ordering
      auto &C = G->clusters[cluster_idx];
      for(auto nbr: C.nbrs){
        if(nbr == C.up_nbr) //do not send to parent again!!
          continue;
        FactorTable del = calcMessage(cluster_idx, nbr);
        msgs[make_pair(cluster_idx, nbr)] = del;
      }
      beliefs[cluster_idx] = calcBelief(cluster_idx);

      set<RandomVariable> eliminated_vars;
      if(C.up_nbr == -1){
        eliminated_vars = C.vars;
      }else{
        eliminated_vars = RandomVariable::set_difference(C.vars, G->clusters[C.up_nbr].vars);
      }
      assert(eliminated_vars.size() == 1);
      for(auto &var: eliminated_vars){
        var_beliefs[var] = cluster_idx;
      }
    }
    assert(var_beliefs.size() == G->clusters.size());
    psis.clear(); //no longer need psis!!!
  }
  FactorTable calcBelief(int i){
    FactorTable b = psis[i];
    for(auto nbr: G->clusters[i].nbrs){
      b = b.product(msgs.at(make_pair(nbr,i)));
      msgs.erase(make_pair(nbr,i)); //should no longer need this msg!
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
    for(const auto &C: G->clusters){
      FactorTable psi;
      psi.init(C.vars);
      for(const auto &phi: C.alpha){
        psi = psi.product(phi);
      }
      psis.push_back(psi);
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