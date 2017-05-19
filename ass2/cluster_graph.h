#pragma once

#include <vector>
#include <set>
#include "rand_var.h"
#include "markov_net.h"
#include "factor_table.h"

using namespace std;

class ClusterGraph;

class ClusterGraphInferer{
public:
  //runs the algo.
  void run(const ClusterGraph &cgraph){}

  //returns log-likelihood
  double infer(vector<int> &res){return 0;}
};

class ClusterGraph{
  friend class ClusterGraphInferer;
  struct Cluster{
    set<RandomVariable> vars;
    vector<int> nbrs; //includes up_nbr also!
    vector<FactorTable> alpha; //factors in this cluster
    int up_nbr;
    Cluster(){
      up_nbr = -1;
    }
  };
public:
  vector<Cluster> clusters; //clusters[0] is fixed as root.
  void to_clique_tree(const MarkovNet &mnet, const vector<RandomVariable> &ordering){
    map< int, set<int> > taus; //map of cluster_idx -> tau's
    map< int, set<int> > induced;
    map< int, set<int> > phi_map; //map from var_id to phi's containing it.
    clusters.clear();
    int factor_idx = 0;
    for(const auto &factor: mnet.factors){ //initialize with existing cliques
      for(auto it1 = factor.vars.begin(); it1 != factor.vars.end(); ++it1){
        phi_map[it1->id].insert(factor_idx);
        auto it2 = it1;
        ++it2;
        for(; it2 != factor.vars.end(); ++it2){
          induced[it1->id].insert(it2->id);
          induced[it2->id].insert(it1->id);
        }
      }
      factor_idx++;
    }

    int cluster_idx=0;
    for(const auto &var: ordering){
      Cluster C;
      C.vars.insert(var);

      for(auto phi: phi_map[var.id]){
        C.alpha.push_back(mnet.factors[phi]);
        for(auto &var_phis: phi_map){
          if(var_phis.first == var.id) //do not delete from present var_phis, as it will invalidate outer loop iterator
            continue;
          var_phis.second.erase(phi);
        }
      }
      phi_map.erase(var.id); //no longer needed

      for(auto &nbr: induced[var.id]){
        auto nbr_var_iter = mnet.vars.find(RandomVariable(nbr));
        assert(nbr_var_iter != mnet.vars.end());
        C.vars.insert(*nbr_var_iter);
        taus[cluster_idx].insert(nbr);
        induced[nbr].erase(var.id);
        for(auto &nbr2: induced[var.id]){
          if(nbr==nbr2)
            continue;
          induced[nbr].insert(nbr2);
        }
      }
      induced.erase(var.id);
      for(auto &idx_tau: taus){
        auto down_idx = idx_tau.first;
        auto &tau = idx_tau.second;
        if(tau.find(var.id) != tau.end())
          C.nbrs.push_back(down_idx);
      }
      //C.nbrs contains only down_nbrs at this point
      //up_nbr is added by this C's parent.
      for(auto down_idx: C.nbrs){
        clusters[down_idx].up_nbr = cluster_idx;
        clusters[down_idx].nbrs.push_back(cluster_idx);
        taus.erase(down_idx);
      }
      clusters.push_back(C);
      cluster_idx++;
    }
    assert(phi_map.size() == 0); //no more phi should be left!
  }

  void to_bethe_graph(const MarkovNet &mnet){
    clusters.clear();
    map<int, int> var_map; //var_id to cluster_idx map.
    int cluster_idx = 0;
    for(const auto &var: mnet.vars){
      Cluster C;
      C.vars.insert(var);
      clusters.push_back(C);
      var_map[var.id] = cluster_idx;
      cluster_idx++;
    }
    for(const auto &factor: mnet.factors){ //initialize with existing cliques
      Cluster C;
      C.alpha.push_back(factor);
      for(const auto &var: factor.vars){
        C.vars.insert(var);
        C.nbrs.push_back(var_map[var.id]);
        clusters[var_map[var.id]].nbrs.push_back(cluster_idx);
      }
      clusters.push_back(C);
      cluster_idx++;
    }
  }

  void print(){
    int idx=0;
    printf("---cluster_graph---\n");
    for(auto &C: clusters){
      printf("%d: vars: < ",idx);
      for(auto &var: C.vars)
        printf("%d ",var.id);
      printf(">\n   phis: [ ");
      for(auto &phi: C.alpha)
      {
        printf("< ");
        for(auto &var: phi.vars)
          printf("%d ",var.id);
        printf("> ");
      }
      printf("]\n   up_nbr:");
      printf("%d nbrs:[ ",C.up_nbr);
      for(auto &nbr: C.nbrs)
        printf("%d ",nbr);
      printf("]\n");
      idx++;
    }
    printf("\n");
  }
};