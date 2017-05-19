#pragma once

#include <set>
#include <vector>
#include <map>
#include "rand_var.h"
#include "factor_table.h"
using namespace std;
class MarkovNet{
public:
  set<RandomVariable> vars;
  vector<FactorTable> factors;
  void add_var(const RandomVariable &var){
    vars.insert(var);
  }
  void add_factor(const FactorTable &factor){
    factors.push_back(factor);
  }
  vector<RandomVariable> get_ordering() const{
    //returns variable elimination ordering by min-fill
    vector<RandomVariable> ordering;
    set<int> eliminated;
    map<int, set<int> > induced; //induced graph adjacency list
    for(const auto &var: vars){ //initialize all vars to empty list.
      induced[var.id] = {};
    }
    for(const auto &factor: factors){ //initialize with existing cliques
      for(auto it1 = factor.vars.begin(); it1 != factor.vars.end(); ++it1){
        auto it2 = it1;
        ++it2;
        for(; it2 != factor.vars.end(); ++it2){
          induced[it1->id].insert(it2->id);
          induced[it2->id].insert(it1->id);
        }
      }
    }

    for(int i=0;i<vars.size();i++){
      int arg_min=-1, min_fill=vars.size();
      for(auto &node_nbrs: induced){
        int fill = 0;
        auto &nbrs = node_nbrs.second;
        for(auto &nb1: nbrs){
          for(auto &nb2: nbrs){
            if(nb1==nb2)
              continue;
            if(induced[nb1].find(nb2) == induced[nb1].end())
              fill++;
          }
        }
        fill = fill/2;
        if(fill < min_fill){
          min_fill = fill;
          arg_min = node_nbrs.first;
        }
      }
      assert(arg_min != -1);
      auto &nbrs = induced[arg_min];
      for(auto &nb1: nbrs){
        induced[nb1].erase(arg_min);
        for(auto &nb2: nbrs){
          if(nb1==nb2)
            continue;
          induced[nb1].insert(nb2);
        }
      }
      induced.erase(arg_min);
      // eliminated.insert(arg_min);
      auto var_iter = vars.find(RandomVariable(arg_min));
      assert(var_iter != vars.end());
      ordering.push_back(*var_iter);
    }
    assert(ordering.size() == vars.size());
    return ordering;
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