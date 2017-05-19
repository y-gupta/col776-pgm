#pragma once

#include <set>
#include <vector>
#include <map>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "rand_var.h"
#include "markov_net.h"
#include "factor_table.h"
#include "sample.h"

using namespace std;
class BayesNet{
  map<int, int> var_map; //var.id to idx in vars
public:
  class Node: public RandomVariable{
  public:
    FactorTable cpt;
    vector<int> parents, children; //var ids stored here
  };
  vector<Node> vars;
  vector<int> unobserved_vars; //var.id's
  map<int, int> val_map; //var.id to val for observed vars

  map<string, int> name_id; //var names to var id
  map<int, map<string, int> > id_name_value; //var.id to var_name->idx map
  map<int, string> id_name; //var.id to name
  map<int, map<int, string> > id_value_name; //var.id to idx->name map

  void add_var(const Node &node){
    var_map[node.id] = vars.size();
    vars.push_back(node);
    unobserved_vars.push_back(node.id);
  }
  Node& get_var(int id){
    assert(var_map.at(id) < vars.size());
    return vars[var_map.at(id)];
  }

  MarkovNet moralize(){
    MarkovNet mnet;
    for(auto &var: vars){
      mnet.add_var(RandomVariable(var.id, var.card));
      set<RandomVariable> phi_vars({RandomVariable(var.id, var.card)});
      for(int p1: var.parents){
        phi_vars.insert(RandomVariable(p1,get_var(p1).card));
      }
      FactorTable phi;
      phi.init(phi_vars);
      mnet.add_factor(phi);
    }
    return mnet;
  }

  size_t learn_count;
  void init(){
    assert(val_map.size() == 0); //only works for fully unobserved nets :(
    learn_count = 0;
    for(auto &var: vars){
      set<RandomVariable> cpt_vars({RandomVariable(var.id, var.card)});
      for(int id: var.parents){
        cpt_vars.insert(RandomVariable(id, get_var(id).card));
      }
      var.cpt.init(cpt_vars);
    }
  }
  void learn(Sample &samp){
    learn_count++;
    for(auto &var: vars){
      auto idx = var.cpt.get_idx(samp.var_vals);
      var.cpt.table[idx]++;
    }
  }
  void normalize(){
    for(auto &var: vars){
      map<int,int> var_vals;
      size_t sz = var.cpt.table.size()/var.card;

      for(auto p:var.parents){
        var_vals[p] = 0;
      }
      var_vals[var.id] = 0;

      for(size_t i = 0;i<sz;i++){
        double sum = 0;
        for(int j=0;j<var.card;j++){
          var_vals[var.id] = j;
          auto idx = var.cpt.get_idx(var_vals);
          var.cpt.table[idx] = var.cpt.table[idx]+1;
          sum += var.cpt.table[idx];
        }
        for(int j=0;j<var.card;j++){
          var_vals[var.id] = j;
          auto idx = var.cpt.get_idx(var_vals);
          var.cpt.table[idx] = log(var.cpt.table[idx]) - log(sum);
        }
        for(auto p:var.parents){
          var_vals[p]++;
          if(var_vals[p] == get_var(p).card)
            var_vals[p] = 0;
          else
            break;
        }
      }
    }
  }

  BayesNet reduce(const map<int, int> &var_vals) const{
    BayesNet res;
    for(auto &var: vars){
      res.add_var(var);
      auto &cpt = res.vars.back().cpt;
      cpt = cpt.reduce(var_vals);
    }
    res.unobserved_vars.clear();
    for(auto id: unobserved_vars){
      if(var_vals.find(id) == var_vals.end())
        res.unobserved_vars.push_back(id);
    }
    res.val_map = val_map;
    for(const auto &var_val: var_vals)
      res.val_map[var_val.first] = var_val.second;
    return res;
  }

  Sample random_sample()const{
    Sample samp;
    for(const auto &var: vars){
      auto it = val_map.find(var.id);
      if(it == val_map.end())
        samp.var_vals[var.id] = rand() % var.card;
      else
        samp.var_vals[var.id] = it->second;
    }
    return samp;
  }

  void next_gibbs_sample(Sample &samp){
    // int idx = rand() % vars.size();
    static int idx = 0;
    idx = (idx+1) % unobserved_vars.size();
    auto &var = get_var(unobserved_vars[idx]);
    samp.var_vals.erase(var.id);
    auto dist = var.cpt.reduce(samp.var_vals); //value of parents taken from here

    for(int child: var.children){
      // auto tmp = get_var(child).cpt.reduce(samp.var_vals).table;
      auto tmp = get_var(child).cpt.reduce_to_one(var.id, samp.var_vals);
      assert(tmp.size() == var.card);
      for(int i=0;i<var.card;i++){
        dist.table[i] += tmp[i];
      }
    }

    double sum = 0;
    for(int i=0;i < var.card; i++){
      sum += exp(dist.table[i]);
      dist.table[i] = sum; //dist.table is now cumulative p's
    }
    auto f = sum * double(rand())/RAND_MAX;
    auto it = lower_bound(dist.table.begin(), dist.table.end(), f); //binary search ftw!
    samp.var_vals[var.id] = it - dist.table.begin();
  }

  void print(){
    int idx=0;
    printf("---bayes_net---\n");
    for(auto &var: vars){
        printf("%d: < ",var.id);
        for(auto &child: var.children){
          printf("%d ", child);
        }
        printf("> cpt sz: %zu\n", var.cpt.table.size());
    }
    printf("\n");
  }

  static BayesNet load(const char* fname){
    string line;
    int find=0;
    ifstream myfile(fname);
    string temp;
    string name;
    map<string, int> val_map;
    map<int, string> inv_val_map;

    BayesNet res;

    auto &var_map = res.name_id; //name -> var.id
    auto &vals_map = res.id_name_value;
    auto &inv_var_map = res.id_name;
    auto &inv_vals_map = res.id_value_name;

    assert(myfile.is_open());
    while(!myfile.eof())
    {
      stringstream ss;
      getline(myfile,line);
      ss.str(line);
      ss>>temp;
      if(temp.compare("variable")==0)
      {
        ss>>name;
        getline(myfile,line);
        stringstream ss2;
        ss2.str(line);
        for(int i=0;i<7;i++)
        {
          ss2>>temp;
        }
        val_map.clear();
        inv_val_map.clear();
        // cout<<name<<": ";
        while(temp.compare("};")!=0)
        {
          if(temp.back()==',')
            temp.erase(temp.end()-1);
          // cout<<temp<<" ";
          val_map.emplace(temp, val_map.size());
          inv_val_map.emplace(val_map.size()-1, temp);
          ss2>>temp;
        }
        Node var;
        var.id = var_map.size();
        var.card = val_map.size();
        var_map.emplace(name, var.id);
        inv_var_map.emplace(var.id, name);

        vals_map.emplace(var.id, val_map);
        inv_vals_map.emplace(var.id, inv_val_map);

        res.add_var(var);
      }else if(temp.compare("probability")==0){
        ss>>temp; //(
        ss>>name;
        // cout<<name<<": ";
        auto &main_var = res.get_var(var_map.at(name));
        ss>>name; // | or )

        while(name.compare(")")!=0){
          ss>>name;
          if(name.back() == ',')
            name.erase(name.end()-1);
          if(name == ")")
            break;
          // cout<<name<<" ";
          auto &parent = res.get_var(var_map.at(name));
          parent.children.push_back(main_var.id);
          main_var.parents.push_back(parent.id);
        }
        // cout<<endl;
      }
    }

    myfile.close();
    return res;
  }

  void output(const char *fname){
    FILE *ft;
    ft = fopen(fname, "w");
    fprintf(ft, "network unknown {\n}\n");
    for(auto &var: vars){
      fprintf(ft, "variable %s {\n  type discrete [ %d ] {", id_name[var.id].c_str(), var.card);
      for(int i=0;i<var.card;i++){
        fprintf(ft, "%s%s",(i==0?" ":", "),id_value_name[var.id][i].c_str());
      }
      fprintf(ft, " };\n}\n");
    }

    for(auto &var: vars){
      fprintf(ft, "probability ( %s ", id_name[var.id].c_str());
      if(var.parents.size() == 0){
        fprintf(ft, ") {\n  table");
        for(int i=0;i<var.card;i++)
          fprintf(ft, "%s%lf",i==0?" ":", ",exp(var.cpt.table[i]));
        fprintf(ft, ";\n}\n");
      }else{
        fprintf(ft, "|");
        size_t sz = 1;
        for(int i=0;i<var.parents.size();i++){
          sz *= get_var(var.parents[i]).card;
          fprintf(ft, "%s%s",(i==0?" ":", "),id_name[var.parents[i]].c_str());
        }
        fprintf(ft, " ) {\n");
        map<int, int> var_vals;
        for(int i=0;i<sz;i++){
          for(int j=0;j<var.parents.size();j++){
            int p = var.parents[j];
            fprintf(ft,"%s%s",j==0?"  (":", ",id_value_name[p][var_vals[p]].c_str());
          }
          fprintf(ft, ")");
          for(int j=0;j<var.card;j++){
            var_vals[var.id] = j;
            fprintf(ft,"%s%f",j==0?" ":", ",exp(var.cpt.get(var_vals)));
          }
          fprintf(ft, ";\n");
          for(int j=0;j<var.parents.size();j++){
            int p = var.parents[j];
            var_vals[p]++;
            if(var_vals[p] == get_var(p).card){
              var_vals[p] = 0;
            }else
              break;
          }
        }
        fprintf(ft, " }\n");
      }
    }
    fclose(ft);
  }

  void output_markov(const char *fname, MarkovNet &mnet){
    FILE *ft;
    ft = fopen(fname, "w");
    fprintf(ft, "network unknown {\n}\n");
    for(auto &var: vars){
      fprintf(ft, "variable %s {\n  type discrete [ %d ] {", id_name[var.id].c_str(), var.card);
      for(int i=0;i<var.card;i++){
        fprintf(ft, "%s%s",(i==0?" ":", "),id_value_name[var.id][i].c_str());
      }
      fprintf(ft, " };\n}\n");
    }

    for(auto &factor: mnet.factors){
      fprintf(ft, "probability (");
      {
        size_t sz = factor.table.size();
        int i=0;
        for(auto &var: factor.vars){
          fprintf(ft, "%s%s",(i==0?" ":", "),id_name[var.id].c_str());
          i++;
        }
        fprintf(ft, " ) {\n");
        map<int, int> var_vals;
        for(int i=0;i<sz;i++){
          int j=0;
          for(auto &var: factor.vars){
            fprintf(ft,"%s%s",j==0?"  (":", ",id_value_name[var.id][var_vals[var.id]].c_str());
            j++;
          }
          fprintf(ft,") %f;\n",exp(factor.get(var_vals)));
          for(auto &var: factor.vars){
            var_vals[var.id]++;
            if(var_vals[var.id] == var.card){
              var_vals[var.id] = 0;
            }else
              break;
          }
        }
        fprintf(ft, " }\n");
      }
    }
    fclose(ft);
  }
};

