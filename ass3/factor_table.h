#pragma once

#include <vector>
#include <cmath>
#include <set>
#include <map>
#include <cassert>
#include <float.h>
#include "rand_var.h"

using namespace std;

class FactorTable{
  vector<size_t> strides;
  map<int, int> var_mapping; //maps var.id to column index in table.
public:
  vector<double> table;
  std::set<RandomVariable> vars;
  FactorTable(){
    table.push_back(0); //for empty table, have atleast one entry.
  }
  void init(const std::set<RandomVariable> &_vars){
    vars = _vars;
    strides.clear();
    strides.reserve(vars.size());
    size_t stride = 1;
    int i=0;
    var_mapping.clear();
    for(const auto &var: vars){
      strides.push_back(stride);
      var_mapping[var.id] = i;
      i++;
      stride *= (size_t)var.card;
    }
    table.clear();
    table.resize(stride,0);
  }

  void remap_vars(const map<RandomVariable, RandomVariable> &vars_vars){ //change vars to given, keeping the table/strides intact.
    map<int, int> new_mapping;
    assert(vars_vars.size() == vars.size());
    decltype(vars) new_vars;
    for(const auto &var_var: vars_vars){
      assert(vars.find(var_var.first) != vars.end() && vars.find(var_var.first)->card == var_var.second.card);
      new_mapping[var_var.second.id] = var_mapping[var_var.first.id];
      new_vars.insert(var_var.second);
    }
    vars = new_vars;
    var_mapping = new_mapping;
  }

  size_t get_idx(const vector<int> &vals) const{
    size_t idx=0;
    for(int i=0;i<vals.size();i++){
      idx += strides[i] * (size_t) vals[i];
    }
    return idx;
  }
  size_t get_idx(const map<int, int> &var_vals) const{ //must be superset of vars in this table
    size_t idx=0;
    for(const auto &var_val: var_vals){
      auto it = var_mapping.find(var_val.first);
      if(it == var_mapping.end())
        continue;
      idx += strides.at(it->second) * (size_t) var_val.second;
    }
    return idx;
  }
  double get(const vector<int> &vals) const{ //takes values of random variables
    assert(vars.size() == vals.size());
    return table[get_idx(vals)];
  }
  double get(const map<int, int> &var_vals) const{ //must be superset of vars in this table
    return table[get_idx(var_vals)];
  }
  void set(const vector<int> &vals, double factor){
    assert(vars.size() == vals.size());
    table[get_idx(vals)] = factor;
  }
  double& operator [](const vector<int> &vals){
    assert(vars.size() == vals.size());
    return table[get_idx(vals)];
  }
  FactorTable marginalize(const std::set<RandomVariable> &vars_to_maginalize) const{
#if MAP_INFERENCE
    return maximize(vars_to_maginalize);
#else
    return sum(vars_to_maginalize);
#endif
  }
  double delta(const FactorTable &b) const{
    double del=0;
    if(table.size() != b.table.size())
      return 1;
    for(int i=0;i<table.size();i++){
      del += abs(table[i] - b.table[i]);
    }
    return del;
  }
  //sum/marginalize over vars_to_sum
  FactorTable sum(const std::set<RandomVariable> &vars_to_sum) const{
    FactorTable res;
    auto final_vars = RandomVariable::set_difference(vars, vars_to_sum);

    res.init(final_vars);

    vector<int> vals(vars.size(),0);

    double normalization_const = table[0];
    for(size_t idx=0;idx<table.size();idx++)
      if(normalization_const < table[idx])
        normalization_const = table[idx];

    for(size_t idx=0, idx2=0;idx < table.size(); idx++){
      res.table[idx2] += exp(table[idx]-normalization_const); //assuming log_phi in table.

      for(auto &var: vars){
      // for(size_t var_idx=0; var_idx < vars.size(); var_idx++){
        auto var_idx = var_mapping.at(var.id);
        vals[var_idx] += 1;
        size_t stride = 0; //default stride to 0 for variables that will be summed-out
        auto var_idx2_iter = res.var_mapping.find(var.id); //find idx of currently changed var in res
        if(var_idx2_iter != res.var_mapping.end()) //not a summed-out var. should update idx2 with valid stride!
          stride = res.strides[var_idx2_iter->second];

        if(vals[var_idx] == var.card){ //completed cycle for present var. switch to next.
          vals[var_idx] = 0;
          idx2 -= (var.card-1)*stride; //revert the idx increment from this var
        }else{
          idx2 += stride; //regular update, continue with next idx.
          break;
        }
      }
    }

    for(auto &val: res.table){
      val = log(val);
    }
    return res;
  }

  //max-marginalize over given vars
  FactorTable maximize(const std::set<RandomVariable> &vars_to_max) const{
    FactorTable res;
    auto final_vars = RandomVariable::set_difference(vars, vars_to_max);

    res.init(final_vars);

    vector<int> vals(vars.size(),0);

    double normalization_const = table[0];
    for(size_t idx=0;idx<table.size();idx++)
      if(normalization_const < table[idx])
        normalization_const = table[idx];

    for(size_t idx2=0;idx2<res.table.size();idx2++)
      res.table[idx2] = -DBL_MAX;

    for(size_t idx=0, idx2=0;idx < table.size(); idx++){
      auto row = table[idx]-normalization_const;
      if(row > res.table[idx2])
        res.table[idx2] = row;

      for(auto &var: vars){
      // for(size_t var_idx=0; var_idx < vars.size(); var_idx++){
        auto var_idx = var_mapping.at(var.id);
        vals[var_idx] += 1;
        size_t stride = 0; //default stride to 0 for variables that will be summed-out
        auto var_idx2_iter = res.var_mapping.find(var.id); //find idx of currently changed var in res
        if(var_idx2_iter != res.var_mapping.end()) //not a summed-out var. should update idx2 with valid stride!
          stride = res.strides[var_idx2_iter->second];

        if(vals[var_idx] == var.card){ //completed cycle for present var. switch to next.
          vals[var_idx] = 0;
          idx2 -= (var.card-1)*stride; //revert the idx increment from this var
        }else{
          idx2 += stride; //regular update, continue with next idx.
          break;
        }
      }
    }

    return res;
  }

  FactorTable product(const FactorTable &b) const{
    FactorTable res;
    auto final_vars = RandomVariable::set_union(vars, b.vars);

    res.init(final_vars);

    vector<int> vals(res.vars.size(),0);
    for(size_t idx=0, idx1=0, idx2=0;idx < res.table.size(); idx++){
      res.table[idx] = table[idx1] + b.table[idx2]; //assuming log_phi in table.

      for(auto &var: res.vars){
        auto var_idx = res.var_mapping[var.id];
        vals[var_idx] += 1;
        size_t stride1 = 0;
        size_t stride2 = 0;
        const auto var_idx1_iter=var_mapping.find(var.id);
        if(var_idx1_iter != var_mapping.end())
         stride1 = strides[var_idx1_iter->second];
        const auto var_idx2_iter=b.var_mapping.find(var.id);
        if(var_idx2_iter != b.var_mapping.end())
         stride2 = b.strides[var_idx2_iter->second];

        if(vals[var_idx] == var.card){ //completed cycle for present var. switch to next.
          vals[var_idx] = 0;
          idx1 -= (var.card-1)*stride1; //revert the idx increment from this var
          idx2 -= (var.card-1)*stride2;
        }else{
          idx1 += stride1; //regular update, continue with next idx.
          idx2 += stride2;
          break;
        }
      }
    }

    return res;
  }

  FactorTable reduce(const map<int, int> var_vals) const{ // reduced table, given vals for some vars. var_vals need not be a subset of vars.
    FactorTable res;
    auto final_vars = vars;
    vector<int> vals(vars.size(),0);

    for(const auto &var_val: var_vals){ //remove vars to be reduced from final
      if(final_vars.erase(RandomVariable(var_val.first))) //if var doesn't exist in vars, ignore it, otherwise assign value.
        vals[var_mapping.at(var_val.first)] = var_val.second; //populate val for reducted vars
    }

    res.init(final_vars);

    size_t idx = get_idx(vals); // non reduced will iterate, reduced fixed.


    for(size_t idx2=0;idx2 < res.table.size(); idx2++){
      res.table[idx2] = table[idx];

      for(auto &var: res.vars){
        auto var_idx = var_mapping.at(var.id);
        vals[var_idx] += 1;
        size_t stride = strides[var_idx]; //default stride to 0 for variables that will be summed-out

        if(vals[var_idx] == var.card){ //completed cycle for present var. switch to next.
          vals[var_idx] = 0;
          idx -= (var.card-1)*stride; //revert the idx increment from this var
        }else{
          idx += stride; //regular update, continue with next idx.
          break;
        }
      }
    }

    return res;
  }


  vector<double> reduce_to_one(int var_id, const map<int, int> var_vals) const{
    auto it = vars.find(RandomVariable(var_id));
    assert(it != vars.end());
    const auto &final_var = *it;
    vector<double> res(final_var.card);

    vector<int> vals(vars.size(),0);

    for(const auto &var: vars){ //remove vars to be reduced from final
      if(var.id != final_var.id) //if var doesn't exist in vars, ignore it, otherwise assign value.
        vals[var_mapping.at(var.id)] = var_vals.at(var.id); //populate val for reducted vars
    }

    size_t idx = get_idx(vals); // non reduced will iterate, reduced fixed.
    auto stride = strides[var_mapping.at(final_var.id)];

    for(size_t idx2=0;idx2 < final_var.card; idx2++){
      res[idx2] = table[idx];
      idx += stride; //regular update, continue with next idx.
    }

    return res;
  }

  //maximize table over a given var, return (logp, val)
  pair<double, int> maximal_val(int var_id) const {
    assert(vars.find(var_id) != vars.end());
    auto to_sum = RandomVariable::set_difference(vars, {RandomVariable(var_id)});
    auto marginal = sum(to_sum);
    double max=0, sum=0; int arg_max=-1;
    assert(marginal.table.size() == vars.find(var_id)->card);
    for(int i=0;i<marginal.table.size();i++){
      sum += exp(marginal.table[i]);
      if(arg_max==-1 || marginal.table[i] > max){
        arg_max = i;
        max = marginal.table[i];
      }
    }
    auto logp = max - log(sum);
    return make_pair(logp, arg_max);
  }

  //returns logp for a given var,val
  double logp(int var_id, int val) const {
    assert(vars.find(var_id) != vars.end());
    assert(val < vars.find(var_id)->card);
    auto to_sum = RandomVariable::set_difference(vars, {RandomVariable(var_id)});
    auto marginal = sum(to_sum);
    double sum=0;
    assert(marginal.table.size() == vars.find(var_id)->card);
    for(int i=0;i<marginal.table.size();i++){
      sum += exp(marginal.table[i]);
    }
    auto logp = marginal.table[val] - log(sum);
    return logp;
  }
};