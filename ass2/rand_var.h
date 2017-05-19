#pragma once

#include <set>
#include <iterator>
#include <algorithm>

using namespace std;

class RandomVariable{
public:
  int id, card;
  RandomVariable(int _id=0,int _card=0):id(_id),card(_card){
  }

  bool operator < (const RandomVariable &b) const{
    return id < b.id;
  }

  static set<RandomVariable> set_difference(const set<RandomVariable> &a, const set<RandomVariable> &b){
    set<RandomVariable> res;
    std::set_difference(a.begin(), a.end(), b.begin(), b.end(), std::inserter(res, res.end()));
    return res;
  }
  static set<RandomVariable> set_union(const set<RandomVariable> &a, const set<RandomVariable> &b){
    set<RandomVariable> res;
    std::set_union(a.begin(), a.end(), b.begin(), b.end(), std::inserter(res, res.end()));
    return res;
  }
  static set<RandomVariable> set_intersection(const set<RandomVariable> &a, const set<RandomVariable> &b){
    set<RandomVariable> res;
    std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::inserter(res, res.end()));
    return res;
  }
};