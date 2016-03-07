#ifndef MATRIX_H
#define MATRIX_H

/*

This file is part of VROOM.

Copyright (c) 2015-2016, Julien Coupey.
All rights reserved (see LICENSE).

*/

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include "./typedefs.h"
#include "../utils/exceptions.h"

template <class T>
class line{
 private:
  const index_t _row;
  const std::vector<std::pair<double, double>>* const _locations_ptr;

 public:
  line(index_t row,
       const std::vector<std::pair<double, double>>& locations):
    _row(row),
    _locations_ptr(&locations){
  }

  T operator[](index_t index) const {
    return (_row == index) ?
      3 * (std::numeric_limits<distance_t>::max() / 4):
      (int) (sqrt(std::pow((*_locations_ptr)[_row].first - (*_locations_ptr)[index].first, 2)
                  + std::pow((*_locations_ptr)[_row].second - (*_locations_ptr)[index].second, 2)) + 0.5);
  }
};

template <class T>
class matrix{
private:
  std::vector<std::pair<double, double>> _locations;
  std::vector<line<T>> _lines;

public:
  const line<T>& operator[](index_t index) const{
    return _lines[index];
  }

  matrix(std::vector<std::pair<double, double>> locations):
    _locations(locations){
    for(index_t i = 0; i < _locations.size(); ++i){
      _lines.emplace_back(i, _locations);
    }
  }

  matrix(){}

  std::size_t size() const{
    return _locations.size();
  }

  matrix<T> get_sub_matrix(const std::vector<index_t>& indices) const{
    std::vector<std::pair<double, double>> new_locations;
    std::transform(indices.begin(), indices.cend(),
                   std::back_inserter(new_locations),
                   [this](index_t i) {return _locations[i];});
    return matrix<T>(new_locations);
  }
};

#endif
