/*

This file is part of VROOM.

Copyright (c) 2015-2016, Julien Coupey.
All rights reserved (see LICENSE).

*/

#include "tsp.h"

tsp::tsp(const cl_args_t& cl_args): 
  _loader(cl_args.input),
  _matrix(_loader.get_locations()),
  _symmetrized_graph(_matrix),
  _cl_args(cl_args) {}

const matrix<distance_t>& tsp::get_matrix() const{
  return _matrix;
}

const undirected_graph<distance_t>& tsp::get_symmetrized_graph() const{
  return _symmetrized_graph;
}

std::size_t tsp::size() const{
  return _matrix.size();
}

distance_t tsp::cost(const std::list<index_t>& tour) const{
  distance_t cost = 0;
  index_t init_step = 0;        // Initialization actually never used.

  auto step = tour.cbegin();
  if(tour.size() > 0){
    init_step = *step;
  }

  index_t previous_step = init_step;
  ++step;
  for(; step != tour.cend(); ++step){
    cost += _matrix[previous_step][*step];
    previous_step = *step;
  }
  if(tour.size() > 0){
    cost += _matrix[previous_step][init_step];
  }
  return cost;
}

void tsp::get_route(const std::list<index_t>& tour,
                    rapidjson::Value& value,
                    rapidjson::Document::AllocatorType& allocator) const{
  assert(tour.size() == _matrix.size());
  _loader.get_route(tour, value, allocator);
}

void tsp::get_tour(const std::list<index_t>& tour,
                   rapidjson::Value& value,
                   rapidjson::Document::AllocatorType& allocator) const{
  assert(tour.size() == _matrix.size());
  _loader.get_tour(tour, value, allocator);
}

void tsp::get_route_infos(const std::list<index_t>& tour,
                          rapidjson::Document& output) const{
  assert(tour.size() == _matrix.size());

  // Back to the starting location when the trip is a loop.
  std::list<index_t> actual_trip (tour);
  actual_trip.push_back(actual_trip.front());
  _loader.get_route_infos(actual_trip, output);
}
