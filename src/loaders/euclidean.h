/*
  VROOM (Vehicle Routing Open-source Optimization Machine)
  Copyright (C) 2015, Julien Coupey

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or (at
  your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef EUCLIDEAN_H
#define EUCLIDEAN_H
#include <vector>
#include <regex>
#include <boost/regex.hpp>
#include "../../include/rapidjson/document.h"
#include "../structures/matrix.h"
#include "../utils/exceptions.h"

class euclidean{

private:
  struct Node {std::size_t index; double x; double y;};
  
  std::vector<Node> _locations;

  void add_lat_lon_location(const std::string location){
    // Regex check for valid location.
    std::regex valid_loc ("loc=-?[0-9]+\\.?[0-9]*,-?[0-9]+\\.?[0-9]*[[:space:]]*");
    if(!std::regex_match(location, valid_loc)){
      throw custom_exception("invalid syntax for location "
                             + std::to_string(_locations.size() + 1)
                             + ", see vroom -h for usage display."
                             );
    }

    // Parsing the location is now safe.
    std::size_t separator_rank = location.find(",");
    std::string lat = location.substr(4, separator_rank);
    std::string lon = location.substr(separator_rank + 1, location.length() -1);
    _locations.push_back({_locations.size(), std::stod(lat, nullptr), std::stod(lon, nullptr)});
  }
  
public:
  euclidean(std::string input){
    bool use_lat_lon_query
      = (input.find("DIMENSION") == std::string::npos);

    if(use_lat_lon_query){
      // Parsing input in locations from loc=lon,lat&... format.
      std::size_t start = 0;
      std::size_t end = input.find("&", start);
      while(end != std::string::npos){
        this->add_lat_lon_location(input.substr(start, end - start));
        start = end + 1;
        end = input.find("&", start);
      }
      // Adding last element, after last "&".
      end = input.length();
      this->add_lat_lon_location(input.substr(start, end - start));
    }
    else{
      // Use nodes from TSPLIB format (mostly duplicate code from
      // tsplib_loader).

      // 1. Get problem dimension.
      boost::regex dim_rgx ("DIMENSION[[:space:]]*:[[:space:]]*([0-9]+)[[:space:]]");
      boost::smatch dim_match;
      if(!boost::regex_search(input, dim_match, dim_rgx)){
        throw custom_exception("Incorrect \"DIMENSION\" key.");
      }
      std::size_t dimension = std::stoul(dim_match[1].str());
    
      // Looking for a node coord section.
      boost::regex ews_rgx ("NODE_COORD_SECTION[[:space:]]*(.+)[[:space:]]*(EOF)?");
      boost::smatch ews_match;
      if(!boost::regex_search(input, ews_match, ews_rgx)){
        throw custom_exception("Incorrect \"NODE_COORD_SECTION\".");
      }

      // Parsing nodes.
      std::istringstream data (ews_match[1].str());
      for(std::size_t i = 0; i < dimension; ++i){
        index_t index;
        double x,y;
        data >> index >> x >> y;
        _locations.push_back({index, x, y});
      }
    }

    if(_locations.size() <= 1){
      throw custom_exception("at least two locations required!");
    }
  }

  std::vector<std::pair<double, double>> get_locations() const{
    std::vector<std::pair<double, double>> locs;
    std::transform(_locations.begin(), _locations.end(),
                   std::back_inserter(locs),
                   [](const auto& node){return std::make_pair(node.x, node.y);});
    return locs;
  }

  void get_route(const std::list<index_t>& tour,
                 rapidjson::Value& value,
                 rapidjson::Document::AllocatorType& allocator) const{
    rapidjson::Value route_array(rapidjson::kArrayType);
    for(auto const& step: tour){
      route_array
        .PushBack(rapidjson::Value(rapidjson::kArrayType)
                  .PushBack(_locations[step].x, allocator)
                  .PushBack(_locations[step].y, allocator),
                  allocator);
    }
    value.Swap(route_array);
  }

  void get_tour(const std::list<index_t>& tour,
                rapidjson::Value& value,
                rapidjson::Document::AllocatorType& allocator) const{
    rapidjson::Value tour_array(rapidjson::kArrayType);
    for(auto const& step: tour){
      // Using input index to describe locations for lon_lat syntax
      // and index provided in tsplib file.
      tour_array.PushBack(_locations[step].index, allocator);
    }
    value.Swap(tour_array);
  }

  void get_route_infos(const std::list<index_t>& tour,
                       rapidjson::Document& output) const{
  }
};

#endif
