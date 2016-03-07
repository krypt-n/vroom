/*

This file is part of VROOM.

Copyright (c) 2015-2016, Julien Coupey.
All rights reserved (see LICENSE).

*/

#include "tsp_strategy.h"

void solve_symmetric_tsp(const cl_args_t& cl_args){
  // Store timings.
  timing_t computing_times;

  // Building problem object with embedded table request.
  auto start_problem_build = std::chrono::high_resolution_clock::now();

  BOOST_LOG_TRIVIAL(info) 
    << "[Matrix] Start matrix computing and problem loading.";

  tsp symmetric_tsp (cl_args);

  auto end_problem_build = std::chrono::high_resolution_clock::now();

  computing_times.matrix_loading =
    std::chrono::duration_cast<std::chrono::milliseconds>
    (end_problem_build - start_problem_build).count();

  BOOST_LOG_TRIVIAL(info) << "[Matrix] Done, took "
                          << computing_times.matrix_loading << " ms.";

  // Applying heuristic.
  auto start_heuristic = std::chrono::high_resolution_clock::now();
  BOOST_LOG_TRIVIAL(info) 
    << "[Heuristic] Start heuristic on symmetrized problem.";

  std::unique_ptr<heuristic> christo_h = std::make_unique<christo_heuristic>();
  std::list<index_t> christo_sol
    = christo_h->build_solution(symmetric_tsp);
  distance_t christo_cost = symmetric_tsp.cost(christo_sol);

  auto end_heuristic = std::chrono::high_resolution_clock::now();

  computing_times.heuristic = 
    std::chrono::duration_cast<std::chrono::milliseconds>
    (end_heuristic - start_heuristic).count();

  BOOST_LOG_TRIVIAL(info) << "[Heuristic] Done, took "
                          << computing_times.heuristic << " ms.";

  BOOST_LOG_TRIVIAL(info) << "[Heuristic] Symmetric solution cost is "
                          << christo_cost << ".";


  // Local search on symmetric problem.
  // Applying deterministic, fast local search to improve the current
  // solution in a small amount of time. All possible moves for the
  // different neighbourhoods are performed, stopping when reaching a
  // local minima.
  auto start_sym_local_search = std::chrono::high_resolution_clock::now();
  BOOST_LOG_TRIVIAL(info) 
    << "[Local search] Start local search on symmetric problem.";
  BOOST_LOG_TRIVIAL(info) 
    << "[Local search] Using " << cl_args.nb_threads << " thread(s).";

  local_search sym_ls (symmetric_tsp.get_matrix(),
                       christo_sol,
                       cl_args.nb_threads);

  distance_t sym_two_opt_gain = 0;
  distance_t sym_relocate_gain = 0;
  distance_t sym_or_opt_gain = 0;

  do{
    // All possible 2-opt moves.
    sym_two_opt_gain = sym_ls.perform_all_two_opt_steps();

    // All relocate moves.
    sym_relocate_gain = sym_ls.perform_all_relocate_steps();

    // All or-opt moves.
    sym_or_opt_gain = sym_ls.perform_all_or_opt_steps();
  }while((sym_two_opt_gain > 0) 
         or (sym_relocate_gain > 0) 
         or (sym_or_opt_gain > 0));

  // Default for first input location.
  index_t first_loc_index = 0;

  std::list<index_t> current_sol = sym_ls.get_tour(first_loc_index);
  auto current_cost = symmetric_tsp.cost(current_sol);

  auto end_sym_local_search = std::chrono::high_resolution_clock::now();

  auto sym_local_search_duration 
    = std::chrono::duration_cast<std::chrono::milliseconds>
    (end_sym_local_search - start_sym_local_search).count();
  BOOST_LOG_TRIVIAL(info) << "[Local search] Done, took "
                          << sym_local_search_duration << " ms.";

  BOOST_LOG_TRIVIAL(info) << "[Local search] Symmetric solution cost is now "
                          << current_cost
                          << " (" 
                          << std::fixed << std::setprecision(2)
                          << 100 *(((double) current_cost) / christo_cost - 1)
                          << "%).";

  computing_times.local_search 
    = sym_local_search_duration;

  logger log (cl_args);
  log.write_solution(symmetric_tsp, current_sol, computing_times);
}
