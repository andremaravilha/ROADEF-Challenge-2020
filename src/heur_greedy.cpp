#include "heur_greedy.h"
#include <vector>
#include <list>
#include <algorithm>


std::tuple<orcs::roadef2020::Schedule, double, double, double> orcs::roadef2020::greedy(const Problem& problem) {

    constexpr double EPS = 1E-6;

    // Create a set of possible starting times for each intervention
    std::list<greedy_entry_> entries;
    std::vector<int> slots;
    slots.reserve(problem.count_interventions());

    for (int i = 0; i < problem.count_interventions(); ++i) {
        slots.push_back(problem.t_max(i) + 1);

        for (int start = 0; start <= problem.t_max(i); ++start) {
            greedy_entry_ entry{i, start, 0.0};
            int end = start + problem.delta(i, start) - 1;

            // Compute a heuristic weight for the entry
            for (int t = 0; t < problem.T(); ++t) {
                double sum = 0.0;
                double max = 0.0;
                for (int s = 0; s < problem.count_scenarios(t); ++s) {
                    sum += problem.risk(i, t, start, s);
                    max = std::max(max, problem.risk(i, t, start, s));
                }

                entry.weight += (problem.alpha() * (sum / problem.count_scenarios(t))) + (1.0 - problem.alpha()) * max;
            }

            entries.push_back(std::move(entry));
        }
    }

    // Build a schedule
    Schedule schedule(problem);
    int last_i = -1;

    // Sort entries according to <slots, weight>
    auto f_sort = [&slots](const greedy_entry_& first, const greedy_entry_& second) -> bool {
        return ((slots[first.i] < slots[second.i]) || (slots[first.i] == slots[second.i] && first.weight < second.weight));
    };

    // Filter entries with the last used intervention and that violate constraints
    auto f_filter = [&EPS, &problem, &slots, &schedule, &last_i](const decltype(entries)::value_type& entry) -> bool {

        // Check if intervention has already been set
        if (entry.i == last_i) {
            slots[entry.i] -= 1;
            return true;
        }

        // Check for exclusions feasibility
        if (last_i != -1) {

            int i = entry.i;
            int i_start = entry.ts;
            int i_end = i_start + problem.delta(i, i_start) - 1;

            int j = last_i;
            int j_start = schedule.get(j);
            int j_end = j_start + problem.delta(j, j_start) - 1;

            // If intervals [i_start, i_end] and [j_start, j_end] overlaps
            if (i_start <= j_end && j_start <= i_end) {
                int t_start = std::max(i_start, j_start);
                int t_end = std::min(i_end, j_end);
                for (int t = t_start; t <= t_end; ++t) {
                    if (problem.is_simultaneity_allowed(i, j, problem.which_season(t))) {
                        slots[entry.i] -= 1;
                        return true;
                    }
                }
            }
        }

//        // Check for workload feasibility
//        for (int r = 0; r < problem.count_resources(); ++r) {
//             for (int t = 0; t < problem.T(); ++t) {
//             //for (int t = entry.ts; t <= entry.ts + problem.delta(entry.i, entry.ts) - 1; ++t) {
//                if (schedule.workload(t, r) + problem.workload(entry.i, r, t, entry.ts) > problem.resource_max(r, t) + EPS) {
//                    slots[entry.i] -= 1;
//                    return true;
//                }
//            }
//        }

        // Entry does not violate constraints
        return false;
    };

    // Build the schedule
    bool failed = false;
    do {

        // Filter entries
        entries.remove_if(f_filter);

        if (!entries.empty()) {

            // Sort entries
            entries.sort(f_sort);

            failed = true;
            for (const auto& entry : entries) {

                // Check for workload feasibility
                bool feasible_workload = true;
                for (int r = 0; r < problem.count_resources() && feasible_workload; ++r) {
                     //for (int t = 0; t < problem.T(); ++t) {
                     for (int t = entry.ts; t <= entry.ts + problem.delta(entry.i, entry.ts) - 1 && feasible_workload; ++t) {
                        if (schedule.workload(t, r) + problem.workload(entry.i, r, t, entry.ts) > problem.resource_max(r, t) + EPS) {
                            feasible_workload = false;
                        }
                    }
                }

                // If a feasible entry is found, use it
                if (feasible_workload) {
                    schedule.set(entry.i, entry.ts);
                    last_i = entry.i;
                    failed = false;
                    break;
                }

            }
        }

    } while (!failed && !entries.empty());

    auto [objective, mean_risk, expected_excess] = schedule.evaluation();
    return {schedule, objective, mean_risk, expected_excess};
}