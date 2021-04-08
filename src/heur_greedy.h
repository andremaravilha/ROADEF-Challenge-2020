#ifndef ROADEF2020_HEUR_GREEDY_H
#define ROADEF2020_HEUR_GREEDY_H

#include "roadef2020.h"

namespace orcs::roadef2020 {

    struct greedy_entry_ {
        int i;
        int ts;
        double weight;
    };

    /**
     * Greedy heuristic for solving the grid operation-based outage maintenance planning.
     * @param problem The instance of the problem to solve.
     * @return A tuple with four elements in the following order: the schedule, the value of objective function, the
     * mean fisk, and the expected excess.
     */
    std::tuple<Schedule, double, double, double> greedy(const Problem &problem);

}

#endif