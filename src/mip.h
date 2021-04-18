#ifndef ROADEF2020_MIP_H
#define ROADEF2020_MIP_H

#include <cxxtimer.hpp>
#include "roadef2020.h"

namespace orcs::roadef2020 {

    /**
     * Implements and solve the grid operation-based outage maintenance problem through a MIP formulation using Gurobi.
     * @param problem The instance of the problem to solve.
     * @return A tuple with four elements in the following order: the schedule, the value of objective function, the
     * mean fisk, and the expected excess.
     */
    std::tuple<Schedule, double, double, double> mip(const Problem &problem, double time_limit,
            const cxxtimer::Timer& timer, int seed, bool verbose);

}

#endif