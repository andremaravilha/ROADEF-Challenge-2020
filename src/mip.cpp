#include "mip.h."
#include <cmath>
#include <algorithm>
#include <gurobi_c++.h>


std::tuple<orcs::roadef2020::Schedule, double, double, double> orcs::roadef2020::mip(const Problem &problem,
        double time_limit, const cxxtimer::Timer& timer, int seed, bool verbose) {

    // Data about the best schedule found
    int* assignment = new int[problem.count_interventions()];
    bool found_solution = false;

    // Compute Big-M valued used by the MIP model
    double** M = new double*[problem.T()];
    for (int t = 0; t < problem.T(); ++t) {
        M[t] = new double[problem.count_scenarios(t)];
        for (int s = 0; s < problem.count_scenarios(t); ++s) {
            M[t][s] = 0.0;
            for (int i = 0; i < problem.count_interventions(); ++i) {
                double aux = 0.0;
                for (int ts = 0; ts <= problem.t_max(i); ++ts) {
                    aux = std::max(problem.risk(i, t, ts, s), aux);
                }
                M[t][s] += aux;
            }
        }
    }

    // Solve the problem with Gurobi solver
    GRBEnv* env = nullptr;
    try {

        // Gurobi environment
        env = new GRBEnv();

        // Master problem
        GRBModel model(*env);

        // Gurobi parameters
        model.getEnv().set(GRB_IntParam_LogToConsole, (verbose ? 1 : 0));
        model.getEnv().set(GRB_IntParam_OutputFlag, (verbose ? 1 : 0));
        model.getEnv().set(GRB_IntParam_Threads, 1);
        model.getEnv().set(GRB_IntParam_Seed, seed);
        model.getEnv().set(GRB_DoubleParam_TimeLimit, GRB_INFINITY);
        model.getEnv().set(GRB_DoubleParam_NodeLimit, GRB_INFINITY);

        // Allocate memory for decision variables and instantiate them
        GRBVar **x = new GRBVar *[problem.count_interventions()];
        for (int i = 0; i < problem.count_interventions(); ++i) {
            x[i] = new GRBVar[problem.t_max(i) + 1];
            for (int ts = 0; ts <= problem.t_max(i); ++ts) {
                x[i][ts] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
            }
        }

        GRBVar *z = new GRBVar[problem.T()];
        GRBVar *f = new GRBVar[problem.T()];
        GRBVar **y = new GRBVar *[problem.T()];
        for (int t = 0; t < problem.T(); ++t) {
            z[t] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
            f[t] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
            y[t] = new GRBVar[problem.count_scenarios(t)];
            for (int s = 0; s < problem.count_scenarios(t); ++s) {
                y[t][s] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
            }
        }

        model.update();

        // Objective function (Eq. 1)
        GRBLinExpr objective1 = 0.0;
        GRBLinExpr objective2 = 0.0;
        for (int t = 0; t < problem.T(); ++t) {
            objective2 += (1.0 / problem.T()) * f[t];
            for (int s = 0; s < problem.count_scenarios(t); ++s) {
                for (int i = 0; i < problem.count_interventions(); ++i) {
                    for (int ts = 0; ts <= problem.t_max(i); ++ts) {
                        objective1 +=  (1.0 / (problem.T() * problem.count_scenarios(t))) * problem.risk(i, t, ts, s) * x[i][ts];
                    }
                }
            }
        }

        GRBLinExpr objective = (problem.alpha() * objective1) + ((1.0 - problem.alpha()) * objective2);
        model.setObjective(objective, GRB_MINIMIZE);

        // Constraints (Eq. 2)
        for (int i = 0; i < problem.count_interventions(); ++i) {
            GRBLinExpr expr = 0.0;
            for (int ts = 0; ts <= problem.t_max(i); ++ts) {
                expr += x[i][ts];
            }
            model.addConstr(expr == 1.0);
        }

        // Constraints (Eq. 3-4)
        for (int c = 0; c < problem.count_resources(); ++c) {
            for (int t = 0; t < problem.T(); ++t) {
                GRBLinExpr expr = 0.0;
                for (int i = 0; i < problem.count_interventions(); ++i) {
                    for (int ts = 0; ts <= problem.t_max(i); ++ts) {
                        expr += problem.workload(i, c, t, ts) * x[i][ts];
                    }
                }
                model.addConstr(expr <= problem.resource_max(c, t));
                model.addConstr(expr >= problem.resource_min(c, t));
            }
        }

        // Constraints (Eq. 5)
        for (const auto&[i, j, season, name] : problem.exclusions()) {
            for (auto t : problem.season(season)) {

                GRBLinExpr expr1 = 0.0;
                for (int ts = 0; ts <= problem.t_max(i); ++ts) {
                    if (ts <= t && t <= ts + problem.delta(i, ts) - 1) {
                        expr1 += x[i][ts];
                    }
                }

                GRBLinExpr expr2 = 0.0;
                for (int ts = 0; ts <= problem.t_max(j); ++ts) {
                    if (ts <= t && t <= ts + problem.delta(j, ts) - 1) {
                        expr2 += x[j][ts];
                    }
                }

                model.addConstr(expr1 + expr2 <= 1.0);
            }
        }

        // Constraints (Eq. 6)
        for (int t = 0; t < problem.T(); ++t) {
            GRBLinExpr expr = 0.0;
            for (int s = 0; s < problem.count_scenarios(t); ++s) {
                expr += y[t][s];
            }
            model.addConstr(expr >= problem.quantile_index(t) + 1);
        }

        // Constraints (Eq. 7)
        for (int t = 0; t < problem.T(); ++t) {
            for (int s = 0; s < problem.count_scenarios(t); ++s) {
                GRBLinExpr expr = 0.0;
                for (int i = 0; i < problem.count_interventions(); ++i) {
                    for (int ts = 0; ts <= problem.t_max(i); ++ts) {
                        expr += problem.risk(i, t, ts, s) * x[i][ts];
                    }
                }
                model.addConstr(z[t] >= expr - M[t][s] * (1.0 - y[t][s]));
            }
        }

        // Constraints (Eq. 8)
        for (int t = 0; t < problem.T(); ++t) {
            GRBLinExpr expr = 0.0;
            for (int s0 = 0; s0 < problem.count_scenarios(t); ++s0) {
                for (int i = 0; i < problem.count_interventions(); ++i) {
                    for (int ts = 0; ts <= problem.t_max(i); ++ts) {
                        expr += problem.risk(i, t, ts, s0) * x[i][ts];
                    }
                }
            }
            model.addConstr(f[t] >= z[t] - (1.0 / problem.count_scenarios(t)) * expr);
        }

        // Solve model model
        model.getEnv().set(GRB_DoubleParam_TimeLimit, time_limit - (timer.count<std::chrono::milliseconds>() / 1000.0));
        model.optimize();

        // Update incumbent schedule
        if (model.get(GRB_IntAttr_SolCount) > 0) {
            for (int i = 0; i < problem.count_interventions(); ++i) {
                for (int ts = 0; ts <= problem.t_max(i); ++ts) {
                    if (x[i][ts].get(GRB_DoubleAttr_X) > 0.5) {
                        assignment[i] = ts;
                    }
                }
            }
        }

        // Deallocate resources
        for (int t = 0; t < problem.T(); ++t) {
            delete[] y[t];
        }

        for (int i = 0; i < problem.count_interventions(); ++i) {
            delete[] x[i];
        }

        delete[] y;
        delete[] f;
        delete[] z;
        delete[] x;

    } catch (...) {

        // Deallocate resources
        if (env != nullptr) {
            delete env;
            env = nullptr;
        }

        for (int t = 0; t < problem.T(); ++t) {
            delete[] M[t];
        }

        delete[] M;
        delete[] assignment;

        // Re-throw the exception
        throw;
    }

    // Get the best schedule found
    Schedule schedule(problem);
    for (int i = 0; i < problem.count_interventions(); ++i) {
        schedule.set(i, assignment[i]);
    }

    // Deallocate resources
    if (env != nullptr) {
        delete env;
        env = nullptr;
    }

    for (int t = 0; t < problem.T(); ++t) {
        delete[] M[t];
    }

    delete[] M;
    delete[] assignment;

    // Return the best schedule found
    auto [objective, mean_risk, expected_excess] = schedule.evaluation();
    return {schedule, objective, mean_risk, expected_excess};

}