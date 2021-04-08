#include "heur_benders.h."
#include <cmath>
#include <algorithm>
#include <gurobi_c++.h>


std::tuple<orcs::roadef2020::Schedule, double, double, double> orcs::roadef2020::benders_heuristic(const Problem &problem,
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

    // Auxiliary data used to sort risks when solving sub-problems
    int* risk_idx = new int[problem.max_scenarios()];
    double* risk_val = new double[problem.max_scenarios()];
    const auto sort_risk = [&risk_val](int first, int second) { return risk_val[first] < risk_val[second]; };

    // Solve the problem with Gurobi solver
    GRBEnv* env = nullptr;
    try {

        // Gurobi environment
        env = new GRBEnv();

        // Master problem
        GRBModel master(*env);

        // Gurobi parameters
        master.getEnv().set(GRB_IntParam_LogToConsole, 0);
        master.getEnv().set(GRB_IntParam_OutputFlag, 0);
        master.getEnv().set(GRB_IntParam_Threads, 1);
        master.getEnv().set(GRB_IntParam_Seed, seed);
        master.getEnv().set(GRB_DoubleParam_TimeLimit, GRB_INFINITY);
        master.getEnv().set(GRB_DoubleParam_NodeLimit, GRB_INFINITY);

        // Allocate memory for decision variables and instantiate them
        GRBVar **x = new GRBVar *[problem.count_interventions()];
        for (int i = 0; i < problem.count_interventions(); ++i) {
            x[i] = new GRBVar[problem.t_max(i) + 1];
            for (int ts = 0; ts <= problem.t_max(i); ++ts) {
                x[i][ts] = master.addVar(0.0, 1.0, 0.0, GRB_BINARY);
            }
        }

        GRBVar *eta = new GRBVar[problem.T()];
        for (int t = 0; t < problem.T(); ++t) {
            eta[t] = master.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS);
        }

        master.update();

        // Objective function (Eq. 1)
        GRBLinExpr objective1 = 0.0;
        GRBLinExpr objective2 = 0.0;
        for (int t = 0; t < problem.T(); ++t) {
            objective2 += (1.0 / problem.T()) * eta[t];
            for (int s = 0; s < problem.count_scenarios(t); ++s) {
                for (int i = 0; i < problem.count_interventions(); ++i) {
                    for (int ts = 0; ts <= problem.t_max(i); ++ts) {
                        objective1 +=  (1.0 / (problem.T() * problem.count_scenarios(t))) * problem.risk(i, t, ts, s) * x[i][ts];
                    }
                }
            }
        }

        GRBLinExpr objective = (problem.alpha() * objective1) + ((1.0 - problem.alpha()) * objective2);
        master.setObjective(objective, GRB_MINIMIZE);

        // Constraints (Eq. 2)
        for (int i = 0; i < problem.count_interventions(); ++i) {
            GRBLinExpr expr = 0.0;
            for (int ts = 0; ts <= problem.t_max(i); ++ts) {
                expr += x[i][ts];
            }
            master.addConstr(expr == 1.0);
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
                master.addConstr(expr <= problem.resource_max(c, t));
                master.addConstr(expr >= problem.resource_min(c, t));
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

                master.addConstr(expr1 + expr2 <= 1.0);
            }
        }

        // Benders iterative process
        int iter = 0;
        double ub = 1e10;
        double lb = -1e10;
        bool stop = false;

        if (verbose) {
            std::printf("---------------------------------------------------------------------------------\n");
            std::printf("|   Iter. |      Current |  Upper bound |  Lower bound | Gap (%%) |  Runtime (s) |\n");
            std::printf("---------------------------------------------------------------------------------\n");
        }

        while (!stop) {

            // Increment iterations counter
            ++iter;

            // Solve master model
            master.getEnv().set(GRB_DoubleParam_TimeLimit, time_limit - (timer.count<std::chrono::milliseconds>() / 1000.0));
            master.optimize();

            // Update lower bound
            lb = master.get(GRB_DoubleAttr_ObjVal);

            // Solve sub-problems
            double obj_subproblem = 0.0;
            for (int t = 0; t < problem.T(); ++t) {

                // Find the quantile
                for (int s = 0; s < problem.count_scenarios(t); ++s) {
                    risk_idx[s] = s;
                    risk_val[s] = 0.0;
                    for (int i = 0; i < problem.count_interventions(); ++i) {
                        for (int ts = 0; ts <= problem.t_max(i); ++ts) {
                            if (x[i][ts].get(GRB_DoubleAttr_X) > 0.5) {
                                risk_val[s] += problem.risk(i, t, ts, s);
                            }
                        }
                    }
                }

                std::nth_element(risk_idx, risk_idx + problem.quantile_index(t), risk_idx + problem.count_scenarios(t), sort_risk);
                int s_quantile = risk_idx[problem.quantile_index(t)];

                // Add cut
                GRBLinExpr cut = 0.0;
                for (int i = 0; i < problem.count_interventions(); ++i) {
                    for (int ts = 0; ts <= problem.t_max(i); ++ts) {
                        cut += problem.risk(i, t, ts, s_quantile) * x[i][ts];
                        for (int s = 0; s < problem.count_scenarios(t); ++s) {
                            cut -= (1.0 / problem.count_scenarios(t)) * problem.risk(i, t, ts, s) * x[i][ts];
                        }
                    }
                }

                double obj_subproblem_t = cut.getValue();
                if (obj_subproblem_t > 1e-6) {
                    obj_subproblem += obj_subproblem_t;
                    master.addConstr(eta[t] >= cut);
                }
            }

            // Compute new upper bound
            bool improved = false;
            double new_ub = problem.alpha() * objective1.getValue() + ((1.0 * problem.alpha()) / problem.T()) * obj_subproblem;
            if (new_ub < ub - 1e-6) {

                // Update upper bound
                ub = new_ub;
                improved = true;

                // Update incumbent schedule
                if (master.get(GRB_IntAttr_SolCount) > 0) {
                    found_solution = true ;
                    for (int i = 0; i < problem.count_interventions(); ++i) {
                        for (int ts = 0; ts <= problem.t_max(i); ++ts) {
                            if (x[i][ts].get(GRB_DoubleAttr_X) > 0.5) {
                                assignment[i] = ts;
                            }
                        }
                    }
                }
            }

            if (verbose) {
                std::printf("| %c %5d | %12.5lf | %12.5lf | %12.5lf | %7.2lf | %12.2lf |\n",
                            (improved ? '*' : ' '), iter, new_ub, ub, lb, 100.0 * ((ub - lb) / lb),
                            timer.count<std::chrono::milliseconds>() / 1000.0);
            }

            // Stopping criteria
            if (ub - lb < 1e-6 || (timer.count<std::chrono::milliseconds>() / 1000.0) >= time_limit) {
                stop = true;
            }

        }

        if (verbose) {
            std::printf("---------------------------------------------------------------------------------\n\n");
        }

        // Deallocate resources
        delete[] eta;
        for (int i = 0; i < problem.count_interventions(); ++i) {
            delete[] x[i];
        }
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

        delete[] risk_idx;
        delete[] risk_val;

        delete[] assignment;

        // Re-throw the exception
        throw;
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

    delete[] risk_idx;
    delete[] risk_val;

    // Get the best schedule found
    Schedule schedule(problem);
    for (int i = 0; i < problem.count_interventions(); ++i) {
        schedule.set(i, assignment[i]);
    }

    delete[] assignment;

    // Return the best schedule found
    auto [objective, mean_risk, expected_excess] = schedule.evaluation();
    return {schedule, objective, mean_risk, expected_excess};

}