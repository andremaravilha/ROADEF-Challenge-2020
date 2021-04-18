#include <cstddef>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <cxxtimer.hpp>
#include "utils.h"
#include "roadef2020.h"
#include "heur_greedy.h"
#include "heur_benders.h"
#include "mip.h"


/**
 * Team identification in the competition.
 */
const char* TEAM_ID = "S36";


/**
 * Utility class to parse command line arguments.
 */
class Options {

public:
    std::string executable;
    std::string instance_file;
    std::string output_file;
    std::string algorithm;
    int seed;
    double time_limit;
    int verbosity;
    bool show_help;
    bool show_team_id;

    Options(const char* executable_name) {
        executable = executable_name;
        show_help = false;
        show_team_id = false;
        verbosity = 0;
        time_limit = std::numeric_limits<double>::max();
        instance_file = "";
        output_file = "";
        algorithm = "benders";
        seed = 0;
    }

    bool parse(int argc, char** argv) {
        int idx = 1;
        while (idx < argc) {

            if (std::strcmp("-h", argv[idx]) == 0 || std::strcmp("--help", argv[idx]) == 0) {
                show_help = true;
                idx = idx + 1;

            } else if (std::strcmp("-name", argv[idx]) == 0 || std::strcmp("--name", argv[idx]) == 0) {
                show_team_id = true;
                idx = idx + 1;

            } else if (std::strcmp("-t", argv[idx]) == 0 || std::strcmp("--timelimit", argv[idx]) == 0) {
                time_limit = std::atof(argv[idx + 1]);
                idx = idx + 2;

            } else if (std::strcmp("-a", argv[idx]) == 0 || std::strcmp("--algorithm", argv[idx]) == 0) {
                algorithm = std::string(argv[idx + 1]);
                idx = idx + 2;

            } else if (std::strcmp("-p", argv[idx]) == 0 || std::strcmp("--instance", argv[idx]) == 0) {
                instance_file = std::string(argv[idx + 1]);
                idx = idx + 2;

            } else if (std::strcmp("-o", argv[idx]) == 0 || std::strcmp("--output", argv[idx]) == 0) {
                output_file = std::string(argv[idx + 1]);
                idx = idx + 2;

            } else if (std::strcmp("-s", argv[idx]) == 0 || std::strcmp("--seed", argv[idx]) == 0) {
                seed = std::atoi(argv[idx + 1]);
                idx = idx + 2;

            } else if (std::strcmp("-v", argv[idx]) == 0 || std::strcmp("--verbosity", argv[idx]) == 0) {
                verbosity = std::atoi(argv[idx + 1]);
                verbosity = std::min(verbosity, 3);
                verbosity = std::max(verbosity, 0);
                idx = idx + 2;
            } else {
                return false;
            }
        }
        return true;
    }

    void print_help() {
        std::printf("ROADEF/EURO Challenge 2020: Maintenance Planning Problem\n");
        std::printf("Usage:\n");
        std::printf("  %s [OPTIONS...]\n\n", executable.c_str());
        std::printf("  -h, --help              Show this help message and exit. Other options are\n"
                    "                          ignored if it is used.\n");
        std::printf("  -v, --verbosity VALUE   Display details about the optimization progress. Valid\n"
                    "                          values range from 0 to 3. The higher the value, the more\n"
                    "                          details are displayed on the screen. (Default is 0)\n");
        std::printf("  -name, --name           Display the team ID. If this is the only option, the\n"
                    "                          executable returns the team ID and quits.\n");
        std::printf("  -p, --instance FILE     Path to the file with data of an instance of the problem\n"
                    "                          to be loaded\n");
        std::printf("  -o, --output FILE       Path to the output file to write the result.\n");
        std::printf("  -s, --seed VALUE        Set the seed used to initialize the random number\n"
                    "                          generator. (Default is 0)\n");
        std::printf("  -t, --timelimit VALUE   Limit the maximum runtime (in seconds) of the program.\n"
                    "                          (Default is infinity)\n");
        std::printf("\n");
    }

};


/*
 * Main function.
 */
int main(int argc, char** argv) {

    // Start timer
    cxxtimer::Timer timer(true);

    try {

        // Parse command line arguments
        Options options(argv[0]);
        if (!options.parse(argc, argv)) {
            throw "Failed to parse command line arguments.";
        }

        if (options.show_help) {
            options.print_help();
            return EXIT_SUCCESS;
        }

        if (options.show_team_id) {
            std::cout << TEAM_ID << std::endl;
            if (argc <= 2) {
                return EXIT_SUCCESS;
            }
        }

        // Algorithm's input parameters
        bool verbose = (options.verbosity == 3);
        int seed = options.seed;
        double time_limit = options.time_limit;

        // Load problem data
        orcs::roadef2020::Problem problem(options.instance_file);

        // Get a starting solution
        //auto [schedule, objective, mean_risk, expected_excess] = orcs::roadef2020::greedy(problem);

//        auto [schedule, objective, mean_risk, expected_excess] = orcs::roadef2020::mip(problem,
//                time_limit, timer, seed, verbose);

        auto [schedule, objective, mean_risk, expected_excess] = orcs::roadef2020::benders_heuristic(problem,
                time_limit, timer, seed, verbose);

        // Stop timer
        timer.stop();

        // If verbose mode is enable, print a summary
        if (options.verbosity > 0) {

            // Check feasibility of the schedule
            int infeasibilities = schedule.count_infeasibilities();

            // Get elapsed time
            double elapsed_time = timer.count<std::chrono::milliseconds>() / 1000.0;

            // Print output for verbose (1)
            if (options.verbosity == 1) {
                if (infeasibilities > 0) {
                    std::cout << infeasibilities << " violated constraints." << std::endl;
                }
                std::printf("%s %.6lf %.6lf\n", ((infeasibilities == 0) ? "FEASIBLE" : "INFEASIBLE"), objective, elapsed_time);
            }

            if (options.verbosity == 2 || options.verbosity == 3) {
                std::printf("FEASIBILITY: %s%s\n"
                            "OBJECTIVE: %.6lf\n"
                            "MEAN.RISK: %.6lf\n"
                            "EXPECTED.EXCESS: %.6lf\n"
                            "RUNTIME.S: %.6lf\n",
                            ((infeasibilities == 0) ? "Feasible" : "Infeasible"),
                            ((infeasibilities == 0) ? "" : orcs::utils::format(" (%d violated constraints)", infeasibilities).c_str()),
                            objective, mean_risk, expected_excess, elapsed_time);
            }
        }

        // Export solution to file
        std::ofstream ofs(options.output_file);
        schedule.export_ROADEF2020(ofs);
        ofs.close();

    } catch (...) {
        std::cerr << "Unexpected error." << std::endl;
        std::cerr << "Type the following command for a correct usage." << std::endl;
        std::cerr << argv[0] << " --help" << std::endl << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}