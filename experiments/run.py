# MIT License
#
# Copyright (c) 2020 André L. Maravilha
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

'''
This script runs the experiments designed for performance evaluation of the solution
strategies proposed for the ROADEF/EURO Challenge 2020: Maintenance Planning Problem.
'''

import argparse
import os
import subprocess
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from multiprocessing import RLock
from solution_checker import check_solution


###################################################################################################
# Re-entrant lock object used to control the access to the output file.
#
lock = RLock()


###################################################################################################
# Initializes the CSV output file and writes the header.
#
def initialize_results_file(filename):
    with open(filename, "w") as f:
        f.write("INSTANCE,STRATEGY,SEED,STATUS,MEAN.RISK,EXPECTED.EXCESS,OBJECTIVE,OPT.GAP,ITERATIONS,RUNTIME.S," +
                "CHECK.MEAN.RISK,CHECK.EXPECTED.EXCESS,CHECK.OBJECTIVE,CHECK.FEASIBILITY\n")


###################################################################################################
# Writes the results of an entry in the output file.
#
def write_results(filename, results):
    content = ",".join([results["instance"], results["strategy"], results["seed"],
                        results["status"], results["mean-risk"], results["expected-excess"],
                        results["objective"], results["opt-gap"], results["iterations"],
                        results["runtime"], results["check-mean-risk"], results["check-expected-excess"],
                        results["check-objective"], results["check-feasibility"]])

    f = open(filename, "a")
    f.write(content + "\n")
    f.close()


###################################################################################################
# Runs an entry.
#
def run_entry(instance, strategy, data, seed):

    # Keep the result
    results = dict()
    results["instance"] = instance
    results["strategy"] = strategy
    results["seed"] = str(seed)
    results["status"] = "ERROR"
    results["mean-risk"] = ""
    results["expected-excess"] = ""
    results["objective"] = ""
    results["opt-gap"] = ""
    results["iterations"] = ""
    results["runtime"] = ""

    try:

        # Create the command to run
        command = ([data["command"], "--seed", str(seed)] + data["strategies"][strategy] +
                   ["--instance", data["instances"][instance]] +
                   ["--output", data["solutions-path"] + strategy + "-" + instance + "-{0:0>3d}.txt".format(seed)])

        # Run the command
        output = subprocess.run(command, stdout=subprocess.PIPE, universal_newlines=True)

        # Check if the command run without errors
        if output.returncode == 0:

            # Get the output as a string
            str_output = [x for x in str(output.stdout).strip('\n\t ').split('\n')]

            status = "ERROR"
            objective = ""
            optgap = ""
            meanrisk = ""
            expectedexcess = ""
            iterations = ""
            runtime = ""
            check_objective = ""
            check_meanrisk = ""
            check_expectedexcess = ""
            check_feasibility = ""

            for line in str_output:
                values = line.strip('\n\t ').split(':')

                if values[0].strip('\n\t ') == "FEASIBILITY":
                    status = values[1].strip('\n\t ')

                if values[0].strip('\n\t ') == "OBJECTIVE":
                    objective = values[1].strip('\n\t ')

                if values[0].strip('\n\t ') == "OPT.GAP":
                    optgap = values[1].strip('\n\t ')

                if values[0].strip('\n\t ') == "MEAN.RISK":
                    meanrisk = values[1].strip('\n\t ')

                if values[0].strip('\n\t ') == "EXPECTED.EXCESS":
                    expectedexcess = values[1].strip('\n\t ')

                if values[0].strip('\n\t ') == "ITERATIONS":
                    iterations = values[1].strip('\n\t ')

                if values[0].strip('\n\t ') == "RUNTIME.S":
                    runtime = values[1].strip('\n\t ')

            if status == "Feasible":
                results["status"] = "FEASIBLE"
                obj_total, obj_1, obj_2, violations = check_solution(data["instances"][instance],
                                                                     data["solutions-path"] + strategy + "-" + instance + "-{0:0>3d}.txt".format(seed))
                check_objective = "{:.6f}".format(obj_total)
                check_meanrisk = "{:.6f}".format(obj_1)
                check_expectedexcess = "{:.6f}".format(obj_2)
                check_feasibility = ("FEASIBLE" if len(violations) == 0 else "INFEASIBLE")

            else:
                results["status"] = "INFEASIBLE"

            results["objective"] = objective
            results["opt-gap"] = optgap
            results["mean-risk"] = meanrisk
            results["expected-excess"] = expectedexcess
            results["iterations"] = iterations
            results["runtime"] = runtime
            results["check-objective"] = check_objective
            results["check-mean-risk"] = check_meanrisk
            results["check-expected-excess"] = check_expectedexcess
            results["check-feasibility"] = check_feasibility

    except:
        results["status"] = "ERROR"
        results["objective"] = ""
        results["opt-gap"] = ""
        results["mean-risk"] = ""
        results["expected-excess"] = ""
        results["iterations"] = ""
        results["runtime"] = ""
        results["check-objective"] = ""
        results["check-mean-risk"] = ""
        results["check-expected-excess"] = ""
        results["check-feasibility"] = ""

    # Write results in the output file and log the progress
    with lock:

        # Write results
        write_results(data["output-file"], results)

        # Progress log
        data["progress"] += 1
        total_entries = len(data["instances"]) * len(data["strategies"]) * len(data["seeds"])
        progress = (data["progress"] / total_entries) * 100

        print("[{:3} of {:3} ({:6.2f}%) completed] {:15} -> {:16} -> {:4} -> {:8}"
              .format(data["progress"], total_entries, progress, strategy, instance, seed, results["status"]))

        sys.stdout.flush()


###################################################################################################
# Main function
#
def main():

    # Configure argument parser
    parser = argparse.ArgumentParser(description="Perform the computational experiments.")
    parser.add_argument("-c", "--continue", dest="continue_previous", action="store_true",
                        help="Resume the experiment from a previous state.")

    # Parse input arguments
    args = parser.parse_args()

    # Map used to store data necessary to run the experiments
    data = dict()
    data["command"] = "../build/challengeRTE"   # program
    data["threads"] = 1                         # entries to solve simultaneously
    data["solutions-path"] = "./solutions/"     # file to write the results
    data["output-file"] = "results.csv"         # file to write the results
    data["instances-path"] = "./instances"      # base path to instances
    data["instances"] = dict()                  # path to each instance
    data["strategies"] = dict()                 # settings of each solution strategy
    data["seeds"] = None                        # seeds
    data["progress"] = 0                        # num. of entries finished

    # Set the seeds used at each repetition of (instance, algorithm)
    data["seeds"] = [0]

    # Load the list of instances
    instances = [file.replace(".json", "") for file in os.listdir(data["instances-path"]) if file.endswith(".json")]

    for instance in instances:
        data["instances"][instance] = data["instances-path"] + "/" + instance + ".json"

    # Algorithms' settings
    params_common_15 = "--verbosity 2 --timelimit 900".split(" ")   # 15 minutes
    params_common_90 = "--verbosity 2 --timelimit 5400".split(" ")  # 90 minutes

    #params_greedy     = params_common_90 + "--algorithm greedy".split(" ")
    #params_mip_15     = params_common_15 + "--algorithm mip".split(" ")
    #params_mip_90     = params_common_90 + "--algorithm mip".split(" ")
    #params_benders_15 = params_common_15 + "--algorithm benders".split(" ")
    #params_benders_90 = params_common_90 + "--algorithm benders".split(" ")
    params_benders_15 = params_common_15
    params_benders_90 = params_common_90

    #data["strategies"]["greedy"] = params_greedy
    #data["strategies"]["mip-15"] = params_mip_15
    #data["strategies"]["mip-90"] = params_mip_90
    data["strategies"]["benders-15"] = params_benders_15
    data["strategies"]["benders-90"] = params_benders_90

    # Create and initialize the output file
    Path(data["solutions-path"]).mkdir(parents=True, exist_ok=True)
    if not os.path.exists(data["output-file"]) or not args.continue_previous:
        initialize_results_file(data["output-file"])

    # Check entries already solved (instance, strategy, seed)
    entries_completed = []
    if args.continue_previous:
        if os.path.exists(data["output-file"]):
            with open(data["output-file"], "r") as f:
                lines = f.readlines()[1:]
                for line in lines:
                    values = line.strip().split(",")
                    print(values)
                    if len(values) > 0:
                        entries_completed.append((values[0], values[1], int(values[2])))

    data["progress"] = len(entries_completed)

    # Run each entry (instance, strategy, seed)
    with ThreadPoolExecutor(max_workers=data["threads"]) as executor:
        for seed in data["seeds"]:
            for instance in list(data["instances"].keys()):
                for strategy in list(data["strategies"].keys()):
                    if entries_completed.count((instance, strategy, seed)) == 0:
                        executor.submit(run_entry, instance, strategy, data, seed)


###################################################################################################
# Main statements
#
if __name__ == "__main__":
    main()

