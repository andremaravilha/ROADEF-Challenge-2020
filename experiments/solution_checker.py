"""
This is an adaptation from the "RTE_ChallengeROADEF2020_checker" file, provided on the competition's official repository
and available on https://github.com/rte-france/challenge-roadef-2020/. In this adaptation, instead of displaying the
information on the screen, the data are returned to be used elsewhere.
"""
import numpy as np
import json


# Global variables
RESOURCES_STR = 'Resources'
SEASONS_STR = 'Seasons'
INTERVENTIONS_STR = 'Interventions'
EXCLUSIONS_STR = 'Exclusions'
T_STR = 'T'
SCENARIO_NUMBER = 'Scenarios_number'
RESOURCE_CHARGE_STR = 'workload'
TMAX_STR = 'tmax'
DELTA_STR = 'Delta'
MAX_STR = 'max'
MIN_STR = 'min'
RISK_STR = 'risk'
START_STR = 'start'
QUANTILE_STR = 'Quantile'
ALPHA_STR = 'Alpha'
FEASIBILITY = 'Feasibility'
VIOLATIONS = 'Violations'


def read_instance(filename: str):
    """Read a json file and return data as a dict object"""

    with open(filename, 'r') as f:
        instance = json.load(f)

    return instance


def read_solution(instance: dict, solution_filename: str):
    """Read a txt formatted Solution file, and store the solution information in Instance"""

    # Load interventions
    instance[FEASIBILITY] = True
    instance[VIOLATIONS] = []
    interventions = instance[INTERVENTIONS_STR]

    # Read file line by line, and store starting time value (no checks yet, except format and duplicate)
    with open(solution_filename, 'r') as solution_file:
        for line in solution_file:

            # Split line to retrieve infos: intervention name and decided starting date
            tmp = line.split(' ')
            intervention_name = tmp[0]
            start_time_str = tmp[1].split('\n')[0]

            # Assert intervention exists
            if not intervention_name in interventions:
                instance[FEASIBILITY] = False
                instance[VIOLATIONS] += ['ERROR: Unexpected Intervention ' + intervention_name + ' in solution file ' + solution_filename + '.']
                continue

            # Assert starting date is an integer
            start_time: int
            try:
                start_time = int(start_time_str)
            except ValueError:
                instance[FEASIBILITY] = False
                instance[VIOLATIONS] += ['ERROR: Unexpected starting time ' + start_time_str + ' for Intervention ' + intervention_name + '. Expect integer value.']
                continue

            # Assert no duplicate
            if START_STR in interventions[intervention_name]:
                instance[FEASIBILITY] = False
                instance[VIOLATIONS] += ['ERROR: Duplicate entry for Intervention ' + intervention_name + '. Only first read value is being considered.']
                continue

            # Store starting time
            interventions[intervention_name][START_STR] = start_time


def compute_resources(instance: dict):
    """Compute effective workload (i.e. resources consumption values) for every resource and every time step"""

    # Retrieve useful infos
    interventions = instance[INTERVENTIONS_STR]
    T_max = instance[T_STR]
    resources = instance[RESOURCES_STR]

    # Init resource usage dictionary for each resource and time
    resources_usage = {}
    for resource_name in resources.keys():
        resources_usage[resource_name] = np.zeros(T_max)

    # Compute value for each resource and time step
    for intervention_name, intervention in interventions.items():

        # Start time should be defined (already checked in scheduled constraint checker)
        if not START_STR in intervention:
            continue

        start_time = intervention[START_STR]
        start_time_idx = start_time - 1  # index of list starts at 0
        intervention_workload = intervention[RESOURCE_CHARGE_STR]
        intervention_delta = int(intervention[DELTA_STR][start_time_idx])

        # Compute effective workload
        for resource_name, intervention_resource_worload in intervention_workload.items():
            for time in range(start_time_idx, start_time_idx + intervention_delta):

                # Null values are not available
                if str(time+1) in intervention_resource_worload and str(start_time) in intervention_resource_worload[str(time+1)]:
                    resources_usage[resource_name][time] += intervention_resource_worload[str(time+1)][str(start_time)]

    return resources_usage


def compute_risk_distribution(interventions: dict, T_max: int, scenario_numbers):
    """Compute risk distributions for all time steps, given the interventions starting times"""

    # Init risk table
    risk = [scenario_numbers[t] * [0] for t in range(T_max)]

    # Compute for each intervention independently
    for intervention in interventions.values():

        # Retrieve Intervention's useful infos
        intervention_risk = intervention[RISK_STR]

        # Start time should be defined (already checked in scheduled constraint checker)
        if not START_STR in intervention:
            continue

        start_time = intervention[START_STR]
        start_time_idx = int(start_time) - 1  # index for list getter
        delta = int(intervention[DELTA_STR][start_time_idx])
        for time in range(start_time_idx, start_time_idx + delta):
            for i, additional_risk in enumerate(intervention_risk[str(time + 1)][str(start_time)]):
                risk[time][i] += additional_risk

    return risk


def compute_mean_risk(risk, T_max: int, scenario_numbers):
    """Compute mean risk values over each time period"""

    mean_risk = np.zeros(T_max)
    for t in range(T_max):
        mean_risk[t] = sum(risk[t]) / scenario_numbers[t]

    return mean_risk


def compute_quantile(risk, T_max: int, scenario_numbers, quantile):
    """Compute Quantile values over each time period"""

    q = np.zeros(T_max)
    for t in range(T_max):
        risk[t].sort()
        q[t] = risk[t][int(np.ceil(scenario_numbers[t] * quantile))-1]

    return q


def compute_costs(instance: dict):
    """Compute costs (mean and expected excess)"""

    # Retrieve useful infos
    T_max = instance[T_STR]
    scenario_numbers = instance[SCENARIO_NUMBER]
    interventions = instance[INTERVENTIONS_STR]
    quantile = instance[QUANTILE_STR]

    # Retrieve risk final distribution
    risk = compute_risk_distribution(interventions, T_max, scenario_numbers)

    # Compute mean risk
    mean_risk = compute_mean_risk(risk, T_max, scenario_numbers)

    # Compute quantile
    q = compute_quantile(risk, T_max, scenario_numbers, quantile)

    return mean_risk, q


def compute_objective(instance: dict, mean_risk, quantile):
    """Compute objective function"""

    # Useful infos
    alpha = instance[ALPHA_STR]
    tmp = np.zeros(len(quantile))

    obj_1 = np.mean(mean_risk)
    obj_2 = np.mean(np.max(np.vstack((quantile - mean_risk, tmp)), axis=0))
    obj_total = alpha * obj_1 + (1-alpha)*obj_2

    return obj_total, obj_1, obj_2


def check_all_constraints(instance: dict):
    """Run all constraint checks"""

    check_schedule(instance)
    check_resources(instance)
    check_exclusions(instance)


def check_schedule(instance: dict):
    """Check schedule constraints"""

    # Continuous interventions: §4.1.1
    #   This constraint is implicitly checked by the resource computation:
    #   computation is done under continuity hypothesis, and resource bounds will ensure the feasibility

    # Checks is done on each Intervention
    interventions = instance[INTERVENTIONS_STR]
    for intervention_name, intervention in interventions.items():

        # Interventions are planned once: §4.1.2
        #   assert a starting time has been assigned to the intervention
        if not START_STR in intervention:
            instance[FEASIBILITY] = False
            instance[VIOLATIONS] += ['ERROR: Schedule constraint 4.1.2: Intervention ' + intervention_name + ' has not been scheduled.']
            continue

        # Starting time validity: no explicit constraint
        start_time = intervention[START_STR]
        horizon_end = instance[T_STR]
        if not (1 <= start_time <= horizon_end):
            instance[FEASIBILITY] = False
            instance[VIOLATIONS] += ['ERROR: Schedule constraint 4.1 time validity: Intervention ' + intervention_name + ' starting time ' + str(start_time)
                                     + ' is not a valid starting date. Expected value between 1 and ' + str(horizon_end) + '.']

            # Remove start time to avoid later access errors
            del intervention[START_STR]
            continue

        # No work left: §4.1.3
        #   assert intervention is not ongoing after time limit or end of horizon
        time_limit = int(intervention[TMAX_STR])
        if time_limit < start_time:
            instance[FEASIBILITY] = False
            instance[VIOLATIONS] += ['ERROR: Schedule constraint 4.1.3: Intervention ' + intervention_name + ' realization exceeds time limit.'
                                     + ' It starts at ' + str(start_time) + ' while time limit is ' + str(time_limit) + '.']

            # Remove start time to avoid later access errors
            del intervention[START_STR]
            continue


def check_resources(instance: dict):
    """Check resources constraints"""

    T_max = instance[T_STR]
    Resources = instance[RESOURCES_STR]

    # Bounds are checked with a tolerance value
    tolerance = 1e-5

    # Compute resource usage
    resource_usage = compute_resources(instance)  # dict on resources and time

    # Compare bounds to usage
    for resource_name, resource in Resources.items():
        for time in range(T_max):

            # retrieve bounds values
            upper_bound = resource[MAX_STR][time]
            lower_bound = resource[MIN_STR][time]

            # Consumed value
            worload = resource_usage[resource_name][time]

            # Check max
            if worload > upper_bound + tolerance:
                instance[FEASIBILITY] = False
                instance[VIOLATIONS] += ['ERROR: Resources constraint 4.2 upper bound: Worload on Resource ' + resource_name + ' at time ' + str(time+1) + ' exceeds upper bound.'
                                         + ' Value ' + str(worload) + ' is greater than bound ' + str(upper_bound) + ' plus tolerance ' + str(tolerance) + '.']

            # Check min
            if worload < lower_bound - tolerance:
                instance[FEASIBILITY] = False
                instance[VIOLATIONS] += ['ERROR: Resources constraint 4.2 lower bound: Worload on Resource ' + resource_name + ' at time ' + str(time+1) + ' does not match lower bound.'
                                         + ' Value ' + str(worload) + ' is lower than bound ' + str(lower_bound) + ' minus tolerance ' + str(tolerance) + '.']


def check_exclusions(instance: dict):
    """Check exclusions constraints"""

    # Retrieve Interventions and Exclusions
    interventions = instance[INTERVENTIONS_STR]
    exclusions = instance[EXCLUSIONS_STR]

    # Assert every exclusion holds
    for exclusion in exclusions.values():

        # Retrieve exclusion infos
        [intervention_1_name, intervention_2_name, season] = exclusion

        # Retrieve concerned interventions...
        intervention_1 = interventions[intervention_1_name]
        intervention_2 = interventions[intervention_2_name]

        # Start time should be defined (already checked in scheduled constraint checker)
        if (not START_STR in intervention_1) or (not START_STR in intervention_2):
            continue

        # ... their respective starting times...
        intervention_1_start_time = intervention_1[START_STR]
        intervention_2_start_time = intervention_2[START_STR]

        # ... and their respective deltas (duration)
        intervention_1_delta = int(intervention_1[DELTA_STR][intervention_1_start_time - 1])  # get index in list
        intervention_2_delta = int(intervention_2[DELTA_STR][intervention_2_start_time - 1])  # get index in list

        # Check overlaps for each time step of the season
        for time_str in instance[SEASONS_STR][season]:
            time = int(time_str)
            if (intervention_1_start_time <= time < intervention_1_start_time + intervention_1_delta) and (intervention_2_start_time <= time < intervention_2_start_time + intervention_2_delta):
                instance[FEASIBILITY] = False
                instance[VIOLATIONS] += ['ERROR: Exclusions constraint 4.3: Interventions ' + intervention_1_name + ' and ' + intervention_2_name
                                         + ' are both ongoing at time ' + str(time) + '.']


def check_solution(instance_file, solution_file):
    """Control checker actions"""

    # Read Instance
    instance = read_instance(instance_file)

    # Read Solution
    read_solution(instance, solution_file)

    # Check all constraints
    check_all_constraints(instance)

    # Compute value of objective function
    mean_risk, quantile = compute_costs(instance)
    obj_total, obj_1, obj_2 = compute_objective(instance, mean_risk, quantile)

    return obj_total, obj_1, obj_2, instance[VIOLATIONS]

