#ifndef ROADEF2020_ROADEF2020_H
#define ROADEF2020_ROADEF2020_H

#include <cstddef>
#include <map>
#include <tuple>
#include <vector>
#include <set>
#include <string>


namespace orcs::roadef2020 {

    // =================================================================================================================
    // Enums and constants to handle seasons.
    // =================================================================================================================

    /**
     * Seasons type.
     */
    enum class Season : int {
        FULL = 0,
        WINTER = 1,
        SUMMER = 2,
        INTER_SEASON = 3
    };

    /**
     * Mapping of season indexes to names.
     */
    const std::map<Season, std::string> SEASON_NAMES = {
            {Season::FULL, "full"},
            {Season::WINTER, "winter"},
            {Season::SUMMER, "summer"},
            {Season::INTER_SEASON, "is"}
    };

    /**
     * Mapping of season names to indexes.
     */
    const std::map<std::string, Season> SEASON_INDEXES = {
            {"full", Season::FULL},
            {"winter", Season::WINTER},
            {"summer", Season::SUMMER},
            {"is", Season::INTER_SEASON}
    };


    // =================================================================================================================
    // Structs and classes used to store the problem's data.
    // =================================================================================================================

    struct Resource;
    struct Intervention;
    struct ExclusionComparator;
    class Problem;
    class Schedule;


    /**
     * Resource's data.
     */
    struct Resource {
        std::string name;
        std::vector<double> min;
        std::vector<double> max;
    };

    /**
     * Intervention's data
     */
    struct Intervention {
        std::string name;
        int t_max;
        std::vector<int> delta;
        std::map<int, std::map<int, std::map<int, double> > > workload;
        std::map<int, std::map<int, std::vector<double> > > risk;
    };

    /**
     * Alias to type that defines an exclusion
     */
    using Exclusion = std::tuple<int, int, Season, std::string>;

    /**
     * Comparator used in the maintenance of sets of Exclusions.
     */
    struct ExclusionComparator {
        bool operator() (const Exclusion& lhs, const Exclusion& rhs) const;
    };

    /**
     * The grid operation-based outage maintenance planning problem.
     */
    class Problem {
    public:

        /**
         * Constructor.
         * @param filename Path to the file containing the problem data.
         */
        Problem(const std::string& filename);

        /**
         * Destructor.
         */
        ~Problem();

        /**
         * Return the instance name.
         */
        const std::string& name() const;

        /**
         * The number of time steps in the time horizon.
         * @return The number of time steps in the time horizon.
         */
        int T() const;

        /**
         * The alpha parameter used in the objective function.
         * @return The value of alpha.
         */
        double alpha() const;

        /**
         * The quantile parameter used to define which scenario to use in the computation of excess.
         * @return The quantile value.
         */
        double quantile() const;

        /**
         * Return the index of the scenario used to compute the excess for a given time slot. The index returned
         * considers the scenarios are ordered by risk.
         * @param t A time step.
         * @return Index of the scenario used to compute the excess.
         */
        int quantile_index(int t) const;

        /**
         * Return a set containing the time steps in a given season.
         * @param season A season,
         * @return A set of time steps in the season.
         */
        const std::set<int>& season(Season season) const;

        /**
         * Return the season in which a given time step belongs.
         * @param t A time step
         * @return A season.
          */
        Season which_season(int t) const;

        /**
         * The number of scenarios available for a given time step.
         * @param t A time step.
         * @return The number of scenarios available for the time step.
         */
        int count_scenarios(int t) const;

        /**
         * Return the number of scenarios in the time step that contains more scenarios.
         * @return The number of scenarios in the time step that contains more scenarios.
         */
        int max_scenarios() const;

        /**
         * Return the list of resources.
         * @return The list of resources.
         */
        const std::vector<Resource>& resources() const;

        /**
         * The number of resources.
         * @return The number of resources.
         */
        int count_resources() const;

        /**
         * Return the name of a given resource.
         * @param r A resource.
         * @return The resource's name.
         */
        const std::string& resource_name(int r) const;

        /**
         * The lower bound on consumption of resource r in time step t.
         * @param r A resource.
         * @param t A time slot.
         * @return The lower bound on resource consumption.
         */
        double resource_min(int r, int t) const;

        /**
         * The upper bound on consumption of resource r in time step t.
         * @param r A resource.
         * @param t A time step.
         * @return The upper bound on resource consumption.
         */
        double resource_max(int r, int t) const;

        /**
         * Return the list of interventions.
         * @return The list of interventions.
         */
        const std::vector<Intervention>& interventions() const;

        /**
         * The number of interventions.
         * @return The number of interventions.
         */
        int count_interventions() const;

        /**
         * Return the name of a given intervention.
         * @param i An intervention.
         * @return The intervention's name.
         */
        const std::string& intervention_name(int i) const;

        /**
         * The last time step in in the time horizone in which the intervention can be assigned to start. That is,
         * the intervention cannot be started after the returned value.
         * @param i An intervention.
         * @return The last time step allowed for starting the intervention.
         */
        int t_max(int i) const;

        /**
         * The duration (number of time steps in the time horizon) of an intervention when started at a given time
         * step.
         * @param i An intervention.
         * @param ts A time step.
         * @return The number of time slots required to perform the intervention.
         */
        int delta(int i, int ts) const;

        /**
         * Workload from resource r required by an intervention i in time step t if the intervention is assigned to
         * start at time step ts.
         * @param i Intervention.
         * @param r Resource.
         * @param t Time step.
         * @param ts Starting time step.
         * @return The workload required.
         */
        double workload(int i, int r, int t, int ts) const;

        /**
         * The previously calculated risk on the s-th scenario of performing an intervention i in time step t if the
         * intervention is assigned to start at time step ts.
         * @param i Intervention.
         * @param t Time step.
         * @param ts Starting time step.
         * @param s Scenario.
         * @return The computed risk.
         */
        double risk(int i, int t, int ts, int s) const;

        /**
         * Return the set of exclusions (disjunctive constraints).
         * Each exclusion is defined by a tuple <i1, i2, season, name> in which i1 and i2 are interventions, season is
         * some season and name is the name of the exclusion. An exclusion <i1, i2, season, name> means that i1 and i2
         * cannot share time steps in season.
         * @return The set of exclusions.
         */
        const std::set<Exclusion, ExclusionComparator>& exclusions() const;

        /**
         * Check if two interventions i1 and i2 are allowed to be performed at same time steps in a given season.
         * @param i1 First intervention.
         * @param i2 Second intervention.
         * @param season A season
         * @return True if interventions are allowed to be performed at same time steps of the given season, or false
         * otherwise.
         */
        bool is_simultaneity_allowed(int i1, int i2, Season season) const;


    private:
        std::string name_;
        int T_;
        double alpha_;
        double quantile_;
        int max_scenarios_;
        std::vector<int> scenarios_;
        std::map< Season, std::set<int> > seasons_;
        std::vector<Season> step_season_;
        std::vector<Resource> resources_;
        std::map<std::string, int> resources_indexes_;
        std::vector<Intervention> interventions_;
        std::map<std::string, int> interventions_indexes_;
        std::set< std::tuple<int, int, Season, std::string>, ExclusionComparator> exclusions_;

    };


    // =================================================================================================================
    // Class that represents a solution to the problem.
    // =================================================================================================================

    /**
     * Schedule is the class that represents a solution for the grid opreation-based outage maintenance planning
     * problem.
     */
    class Schedule {

    public:

        /**
         * Constructor.
         * @param problem An instance of the maintenance planning problem to which this schedule refers.
         */
        Schedule(const Problem& problem);

        /**
         * Constructor.
         * @param problem An instance of the maintenance planning problem to which this schedule refers.
         * @param filename Path to solution file.
         */
        Schedule(const Problem& problem, const std::string& filename);

        /**
         * Copy constructor.
         * @param other Schedule to be copied.
         */
        Schedule(const Schedule& other);

        /**
         * Move constructor.
         * @param other Schedule to be moved.
         */
        Schedule(Schedule&& other);

        /**
         * Destructor.
         */
        ~Schedule();

        /**
         * Assignment operator.
         * @param other Schedule to be copied.
         * @return A reference to this.
         */
        Schedule& operator = (const Schedule& other);

        /**
         * Move operator.
         * @param other Schedule to be moved.
         * @return A reference to this.
         */
        Schedule& operator = (Schedule&& other);

        /**
         * Check if an intervention has the starting time set.
         * @param i An intervention.
         * @return True if the intervention has the starting time set, false otherwise.
         */
        bool is_set(int i) const;

        /**
         * Assign a starting time to an intervention.
         * @param i Intervention.
         * @param start Starting time.
         */
        void set(int i, int ts);

        /**
         * Unset the starting time previously assigned (if any) to an intervention i.
         * @param i Intervention.
         */
        void unset(int i);

        /**
         * Return the starting time assigned to intervention i. If the intervention does not have any starting time
         * assigned, -1 is returned.
         * @param i Intervention.
         * @return The starting time assigned to the intervention.
         */
        int get(int i) const;

        /**
         * Return a tuple with the following values: objective, mean risk, and expected excess.
         * @return A tuple with the following values: objective, mean risk, and expected excess.
         */
        std::tuple<double, double, double> evaluation() const;

        /**
         * Check if this schedule is feasible, i.e., it satisfy all constraints imposed by the problem.
         * @return True if this schedule is feasible, false otherwise.
         */
        bool is_feasible() const;

        /**
         * Count the number of violated constraints.
         * @return The number of violated constraints.
         */
        int count_infeasibilities() const;

        /**
         * Return the current risk in time step t and scenario s.
         * @param t Time step.
         * @param s Scenario.
         * @return The current risk.
         */
        double risk(int t, int s) const;

        /**
         * Return the current workload in time step t for a resource r.
         * @param t Time step.
         * @param s Resource
         * @return The current workload.
         */
        double workload(int t, int r) const;

        /**
         * Export this solution to a file with the format specified by the ROADEF 2020 Challenge.
         * @param os output stream.
         */
        void export_ROADEF2020(std::ostream& os) const;

    private:
        const Problem& problem_;
        int* start_;
        double** risk_;
        double** workload_;

    };

}

#endif
