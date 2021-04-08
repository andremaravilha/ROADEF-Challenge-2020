#include "roadef2020.h"

#include <cassert>
#include <cmath>
#include <algorithm>
#include <utility>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <json.hpp>


// =====================================================================================================================
// Implementation of member functions from ExclusionComparator struct.
// =====================================================================================================================

bool orcs::roadef2020::ExclusionComparator::operator() (const Exclusion& lhs, const Exclusion& rhs) const {
    const auto& [lhs_i1, lhs_i2, lhs_season, lhs_name] = lhs;
    const auto& [rhs_i1, rhs_i2, rhs_season, rhs_name] = rhs;
    if (lhs_i1 < rhs_i1) return true;
    if (lhs_i1 == rhs_i1 && lhs_i2 < rhs_i2) return true;
    if (lhs_i1 == rhs_i1 && lhs_i2 == rhs_i2 && lhs_season < rhs_season) return true;
    return false;
}


// =====================================================================================================================
// Implementation of member functions from Problem class.
// =====================================================================================================================

orcs::roadef2020::Problem::Problem(const std::string& filename) {

    // Parse JSON data
    std::ifstream ifs(filename, std::ifstream::in);
    nlohmann::json json = nlohmann::json::parse(ifs);
    ifs.close();

    // Get instance name
    name_ = std::filesystem::path(filename).replace_extension("").filename().string();

    // Read some data
    json.at("T").get_to(T_);
    json.at("Alpha").get_to(alpha_);
    json.at("Quantile").get_to(quantile_);
    json.at("Scenarios_number").get_to(scenarios_);
    max_scenarios_ = *(std::max_element(scenarios_.begin(), scenarios_.end()));

    // Read seasons' data
    const auto& seasons = json.at("Seasons");
    step_season_.insert(step_season_.begin(), T_, Season::FULL);
    for (const auto& [name, index] : SEASON_INDEXES) {
        seasons_.insert(std::pair(index, std::set<int>()));
        const auto steps = seasons.find(name);
        if (steps != seasons.end()) {
            for (const auto& step : *steps) {
                int t = (step.is_number() ? step.get<int>() : std::stoi(step.get<std::string>())) - 1;
                seasons_.at(index).insert(t);
                step_season_[t] = index;
            }
        }
    }

    // Read resources' data
    const auto& resources = json.at("Resources");
    resources_.reserve(resources.size());
    int resource_idx = 0;
    for (const auto& [name, data] : resources.items()) {
        Resource resource {name, {}, {}};
        data.at("min").get_to(resource.min);
        data.at("max").get_to(resource.max);
        resources_.push_back(std::move(resource));
        resources_indexes_.insert(std::pair(name, resource_idx));
        ++resource_idx;
    }

    // Read interventions' data
    const auto& interventions = json.at("Interventions");
    interventions_.reserve(interventions.size());
    int intervention_idx = 0;
    for (const auto& [name, data] : interventions.items()) {
        Intervention intervention {name, 0, {}, {}, {}};
        intervention.t_max = (data.at("tmax").is_number() ? data.at("tmax").get<int>() : std::stoi(data.at("tmax").get<std::string>())) - 1;
        data.at("Delta").get_to(intervention.delta);

        // Workload
        const auto& workload = data.at("workload");
        for (const auto& [resource_name, workload_data] : workload.items()) {
            int r = resources_indexes_.at(resource_name);
            intervention.workload.insert(std::pair(r, std::map<int, std::map<int, double > >()));
            for (const auto& [t_str, t_data] : workload_data.items()) {
                int t = std::stoi(t_str) - 1;
                intervention.workload.at(r).insert(std::pair(t, std::map<int, double >()));
                for (const auto& [ts_str, ts_data] : t_data.items()) {
                    int ts = std::stoi(ts_str) - 1;
                    intervention.workload.at(r).at(t).insert(std::pair(ts, ts_data.get<double>()));
                }
            }
        }

        // Risk
        const auto& risk = data.at("risk");
        for (const auto& [t_str, t_data] : risk.items()) {
            int t = std::stoi(t_str) - 1;
            intervention.risk.insert(std::pair(t, std::map<int, std::vector<double> >()));
            for (const auto& [ts_str, ts_data] : t_data.items()) {
                int ts = std::stoi(ts_str) - 1;
                intervention.risk.at(t).insert(std::pair(ts, std::vector<double>()));
                intervention.risk.at(t).at(ts).reserve(ts_data.size());
                ts_data.get_to(intervention.risk.at(t).at(ts));
            }
        }

        interventions_.push_back(std::move(intervention));
        interventions_indexes_.insert(std::pair(name, intervention_idx));
        ++intervention_idx;
    }

    // Read exclusions' data
    const auto& exclusions = json.at("Exclusions");
    for (const auto& [name, data] : exclusions.items()) {
        int i1 = interventions_indexes_.at(data[0].get<std::string>());
        int i2 = interventions_indexes_.at(data[1].get<std::string>());
        Season season = SEASON_INDEXES.at(data[2].get<std::string>());

        if (i1 > i2) std::swap(i1, i2);

        exclusions_.insert(std::make_tuple(i1, i2, season, name));
    }

}


orcs::roadef2020::Problem::~Problem() {
    // It does nothing here.
}


const std::string& orcs::roadef2020::Problem::name() const {
    return name_;
}


int orcs::roadef2020::Problem::T() const {
    return T_;
}


double orcs::roadef2020::Problem::alpha() const {
    return alpha_;
}


double orcs::roadef2020::Problem::quantile() const {
    return quantile_;
}


int orcs::roadef2020::Problem::quantile_index(int t) const {
    assert(("Invalid time step", t >= 0 && t < T()));
    return ((int) (std::ceil(quantile_ * scenarios_[t]) + 0.5)) - 1;
}


const std::set<int>& orcs::roadef2020::Problem::season(Season season) const {
    return seasons_.at(season);
}


orcs::roadef2020::Season orcs::roadef2020::Problem::which_season(int t) const {
    assert(("Invalid time step", t >= 0 && t < T()));
    return step_season_[t];
}


int orcs::roadef2020::Problem::count_scenarios(int t) const {
    assert(("Invalid time step", t >= 0 && t < T()));
    return scenarios_[t];
}


int orcs::roadef2020::Problem::max_scenarios() const {
    return max_scenarios_;
}


const std::vector<orcs::roadef2020::Resource>& orcs::roadef2020::Problem::resources() const {
    return resources_;
}


int orcs::roadef2020::Problem::count_resources() const {
    return resources_.size();
}


const std::string& orcs::roadef2020::Problem::resource_name(int r) const {
    assert(("Invalid resource", r >= 0 && r < resources_.size()));
    return resources_[r].name;
}


double orcs::roadef2020::Problem::resource_min(int r, int t) const {
    assert(("Invalid resource", r >= 0 && r < resources_.size()));
    assert(("Invalid time step", t >= 0 && t < T()));
    return resources_[r].min[t];
}


double orcs::roadef2020::Problem::resource_max(int r, int t) const {
    assert(("Invalid resource", r >= 0 && r < resources_.size()));
    assert(("Invalid time step", t >= 0 && t < T()));
    return resources_[r].max[t];
}


const std::vector<orcs::roadef2020::Intervention>& orcs::roadef2020::Problem::interventions() const {
    return interventions_;
}


int orcs::roadef2020::Problem::count_interventions() const {
    return interventions_.size();
}


const std::string& orcs::roadef2020::Problem::intervention_name(int i) const {
    assert(("Invalid intervention", i >= 0 && i < interventions_.size()));
    return interventions_[i].name;
}

int orcs::roadef2020::Problem::t_max(int i) const {
    assert(("Invalid intervention", i >= 0 && i < interventions_.size()));
    return interventions_[i].t_max;
}


int orcs::roadef2020::Problem::delta(int i, int ts) const {
    assert(("Invalid intervention", i >= 0 && i < interventions_.size()));
    assert(("Invalid starting step time", ts >= 0 && ts <= interventions_[i].t_max));
    return interventions_[i].delta[ts];
}


double orcs::roadef2020::Problem::workload(int i, int r, int t, int ts) const {
    assert(("Invalid intervention", i >= 0 && i < interventions_.size()));
    assert(("Invalid resource", r >= 0 && r < resources_.size()));
    assert(("Invalid time step", t >= 0 && t < T()));
    assert(("Invalid starting step time", ts >= 0 && ts <= interventions_[i].t_max));

    auto iter_r = interventions_[i].workload.find(r);
    if (iter_r != interventions_[i].workload.end()) {
        auto iter_r_t = iter_r->second.find(t);
        if (iter_r_t != iter_r->second.end()) {
            auto iter_r_t_ts = iter_r_t->second.find(ts);
            if (iter_r_t_ts != iter_r_t->second.end()) {
                return iter_r_t_ts->second;
            }
        }
    }

    return 0.0;
}


double orcs::roadef2020::Problem::risk(int i, int t, int ts, int s) const {
    assert(("Invalid intervention", i >= 0 && i < interventions_.size()));
    assert(("Invalid time step", t >= 0 && t < T()));
    assert(("Invalid starting step time", ts >= 0 && ts <= interventions_[i].t_max));
    assert(("Invalid scenario", s >= 0 && s < scenarios_[t]));

    auto iter_t = interventions_[i].risk.find(t);
    if (iter_t != interventions_[i].risk.end()) {
        auto iter_t_ts = iter_t->second.find(ts);
        if (iter_t_ts != iter_t->second.end()) {
            return iter_t_ts->second[s];
        }
    }

    return 0.0;
}


const std::set<orcs::roadef2020::Exclusion,  orcs::roadef2020::ExclusionComparator>& orcs::roadef2020::Problem::exclusions() const {
    return exclusions_;
}


bool orcs::roadef2020::Problem::is_simultaneity_allowed(int i1, int i2, Season season) const {
    assert(("Invalid intervention", i1 >= 0 && i1 < interventions_.size()));
    assert(("Invalid intervention", i2 >= 0 && i2 < interventions_.size()));
    if (i1 > i2) std::swap(i1, i2);
    return (exclusions_.find(std::make_tuple(i1, i2, season, "")) != exclusions_.end());
}


// =====================================================================================================================
// Implementation of member functions from Schedule class.
// =====================================================================================================================

orcs::roadef2020::Schedule::Schedule(const Problem& problem) : problem_(problem) {
    start_ = new int[problem_.count_interventions()];
    for (int i = 0; i < problem_.count_interventions(); ++i) {
        start_[i] = -1;
    }

    risk_ = new double*[problem_.T()];
    workload_ = new double*[problem_.T()];
    for (int t = 0; t < problem_.T(); ++t) {

        risk_[t] = new double[problem_.count_scenarios(t)];
        for (int s = 0; s < problem_.count_scenarios(t); ++s) {
            risk_[t][s] = 0.0;
        }

        workload_[t] = new double[problem_.count_resources()];
        for (int r = 0; r < problem_.count_resources(); ++r) {
            workload_[t][r] = 0.0;
        }
    }
}


orcs::roadef2020::Schedule::Schedule(const Problem& problem, const std::string& filename) : problem_(problem) {
    // TODO
}


orcs::roadef2020::Schedule::Schedule(const Schedule& other): problem_(other.problem_) {
    start_ = new int[problem_.count_interventions()];
    for (int i = 0; i < problem_.count_interventions(); ++i) {
        start_[i] = other.start_[i];
    }

    risk_ = new double*[problem_.T()];
    workload_ = new double*[problem_.T()];
    for (int t = 0; t < problem_.T(); ++t) {

        risk_[t] = new double[problem_.count_scenarios(t)];
        for (int s = 0; s < problem_.count_scenarios(t); ++s) {
            risk_[t][s] = other.risk_[t][s];
        }

        workload_[t] = new double[problem_.count_resources()];
        for (int r = 0; r < problem_.count_resources(); ++r) {
            workload_[t][r] = other.workload_[t][r];
        }
    }
}


orcs::roadef2020::Schedule::Schedule(Schedule&& other): problem_(other.problem_), start_(nullptr),
risk_(nullptr), workload_(nullptr) {
    std::swap(start_, other.start_);
    std::swap(risk_, other.risk_);
    std::swap(workload_, other.workload_);
}


orcs::roadef2020::Schedule::~Schedule() {
    if (start_ != nullptr) {
        delete[] start_;
        start_ = nullptr;
    }

    if (workload_ != nullptr) {
        delete[] workload_;
        workload_ = nullptr;
    }

    if (risk_ != nullptr) {
        for (int t = 0; t < problem_.T(); ++t) {
            if (risk_[t] != nullptr) {
                delete[] risk_[t];
                risk_[t] = nullptr;
            }
        }
        delete[] risk_;
        risk_ = nullptr;
    }

    if (workload_ != nullptr) {
        for (int t = 0; t < problem_.T(); ++t) {
            if (workload_[t] != nullptr) {
                delete[] workload_[t];
                workload_[t] = nullptr;
            }
        }
        delete[] workload_;
        workload_ = nullptr;
    }
}


orcs::roadef2020::Schedule& orcs::roadef2020::Schedule::operator = (const Schedule& other) {
    if (this != &other) {
        assert(("Different problem instances", &problem_ == &(other.problem_)));

        for (int i = 0; i < problem_.count_interventions(); ++i) {
            start_[i] = other.start_[i];
        }

        for (int t = 0; t < problem_.T(); ++t) {

            for (int s = 0; s < problem_.count_scenarios(t); ++s) {
                risk_[t][s] = other.risk_[t][s];
            }

            for (int r = 0; r < problem_.count_resources(); ++r) {
                workload_[t][r] = other.workload_[t][r];
            }

        }
    }

    return *this;
}


orcs::roadef2020::Schedule& orcs::roadef2020::Schedule::operator = (Schedule&& other) {
    if (this != &other) {
        assert(("Different problem instances", &problem_ == &(other.problem_)));

        if (start_ != nullptr) {
            delete[] start_;
            start_ = nullptr;
        }

        if (risk_ != nullptr) {
            for (int t = 0; t < problem_.T(); ++t) {
                if (risk_[t] != nullptr) {
                    delete[] risk_[t];
                    risk_[t] = nullptr;
                }
            }
            delete[] risk_;
            risk_ = nullptr;
        }

        if (workload_ != nullptr) {
            for (int t = 0; t < problem_.T(); ++t) {
                if (workload_[t] != nullptr) {
                    delete[] workload_[t];
                    workload_[t] = nullptr;
                }
            }
            delete[] workload_;
            workload_ = nullptr;
        }

        std::swap(start_, other.start_);
        std::swap(risk_, other.risk_);
        std::swap(workload_, other.workload_);
    }

    return *this;
}


bool orcs::roadef2020::Schedule::is_set(int i) const {
    assert(("Invalid intervention", i >= 0 && i < problem_.count_interventions()));
    return (start_[i] != -1);
}


void orcs::roadef2020::Schedule::set(int i, int ts) {
    assert(("Invalid intervention", i >= 0 && i < problem_.count_interventions()));
    assert(("Invalid time step", ts >= 0 && ts <= problem_.t_max(i)));

    for (int t = 0; t < problem_.T(); ++t) {

        for (int s = 0; s < problem_.count_scenarios(t); ++s) {
            risk_[t][s] += problem_.risk(i, t, ts, s) - (is_set(i) ? problem_.risk(i, t, start_[i], s) : 0.0);
        }

        for (int r = 0; r < problem_.count_resources(); ++r) {
            workload_[t][r] += problem_.workload(i, r, t, ts) - (is_set(i) ? problem_.workload(i, r, t, start_[i]) : 0.0);
        }

    }

    start_[i] = ts;
}


void orcs::roadef2020::Schedule::unset(int i) {
    assert(("Invalid intervention", i >= 0 && i < problem_.count_interventions()));

    if (is_set(i)) {
        for (int t = 0; t < problem_.T(); ++t) {

            for (int s = 0; s < problem_.count_scenarios(t); ++s) {
                risk_[t][s] -= problem_.risk(i, t, start_[i], s);
            }

            for (int r = 0; r < problem_.count_resources(); ++r) {
                workload_[t][r] -= problem_.workload(i, r, t, start_[i]);
            }
        }

        start_[i] = -1;
    }

}


int orcs::roadef2020::Schedule::get(int i) const {
    assert(("Invalid intervention", i >= 0 && i < problem_.count_interventions()));
    return start_[i];
}


std::tuple<double, double, double> orcs::roadef2020::Schedule::evaluation() const {
    double mean_risk = 0.0;
    double expected_excess = 0.0;
    double** risk_by_scenario = new double*[problem_.max_scenarios()];
    const auto sort_risks = [](double* first, double* second) -> bool { return *first < *second; };

    for (int t = 0; t < problem_.T(); ++t) {
        double risk_t = 0.0;
        for (int s = 0; s < problem_.count_scenarios(t); ++s) {
            risk_t += risk_[t][s];
            risk_by_scenario[s] = &(risk_[t][s]);
        }

        std::nth_element(risk_by_scenario, risk_by_scenario + problem_.quantile_index(t),
                         risk_by_scenario + problem_.count_scenarios(t), sort_risks);

        risk_t = (risk_t / problem_.count_scenarios(t));
        double excess_t = std::max(0.0, *(risk_by_scenario[problem_.quantile_index(t)]) - risk_t);

        mean_risk += risk_t;
        expected_excess += excess_t;
    }

    delete[] risk_by_scenario;
    mean_risk /= problem_.T();
    expected_excess /= problem_.T();
    double objective = ((problem_.alpha() * mean_risk) + ((1.0 - problem_.alpha()) * expected_excess));

    return {objective, mean_risk, expected_excess};
}


bool orcs::roadef2020::Schedule::is_feasible() const {

    constexpr double eps = 1e-6;

    // Check if all interventions has a valid stating time assigned
    for (int i = 0; i < problem_.count_interventions(); ++i) {
        if (start_[i] < 0 || start_[i] > problem_.t_max(i)) {
            return false;
        }
    }

    // Check workload
    for (int t = 0; t < problem_.T(); ++t) {
        for (int r = 0; r < problem_.count_resources(); ++r) {
            if (workload_[t][r] < problem_.resource_min(r, t) - eps || workload_[t][r] > problem_.resource_max(r, t) + eps) {
                return false;
            }
        }
    }

    // Check exclusions
    for (const auto& [i1, i2, season, name] : problem_.exclusions()) {

        int i1_start = start_[i1];
        int i1_end = i1_start + problem_.delta(i1, i1_start) - 1;

        int i2_start = start_[i2];
        int i2_end = i2_start + problem_.delta(i2, i2_start) - 1;

        // If intervals [i1_start, i1_end] and [i2_start, i2_end] overlaps
        if (i1_start <= i2_end && i2_start <= i1_end) {
            int t_start = std::max(i1_start, i2_start);
            int t_end = std::min(i1_end, i2_end);
            for (int t = t_start; t <= t_end; ++t) {
                if (problem_.which_season(t) == season) {
                    return false;
                }
            }
        }
    }

    // All constraints are satisfied.
    return true;
}


int orcs::roadef2020::Schedule::count_infeasibilities() const {

    int infeasibilities = 0;
    constexpr double eps = 1e-6;

    // Check if all interventions has a valid stating time assigned
    for (int i = 0; i < problem_.count_interventions(); ++i) {
        if (start_[i] < 0 || start_[i] > problem_.t_max(i)) {
            ++infeasibilities;
            //std::cout << "Intervention " << i + 1 << " has invalid starting time (" << start_[i] << ")." << std::endl;
        }
    }

    // Check workload
    for (int t = 0; t < problem_.T(); ++t) {
        for (int r = 0; r < problem_.count_resources(); ++r) {
            if (workload_[t][r] < problem_.resource_min(r, t) - eps || workload_[t][r] > problem_.resource_max(r, t) + eps) {
                ++infeasibilities;
                //std::cout << "Resource " << r + 1 << " out of bounds (" << workload_[t][r] <<
                //", [" << problem_.resource_min(r, t) << ", " << problem_.resource_max(r, t) << "])." << std::endl;
            }
        }
    }

    // Check exclusions
    for (const auto& [i1, i2, season, name] : problem_.exclusions()) {
        if (is_set(i1) && is_set(i2)) {

            int i1_start = start_[i1];
            int i1_end = i1_start + problem_.delta(i1, i1_start) - 1;

            int i2_start = start_[i2];
            int i2_end = i2_start + problem_.delta(i2, i2_start) - 1;

            // If intervals [i1_start, i1_end] and [i2_start, i2_end] overlaps
            if (i1_start <= i2_end && i2_start <= i1_end) {
                int t_start = std::max(i1_start, i2_start);
                int t_end = std::min(i1_end, i2_end);
                for (int t = t_start; t <= t_end; ++t) {
                    if (problem_.which_season(t) == season) {
                        //std::cout << "Disjunction (" << i1 + 1 << ", " << i2 + 1 << ") in " << SEASON_NAMES.at(season) << " is violated." << std::endl;
                        ++infeasibilities;
                    }
                }
            }
        }
    }

    return infeasibilities;
}


double orcs::roadef2020::Schedule::risk(int t, int s) const {
    assert(("Invalid time step", t >= 0 && t < problem_.T()));
    assert(("Invalid scenario", s >= 0 && s < problem_.count_scenarios(t)));
    return risk_[t][s];
}


double orcs::roadef2020::Schedule::workload(int t, int r) const {
    assert(("Invalid time step", t >= 0 && t < problem_.T()));
    assert(("Invalid resource", r >= 0 && r < problem_.count_resources()));
    return workload_[t][r];
}


void orcs::roadef2020::Schedule::export_ROADEF2020(std::ostream& os) const {
    for (int i = 0; i < problem_.count_interventions(); ++i) {
        os << problem_.intervention_name(i) << " " << start_[i] + 1 << std::endl;
    }
}
