#include <iostream>
#include <unistd.h>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <random>
#include <sstream>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <optional>
#include <map>

#include <chrono>

#include <stdio.h>
#include <omp.h>

#include <array>
#include <cstdint>
#include <zlib.h>

#include "model_setup.cpp"

int in_which_state(int state_nr_start, int state_nr_end, int t_start, int t_end, int t) {
    int state_nr;
    bool is_non_abm_start = (model_type_of_domain[state_nr_start] != abm_idx); 
    bool is_non_abm_end = (model_type_of_domain[state_nr_end] != abm_idx);
    if (is_non_abm_start) { // state_nr_start in stateIndices_ODE or stateIndices_PDE -> started commuting in PDE or ODE domain -> add person to ABM state as soon as possible if agent is commuting to ABM domain
        if (is_non_abm_end) { // state_nr_end in stateIndices_ODE or stateIndices_PDE -> ended commuting in PDE or ODE domain again -> add person to state it is closest to 
            if (std::abs(t_end-t) < std::abs(t_start-t)) {
                state_nr = state_nr_end;
            }
            else {
                state_nr = state_nr_start;
            }
        }
        else {
            state_nr = state_nr_end;
        }
    }
    else { // started in ABM domain
        // ended in ??
        if (is_non_abm_end) { // state_nr in stateIndices_ODE or stateIndices_PDE -> ended commuting in PDE or ODE domain -> add person to ABM state
            state_nr = state_nr_start;
        }
        else { // figure out where it is closest to
            if (std::abs(t_end-t) < std::abs(t_start-t)) {
                state_nr = state_nr_end;
            }
            else {
                state_nr = state_nr_start;
            }
        }
    }
    return state_nr;
}

class ABModel {
private:
    double population_scale;
    bool with_restriction;
    int dt_inv;
    
public:
    ABModel(double population_scale_, bool with_restriction_, int dt_inv_): population_scale(population_scale_), with_restriction(with_restriction_), dt_inv(dt_inv_) {
        day_shift_mask_wearing_all.resize(nbr_federal_states);
        if (with_restriction) {
            for (int state=0; state<nbr_federal_states; state++) {
                day_shift_mask_wearing_all[state] = 56;
            }
        }
        else {
            for (int state=0; state<nbr_federal_states; state++) {
                day_shift_mask_wearing_all[state] = 1e+6; // no restriction 
            }
        }
    }
    std::vector<int> day_shift_mask_wearing_all;

    std::vector<std::vector<double>> activityChangePercentage;
    std::vector<std::vector<std::vector<std::vector<double>>>> facilityCoordinates_allCategories;
    std::vector<std::vector<std::vector<int>>> facilityLabels_allCategories; 
    std::vector<std::unordered_map<int, int>> facilityCoordinates_home;

    std::vector<std::vector<std::vector<int>>> eventData; 
    std::vector<std::vector<std::vector<std::vector<int>>>> jump_matrix;
    std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>> germany_jump_matrix_agent_IDs;

    std::vector<std::vector<int>> nbr_agents;

    int nbr_trajectories = 0; // total number of given trajectories
    int nbr_categories;
    std::vector<std::vector<std::vector<int>>> all_initial_agent_ids;
    std::vector<std::map<std::pair<int, int>, int>> max_agents_per_facility_category;

    std::vector<std::string> facility_categories = {
        "home", "leisure", "shop_daily", "errands", "business", "shop_other", "work", "educ_higher", 
        "educ_other", "educ_kiga", "visit", "educ_secondary", "educ_primary", "educ_tertiary"
    };

    std::vector<std::string> days_str = {"wt", "sa", "so"};

    void setLeavingNumberOfPeople() {
        std::string population_scale_str;
        if (int(population_scale+0.1) == 1) {
            population_scale_str = "1";
        }
        else {
            population_scale_str = "4";
        }
        
        jump_matrix.resize(3);
        germany_jump_matrix_agent_IDs.resize(3);
        for (int nr_day=0; nr_day<3; nr_day++) {
            std::string filename = "../../work/input_data/global/germany_jump_matrix_agent_IDs_" + days_str[nr_day] + addOn + "_dt_inv_" + std::to_string(dt_inv) + ".bin";
            std::cout << filename << std::endl;
            std::ifstream file(filename, std::ios::binary);
            if (!file) {
                std::cerr << "Fehler beim Öffnen der Datei " << filename << std::endl;
                return; 
            }
            
            jump_matrix[nr_day].resize(dt_inv, std::vector<std::vector<int>>(nbr_federal_states, std::vector<int>(nbr_federal_states))); 
            germany_jump_matrix_agent_IDs[nr_day].resize(dt_inv, std::vector<std::vector<std::vector<int>>>(nbr_federal_states, std::vector<std::vector<int>>(nbr_federal_states))); 
            int nbr_agent_IDs, agent_ID;
            for (std::size_t t = 0; t < dt_inv; t++) {
                if (t<=1) {
                    std::cout << "time step " << t << std::endl;
                }
                for (int state_1 = 0; state_1 < nbr_federal_states; state_1++) {
                    for (int state_2 = 0; state_2 < nbr_federal_states; state_2++) {
                        file.read(reinterpret_cast<char*>(&nbr_agent_IDs), sizeof(int));
                        jump_matrix[nr_day][t][state_1][state_2] = nbr_agent_IDs;
                        germany_jump_matrix_agent_IDs[nr_day][t][state_1][state_2].resize(nbr_agent_IDs);

                        for (int nbr_agent_ID=0; nbr_agent_ID<nbr_agent_IDs; nbr_agent_ID++) {
                            file.read(reinterpret_cast<char*>(&agent_ID), sizeof(int)); 
                            germany_jump_matrix_agent_IDs[nr_day][t][state_1][state_2][nbr_agent_ID] = agent_ID;
                        }
                        if (t<=1) {
                            std::cout << nbr_agent_IDs << " ";
                        }
                    }
                    if (t<=1) {
                        std::cout << std::endl;
                    }
                }
            }
            file.close();
        }
    }

    void readEventData(int nr_day_specific=-1) {
        std::string population_scale_str;
        if (int(population_scale+0.1) == 1) {
            population_scale_str = "1";
        }
        else {
            population_scale_str = "4";
        }

        int total_nr_days;
        int init_nr_day;
        if (nr_day_specific == -1) {
            total_nr_days = 3;
            init_nr_day = 0;
        }
        else {
            total_nr_days = 1;
            init_nr_day = nr_day_specific;
        }
        eventData.resize(3);
        const int cols = 5; // time in seconds, activity type, person id, facility index (corresponding to facility coordinates of that category), category index (corresponding to facility_categories)
        for (int nr_day=init_nr_day; nr_day<init_nr_day+total_nr_days; nr_day++) {
            std::cout << "trying to read file with nr_day " << nr_day << std::endl;
            std::string filename = "../../work/input_data/global/germany_event_data_" + std::to_string(nr_day) + "_" + population_scale_str + ".bin.gz";
            gzFile file = gzopen(filename.c_str(), "rb");
            if (!file) {
                std::cerr << "Fehler beim Öffnen der Datei " << filename << std::endl;
                return; 
            }

            int cnt_trajectories = 0;
            int previous_agent_id = -1;
            int64_t row_idx = 0;

            bool reached_eof = false;

            while (!reached_eof) {
                std::vector<int> row(cols);

                for (int j = 0; j < cols; j++) {
                    int64_t value = 0;
                    int bytes_read = gzread(file, &value, sizeof(int64_t));

                    if (bytes_read == 0) {
                        if (j != 0) {
                            std::cerr << "Fehler: unvollständiger Datensatz in Datei "
                                    << filename << " bei Zeile " << row_idx << std::endl;
                            gzclose(file);
                            return;
                        }

                        reached_eof = true;
                        break;
                    }

                    if (bytes_read != static_cast<int>(sizeof(int64_t))) {
                        std::cerr << "Fehler beim Lesen der Datei " << filename
                                << " bei Zeile " << row_idx
                                << ", Spalte " << j << std::endl;
                        gzclose(file);
                        return;
                    }
                    row[j] = static_cast<int>(value);
                }

                if (reached_eof) {
                    break;
                }

                eventData[nr_day].push_back(row);

                if (row[2] != previous_agent_id) {
                    cnt_trajectories++;
                    previous_agent_id = row[2];
                }

                if (row_idx < 100) {
                    for (int j = 0; j < cols; j++) {
                        std::cout << row[j] << " ";
                    }
                    std::cout << std::endl;
                }

                row_idx++;
            } 
            if (nr_day == initial_nr_day) { 
                nbr_trajectories = cnt_trajectories;
                std::cout << "nr_day " << nr_day << ", cnt_trajectories " << cnt_trajectories << std::endl;
            }
            gzclose(file);
        }
    }
    
    void readFacilityCoordinates_home() {
        std::string population_scale_str;
        if (int(population_scale+0.1) == 1) {
            population_scale_str = "1";
        }
        else {
            population_scale_str = "4";
        }

        facilityCoordinates_home.resize(3);
        for (int nr_day=0; nr_day<3; nr_day++) {
            std::string filename = "../../work/input_data/global/germany_agent_ids_with_home_facility_coordinates_index_" + std::to_string(nr_day) + "_" + population_scale_str + ".bin"; 
            std::ifstream infile(filename, std::ios::binary);

            if (!infile) {
                std::cerr << "Datei " << filename << " konnte nicht geöffnet werden" << std::endl;
                return; 
            }
            
            // get filesize
            infile.seekg(0, std::ios::end);
            int64_t fileSize = infile.tellg();
            infile.seekg(0, std::ios::beg);
            int cols = 2; // person id, facility index (corresponding to facility coordinates of that category home)
            int64_t num_pairs = fileSize / (cols * sizeof(double));
            std::cout << "category home, num_pairs " << num_pairs << std::endl;
            
            for (unsigned int i = 0; i < num_pairs; i++) {
                int64_t agent_id, facility_id;
                infile.read(reinterpret_cast<char*>(&agent_id), sizeof(agent_id));
                infile.read(reinterpret_cast<char*>(&facility_id), sizeof(facility_id));
                if (!infile) {
                    std::cerr << "Fehler beim Lesen von agent_id,facility_id" << std::endl;
                    break;
                }
                facilityCoordinates_home[nr_day][agent_id] = facility_id;
            }
            infile.close();
        }
    }

    void readFacilityCoordinates_allCategories() {
        std::string population_scale_str;
        if (int(population_scale+0.1) == 1) {
            population_scale_str = "1";
        }
        else {
            population_scale_str = "4";
        }

        facilityCoordinates_allCategories.resize(3);
        facilityLabels_allCategories.resize(3);
        for (int nr_day=0; nr_day<3; nr_day++) {
            std::string filename = "../../work/input_data/global/germany_facility_coordinates_" + std::to_string(nr_day) + "_" + population_scale_str + "_labeled.bin";
            std::ifstream infile(filename, std::ios::binary);
            if (!infile) {
                std::cerr << "Datei " << filename << " konnte nicht geöffnet werden" << std::endl;
                return; 
            }
            
            std::cout << "facility_categories.size() " << facility_categories.size() << std::endl;
            facilityCoordinates_allCategories[nr_day].resize(facility_categories.size());
            facilityLabels_allCategories[nr_day].resize(facility_categories.size());

            for (int i = 0; i < facility_categories.size(); i++) {
                unsigned int num_triples;
                infile.read(reinterpret_cast<char*>(&num_triples), sizeof(num_triples));
                std::cout << "category " << i << ", num_triples " << num_triples << std::endl;

                if (!infile) {
                    std::cerr << "Fehler beim Lesen der Anzahl der Triples" << std::endl;
                    break;
                }
                
                facilityCoordinates_allCategories[nr_day][i].resize(num_triples, std::vector<double>(2));
                facilityLabels_allCategories[nr_day][i].resize(num_triples); 

                for (unsigned int j = 0; j < num_triples; j++) {
                    double x, y;
                    infile.read(reinterpret_cast<char*>(&x), sizeof(x));
                    infile.read(reinterpret_cast<char*>(&y), sizeof(y));
                    if (!infile) {
                        std::cerr << "Fehler beim Lesen der Koordinaten" << std::endl;
                        break;
                    }

                    uint32_t label_idx;
                    infile.read(reinterpret_cast<char*>(&label_idx), sizeof(label_idx));
                    if (!infile) {
                        std::cerr << "Fehler beim Lesen des Labels!" << std::endl;
                        break;
                    }

                    facilityCoordinates_allCategories[nr_day][i][j][0] = x;
                    facilityCoordinates_allCategories[nr_day][i][j][1] = y;
                    facilityLabels_allCategories[nr_day][i][j] = label_idx; 
                }
            }
            infile.close();
        }
        nbr_categories = facilityCoordinates_allCategories[0].size();
    }

    void set_max_number_of_agents_per_facility() {
        std::string population_scale_str;
        if (int(population_scale+0.1) == 1) {
            population_scale_str = "1";
        }
        else {
            population_scale_str = "4";
        }

        max_agents_per_facility_category.resize(3);
        for (int nr_day=0; nr_day<3; nr_day++) {
            std::string filename = "../../work/input_data/global/germany_max_agents_per_facility_" + std::to_string(nr_day) + "_" + population_scale_str + ".bin";
            std::ifstream file(filename, std::ios::binary);
            if (!file) {
                std::cerr << "Fehler beim Öffnen der Datei " << filename << std::endl;
                return; 
            }

            std::vector<int64_t> data; 
            
            file.seekg(0, std::ios::end);
            size_t filesize = file.tellg();
            file.seekg(0, std::ios::beg);

            data.resize(filesize / sizeof(double));
            file.read(reinterpret_cast<char*>(data.data()), filesize);

            file.close();

            for (size_t i = 0; i < data.size(); i += 3) {
                int64_t facility_id = data[i]; 
                int64_t category = data[i + 1];
                int64_t max_agents = data[i + 2];
                max_agents_per_facility_category[nr_day][{facility_id, category}] = max_agents;
            }
            std::cout << "data.size()/3.0 " << (int(data.size()/3.0)) << std::endl;
        }
    }

    void set_initial_agent_ids() {
        std::string population_scale_str;
        if (int(population_scale+0.1) == 1) {
            population_scale_str = "1";
        }
        else {
            population_scale_str = "4";
        }
        
        nbr_agents.resize(3, std::vector<int>(nr_ABM_states));
        all_initial_agent_ids.resize(3, std::vector<std::vector<int>>(nr_ABM_states));
        for (int nr_day = 0; nr_day < 3; nr_day++) {
            std::cout << "\n nr_day " << nr_day << std::endl;
            std::string filenameTotalNumbers = "../../work/input_data/global/germany_nbr_individuals_id_t_0_in_states_" + std::to_string(nr_day) + "_" + population_scale_str + ".bin"; 
            std::ifstream fileTotalNumbers(filenameTotalNumbers, std::ios::binary); 
            if ((fileTotalNumbers.is_open())) { 
                int64_t total_nbr_individuals;
                for (int state_idx=0; state_idx<nbr_federal_states; state_idx++) { 
                    fileTotalNumbers.read(reinterpret_cast<char*>(&total_nbr_individuals), sizeof(int64_t));    
                    
                    if (model_type_of_domain[state_idx] == -1) { // federal state is not in any of the models (ABM, PDE, ODE) -> skip it
                        continue;
                    }

                    bool is_abm = (model_type_of_domain[state_idx] == abm_idx); 
                    if (is_abm) { //state_idx in stateIndices_ABM
                        int idx_stateIndices_ABM = local_index[state_idx];
                        nbr_agents[nr_day][idx_stateIndices_ABM] = total_nbr_individuals;
                    }

                    if (nr_day == initial_nr_day) {
                        bool is_ode = (model_type_of_domain[state_idx] == ode_idx);
                        if (is_ode) { //state_idx in stateIndices_ODE
                            int idx_stateIndices_ODE = local_index[state_idx];
                            population_ode_t_0[idx_stateIndices_ODE] = total_nbr_individuals;
                        }

                        bool is_pde = (model_type_of_domain[state_idx] == pde_idx);
                        if (is_pde) { //state_idx in stateIndices_PDE
                            int idx_stateIndices_PDE = local_index[state_idx];
                            population_pde_t_0[idx_stateIndices_PDE] = total_nbr_individuals;
                        }
                    }
                    std::cout << "total_nbr_individuals " << total_nbr_individuals << std::endl;
                    
                }
                fileTotalNumbers.close();
            } else {
                std::cerr << "Error opening the file " << filenameTotalNumbers << std::endl;
            }

            std::string filename = "../../work/input_data/global/germany_agent_id_t_0_in_ABM_" + std::to_string(nr_day) + "_" + population_scale_str + "_all_fed_states.bin";
            std::ifstream file(filename, std::ios::binary); 
            if ((file.is_open())) {
                int64_t bundesland_nr;
                for (int agent_id=0; agent_id<nbr_trajectories; agent_id++) {
                    file.read(reinterpret_cast<char*>(&bundesland_nr), sizeof(int64_t));
                    // if (ABMstate_nr >= 0 && nr_ABM_states > ABMstate_nr) { // if ABMstate_nr == -1, then it's not in an ABM state but PDE or ODE state
                    if (model_type_of_domain[bundesland_nr] == abm_idx) {
                        all_initial_agent_ids[nr_day][local_index[bundesland_nr]].push_back(agent_id);
                    }
                }
                file.close();
            } else {
                std::cerr << "Error opening the file " << filename << std::endl;
            }
        }
    }


    void readActivityChangeData(){        
        std::string filename = "../../work/input_data/global/activity_data.txt";

        std::ifstream infile(filename, std::ios::binary);

        if (!infile) {
            std::cerr << "Datei " << filename << " konnte nicht geöffnet werden" << std::endl;
            return;
        }
        // get filesize
        infile.seekg(0, std::ios::end);
        std::streampos fileSize = infile.tellg();
        infile.seekg(0, std::ios::beg);

        std::size_t numPoints = fileSize / (2 * sizeof(double));
        activityChangePercentage.resize(numPoints, std::vector<double>(2,0.0)); // first column: home, second column: notAtHome
        if (not with_restriction) { // return activityChangePercentage with zeros  
            return;
        }

        for (int i = 0; i < numPoints; i++) {
            double value1, value2; 
            infile.read(reinterpret_cast<char*>(&value1), sizeof(value1));
            infile.read(reinterpret_cast<char*>(&value2), sizeof(value2));

            if (!infile) {
                std::cerr << "Fehler beim Lesen." << std::endl;
                break;
            }

            activityChangePercentage[i][0] = value1;
            activityChangePercentage[i][1] = value2;

            std::cout << "point " << i << ", value1 " << value1 << ", value2 " << value2 << std::endl;
        }
        infile.close();
        return;
    }    

    std::tuple< int, std::unordered_map<int, int>, std::unordered_map<int, int> > filter_event_data(const std::vector<int>& initial_agent_ids, std::vector<std::vector<int>>& filtered_event_data,
                                                                                std::ranlux24_base& gen_local, int step, int nr_day, int verbosity) {
        double dt = 1.0/dt_inv;

        std::vector<int> day_shift_edu_lockdown_all(nbr_federal_states, dt * step + 1); // no lockdown. if with_restriction -> overwrite!
        if (with_restriction) {
            for (int state=0; state<nbr_federal_states; state++) {
                if (state == 15) { // Thueringen
                    day_shift_edu_lockdown_all[state] = 17; // march 17
                }
                else if (state == 4 || state == 11) { // Nordrhein-Westfalen or Brandenburg
                    day_shift_edu_lockdown_all[state] = 18; // march 18
                }
                else {
                    day_shift_edu_lockdown_all[state] = 16; // march 16
                }
            }
        }

        std::vector<std::vector<std::vector<double>>> facilityCoordinates_allCategories_nr_day = facilityCoordinates_allCategories.at(nr_day);
        std::unordered_map<int, int> facilityCoordinates_home_nr_day = facilityCoordinates_home.at(nr_day);
        auto& dayData = eventData.at(nr_day);
        std::vector<std::vector<int>> facilityLabels_allCategories_nr_day = facilityLabels_allCategories.at(nr_day);

        std::unordered_map<int, int> agentID_location_commuting_start, agentID_location_commuting_end;
        std::unordered_map<int, int> map_agent_id_idx_first_event_data;
            
        int previous_person_id = -1;
        int previous_person_id_map_first_event = -1;
        int row_idx = 0;
        bool added_extra_event = false;

        int corrected_step = step % dt_inv;

        int timestep_in_seconds = 60*30*corrected_step;
        int next_timestep_in_seconds = 60*30*(corrected_step+1);
        if (dt_inv != 48) {
            int whole_day_in_seconds = 60*60*24;
            timestep_in_seconds = whole_day_in_seconds*dt*corrected_step;
            next_timestep_in_seconds = whole_day_in_seconds*dt*(corrected_step+1);
        } 

        int dayData_size = dayData.size();
        if (step == 0) {
            filtered_event_data.resize(dayData_size, std::vector<int>(5,0));
        }

        std::uniform_real_distribution<> dis(0, 1); 

        bool previous_event_ignored = false;
        bool last_event = false;
        bool event_of_last_agent = false;

        int total_nbr_agents = initial_agent_ids.size();
        int last_agent_id = initial_agent_ids.at(total_nbr_agents-1); 
        int previous_event_ignored_agent_id = -1;

        int day_shift_edu_lockdown;
        int idx_initial_agents = 0;
        int state_nr;
        for (int i = 0; i < dayData_size; i++) {
            const auto& row = dayData.at(i);
            int time = row.at(0);
            int activity_type = row.at(1);
            int agent_id = row.at(2);
            int facility_id = row.at(3);
            int category_id = row.at(4);

            if (i == dayData_size-1) {
                last_event = true; 
            }

            if (agent_id > last_agent_id) {
                if (added_extra_event) { // make sure that extra event for previous agent was added before exiting loop
                    break;
                }
                event_of_last_agent = true;
            }

            if (previous_event_ignored) {
                previous_event_ignored = false;
                if (previous_event_ignored_agent_id == agent_id && activity_type == 1) { // activity end of same agent as before
                    continue; // make sure to ignore start and end of events
                }
            }

            if (agent_id != previous_person_id_map_first_event) { // new agent
                map_agent_id_idx_first_event_data[agent_id] = i;
                previous_person_id_map_first_event = agent_id;
            }

            if (idx_initial_agents == total_nbr_agents) { 
                break;
            }

            while ((!event_of_last_agent) && (agent_id > initial_agent_ids.at(idx_initial_agents))) { // do not skip agent without event data (always at home) 
                int skipping_agent_id = initial_agent_ids.at(idx_initial_agents);
                if (row_idx >= 1 && filtered_event_data.at(row_idx-1).at(2) != skipping_agent_id) { // no event added yet for this agent
                    auto it = facilityCoordinates_home_nr_day.find(skipping_agent_id);
                    int facility_id_skipping_agent = it->second;

                    filtered_event_data.at(row_idx).at(0) = timestep_in_seconds; //time
                    filtered_event_data.at(row_idx).at(1) = 0; // activity start
                    filtered_event_data.at(row_idx).at(2) = skipping_agent_id;
                    filtered_event_data.at(row_idx).at(3) = facility_id_skipping_agent;
                    filtered_event_data.at(row_idx).at(4) = 0; //category index home
                    row_idx++;

                    filtered_event_data.at(row_idx).at(0) = next_timestep_in_seconds; //time
                    filtered_event_data.at(row_idx).at(1) = 1; // activity end
                    filtered_event_data.at(row_idx).at(2) = skipping_agent_id;
                    filtered_event_data.at(row_idx).at(3) = facility_id_skipping_agent; //facility index
                    filtered_event_data.at(row_idx).at(4) = 0; //category index home
                    row_idx++;
                }
                else if (row_idx >= 1 && filtered_event_data.at(row_idx-1).at(2) == skipping_agent_id && filtered_event_data.at(row_idx-1).at(0) != next_timestep_in_seconds) { // there is at least one event but not at the end of the time step
                    auto it = map_agent_id_idx_first_event_data.find(skipping_agent_id);

                    int start_index_activity = i-1; 
                    int facility_id_start = dayData.at(start_index_activity).at(3);
                    int category_id_start = dayData.at(start_index_activity).at(4);

                    int end_index_agent_id = it->second;
                    int facility_id_end = dayData.at(end_index_agent_id).at(3);
                    int category_id_end = dayData.at(end_index_agent_id).at(4);

                    // commuting or not?
                    filtered_event_data.at(row_idx).at(0) = next_timestep_in_seconds; //time
                    filtered_event_data.at(row_idx).at(2) = skipping_agent_id; //person id
                    if (dayData.at(start_index_activity).at(1) == 0) { // last event is activity start -> stays at same location
                        filtered_event_data.at(row_idx).at(1) = 0; //activity type start
                        filtered_event_data.at(row_idx).at(3) = facility_id_start; //facility index
                        filtered_event_data.at(row_idx).at(4) = category_id_start; //category index        
                    }
                    else if (facility_id_start == facility_id_end && category_id_start == category_id_end) { 
                        filtered_event_data.at(row_idx).at(1) = 0; //activity type start
                        filtered_event_data.at(row_idx).at(3) = facility_id_end; //facility index
                        filtered_event_data.at(row_idx).at(4) = category_id_end; //category index
                    }
                    else { // commuting: compute location of agent at next_timestep_in_seconds by interpolating between the two activities at index end_index_agent_id and start_index_activity
                        filtered_event_data.at(row_idx).at(1) = 1; // activity end 
                        filtered_event_data.at(row_idx).at(3) = -1; //facility index unknown
                        filtered_event_data.at(row_idx).at(4) = -1; //category index unknown

                        int t_end = dayData.at(end_index_agent_id).at(0);
                        int t_start = dayData.at(start_index_activity).at(0);
                        int t = next_timestep_in_seconds; 
                        if (t_end <= t_start) { // agent is changing location during the night
                            t_end += 24*60*60;
                        }
                        if (t < t_start) { // agent is changing location during the night
                            t += 24*60*60;
                        }

                        int state_nr_start = facilityLabels_allCategories_nr_day.at(category_id_start).at(facility_id_start); // in which state did the agent start his/her journey?
                        int state_nr_end = facilityLabels_allCategories_nr_day.at(category_id_end).at(facility_id_end); 
                        int state_nr = in_which_state(state_nr_start, state_nr_end, t_start, t_end, t);
                        agentID_location_commuting_end[skipping_agent_id] = state_nr; 
                    }
                    row_idx++; 
                }
                added_extra_event = true;
                if (idx_initial_agents < total_nbr_agents - 1) {
                    idx_initial_agents++; // idx_initial_agents < total_nbr_agents -> ok
                }
                else {
                    break; // idx_initial_agents++ -> idx_initial_agents >= total_nbr_agents -> not ok! break!
                }
            }

            if (agent_id < initial_agent_ids.at(idx_initial_agents)) { // need to look for agent IDs which are in initial_agent_ids -> next agent!
                continue;
            }

            // correct event locations in case it does not make sense
            if (i+2 < dayData_size && agent_id == dayData.at(i+1).at(2) && activity_type == 1 && agent_id == dayData.at(i+2).at(2)) {
                if (facility_id != dayData.at(i+1).at(3) || category_id != dayData.at(i+1).at(4)) { // does not make sense, activity should start and end at same location
                    dayData.at(i+2).at(3) = dayData.at(i+1).at(3);
                    dayData.at(i+2).at(4) = dayData.at(i+1).at(4);
                }
                else {
                    dayData.at(i+1).at(3) = dayData.at(i+2).at(3);
                    dayData.at(i+1).at(4) = dayData.at(i+2).at(4);
                }
            }
            if (i+1 < dayData_size && agent_id == dayData.at(i+1).at(2) && activity_type == 0 && (facility_id != dayData.at(i+1).at(3) || category_id != dayData.at(i+1).at(4)) ) { // does not make sense, activity should start and end at same location
                if (i+2 < dayData_size && agent_id == dayData.at(i+2).at(2)) {
                    if (facility_id != dayData.at(i+2).at(3) || category_id != dayData.at(i+2).at(4)) {
                        dayData.at(i+1).at(3) = facility_id;
                        dayData.at(i+1).at(4) = category_id;
                    }
                    else {
                        facility_id = dayData.at(i+1).at(3);
                        category_id = dayData.at(i+1).at(4);
                    }
                }
                else { // only two consecutive events left
                    dayData.at(i+1).at(3) = facility_id;
                    dayData.at(i+1).at(4) = category_id;
                }
            }

            if ((!event_of_last_agent) && last_event && (previous_person_id != agent_id)) { //previous_person_id != agent_id means no initial event was added for this person (agent_id) because there is only one event for this person (agent_id)!
                filtered_event_data.at(row_idx).at(0) = timestep_in_seconds; //time
                filtered_event_data.at(row_idx).at(1) = 0; // activity start
                filtered_event_data.at(row_idx).at(2) = agent_id;
                filtered_event_data.at(row_idx).at(3) = facility_id;
                filtered_event_data.at(row_idx).at(4) = category_id;
                row_idx++;

                filtered_event_data.at(row_idx).at(0) = next_timestep_in_seconds; //time
                filtered_event_data.at(row_idx).at(1) = 1; // activity end
                filtered_event_data.at(row_idx).at(2) = agent_id;
                filtered_event_data.at(row_idx).at(3) = facility_id;
                filtered_event_data.at(row_idx).at(4) = category_id;
                row_idx++;
                continue;
            }

            day_shift_edu_lockdown = day_shift_edu_lockdown_all.at(facilityLabels_allCategories_nr_day.at(category_id).at(facility_id));
            if ((!event_of_last_agent) && (dt*step >= day_shift_edu_lockdown) && (category_id >= 7) && (category_id != 10)) { // category_id == 7, 8, 9, 11, 12, 13
                // gather all events for this agent_id
                std::vector<int> event_ids_of_agent_id_not_edu;
                bool agent_id_found = false;
                int first_index_agent_id = map_agent_id_idx_first_event_data.at(agent_id);
                for (int event_j = first_index_agent_id; event_j < dayData_size; event_j++) {
                    if (dayData.at(event_j).at(2) == agent_id) {
                        agent_id_found = true;
                        if (dayData.at(event_j).at(4) < 7 || dayData.at(event_j).at(4) == 10) { // event is not a school event
                            event_ids_of_agent_id_not_edu.push_back(event_j);
                        }
                    }
                    if (agent_id_found && dayData.at(event_j).at(2) != agent_id) { // next agent -> break
                        break;
                    }
                }
                if (event_ids_of_agent_id_not_edu.size() == 0) { // no event found for this agent_id besides education
                    int home_facility_id = facilityCoordinates_home_nr_day[agent_id];
                    filtered_event_data.at(row_idx).at(0) = timestep_in_seconds; //time
                    filtered_event_data.at(row_idx).at(1) = 0;
                    filtered_event_data.at(row_idx).at(2) = agent_id;
                    filtered_event_data.at(row_idx).at(3) = home_facility_id;
                    filtered_event_data.at(row_idx).at(4) = 0; //category index home
                    row_idx++;

                    filtered_event_data.at(row_idx).at(0) = next_timestep_in_seconds; //time
                    filtered_event_data.at(row_idx).at(1) = 1; // activity end
                    filtered_event_data.at(row_idx).at(2) = agent_id;
                    filtered_event_data.at(row_idx).at(3) = home_facility_id; //facility index
                    filtered_event_data.at(row_idx).at(4) = 0; //category index home
                    row_idx++;
                }
                continue; // education lockdown, ignore all events of this category
            }
            else if ((!event_of_last_agent) && (category_id != 0) && with_restriction) { // remove activity with a certain probability before the activity starts
                double probability_notAtHome = activityChangePercentage.at(int(dt*step)).at(1)/100.0; // not at home
                if (probability_notAtHome < 0) { // negative means we have to ignore some events
                    probability_notAtHome = - probability_notAtHome; // make sure it is positive (easier to understand)
                    double random_number = dis(gen_local); // pick a random number between 0 and 1
                    if (random_number <= probability_notAtHome) {
                        previous_event_ignored = true; // make sure to ignore start and end of events
                        previous_event_ignored_agent_id = agent_id;

                        std::vector<int> event_ids_of_agent_id;
                        int first_index_agent_id = map_agent_id_idx_first_event_data.at(agent_id);
                        for (int event_j = first_index_agent_id; event_j < dayData_size; event_j++) {
                            if (dayData.at(event_j).at(2) == agent_id) {
                                event_ids_of_agent_id.push_back(event_j);
                            }
                            if (dayData.at(event_j).at(2) != agent_id) { // next agent -> break 
                                break;
                            }
                        }
                        if ( (event_ids_of_agent_id.size() == 2 && dayData.at(event_ids_of_agent_id.at(0)).at(1) == 0)
                            || (event_ids_of_agent_id.size() == 1) ) { // no event found for this agent_id besides current event (and next event if its end of same event)
                            int home_facility_id = facilityCoordinates_home_nr_day[agent_id];
                            
                            filtered_event_data.at(row_idx).at(0) = timestep_in_seconds; //time
                            filtered_event_data.at(row_idx).at(1) = 0;
                            filtered_event_data.at(row_idx).at(2) = agent_id;
                            filtered_event_data.at(row_idx).at(3) = home_facility_id;    
                            filtered_event_data.at(row_idx).at(4) = 0; //category index home
                            row_idx++;

                            filtered_event_data.at(row_idx).at(0) = next_timestep_in_seconds; //time
                            filtered_event_data.at(row_idx).at(1) = 1; // activity or end
                            filtered_event_data.at(row_idx).at(2) = agent_id;
                            filtered_event_data.at(row_idx).at(3) = home_facility_id; //facility index   
                            filtered_event_data.at(row_idx).at(4) = 0; //category index home
                            row_idx++;
                        }
                        if (last_event && !added_extra_event) { // last event can not be skipped because no final event was added for current agent yet
                            agent_id = -1; // artificial "previous_person_id != agent_id" to add event at end of time step
                            event_of_last_agent = true; // if the next block (if ((!event_of_last_agent) && last_event && (previous_person_id != agent_id))) does not exist anymore: delete this row!!
                        }
                        else {
                            continue; // continue with next event
                        }
                    }
                }
            }

            if ((!event_of_last_agent) && last_event && (previous_person_id != agent_id)) { //previous_person_id != agent_id means no initial event was added for this person (agent_id) because there is only one event for this person (agent_id)!
                continue;
            }

            if ((last_event || ((previous_person_id != -1) && (previous_person_id != agent_id)))
                && !added_extra_event) { // new person but no final data was added to previous agent
                added_extra_event = true;
                
                int start_index_activity = i-1; 
                int facility_id_start = dayData.at(start_index_activity).at(3);
                int category_id_start = dayData.at(start_index_activity).at(4);

                int end_index_agent_id = map_agent_id_idx_first_event_data.at(previous_person_id); // first event of the day of previous_person_id
                int facility_id_end = dayData.at(end_index_agent_id).at(3);
                int category_id_end = dayData.at(end_index_agent_id).at(4);

                // commuting or not?
                filtered_event_data.at(row_idx).at(0) = next_timestep_in_seconds; //time
                filtered_event_data.at(row_idx).at(2) = previous_person_id; //person id
                if (dayData.at(start_index_activity).at(1) == 0) { // last event is activity start -> stays at same location
                    filtered_event_data.at(row_idx).at(1) = 0; //activity start
                    filtered_event_data.at(row_idx).at(3) = facility_id_start; //facility index
                    filtered_event_data.at(row_idx).at(4) = category_id_start; //category index
                }
                else if (facility_id_start == facility_id_end && category_id_start == category_id_end) {
                    filtered_event_data.at(row_idx).at(1) = 0; //activity start
                    filtered_event_data.at(row_idx).at(3) = facility_id_end; //facility index
                    filtered_event_data.at(row_idx).at(4) = category_id_end; //category index
                }
                else { // commuting: compute location of agent at next_timestep_in_seconds by interpolating between the two activities at index end_index_agent_id and start_index_activity
                    filtered_event_data.at(row_idx).at(1) = 1; // activity end 
                    filtered_event_data.at(row_idx).at(3) = -1; //facility index unknown
                    filtered_event_data.at(row_idx).at(4) = -1; //category index unknown

                    int t_end = dayData.at(end_index_agent_id).at(0);
                    int t_start = dayData.at(start_index_activity).at(0);
                    int t = next_timestep_in_seconds; 
                    if (t_end <= t_start) { // agent is changing location during the night
                        t_end += 24*60*60;
                    }
                    if (t < t_start) { // agent is changing location during the night
                        t += 24*60*60;
                    }

                    int state_nr_start = facilityLabels_allCategories_nr_day.at(category_id_start).at(facility_id_start); // in which state did the agent start his/her journey?
                    int state_nr_end = facilityLabels_allCategories_nr_day.at(category_id_end).at(facility_id_end); 
                    int state_nr = in_which_state(state_nr_start, state_nr_end, t_start, t_end, t);
                    agentID_location_commuting_end[previous_person_id] = state_nr;
                }
                row_idx++;  
            }

            if (event_of_last_agent || last_event) {
                break;
            }

            if (previous_person_id != agent_id) { // new person -> add initial event data 
                added_extra_event = false;

                if (time == timestep_in_seconds) { // activity is happening exactly at timestep_in_seconds: add event
                    filtered_event_data.at(row_idx).at(0) = timestep_in_seconds; //time
                    filtered_event_data.at(row_idx).at(1) = activity_type; //activity type
                    filtered_event_data.at(row_idx).at(2) = agent_id; //person id
                    filtered_event_data.at(row_idx).at(3) = facility_id; //facility index
                    filtered_event_data.at(row_idx).at(4) = category_id; //category index
                    row_idx++;
                }
                else { // (time != timestep_in_seconds)
                    filtered_event_data.at(row_idx).at(0) = timestep_in_seconds; //time
                    filtered_event_data.at(row_idx).at(2) = agent_id;

                    if ((time > timestep_in_seconds) && (activity_type == 1)) { // activity end -> stays at same location
                        filtered_event_data.at(row_idx).at(1) = 0; // activity started before
                        filtered_event_data.at(row_idx).at(3) = facility_id; //facility index is the same
                        filtered_event_data.at(row_idx).at(4) = category_id; //category index is the same
                    }
                    else {
                        int initial_activity_index = -1;
                        int next_index_activity;
                        if (time > timestep_in_seconds) { // first activity is happening after timestep_in_seconds -> find activity at the end of the day
                            int last_index_activity = -1;
                            for (int j = i; j < dayData_size; j++) {
                                if (dayData.at(j).at(2) != agent_id) {
                                    last_index_activity = j-1;
                                    break;
                                }
                            }
                            if (last_index_activity == -1) { // last agent
                                last_index_activity = dayData_size-1;
                            }
                            initial_activity_index = last_index_activity;
                            next_index_activity = i;
                        }
                        else { // activity is happening too early
                            // find last event before timestep_in_seconds
                            for (int j = i; j < dayData_size; j++) { 
                                if (dayData.at(j).at(2) != agent_id) { // next agent
                                    break;
                                }
                                if (dayData.at(j).at(0) <= timestep_in_seconds) {
                                    initial_activity_index = j;
                                }
                            }
                            next_index_activity = initial_activity_index+1;
                            if (next_index_activity >= dayData_size || dayData.at(next_index_activity).at(2) != agent_id) { // next event does not exist for this agent -> first activity of the day is the next event
                                next_index_activity = i;
                            }
                        }

                        if (dayData.at(initial_activity_index).at(0) == timestep_in_seconds) { 
                            filtered_event_data.at(row_idx).at(1) = 0; // activity start
                            filtered_event_data.at(row_idx).at(3) = dayData.at(initial_activity_index).at(3); //facility index
                            filtered_event_data.at(row_idx).at(4) = dayData.at(initial_activity_index).at(4); //category index
                        }
                        else if (dayData.at(initial_activity_index).at(1) == 0 && dayData.at(initial_activity_index).at(3) == dayData.at(next_index_activity).at(3) && dayData.at(initial_activity_index).at(4) == dayData.at(next_index_activity).at(4)) {
                            filtered_event_data.at(row_idx).at(1) = 0; // activity start
                            filtered_event_data.at(row_idx).at(3) = dayData.at(initial_activity_index).at(3); //facility index
                            filtered_event_data.at(row_idx).at(4) = dayData.at(initial_activity_index).at(4); //category index
                        }
                        else {
                            if (next_index_activity == initial_activity_index) { // there is only one event -> initial_activity_index == i
                                filtered_event_data.at(row_idx).at(1) = 0; // activity start
                                filtered_event_data.at(row_idx).at(3) = facility_id; //facility index
                                filtered_event_data.at(row_idx).at(4) = category_id; //category index
                            }
                            else {
                                // compute location of agent at timestep_in_seconds by interpolating between the two activities at index initial_activity_index and next_index_activity
                                filtered_event_data.at(row_idx).at(1) = 1; // activity end 
                                filtered_event_data.at(row_idx).at(3) = -1; //facility index unknown
                                filtered_event_data.at(row_idx).at(4) = -1; //category index unknown

                                const auto& row_start = dayData.at(initial_activity_index);
                                int facility_id_start = row_start.at(3);
                                int category_id_start = row_start.at(4);
                                int t_start = row_start.at(0);
                            
                                const auto& row_end = dayData.at(next_index_activity);
                                int facility_id_end = row_end.at(3);
                                int category_id_end = row_end.at(4);
                                int t_end = row_end.at(0);

                                int t = timestep_in_seconds; 
                                if (t_end <= t_start) { // agent is changing location during the night
                                  t_end += 24*60*60;
                                }
                                if (t < t_start) { // agent is changing location during the night
                                  t += 24*60*60;
                                }

                                int state_nr_start = facilityLabels_allCategories_nr_day.at(category_id_start).at(facility_id_start); // in which state did the agent start his/her journey?
                                int state_nr_end = facilityLabels_allCategories_nr_day.at(category_id_end).at(facility_id_end); 
                                int state_nr = in_which_state(state_nr_start, state_nr_end, t_start, t_end, t);
                                agentID_location_commuting_start[agent_id] = state_nr;
                            }
                        }
                    }
                    row_idx++;
                }
                previous_person_id = agent_id;
            }

            if (time < next_timestep_in_seconds && time >= timestep_in_seconds) { // event falls into the time window! add event !
                filtered_event_data.at(row_idx).at(0) = time; //time
                filtered_event_data.at(row_idx).at(1) = activity_type; //activity type
                filtered_event_data.at(row_idx).at(2) = agent_id; //person id
                filtered_event_data.at(row_idx).at(3) = facility_id; //facility index
                filtered_event_data.at(row_idx).at(4) = category_id; //category index
                row_idx++;
            }
            else if (time == next_timestep_in_seconds && (!added_extra_event)) { // perfect timing 
                added_extra_event = true;
                filtered_event_data.at(row_idx).at(0) = time; //time
                filtered_event_data.at(row_idx).at(1) = activity_type; //activity type
                filtered_event_data.at(row_idx).at(2) = agent_id; //person id
                filtered_event_data.at(row_idx).at(3) = facility_id; //facility index
                filtered_event_data.at(row_idx).at(4) = category_id; //category index
                row_idx++;
            }
            else if ((time > next_timestep_in_seconds) && (!added_extra_event)) { // no event added yet for this agent besides the artificial event at time timestep_in_seconds and next will be after next_timestep_in_seconds
                added_extra_event = true;

                filtered_event_data.at(row_idx).at(0) = next_timestep_in_seconds; //time
                filtered_event_data.at(row_idx).at(2) = agent_id; //person id
                if (activity_type == 1) { // activity end after next_timestep_in_seconds: agent stays at same location
                    filtered_event_data.at(row_idx).at(1) = 0; //activity type start
                    filtered_event_data.at(row_idx).at(3) = facility_id; //facility index
                    filtered_event_data.at(row_idx).at(4) = category_id; //category index
                }
                else { // commuting: compute final location
                    filtered_event_data.at(row_idx).at(1) = 1; // activity end 
                    filtered_event_data.at(row_idx).at(3) = -1; //facility index unknown
                    filtered_event_data.at(row_idx).at(4) = -1; //category index unknown

                    int facility_id_end = facility_id;
                    int category_id_end = category_id;

                    int idx_previous_event;
                    if (dayData.at(i-1).at(2) == agent_id) { // previous event is from the same agent
                        idx_previous_event = i-1;
                    }
                    else { // no previous event
                        // get last event of the day
                        int last_index_activity = -1;
                        for (int j = i; j < dayData_size; j++) {
                            if (dayData.at(j).at(2) != agent_id) {
                                last_index_activity = j-1;
                                break;
                            }
                        }
                        if (last_index_activity == -1) { // last agent
                            last_index_activity = dayData_size-1; 
                        }
                        idx_previous_event = last_index_activity;
                    }

                    int previous_event_facility_id = dayData.at(idx_previous_event).at(3);
                    int previous_event_category_id = dayData.at(idx_previous_event).at(4);

                    int t_end = time;
                    int t_start = dayData.at(idx_previous_event).at(0);
                    int t = next_timestep_in_seconds;
                    if (t_end <= t_start) { // agent is changing location during the night
                        t_end += 24*60*60;
                    }
                    if (t < t_start) { // agent is changing location during the night
                        t += 24*60*60;
                    }

                    int state_nr_start = facilityLabels_allCategories_nr_day.at(previous_event_category_id).at(previous_event_facility_id); // in which state did the agent start his/her journey?
                    int state_nr_end = facilityLabels_allCategories_nr_day.at(category_id_end).at(facility_id_end); 
                    int state_nr = in_which_state(state_nr_start, state_nr_end, t_start, t_end, t);
                    agentID_location_commuting_end[agent_id] = state_nr;
                }
                row_idx++;
            }
        }
        return std::make_tuple(row_idx, agentID_location_commuting_start, agentID_location_commuting_end);
    }
    
    std::tuple< std::vector<std::unordered_map<int, int>>, std::vector<int>, std::vector<std::vector<int>>, std::vector<std::vector<int>> > step(const std::vector<std::vector<int>>& event_data, int filtered_event_data_size, std::vector<std::vector<bool>> with_coupling_states, std::unordered_map<int, int>& agentID_location_commuting_start, std::unordered_map<int, int>& agentID_location_commuting_end, std::unordered_map<int, int> agentID_healthStatus, int& step, std::vector<double>& infection_rates, int& nr_day, double& population_scale, std::ranlux24_base& gen_local, int& verbosity, int run, int correction_step) { 
        double dt = 1.0/dt_inv;
        int corrected_step = step % dt_inv;

        double sigma = 1.0/3.5, gamma = 1.0/2, eta = 1.0/4, kappa = 1.0, eta_c = 1.0/21;
        double phi_i = 1.0/4, phi_sy = 1.0/8, phi_h = 1.0/14, phi_hc = 1.0/7;

        std::uniform_real_distribution<double> dis_uniform(0.0, 1.0);

        ////////////////////////////////////////
        // define probabilities for changing status
        double prob_exposed_changing_state = 1 - std::exp(- (sigma) * dt);
        double prob_infectious_not_symp_changing_state = 1 - std::exp(- (phi_i+gamma) * dt);
        double prob_infectious_symptomatic_changing_state = 1 - std::exp(- (phi_sy+eta) * dt);
        double prob_req_hospitalization_changing_state = 1 - std::exp(- (phi_h+kappa) * dt);
        double prob_critical_changing_state = 1 - std::exp(- (eta_c) * dt);
        double prob_req_hospitalization_after_critical_changing_state = 1 - std::exp(- (phi_hc) * dt);

        const auto& facilityLabels_allCategories_day = facilityLabels_allCategories[nr_day];
        std::unordered_map<int, int> map_agent_id_idx_first_filtered_event_data; // first event of agent in timestep of filtered_event_data

        int previous_agent_id = -1;
        int end_of_day_in_seconds = 60*60*24; // 24 hours * 60 minutes * 60 seconds
        int end_of_day_in_seconds_dt = 60*30;
        if (dt_inv != 48) {
            end_of_day_in_seconds_dt = end_of_day_in_seconds*dt;
        }
        int start_time_dt = end_of_day_in_seconds_dt * corrected_step;
        int end_time_dt = end_of_day_in_seconds_dt * (corrected_step + 1);
        int start_time_total = end_of_day_in_seconds;
        int end_time_total = end_of_day_in_seconds * (1 + dt);

        for (int i = 0; i < filtered_event_data_size; i++) {
            const auto& row = event_data[i];
            int agent_id = row[2];
            if (agent_id != previous_agent_id) { // new agent
                map_agent_id_idx_first_filtered_event_data[agent_id] = i;
                previous_agent_id = agent_id;
            }
        }

        /////////        
        // divide agents into those that are in a container and those that are changing their location (every event that happens in this timestep! agent ids may be repeated)
        const auto& facilityCoordinates_allCategories_nr_day = facilityCoordinates_allCategories[nr_day]; 
        std::vector<std::vector<std::vector<std::vector<int>>>> containers(nbr_categories); // contains all agent ids in specific container. dimensions: category nr of facilities, facility nr, nbr of agents (id) in category of facility
        for (int i = 0; i < filtered_event_data_size; i++) { 
            const auto& row = event_data[i];
            int agent_id = row[2];
            int healthState = agentID_healthStatus.at(agent_id);
            if (healthState == 1 || healthState > 3) { // not susceptible or infectious
                continue;
            }

            int category_index = row[4];
            if (category_index == -1) { //agent is changing location -> there is no container
                continue; // with next event
            }    

            int facility_index = row[3];
            int time = row[0];

            int time_leaving, time_entering;
            if (event_data[i][1] == 0) { // activity start -> we know where the agent is. add to container
                time_entering = time;            
                int index_event_data = -1;
                int first_index_agent_id = map_agent_id_idx_first_filtered_event_data.at(agent_id);
                for (int j = first_index_agent_id; j < filtered_event_data_size; j++) {
                    if ((event_data[j][2] == agent_id) && (event_data[j][0] == time_entering)) {
                        index_event_data = j;
                        break;
                    }
                }
                time_leaving = 24*60*60; // assuming that the activity ends at the end of the day 
                if ((index_event_data != -1) && (index_event_data+1 < filtered_event_data_size) && (event_data[index_event_data+1][2] == agent_id)) {
                    time_leaving = event_data[index_event_data+1][0]; // time when activity ends is given -> overwrite
                }
                else {
                    time_leaving = event_data[first_index_agent_id][0]; // check beginning of day
                }
            }
            else { // activity end (1)
                time_leaving = time;   
                time_entering = 0;
                if ((i>0) && (i<filtered_event_data_size) && (event_data[i-1][2] == agent_id)) {
                    time_entering = event_data[i-1][0];
                }
                else {
                    int index_event_data = -1;
                    int first_index_agent_id = map_agent_id_idx_first_filtered_event_data.at(agent_id);
                    for (int j = first_index_agent_id; j < filtered_event_data_size; j++) {
                        if ((event_data[j][2] == agent_id) && (event_data[j][0] == time_leaving)) {
                            index_event_data = j-1;
                            break;
                        }
                    }
                    if ((index_event_data != -1) && (event_data[index_event_data][2] == agent_id)) { 
                        time_entering = event_data[index_event_data][0];
                    }
                    else { // agent is in container at midnight
                        // get last event index of agent
                        int last_index_agent_id = first_index_agent_id;
                        for (int j = first_index_agent_id; j < filtered_event_data_size; j++) {
                            if (event_data[j][2] == agent_id) {
                                last_index_agent_id = j;
                            }
                            if (event_data[j][2] != agent_id) {
                                break;
                            }
                        }
                        time_entering = event_data[last_index_agent_id][0];
                    }
                }
            }
            if (containers.at(category_index).empty()) {
                containers.at(category_index).resize(facilityCoordinates_allCategories_nr_day.at(category_index).size());
            }
            containers.at(category_index).at(facility_index).push_back({agent_id, time_entering, time_leaving});
        }
        std::cout << "after containers" << std::endl;

        // find out whether an agent leaves a container, for checking whether the health status changes or not!
        // deep copy health state map
        std::unordered_map<int, int> agentID_healthStatus_next = agentID_healthStatus;
        for (int i = 0; i < filtered_event_data_size; i++) {
            int agent_id = event_data[i][2];
            int healthState = agentID_healthStatus.at(agent_id);

            if (healthState > 0) { // not susceptible 
                continue;
            }

            if (event_data[i][1] == 1) { // activity ends
                int category_index = event_data[i][4];
                int facility_index = event_data[i][3];
                if ((category_index == -1)) { //missing information: agent is changing location
                    continue;
                }

                int nbr_facilities = containers.at(category_index).size();
                int nbr_agents_in_container = containers.at(category_index).at(facility_index).size();
                int index_of_susceptible_in_container = -1;      
                for (int i_container = 0; i_container < nbr_agents_in_container; i_container++) {
                    int agent_id_i = containers.at(category_index).at(facility_index)[i_container][0];
                    if (agent_id == agent_id_i) {
                        index_of_susceptible_in_container = i_container;
                        break;
                    }
                }

                int entering_container = containers.at(category_index).at(facility_index)[index_of_susceptible_in_container][1];
                int leaving_container = containers.at(category_index).at(facility_index)[index_of_susceptible_in_container][2];
                if (leaving_container <= entering_container) {
                    leaving_container += 24*60*60; // activity start at the next day
                }
                int total_duration_of_interaction = 0;
                std::vector<int> durations_of_interaction;
                int cnt_infectious_in_container = 0;  
                std::vector<bool> is_aware_of_being_infected;
                for (int i_container = 0; i_container < nbr_agents_in_container; i_container++) {
                    int agent_id_i = containers.at(category_index).at(facility_index)[i_container][0];
                    healthState = agentID_healthStatus.at(agent_id_i);

                    if (healthState == 2 || healthState == 3) { // infectious
                        cnt_infectious_in_container++;
                        if (healthState == 2) {
                            is_aware_of_being_infected.push_back(false);
                        }
                        else {
                            is_aware_of_being_infected.push_back(true);
                        }
                    }
                    else {
                        continue; // not infectious
                    }
                    int duration_of_interaction = 0;
                    int entering_container_i = containers.at(category_index).at(facility_index)[i_container][1];
                    int leaving_container_i = containers.at(category_index).at(facility_index)[i_container][2];
                    if (leaving_container_i <= entering_container_i) {
                        leaving_container_i += 24*60*60; // activity start at the next day
                    }
                    if (leaving_container >= entering_container_i
                        && entering_container <= leaving_container_i) { // if they intersect
                        int earliest_start = std::max(entering_container, entering_container_i);
                        int latest_end = std::min(leaving_container, leaving_container_i);
                        duration_of_interaction = latest_end - earliest_start;
                    }
                    durations_of_interaction.push_back(duration_of_interaction);
                    total_duration_of_interaction += duration_of_interaction;
                }
                if (cnt_infectious_in_container == 0) {
                    continue; // no infectious in container
                }

                // pick, which masks are worn
                // probability of picking FFP2 with sh=0.6 is 90% and cloth mask with sh=0.15 is 10%
                double sh = 1.0, in = 1.0; // sh shedding, in intake
                double rand_val;
                double exponent = total_duration_of_interaction * sh * in;
                int day_shift_mask_wearing = day_shift_mask_wearing_all.at(facilityLabels_allCategories_day.at(category_index).at(facility_index));
            
                if (end_of_day_in_seconds_dt*step > day_shift_mask_wearing && (category_index == 2 || category_index == 5)) { // masks are only worn in shops and after April 29th
                    exponent = 0.0;
                    for (int idx_infectious = 0; idx_infectious < cnt_infectious_in_container; idx_infectious++) {
                        rand_val = dis_uniform(gen_local);
                        if (rand_val <= 0.9) {
                            sh = 0.6;
                            in = 0.5;
                        }
                        else {
                            sh = 0.15;
                            in = 0.025;
                        }
                        exponent += durations_of_interaction[idx_infectious] * sh * in;
                    } 
                }

                double ci_normalized; // "normalized" contact intensity 
                double max_number_of_persons_in_facility = max_agents_per_facility_category[nr_day][std::make_pair(facility_index, category_index)];
                double number_spaces_per_facility = 20.0; // will be overwritten by 1.0 if home

                // categories: ['home' 'leisure' 'shop_daily' 'errands' 'business' 'shop_other' 'work'
                //              'educ_higher' 'educ_other' 'educ_kiga' 'visit' 'educ_secondary'
                //              'educ_primary' 'educ_tertiary']
                if (category_index == 0) { // home
                    ci_normalized = 1.0;
                    number_spaces_per_facility = 1.0;
                }
                else if (category_index == 1 || category_index == 10) { // leisure, visit
                    ci_normalized = 9.24;
                }
                else if (category_index == 2 || category_index == 5) { // shop_daily, shop_other
                    ci_normalized = 0.88;
                }
                else if (category_index == 3 || category_index == 4 || category_index == 6) { // errands, business, work
                    ci_normalized = 1.47;
                }
                else if (category_index == 7) { // educ_higher 
                    ci_normalized = 5.5;
                }
                else { 
                    ci_normalized = 11.0;
                }
                
                double max_number_of_persons_in_room = max_number_of_persons_in_facility / number_spaces_per_facility;
                double ci = ci_normalized / max_number_of_persons_in_room; // contact intensity
                double random_number = dis_uniform(gen_local);
                int state_nr = facilityLabels_allCategories_day.at(category_index).at(facility_index);
                int ABM_state_nr;
                if (model_type_of_domain[state_nr] == abm_idx) {
                    ABM_state_nr = local_index[state_nr];
                }
                else {
                    continue;
                }
                
                if (random_number <= (1 - std::exp(- (infection_rates[ABM_state_nr]*ci*exponent)))) { // change happens
                    agentID_healthStatus_next[agent_id] = 1;
                }
            }
        }
        std::cout << "after health status change of susceptible" << std::endl;

        std::vector<int> ids;
        ids.reserve(agentID_healthStatus.size());
        for (auto const& kv : agentID_healthStatus)
            ids.push_back(kv.first);
        std::sort(ids.begin(), ids.end());

        // health status change of everyone else
        for (int agent_id : ids) {
            int healthStatus_i = agentID_healthStatus[agent_id];
            if (healthStatus_i == 0 || healthStatus_i == 7) { // susceptible was already checked above, removed does not change
                continue;
            }

            // not susceptible -> iterate over every agent only once
            double random_number = dis_uniform(gen_local);
            if (healthStatus_i == 1) { // exposed
                if (random_number <= prob_exposed_changing_state) {
                    agentID_healthStatus_next[agent_id] = 2;  // new status infectious but not symptomatic
                }
            }
            else if (healthStatus_i == 2) { // infectious but not symptomatic
                if (random_number <= prob_infectious_not_symp_changing_state) { // change happens
                    random_number = dis_uniform(gen_local);
                    if (random_number < gamma/(phi_i+gamma)) {
                        agentID_healthStatus_next[agent_id] = 3; // new status infectious and symptomatic
                    }
                    else {
                        agentID_healthStatus_next[agent_id] = 7; // new status removed
                    }
                }
            }
            else if (healthStatus_i == 3) { // infectious and symptomatic
                if (random_number <= prob_infectious_symptomatic_changing_state) { // change happens
                    random_number = dis_uniform(gen_local);
                    if (random_number < eta/(phi_sy+eta)) {
                        agentID_healthStatus_next[agent_id] = 4; // new status requires hospitalization
                    }
                    else {
                        agentID_healthStatus_next[agent_id] = 7; // new status removed
                    }
                }
            }
            else if (healthStatus_i == 4) { // requires hospitalization
                if (random_number <= prob_req_hospitalization_changing_state) { // change happens
                    random_number = dis_uniform(gen_local);
                    if (random_number < kappa/(phi_h+kappa)) {
                        agentID_healthStatus_next[agent_id] = 5; // new status critical
                    }
                    else {
                        agentID_healthStatus_next[agent_id] = 7; // new status removed
                    }
                }
            }
            else if (healthStatus_i == 5) { // critical
                if (random_number <= prob_critical_changing_state) { // change happens
                    agentID_healthStatus_next[agent_id] = 6; // new status requires hospitalization after being critical
                }
            }
            else if (healthStatus_i == 6) { // requires hospitalization after being critical
                if (random_number <= prob_req_hospitalization_after_critical_changing_state) { // change happens
                    agentID_healthStatus_next[agent_id] = 7; // new status removed
                }
            }
        }
        std::cout << "after health status change of everyone else" << std::endl;
        
        ///////////////
        // compute final location of agents and check if they are in PDE/ODE domain or not
        std::vector<int> agent_IDs_in_ABM;
        std::vector<std::vector<double>> agents_in_ABM_position;
        std::vector<std::unordered_map<int, int>> healthStatus_agents_in_ABM(nr_ABM_states);
        std::vector<std::vector<int>> healthStatus_agents_in_PDE(nr_PDE_states, std::vector<int>(nbr_compartments,0));
        std::vector<std::vector<int>> healthStatus_agents_in_ODE(nr_ODE_states, std::vector<int>(nbr_compartments,0));

        bool in_PDE = false, in_ODE = false;
        int bundesland_nr;
        int stateIdxSpecificModelType = -1;
        int agent_id;

        if ( filtered_event_data_size > 0 &&  // during the night there are sometimes timesteps without events
            ( (nr_ABM_states > 1)
               || 
              ((nr_PDE_states + nr_ODE_states) > 0 && (nr_ABM_states == 1)) ) ) {
            previous_agent_id = event_data[0][2];
            for (int i = 0; i < filtered_event_data_size+1; i++) {
                if (i < filtered_event_data_size) {
                    agent_id = event_data[i][2];
                }
                else { //(i == filtered_event_data_size) {
                    agent_id = -1; // to process the last agent (and because event_data[i][2] for i=filtered_event_data_size does not exist!)
                }

                if (previous_agent_id != agent_id) { // new agent
                    //alten agenten zu Ende bearbeiten
                    
                    // get coordinates of agent
                    // search for previous_agent_id and time in seconds in event_data 
                    int index_event_data = i-1; 
                    int category_index = event_data[index_event_data][4];
                    int facility_index = event_data[index_event_data][3];

                    double x = -1.0;
                    double y = -1.0;
                    if (category_index == -1) { // missing information: agent is changing location. 
                        // artificial events are at .. 
                        if (event_data[index_event_data][0] == end_of_day_in_seconds_dt*(corrected_step+1)) {
                            bundesland_nr = agentID_location_commuting_end[previous_agent_id];
                        }
                        else { // if (event_data[index_event_data][0] == end_of_day_in_seconds_dt*(corrected_step)) {
                            bundesland_nr = agentID_location_commuting_start[previous_agent_id];
                        }

                        int model_type = model_type_of_domain[bundesland_nr];
                        in_PDE = (model_type == pde_idx);
                        in_ODE = (model_type == ode_idx);
                        stateIdxSpecificModelType = local_index[bundesland_nr]; 
                    }
                    else { 
                        if ((verbosity == 2) && (step % (8*3) == 0)) {
                            x = facilityCoordinates_allCategories_nr_day.at(category_index).at(facility_index)[0];
                            y = facilityCoordinates_allCategories_nr_day.at(category_index).at(facility_index)[1];
                        }
                        bundesland_nr = facilityLabels_allCategories_day.at(category_index).at(facility_index);
                        int model_type = model_type_of_domain[bundesland_nr];
                        in_PDE = (model_type == pde_idx);
                        in_ODE = (model_type == ode_idx);
                        stateIdxSpecificModelType = local_index[bundesland_nr]; 
                    }
                    
                    if (in_ODE && with_coupling_states[ode_idx].at(stateIdxSpecificModelType)) { // is inside the ODE domain and transmission into that ODE domain is allowed
                        healthStatus_agents_in_ODE[stateIdxSpecificModelType][agentID_healthStatus_next[previous_agent_id]]++;
                    } 
                    else if (in_PDE && with_coupling_states[pde_idx].at(stateIdxSpecificModelType)) { // is inside the PDE domain and transmission into that PDE domain is allowed
                        healthStatus_agents_in_PDE[stateIdxSpecificModelType][agentID_healthStatus_next[previous_agent_id]]++;
                    } 
                    else { // bleibt im ABM als agent
                        if (stateIdxSpecificModelType < 0) { // federal state is not modeled at all
                            stateIdxSpecificModelType = 0; // add agents outside of modeling domain to one of the ABM states (e.g. state 0) so that they are not lost and can be tracked in the ABM, and also possible infection events are not lost
                        }
                        agent_IDs_in_ABM.push_back(previous_agent_id); 
                        healthStatus_agents_in_ABM[stateIdxSpecificModelType][previous_agent_id] = agentID_healthStatus_next[previous_agent_id];

                        if ((verbosity == 2) && (step % (8*3) == 0)) {
                            if (x>0 && y>0) { 
                                agents_in_ABM_position.push_back({x,y});
                            }
                        } 
                    }  
                }
                previous_agent_id = agent_id;
            }
        }
        else { // all agents stay in ABM domain
            for (const auto& kv : agentID_healthStatus_next) {
                agent_IDs_in_ABM.push_back(kv.first);
                healthStatus_agents_in_ABM[0][kv.first] = kv.second;
            }
        }
        
        // ///////////////
        if ((verbosity == 2) && (step % (8*3) == 0) && (run == 0)) {
            std::string filename_position = "../../work/output/output_data_optim_" + std::to_string(correction_step) + "/agents_position_0.txt";
            std::ofstream outfile_position;
            if (step == 0) { 
                outfile_position.open(filename_position);  // Datei überschreiben im ersten Durchlauf
            }
            else {
                outfile_position.open(filename_position, std::ios::app);  // Datei im Anhängemodus öffnen 
            }
            outfile_position << agents_in_ABM_position.size() << " " << 1 << "\n";
            for (int i = 0; i < agents_in_ABM_position.size(); i++) {
                outfile_position << agents_in_ABM_position[i][0] << " " << agents_in_ABM_position[i][1] << "\n";
            }
            outfile_position.close();


            // save additionally the status of agents
            std::string filename_healthStatus = "../../work/output/output_data_optim_" + std::to_string(correction_step) + "/agents_healthStatus_0.txt";
            std::ofstream outfile_healthStatus;
            if (step == 0) { 
                outfile_healthStatus.open(filename_healthStatus);  // Datei überschreiben im ersten Durchlauf
            }
            else {
                outfile_healthStatus.open(filename_healthStatus, std::ios::app);  // Datei im Anhängemodus öffnen 
            }
            size_t total_entries = 0;
            for (int state_ABM_idx=0; state_ABM_idx<nr_ABM_states; state_ABM_idx++) {
                total_entries += healthStatus_agents_in_ABM.at(state_ABM_idx).size();
            }
            outfile_healthStatus << total_entries << "\n";
            for (int state_ABM_idx=0; state_ABM_idx<nr_ABM_states; state_ABM_idx++) {
                for (const auto& agent_id_healthState : healthStatus_agents_in_ABM.at(state_ABM_idx)) {
                    int healthState = agent_id_healthState.second;
                    outfile_healthStatus << healthState << "\n";
                }
            }
            outfile_healthStatus.close();
        } // end if verbosity

        return std::make_tuple(healthStatus_agents_in_ABM, agent_IDs_in_ABM, healthStatus_agents_in_PDE, healthStatus_agents_in_ODE);
    }
};