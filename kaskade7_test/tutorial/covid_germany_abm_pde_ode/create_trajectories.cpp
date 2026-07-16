#ifndef TRAJECTORIES_HH
#define TRAJECTORIES_HH

#include <iostream>
#include <string> //
#include <cstdlib> // std::exit(EXIT_FAILURE)

#include <chrono>
// #include <thread>
#include <filesystem>

#include <zlib.h>

#if HAVE_LIBAMIRA == 1
#include "io/amira.hh"
#endif

#include "ABM_threadsafe.cpp"

std::tuple<std::vector<std::vector<std::vector<double>>>, std::vector<std::vector<int>>> create_trajectories(int nbr_trajectories, std::vector<std::vector<int>> trajectories_nr_state_previous, 
  std::vector<std::vector<std::vector<int>>> eventData, std::vector<std::vector<std::vector<std::vector<double>>>> facilityCoordinates_allCategories_all_days, std::vector<std::unordered_map<int, int>> facilityCoordinates_home_all_days, 
  std::vector<std::vector<std::vector<int>>> facilityLabels_allCategories_all_days, int dt_inv, int step, int nr_day, std::vector<std::vector<std::vector<int>>>& jump_matrix, std::vector<std::vector<std::vector<std::vector<int>>>>& jump_matrix_agent_IDs) {

  auto check_bounds = [&](int prev_state_nr, int state_nr, int step) {
      if (prev_state_nr < 0 || prev_state_nr >= nbr_domains ||
          state_nr      < 0 || state_nr      >= nbr_domains) {
          std::cout << "INVALID INDEX in jump_matrix: "
                    << "step=" << step
                    << ", prev_state_nr=" << prev_state_nr
                    << ", state_nr=" << state_nr << std::endl;
          std::exit(EXIT_FAILURE);
      }
  };

  double dt = 1.0/dt_inv;

  std::vector<std::vector<std::vector<double>>> facilityCoordinates_allCategories = facilityCoordinates_allCategories_all_days.at(nr_day);
  std::unordered_map<int, int> facilityCoordinates_home = facilityCoordinates_home_all_days.at(nr_day);
  auto& dayData = eventData.at(nr_day); 
  std::vector<std::vector<int>> facilityLabels_allCategories_nr_day = facilityLabels_allCategories_all_days[nr_day];

  std::vector<std::vector<std::vector<double>>> trajectories(nbr_trajectories, std::vector<std::vector<double>>(2, std::vector<double>(2,-1.0))); 
  std::unordered_map<int, int> map_agent_id_idx_first_event_data;

  std::vector<std::vector<int>> trajectories_nr_state_storage;  
  std::vector<std::vector<int>>* trajectories_nr_state = nullptr;
  int timestep_check;
  if (step == 0) {
      // trajectories_nr_state_previous und trajectories_nr_state zeigen auf den selben Speicher, weil für step 0 es kein trajectories_nr_state_previous gibt
      trajectories_nr_state = &trajectories_nr_state_previous; 
      timestep_check = 0;
  } else {
      // eigenes Objekt
      trajectories_nr_state = &trajectories_nr_state_storage;
      timestep_check = 1; // fuer trajectories_nr_state_previous
  }
  trajectories_nr_state->resize(nbr_trajectories, std::vector<int>(2,-1));

    
  int previous_person_id = -1;
  int previous_person_id_map_first_event = -1;
  bool added_extra_event = false;

  int corrected_step = step % dt_inv;
  int whole_day_in_seconds = 60*60*24;

  int timestep_in_seconds = 60*30*corrected_step;
  int next_timestep_in_seconds = 60*30*(corrected_step+1);
  if (dt_inv != 48) {  
    timestep_in_seconds = whole_day_in_seconds*dt*corrected_step;
    next_timestep_in_seconds = whole_day_in_seconds*dt*(corrected_step+1);
  } 

  bool last_event = false;
  bool event_of_last_agent = false;

  int last_agent_id = nbr_trajectories-1;
  std::cout << "last_agent_id " << last_agent_id << std::endl;

  int idx_initial_agents = 0;
  int state_nr;
  double x,y;
  int dayData_size = dayData.size();
  std::cout << "dayData_size " << dayData_size << std::endl;
  for (int i = 0; i < dayData_size; i++) { // iterate over event data
    const auto& row = dayData.at(i);
    int agent_id = row.at(2);

    if (i == dayData_size-1) {
      last_event = true;
    }

    if (agent_id > last_agent_id) {
        if (added_extra_event) { // make sure that extra event for previous agent was added before exiting loop
            break;
        }
        event_of_last_agent = true;
    }

    if (agent_id != previous_person_id_map_first_event) { // new agent
      map_agent_id_idx_first_event_data[agent_id] = i;
      previous_person_id_map_first_event = agent_id;
    }

    if (idx_initial_agents == nbr_trajectories) { 
      break;
    }

    while ((!event_of_last_agent) && (agent_id > idx_initial_agents)) { // do not skip agent without event data (always at home) 
      int skipping_agent_id = idx_initial_agents;
      if (trajectories[skipping_agent_id][0][0] < 0.0) { // no event added yet for this agent
        auto it = facilityCoordinates_home.find(skipping_agent_id);
        if (it == facilityCoordinates_home.end()) {
            std::cout << "[ERR] no home for agent " << skipping_agent_id << std::endl;
            std::exit(EXIT_FAILURE);
        }
        int facility_id_skipping_agent = it->second;
        
        x = facilityCoordinates_allCategories[0][facility_id_skipping_agent][0]; //category_id=0 => home
        y = facilityCoordinates_allCategories[0][facility_id_skipping_agent][1]; //category_id=0 => home

        state_nr = facilityLabels_allCategories_nr_day.at(0).at(facility_id_skipping_agent);
        trajectories[skipping_agent_id][0] = std::vector<double>{x, y}; 
        (*trajectories_nr_state)[skipping_agent_id][0] = state_nr; //number of state in Germany

        trajectories[skipping_agent_id][1] = std::vector<double>{x, y}; 
        (*trajectories_nr_state)[skipping_agent_id][1] = state_nr; //number of state in Germany
      }
      else if (trajectories[skipping_agent_id][1][0] < 0.0) { // final event not added
        int end_index_agent_id = map_agent_id_idx_first_event_data.at(skipping_agent_id); // first event of the day of skipping_agent_id
        int start_index_activity = i-1; // last event of skipping_agent_id

        int facility_id_start = dayData.at(start_index_activity).at(3);
        int category_id_start = dayData.at(start_index_activity).at(4);

        int facility_id_end = dayData.at(end_index_agent_id).at(3);
        int category_id_end = dayData.at(end_index_agent_id).at(4);

        if (dayData.at(start_index_activity).at(1) == 0) { // last event is activity start -> stays at same location
          x = facilityCoordinates_allCategories[category_id_start][facility_id_start][0]; 
          y = facilityCoordinates_allCategories[category_id_start][facility_id_start][1]; 
          state_nr = facilityLabels_allCategories_nr_day.at(category_id_start).at(facility_id_start);

          trajectories[skipping_agent_id][1] = std::vector<double>{x, y}; 
          (*trajectories_nr_state)[skipping_agent_id][1] = state_nr; //number of state in Germany 

          int prev_state_nr = trajectories_nr_state_previous.at(skipping_agent_id).at(timestep_check);
          if (prev_state_nr > -1 && prev_state_nr != state_nr) { 
            check_bounds(prev_state_nr, state_nr, corrected_step);
            jump_matrix[corrected_step][prev_state_nr][state_nr] += 1;
            jump_matrix_agent_IDs[corrected_step][prev_state_nr][state_nr].push_back(skipping_agent_id);
          }
        }
        else if (facility_id_start == facility_id_end && category_id_start == category_id_end) { // same facility and category -> same location during the night
          x = facilityCoordinates_allCategories[category_id_end][facility_id_end][0]; 
          y = facilityCoordinates_allCategories[category_id_end][facility_id_end][1]; 
          state_nr = facilityLabels_allCategories_nr_day.at(category_id_end).at(facility_id_end);
          trajectories[skipping_agent_id][1] = std::vector<double>{x, y}; 
          (*trajectories_nr_state)[skipping_agent_id][1] = state_nr; //number of state in Germany

          int prev_state_nr = trajectories_nr_state_previous.at(skipping_agent_id).at(timestep_check);
          if (prev_state_nr > -1 && prev_state_nr != state_nr) {
            check_bounds(prev_state_nr, state_nr, corrected_step);
            jump_matrix[corrected_step][prev_state_nr][state_nr] += 1;
            jump_matrix_agent_IDs[corrected_step][prev_state_nr][state_nr].push_back(skipping_agent_id);
          }
        }
        else { // commuting: compute location of agent at next_timestep_in_seconds by interpolating between the two activities at index end_index_agent_id and start_index_activity
          const auto& facilities_end = facilityCoordinates_allCategories.at(category_id_end).at(facility_id_end);
          double x_end = facilities_end.at(0);
          double y_end = facilities_end.at(1);

          double x_start, y_start;
          if (category_id_start == -1) { 
            x_start = trajectories[skipping_agent_id][0][0];
            y_start = trajectories[skipping_agent_id][0][1];
          }
          else {
            const auto& facilities_start = facilityCoordinates_allCategories.at(category_id_start).at(facility_id_start);
            x_start = facilities_start.at(0);
            y_start = facilities_start.at(1);  
          }
          int t_end = dayData.at(end_index_agent_id).at(0);
          int t_start = dayData.at(start_index_activity).at(0);
          int t = next_timestep_in_seconds; 
          if (t_end <= t_start) { // agent is changing location during the night
            t_end += 24*60*60;
          }
          if (t < t_start) { // agent is changing location during the night
            t += 24*60*60;
          }

          x = x_start + ((x_end - x_start) * (t - t_start)) / (t_end - t_start); 
          y = y_start + ((y_end - y_start) * (t - t_start)) / (t_end - t_start);

          trajectories[skipping_agent_id][1] = std::vector<double>{x, y};

          int state_nr_start = facilityLabels_allCategories_nr_day.at(category_id_start).at(facility_id_start); // in which state did the agent start his/her journey?
          int state_nr_end = facilityLabels_allCategories_nr_day.at(category_id_end).at(facility_id_end); 
          state_nr = in_which_state(state_nr_start, state_nr_end, t_start, t_end, t);
          (*trajectories_nr_state)[skipping_agent_id][1] = state_nr; //number of state in Germany

          int prev_state_nr = trajectories_nr_state_previous.at(skipping_agent_id).at(timestep_check);
          if (prev_state_nr > -1 && prev_state_nr != state_nr) {
            check_bounds(prev_state_nr, state_nr, corrected_step);
            jump_matrix[corrected_step][prev_state_nr][state_nr] += 1;
            jump_matrix_agent_IDs[corrected_step][prev_state_nr][state_nr].push_back(skipping_agent_id);
          }
        }
      }
      added_extra_event = true;
      if (idx_initial_agents < nbr_trajectories - 1) {
        idx_initial_agents++; // idx_initial_agents < nbr_trajectories -> ok
      }
      else {
        break; // idx_initial_agents++ -> idx_initial_agents >= nbr_trajectories -> not ok! break!
      }
    }

    if (agent_id < idx_initial_agents) { // -> next agent!
      continue;
    }

    int time = row.at(0);
    int activity_type = row.at(1);
    int facility_id = row.at(3);
    int category_id = row.at(4);

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
      x = facilityCoordinates_allCategories[category_id][facility_id][0]; 
      y = facilityCoordinates_allCategories[category_id][facility_id][1]; 
      state_nr = facilityLabels_allCategories_nr_day.at(category_id).at(facility_id);

      trajectories[agent_id][0] = std::vector<double>{x, y};
      (*trajectories_nr_state)[agent_id][0] = state_nr; //number of state in Germany 

      trajectories[agent_id][1] = std::vector<double>{x, y}; 
      (*trajectories_nr_state)[agent_id][1] = state_nr; //number of state in Germany 
      continue;
    }

    if ((last_event || ((previous_person_id != -1) && (previous_person_id != agent_id)))
        && !added_extra_event) { // new person but no final data was added to previous agent
      added_extra_event = true; 

      // add final event data for previous agent! back then, there was no event after next_timestep_in_seconds given, i.e., we need to get the first event of the day and compute the location of the agent at next_timestep_in_seconds 
      int end_index_agent_id = map_agent_id_idx_first_event_data.at(previous_person_id); // first event of the day of previous_person_id
      int start_index_activity = i-1; // last event of previous_person_id

      int facility_id_start = dayData.at(start_index_activity).at(3);
      int category_id_start = dayData.at(start_index_activity).at(4);

      int facility_id_end = dayData.at(end_index_agent_id).at(3);
      int category_id_end = dayData.at(end_index_agent_id).at(4);

      // commuting or not?
      if (dayData.at(start_index_activity).at(1) == 0) { // last event is activity start -> stays at same location
        x = facilityCoordinates_allCategories[category_id_start][facility_id_start][0]; 
        y = facilityCoordinates_allCategories[category_id_start][facility_id_start][1]; 
        state_nr = facilityLabels_allCategories_nr_day.at(category_id_start).at(facility_id_start);
        trajectories[previous_person_id][1] = std::vector<double>{x, y}; 
        (*trajectories_nr_state)[previous_person_id][1] = state_nr; //number of state in Germany 

        int prev_state_nr = trajectories_nr_state_previous.at(previous_person_id).at(timestep_check);
        if (prev_state_nr > -1 && prev_state_nr != state_nr) { 
          check_bounds(prev_state_nr, state_nr, corrected_step);
          jump_matrix[corrected_step][prev_state_nr][state_nr] += 1;
          jump_matrix_agent_IDs[corrected_step][prev_state_nr][state_nr].push_back(previous_person_id);
        }
      }
      else if (facility_id_start == facility_id_end && category_id_start == category_id_end) { // same facility and category -> same location during the night        
        x = facilityCoordinates_allCategories[category_id_end][facility_id_end][0]; 
        y = facilityCoordinates_allCategories[category_id_end][facility_id_end][1]; 
        state_nr = facilityLabels_allCategories_nr_day.at(category_id_end).at(facility_id_end);
        trajectories[previous_person_id][1] = std::vector<double>{x, y}; 
        (*trajectories_nr_state)[previous_person_id][1] = state_nr; //number of state in Germany

        int prev_state_nr = trajectories_nr_state_previous.at(previous_person_id).at(timestep_check);
        if (prev_state_nr > -1 && prev_state_nr != state_nr) {
          check_bounds(prev_state_nr, state_nr, corrected_step);
          jump_matrix[corrected_step][prev_state_nr][state_nr] += 1;
          jump_matrix_agent_IDs[corrected_step][prev_state_nr][state_nr].push_back(previous_person_id);
        }
      }
      else { // commuting: compute location of agent at next_timestep_in_seconds by interpolating between the two activities at index end_index_agent_id and start_index_activity
        const auto& facilities_end = facilityCoordinates_allCategories.at(category_id_end).at(facility_id_end);
        double x_end = facilities_end.at(0);
        double y_end = facilities_end.at(1);

        double x_start, y_start;
        if (category_id_start == -1) { 
          x_start = trajectories[previous_person_id][0][0];
          y_start = trajectories[previous_person_id][0][1];
        }
        else {
          const auto& facilities_start = facilityCoordinates_allCategories.at(category_id_start).at(facility_id_start);
          x_start = facilities_start.at(0);
          y_start = facilities_start.at(1);  
        }
        int t_end = dayData.at(end_index_agent_id).at(0);
        int t_start = dayData.at(start_index_activity).at(0);
        int t = next_timestep_in_seconds; 
        if (t_end <= t_start) { // agent is changing location during the night
          t_end += 24*60*60;
        }
        if (t < t_start) { // agent is changing location during the night
          t += 24*60*60;
        }

        x = x_start + ((x_end - x_start) * (t - t_start)) / (t_end - t_start); 
        y = y_start + ((y_end - y_start) * (t - t_start)) / (t_end - t_start);
 
        trajectories[previous_person_id][1] = std::vector<double>{x, y}; 

        int state_nr_start = facilityLabels_allCategories_nr_day.at(category_id_start).at(facility_id_start); // in which state did the agent start his/her journey?
        int state_nr_end = facilityLabels_allCategories_nr_day.at(category_id_end).at(facility_id_end); 
        state_nr = in_which_state(state_nr_start, state_nr_end, t_start, t_end, t);
        (*trajectories_nr_state)[previous_person_id][1] = state_nr; //number of state in Germany

        int prev_state_nr = trajectories_nr_state_previous.at(previous_person_id).at(timestep_check);
        if (prev_state_nr > -1 && prev_state_nr != state_nr) {
          check_bounds(prev_state_nr, state_nr, corrected_step);
          jump_matrix[corrected_step][prev_state_nr][state_nr] += 1;
          jump_matrix_agent_IDs[corrected_step][prev_state_nr][state_nr].push_back(previous_person_id);
        }
      }
    }

    if (event_of_last_agent || last_event) {
      break;
    }

    if (previous_person_id != agent_id) { // new person -> add initial event data 
      added_extra_event = false;

      if (time == timestep_in_seconds) { // activity is happening exactly at timestep_in_seconds: add event
        x = facilityCoordinates_allCategories[category_id][facility_id][0]; 
        y = facilityCoordinates_allCategories[category_id][facility_id][1]; 
        state_nr = facilityLabels_allCategories_nr_day.at(category_id).at(facility_id);
        trajectories[agent_id][0] = std::vector<double>{x, y}; 
        (*trajectories_nr_state)[agent_id][0] = state_nr; //number of state in Germany
      }
      else { // (time != timestep_in_seconds)
          if ((time > timestep_in_seconds) && (activity_type == 1)) { // activity end -> stays at same location
            x = facilityCoordinates_allCategories[category_id][facility_id][0]; 
            y = facilityCoordinates_allCategories[category_id][facility_id][1]; 
            state_nr = facilityLabels_allCategories_nr_day.at(category_id).at(facility_id);
            trajectories[agent_id][0] = std::vector<double>{x, y}; 
            (*trajectories_nr_state)[agent_id][0] = state_nr; //number of state in Germany
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
              int facility_idx_initial = dayData.at(initial_activity_index).at(3); //facility index
              int category_idx_initial = dayData.at(initial_activity_index).at(4); //category index
              x = facilityCoordinates_allCategories[category_idx_initial][facility_idx_initial][0]; 
              y = facilityCoordinates_allCategories[category_idx_initial][facility_idx_initial][1]; 
              state_nr = facilityLabels_allCategories_nr_day.at(category_idx_initial).at(facility_idx_initial);
              trajectories[agent_id][0] = std::vector<double>{x, y}; 
              (*trajectories_nr_state)[agent_id][0] = state_nr; //number of state in Germany
            }
            else if (dayData.at(initial_activity_index).at(1) == 0 && dayData.at(initial_activity_index).at(3) == dayData.at(next_index_activity).at(3) && dayData.at(initial_activity_index).at(4) == dayData.at(next_index_activity).at(4)) { // activity start and same location -> stays at same location
              int facility_idx_initial = dayData.at(initial_activity_index).at(3); //facility index
              int category_idx_initial = dayData.at(initial_activity_index).at(4); //category index
              x = facilityCoordinates_allCategories[category_idx_initial][facility_idx_initial][0]; 
              y = facilityCoordinates_allCategories[category_idx_initial][facility_idx_initial][1]; 
              state_nr = facilityLabels_allCategories_nr_day.at(category_idx_initial).at(facility_idx_initial);
              trajectories[agent_id][0] = std::vector<double>{x, y}; 
              (*trajectories_nr_state)[agent_id][0] = state_nr; //number of state in Germany
            }
            else {
              if (next_index_activity == initial_activity_index) { // there is only one event -> initial_activity_index == i
                x = facilityCoordinates_allCategories[category_id][facility_id][0]; 
                y = facilityCoordinates_allCategories[category_id][facility_id][1]; 
                state_nr = facilityLabels_allCategories_nr_day.at(category_id).at(facility_id);
                trajectories[agent_id][0] = std::vector<double>{x, y};
                (*trajectories_nr_state)[agent_id][0] = state_nr; //number of state in Germany  
              }
              else {
                // compute location of agent at timestep_in_seconds by interpolating between the two activities at index initial_activity_index and next_index_activity
                const auto& row_start = dayData.at(initial_activity_index);
                int facility_id_start = row_start.at(3);
                int category_id_start = row_start.at(4);
                const auto& facilities_start = facilityCoordinates_allCategories.at(category_id_start).at(facility_id_start);
                double x_start = facilities_start.at(0);
                double y_start = facilities_start.at(1);
                int t_start = row_start.at(0);
          
                const auto& row_end = dayData.at(next_index_activity);
                int facility_id_end = row_end.at(3);
                int category_id_end = row_end.at(4);
                const auto& facilities_end = facilityCoordinates_allCategories.at(category_id_end).at(facility_id_end);
                double x_end = facilities_end.at(0);
                double y_end = facilities_end.at(1);
                int t_end = row_end.at(0);

                int t = timestep_in_seconds; 
                if (t_end <= t_start) { // agent is changing location during the night
                  t_end += 24*60*60;
                }
                if (t < t_start) { // agent is changing location during the night
                  t += 24*60*60;
                }
                x = x_start + ((x_end - x_start) * (t - t_start)) / (t_end - t_start); 
                y = y_start + ((y_end - y_start) * (t - t_start)) / (t_end - t_start);

                trajectories[agent_id][0] = std::vector<double>{x, y};

                int state_nr_start = facilityLabels_allCategories_nr_day.at(category_id_start).at(facility_id_start); // in which state did the agent start his/her journey?
                int state_nr_end = facilityLabels_allCategories_nr_day.at(category_id_end).at(facility_id_end); 
                state_nr = in_which_state(state_nr_start, state_nr_end, t_start, t_end, t);
                (*trajectories_nr_state)[agent_id][0] = state_nr; //number of state in Germany  
              }
            }      
          }
      }
      previous_person_id = agent_id;
    }

    if ((time == next_timestep_in_seconds) && (!added_extra_event)) { // no event added yet for this agent besides the artificial event at time timestep_in_seconds and next will be excatly at next_timestep_in_seconds, so not artificial!
      added_extra_event = true;

      x = facilityCoordinates_allCategories[category_id][facility_id][0]; 
      y = facilityCoordinates_allCategories[category_id][facility_id][1]; 
      state_nr = facilityLabels_allCategories_nr_day.at(category_id).at(facility_id);
      trajectories[agent_id][1] = std::vector<double>{x, y}; 
      (*trajectories_nr_state)[agent_id][1] = state_nr; //number of state in Germany  

      int prev_state_nr = trajectories_nr_state_previous.at(agent_id).at(timestep_check);
      if (prev_state_nr > -1 && prev_state_nr != state_nr) { //(*trajectories_nr_state).at(agent_id).at(0) 
        check_bounds(prev_state_nr, state_nr, corrected_step);
        jump_matrix[corrected_step][prev_state_nr][state_nr] += 1;
        jump_matrix_agent_IDs[corrected_step][prev_state_nr][state_nr].push_back(agent_id);
      }
    }
    else if ((time > next_timestep_in_seconds) && (!added_extra_event)) { // no event added yet for this agent besides the (artificial) event at time timestep_in_seconds (and next will be after next_timestep_in_seconds)
      added_extra_event = true;
      if (activity_type == 1) { // activity end after next_timestep_in_seconds: agent stays at same location
        x = facilityCoordinates_allCategories[category_id][facility_id][0]; 
        y = facilityCoordinates_allCategories[category_id][facility_id][1]; 
        state_nr = facilityLabels_allCategories_nr_day.at(category_id).at(facility_id);

        trajectories[agent_id][1] = std::vector<double>{x, y}; 
        (*trajectories_nr_state)[agent_id][1] = state_nr; //number of state in Germany  

        int prev_state_nr = trajectories_nr_state_previous.at(agent_id).at(timestep_check);
        if (prev_state_nr > -1 && prev_state_nr != state_nr) {
          check_bounds(prev_state_nr, state_nr, corrected_step);
          jump_matrix[corrected_step][prev_state_nr][state_nr] += 1;
          jump_matrix_agent_IDs[corrected_step][prev_state_nr][state_nr].push_back(agent_id);
        }
      }
      else { 
        int facility_id_end = facility_id;
        int category_id_end = category_id;
        const auto& facilities_end = facilityCoordinates_allCategories.at(category_id_end).at(facility_id_end);
        double x_end = facilities_end.at(0);
        double y_end = facilities_end.at(1);  

        int idx_previous_event;
        if (dayData.at(i-1).at(2) == agent_id) { // previous event is from the same agent
          idx_previous_event = i-1;
        }
        else { //no previous event 
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

        double x_start, y_start;
        const auto& facilities_start = facilityCoordinates_allCategories.at(previous_event_category_id).at(previous_event_facility_id);
        x_start = facilities_start.at(0);
        y_start = facilities_start.at(1);
        int t_end = time;
        int t_start = dayData.at(idx_previous_event).at(0); 
        int t = next_timestep_in_seconds;
        if (t_end <= t_start) { // agent is changing location during the night
          t_end += 24*60*60;
        }
        if (t < t_start) { // agent is changing location during the night
          t += 24*60*60;
        }

        x = x_start + ((x_end - x_start) * (t - t_start)) / (t_end - t_start); 
        y = y_start + ((y_end - y_start) * (t - t_start)) / (t_end - t_start);

        trajectories[agent_id][1] = std::vector<double>{x, y}; 

        int state_nr_start = facilityLabels_allCategories_nr_day.at(previous_event_category_id).at(previous_event_facility_id); // in which state did the agent start his/her journey?
        int state_nr_end = facilityLabels_allCategories_nr_day.at(category_id_end).at(facility_id_end); 
        state_nr = in_which_state(state_nr_start, state_nr_end, t_start, t_end, t);
        (*trajectories_nr_state)[agent_id][1] = state_nr; //number of state in Germany 

        int prev_state_nr = trajectories_nr_state_previous.at(agent_id).at(timestep_check); 
        if (prev_state_nr > -1 && prev_state_nr != state_nr) {
          check_bounds(prev_state_nr, state_nr, corrected_step);
          jump_matrix[corrected_step][prev_state_nr][state_nr] += 1;
          jump_matrix_agent_IDs[corrected_step][prev_state_nr][state_nr].push_back(agent_id);
        }
      }
    }
  }
  return std::make_tuple(trajectories, *trajectories_nr_state); 
}


int get_nr_day(int week_day_) {
  int nr_day;
  if (week_day_ % 7 == 5) { // saturday
    nr_day = 1;
  }
  else if (week_day_ % 7 == 6) { // sunday
    nr_day = 2;
  }
  else {
    nr_day = 0;
  }
  return nr_day;
}

int main(int argc, char *argv[]) {
  int week_day = 6; //5; //0; // 
  // 0: monday, 5: saturday, 6: sunday

  int nr_day = get_nr_day(week_day);
  initial_nr_day = nr_day;
  std::cout << "nr_day " << nr_day << std::endl;

  bool with_restriction = true; 

  double population_scale = 4.0;
  int dt_inv = 48; // 30 min

  ABModel abm(population_scale, with_restriction, dt_inv); 
  abm.readEventData(nr_day); 
  abm.readFacilityCoordinates_allCategories();
  abm.readFacilityCoordinates_home();

  double dt = 1.0/dt_inv;
  double T = 1.0;
  double t_0 = 0.0;
  int maxSteps = (T-t_0)*dt_inv; 

  int nbr_trajectories = abm.nbr_trajectories;
  std::cout << "total number of Trajectories: " << (nbr_trajectories) << std::endl;

  std::vector<std::vector<std::vector<int>>> eventData = abm.eventData; //num Events, cols: time in seconds, activity type, person id, facility index (corresponding to all facilities), category index (corresponding to facility_categories)
  std::vector<std::vector<std::vector<std::vector<double>>>> facilityCoordinates_allCategories = abm.facilityCoordinates_allCategories;
  std::vector<std::unordered_map<int, int>> facilityCoordinates_home = abm.facilityCoordinates_home;
  std::vector<std::vector<std::vector<int>>> facilityLabels_allCategories = abm.facilityLabels_allCategories;

  int steps_initial = 0;

  std::vector<std::vector<std::tuple<double, double, int>>> trajectories(nbr_trajectories, std::vector<std::tuple<double, double, int>>(maxSteps + 1, std::make_tuple(0.0, 0.0, 0)));

  std::vector<std::vector<std::vector<int>>> jump_matrix(dt_inv, std::vector<std::vector<int>>(nbr_domains, std::vector<int>(nbr_domains,0)));
  std::vector<std::vector<std::vector<std::vector<int>>>> jump_matrix_agent_IDs(dt_inv, std::vector<std::vector<std::vector<int>>>(nbr_domains, std::vector<std::vector<int>>(nbr_domains)));
  std::vector<std::vector<std::vector<double>>> trajectory_partial;
  std::vector<std::vector<int>> trajectory_partial_nr_state;
  std::tie(trajectory_partial, trajectory_partial_nr_state) = create_trajectories(nbr_trajectories, trajectory_partial_nr_state, eventData, facilityCoordinates_allCategories, facilityCoordinates_home, facilityLabels_allCategories,
    dt_inv, steps_initial, nr_day, jump_matrix, jump_matrix_agent_IDs);

  for (int agent_idx=0; agent_idx<nbr_trajectories; agent_idx++) {
    std::get<0>(trajectories[agent_idx][0]) = trajectory_partial[agent_idx][0][0]; // x
    std::get<1>(trajectories[agent_idx][0]) = trajectory_partial[agent_idx][0][1]; // y
    std::get<2>(trajectories[agent_idx][0]) = trajectory_partial_nr_state[agent_idx][0]; // state

    std::get<0>(trajectories[agent_idx][1]) = trajectory_partial[agent_idx][1][0]; // x
    std::get<1>(trajectories[agent_idx][1]) = trajectory_partial[agent_idx][1][1]; // y
    std::get<2>(trajectories[agent_idx][1]) = trajectory_partial_nr_state[agent_idx][1]; // state
  }

  for (int steps=1; steps<maxSteps; ++steps) {
    std::cout << "step " << steps << std::endl;
    std::tie(trajectory_partial, trajectory_partial_nr_state) = create_trajectories(nbr_trajectories, trajectory_partial_nr_state, eventData, facilityCoordinates_allCategories, facilityCoordinates_home, facilityLabels_allCategories,
      dt_inv, steps, nr_day, jump_matrix, jump_matrix_agent_IDs);

    for (int agent_idx=0; agent_idx<nbr_trajectories; agent_idx++) {
      std::get<0>(trajectories[agent_idx][steps+1]) = trajectory_partial[agent_idx][1][0]; //x
      std::get<1>(trajectories[agent_idx][steps+1]) = trajectory_partial[agent_idx][1][1]; //y
      std::get<2>(trajectories[agent_idx][steps+1]) = trajectory_partial_nr_state[agent_idx][1]; //state number
    }
  }

  std::vector<std::string> days_str = {"wt", "sa", "so"};
  std::string filename = "../../work/output/germany_trajectories_" + days_str[nr_day] + "_dt_inv_" + std::to_string(dt_inv) + ".bin.gz";
  gzFile outfile = gzopen(filename.c_str(), "wb");
  for (const auto& agent : trajectories) {
    for (const auto& entry : agent) {
      float x = static_cast<float>(std::get<0>(entry));
      float y = static_cast<float>(std::get<1>(entry));
      uint8_t state = static_cast<uint8_t>(std::get<2>(entry));
      gzwrite(outfile, &x, sizeof(x));
      gzwrite(outfile, &y, sizeof(y));
      gzwrite(outfile, &state, sizeof(state));
    }
  }
  gzclose(outfile);

  for (int j=0; j<trajectories[0].size(); j++) {
    std::cout << std::get<0>(trajectories[trajectories.size()-2][j]) << ", " << std::get<1>(trajectories[trajectories.size()-2][j]) << std::endl;
  }
  std::cout << "jump_matrix[0][0][0] " << jump_matrix[0][0][0] << std::endl;
  std::cout << "jump_matrix[jump_matrix.size()-1][nbr_domains-1][nbr_domains-1] " << jump_matrix[jump_matrix.size()-1][nbr_domains-1][nbr_domains-1] << std::endl;

  std::string filename_jump_agent_IDs = "../../work/input/global/germany_jump_matrix_agent_IDs_" + days_str[nr_day] + addOn + "_dt_inv_" + std::to_string(dt_inv) + ".bin";
  std::ofstream file(filename_jump_agent_IDs, std::ios::binary);
  for (const auto& mat2d : jump_matrix_agent_IDs) {          // time step
    for (const auto& row2d : mat2d) {                        // from-state
      for (const auto& cell : row2d) {                       // to-state
        int len = static_cast<int>(cell.size());
        file.write(reinterpret_cast<const char*>(&len), sizeof(int));
        if (len > 0) {
            file.write(reinterpret_cast<const char*>(cell.data()), len * sizeof(int));
        }
      }
    }
  }
  file.close();

  std::cout << "done" << std::endl;
}
#endif