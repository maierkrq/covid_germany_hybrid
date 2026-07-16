#ifndef INTEGRATE_HH
#define INTEGRATE_HH

#include <iostream>
#include <string> //
#include <cstdlib> // std::exit(EXIT_FAILURE)

#include <chrono>
#include <thread>
#include <filesystem>

#include "io/vtkwriter.hh"
#include "io/vtkwriter.cpp"

#if HAVE_LIBAMIRA == 1
#include "io/amira.hh"
#endif

#include <unordered_set>

#include "covid.hh"
#include "ABM.cpp"

template<int nbrCompartment, typename T>
std::vector<std::vector<double>> getNumberOfJumps(T& x,
                        std::vector<std::vector<double>>& compartmentJumpingPersons, int state_PDE_idx) {
  auto& coeff = boost::fusion::at_c<nbrCompartment>(x.data).coefficients();
  for (size_t i = 0; i < coeff.N(); i++) {
    if (coeff[i][0] < 0.0) {
      std::cout << "negativ VOR dem Sprung" << std::endl;
      std::cout << "Compartment " << nbrCompartment << std::endl;
      std::cout << "Gitterpunkt " << i << std::endl;
      std::cout << "coeff[i][0] " << coeff[i][0] << std::endl;
      std::cout << "state_PDE_idx " << state_PDE_idx << std::endl;
      std::exit(EXIT_FAILURE);
    }
    double nbr_jumping_persons = coeff[i][0] / minDensityForTransition_PDE.at(state_PDE_idx)[i];
    compartmentJumpingPersons[i][nbrCompartment] = nbr_jumping_persons;
  }
  return compartmentJumpingPersons;
}

template<int N, typename T> 
T processCompartment(T& x, std::vector<std::vector<double>>& compartmentJumpingPersons, int& sign, int state_PDE_idx) {
  for (size_t i = 0; i < boost::fusion::at_c<N>(x.data).coefficients().N(); i++) {
    double nbr_jumping_persons = compartmentJumpingPersons[i][N];
    if ((sign == -1) && (boost::fusion::at_c<N>(x.data).coefficients()[i][0]>0) && (boost::fusion::at_c<N>(x.data).coefficients()[i][0] + sign * nbr_jumping_persons * minDensityForTransition_PDE.at(state_PDE_idx)[i] < 0)) {
      std::cout << "wird negativ nach Sprung" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    boost::fusion::at_c<N>(x.data).coefficients()[i][0] += sign * nbr_jumping_persons * minDensityForTransition_PDE.at(state_PDE_idx)[i];
  }
  return x; 
}

template<typename T>
void processCompartments(T& x, std::vector<std::vector<double>>& compartmentJumpingPersons, int sign, int state_PDE_idx, int max_jumps=0, bool targetIsNonABM=true) {
  if (sign == 1) { // (from ABM/PDE/ODE) to PDE
    //count number of persons jumping per compartment (compartmentJumpingPersons)
    std::vector<double> number_of_jumps_compartment(compartmentJumpingPersons[0].size(), 0.0);  
    for (int j=0; j<compartmentJumpingPersons[0].size(); j++) { //compartments 
        number_of_jumps_compartment[j] += compartmentJumpingPersons[0][j];
        max_jumps += (compartmentJumpingPersons[0][j]+0.1); // round to whole numbers
    }

    for (auto& row : compartmentJumpingPersons) {
      std::fill(row.begin(), row.end(), 0.0);
    }

    // get max number of jumps for all compartments
    compartmentJumpingPersons = getNumberOfJumps<0>(x, compartmentJumpingPersons, state_PDE_idx); 
    compartmentJumpingPersons = getNumberOfJumps<1>(x, compartmentJumpingPersons, state_PDE_idx);
    compartmentJumpingPersons = getNumberOfJumps<2>(x, compartmentJumpingPersons, state_PDE_idx);
    compartmentJumpingPersons = getNumberOfJumps<3>(x, compartmentJumpingPersons, state_PDE_idx);
    compartmentJumpingPersons = getNumberOfJumps<4>(x, compartmentJumpingPersons, state_PDE_idx);
    compartmentJumpingPersons = getNumberOfJumps<5>(x, compartmentJumpingPersons, state_PDE_idx);
    compartmentJumpingPersons = getNumberOfJumps<6>(x, compartmentJumpingPersons, state_PDE_idx);
    compartmentJumpingPersons = getNumberOfJumps<7>(x, compartmentJumpingPersons, state_PDE_idx);

    std::vector<double> sum_per_compartment(compartmentJumpingPersons[0].size(), 0.0);
    for (int j=0; j<compartmentJumpingPersons[0].size(); j++) { //compartments
      for (int i=0; i<compartmentJumpingPersons.size(); i++) { // grid nodes
        sum_per_compartment[j] += compartmentJumpingPersons[i][j];
      }
    }

    double sum_non_susceptible_jumps = 0.0;
    for (int j=1; j<compartmentJumpingPersons[0].size(); j++) { //compartments
      sum_non_susceptible_jumps += number_of_jumps_compartment[j]; 
      if (sum_per_compartment[j] > 0) {
        double ratio_jumping = number_of_jumps_compartment[j] * (1.0/sum_per_compartment[j]);
        for (int i=0; i<compartmentJumpingPersons.size(); i++) { // grid nodes
          compartmentJumpingPersons[i][j] *= ratio_jumping;
        }
      }
      else { // if compartment is still empty, add contribution directly
        for (int i=0; i<compartmentJumpingPersons.size(); i++) { // grid nodes
          compartmentJumpingPersons[i][j] = number_of_jumps_compartment[j] * normalized_inv_V_PDE.at(state_PDE_idx).at(i) / minDensityForTransition_PDE.at(state_PDE_idx)[i]; // add according to landscape like initial values
        }
      }
    }
    // susceptible jump
    if (sum_per_compartment[0] > 0) {
      for (int i=0; i<compartmentJumpingPersons.size(); i++) { // grid nodes
        compartmentJumpingPersons[i][0] *= number_of_jumps_compartment[0] * (1.0/sum_per_compartment[0]);
      }
    }
    else {
      for (int i=0; i<compartmentJumpingPersons.size(); i++) { // grid nodes
        compartmentJumpingPersons[i][0] = number_of_jumps_compartment[0] * normalized_inv_V_PDE.at(state_PDE_idx).at(i) / minDensityForTransition_PDE.at(state_PDE_idx)[i]; // add according to landscape like initial values
      }
    }
  }
  else if (sign == -1) { // from PDE (to ABM/PDE/ODE)
    // set all values of compartmentJumpingPersons to zero
    for (auto& row : compartmentJumpingPersons) {
      std::fill(row.begin(), row.end(), 0.0);
    }

    // get max number of jumps for all compartments
    compartmentJumpingPersons = getNumberOfJumps<0>(x, compartmentJumpingPersons, state_PDE_idx); 
    compartmentJumpingPersons = getNumberOfJumps<1>(x, compartmentJumpingPersons, state_PDE_idx);
    compartmentJumpingPersons = getNumberOfJumps<2>(x, compartmentJumpingPersons, state_PDE_idx);
    compartmentJumpingPersons = getNumberOfJumps<3>(x, compartmentJumpingPersons, state_PDE_idx);
    compartmentJumpingPersons = getNumberOfJumps<4>(x, compartmentJumpingPersons, state_PDE_idx);
    compartmentJumpingPersons = getNumberOfJumps<5>(x, compartmentJumpingPersons, state_PDE_idx);
    compartmentJumpingPersons = getNumberOfJumps<6>(x, compartmentJumpingPersons, state_PDE_idx);
    compartmentJumpingPersons = getNumberOfJumps<7>(x, compartmentJumpingPersons, state_PDE_idx);

    // get probabilities of jumping for each compartment
    double cnt_persons_able_to_jump = 0.0;
    std::vector<double> sum_per_compartment(compartmentJumpingPersons[0].size(), 0.0);
    for (int i=0; i<compartmentJumpingPersons.size(); i++) { // grid nodes
      for (int j=0; j<compartmentJumpingPersons[0].size(); j++) { //compartments
        cnt_persons_able_to_jump += compartmentJumpingPersons[i][j];
        sum_per_compartment[j] += compartmentJumpingPersons[i][j];
      }
    }

    double inverse_cnt_persons_able_to_jump = 1.0 / cnt_persons_able_to_jump;
    double sum_non_susceptible_jumps = 0.0;
    for (int j=1; j<compartmentJumpingPersons[0].size(); j++) { //compartments
      if (sum_per_compartment[j] > 0) {
        double percentage_of_compartment = sum_per_compartment[j] * inverse_cnt_persons_able_to_jump;
        double number_of_jumps_compartment = percentage_of_compartment * max_jumps;
        if (!targetIsNonABM) { // compartmentJumpingPersons needs to work as an input for jumps into ABM domain -> total jumps per compartment!
          number_of_jumps_compartment = int(number_of_jumps_compartment);
        }
        double ratio_jumping = number_of_jumps_compartment * (1.0/sum_per_compartment[j]);
        sum_non_susceptible_jumps += number_of_jumps_compartment; // sum to later account for missing jumps!

        for (int i=0; i<compartmentJumpingPersons.size(); i++) { // grid nodes
          compartmentJumpingPersons[i][j] *= ratio_jumping;
        }
      }
    }

    for (int i=0; i<compartmentJumpingPersons.size(); i++) { // grid nodes
      compartmentJumpingPersons[i][0] *= (max_jumps - sum_non_susceptible_jumps) * (1.0/sum_per_compartment[0]); // account for missing jumps!
    }
  }

  // process jumping persons for each compartment 
  x = processCompartment<0>(x, compartmentJumpingPersons, sign, state_PDE_idx); 
  x = processCompartment<1>(x, compartmentJumpingPersons, sign, state_PDE_idx);
  x = processCompartment<2>(x, compartmentJumpingPersons, sign, state_PDE_idx);
  x = processCompartment<3>(x, compartmentJumpingPersons, sign, state_PDE_idx);
  x = processCompartment<4>(x, compartmentJumpingPersons, sign, state_PDE_idx);
  x = processCompartment<5>(x, compartmentJumpingPersons, sign, state_PDE_idx);
  x = processCompartment<6>(x, compartmentJumpingPersons, sign, state_PDE_idx);
  x = processCompartment<7>(x, compartmentJumpingPersons, sign, state_PDE_idx);
}

void processODECompartments(std::vector<double>& x, std::vector<double>& compartmentJumpingPersons, int& sign, int max_jumps=0, bool targetIsNonABM=true) {
  if (sign == 1) { // (from ABM/PDE/ODE) to ODE
    for (int i=0; i<x.size(); i++) {
      x[i] += compartmentJumpingPersons[i];
    }
  }
  else if (sign == -1) { // from ODE to (ABM/PDE/ODE)
    // get probabilities of jumping for each compartment
    double cnt_persons_able_to_jump = 0.0;
    for (int i=0; i<nbr_compartments; i++) { 
      compartmentJumpingPersons[i] = x[i]; //overwrite
      cnt_persons_able_to_jump += x[i];
    }

    double inverse_cnt_persons_able_to_jump = 1.0 / cnt_persons_able_to_jump;
    double sum_non_susceptible_jumps = 0.0;
    for (int j=1; j<nbr_compartments; j++) { // all but susceptible compartment
      if (compartmentJumpingPersons[j] > 0) {
        double percentage_of_compartment = compartmentJumpingPersons[j] * inverse_cnt_persons_able_to_jump;
        double number_of_jumps_compartment = percentage_of_compartment * max_jumps; 
        if (!targetIsNonABM) { // compartmentJumpingPersons needs to work as an input for jumps into ABM domain -> total jumps per compartment!
          number_of_jumps_compartment = int(number_of_jumps_compartment);
        }
        sum_non_susceptible_jumps += number_of_jumps_compartment; // sum to later account for missing jumps! because of rounding to whole jumps!
        compartmentJumpingPersons[j] *= number_of_jumps_compartment * (1.0/compartmentJumpingPersons[j]);
      }
    }
    compartmentJumpingPersons[0] *= (max_jumps - sum_non_susceptible_jumps) * (1.0/compartmentJumpingPersons[0]); // account for missing jumps! because of rounding to whole jumps!

    for (int i=0; i<x.size(); i++) {
      x[i] -= compartmentJumpingPersons[i];
    }
  }
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
    nr_day = 0; // monday-friday
  }
  return nr_day;
}

template <class Grid, class Equation, class VariableSet>
std::tuple< std::vector<std::unordered_map<int, int>>, std::vector<int> >
integrate(
                                                std::vector<std::unique_ptr<GridManager<Grid>>>& gridManagers_PDE,
                                                std::vector<Equation>& eq, std::vector<VariableSet>const& variableSet, 
                                                int dt_inv, double population_scale,
                                                std::ranlux24_base& gen_local, //& important here
                                                int correction_step,
                                                int week_day, double T, int maxSteps, int extrapolOrder, 
                                                std::vector<typename VariableSet::VariableSet>& x_PDE,
                                                std::vector<std::vector<double>>& x_ODE,
                                                ABModel& abm, std::vector<double> param_ABMs, 
                                                std::vector<double> params_ODEs,
                                                int run, int initial_run, 
                                                std::vector<std::vector<bool>> with_coupling_states,
                                                DirectType directType,
                                                int verbosity)
{ 
  double dt = 1.0/dt_inv;
  int nr_day = get_nr_day(week_day);

  std::vector<std::vector<std::pair<double,double>>> tolX_PDE;
  std::vector<Limex<Equation>> limex_PDE;
  for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) { // iterate over all PDE domains
    tolX_PDE.push_back(std::vector<std::pair<double,double>>(variableSet.at(state_PDE_idx).noOfVariables));
    limex_PDE.push_back(Limex<Equation>(*gridManagers_PDE.at(state_PDE_idx),eq.at(state_PDE_idx),variableSet.at(state_PDE_idx),directType,PrecondType::ILUK,verbosity));
  }

  bool done = false;
  double t_0 = eq[0].time();
  int init_time_step = (int) t_0;
  std::cout << "init_time_step: " << init_time_step << std::endl;

  std::vector<int> nbr_points_PDE(nr_PDE_states);
  for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) { // iterate over all PDE domains
    Kaskade::VTKWriterDetail::VTKGridInfo<typename VariableSet::VariableSet::Descriptions::GridView> gridInfo(x_PDE.at(state_PDE_idx).descriptions.gridView, IoOptions().setPrecision(10).setDataMode(IoOptions::conforming));
    nbr_points_PDE.at(state_PDE_idx) = gridInfo.npoints;
    std::cout << "nbr_points_PDE.at(state_PDE_idx) " << nbr_points_PDE.at(state_PDE_idx) << std::endl;
  }

  std::vector<std::vector<double>> data_F(nr_PDE_states);
  std::vector<std::vector<std::vector<double>>> data_DF(nr_PDE_states);
  for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) { // iterate over all PDE domains
    data_F.at(state_PDE_idx).resize(nbr_points_PDE.at(state_PDE_idx) * (maxSteps+1),0.0);
    data_DF.at(state_PDE_idx).resize(nbr_points_PDE.at(state_PDE_idx) * (maxSteps+1),std::vector<double>(1,0.0));
  }

  std::vector<std::vector<double>> F_PDE(nr_PDE_states,std::vector<double>(maxSteps+1, 0.0));
  std::vector<std::vector<double>> F_ODE(nr_ODE_states,std::vector<double>(maxSteps+1, 0.0));

  std::vector<std::vector<double>> n_PDE(nr_PDE_states,std::vector<double>(maxSteps+1, 0.0));
  std::vector<std::vector<double>> n_ODE(nr_ODE_states,std::vector<double>(maxSteps+1, 0.0));

  std::vector<std::string> variables_symptomatic = {"sy"};   
  std::vector<std::vector<double>> ones;
  for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) {
    ones.push_back(std::vector<double>(nbr_points_PDE.at(state_PDE_idx),1.0));
  }

  std::vector<int> nbr_agents = abm.nbr_agents[nr_day];
  int total_nbr_agents = 0;
  for (int i=0; i<nbr_agents.size(); i++) { // iterate over all ABM domains
    total_nbr_agents += nbr_agents[i];
  }
  int nbr_trajectories = abm.nbr_trajectories;
  std::cout << "total number of agents: " << (total_nbr_agents) << ", nr_day " << nr_day << std::endl;
  std::cout << "total number of Trajectories: " << (nbr_trajectories) << std::endl;

  std::vector<std::vector<int>> filtered_event_data; 
  std::unordered_map<int, int> agentID_location_commuting_start, agentID_location_commuting_end;  
  std::vector<std::vector<std::vector<int>>> initial_agent_ids = abm.initial_agent_ids; 
  std::vector<int> agent_IDs_in_ABM;
  std::unordered_set<int> agent_IDs_in_ABM_set;
  for (int state_ABM_idx=0; state_ABM_idx<nr_ABM_states; state_ABM_idx++) {
    for (int i=0; i<initial_agent_ids[nr_day].at(state_ABM_idx).size(); i++) {
      agent_IDs_in_ABM.push_back(initial_agent_ids[nr_day].at(state_ABM_idx)[i]);
      agent_IDs_in_ABM_set.insert(initial_agent_ids[nr_day].at(state_ABM_idx)[i]);
    }
  }
  std::sort(agent_IDs_in_ABM.begin(), agent_IDs_in_ABM.end()); // sort by id

  std::vector<std::unordered_map<int, int>> agentID_healthStatus_map(nr_ABM_states);
  std::unordered_map<int, int> agentID_healthStatus_map_flatten;
  for (int state_ABM_idx=0; state_ABM_idx<nr_ABM_states; state_ABM_idx++) {
    int number_of_exposed_agents = int(20 * exposed[run-initial_run][0][abm_idx].at(state_ABM_idx)); 
    if (number_of_exposed_agents > nbr_agents.at(state_ABM_idx)) {
      number_of_exposed_agents = nbr_agents.at(state_ABM_idx);
    }
    std::cout << "state_ABM_idx " << state_ABM_idx << ", number of exposed agents " << number_of_exposed_agents << std::endl;

    // mische initial_agent_ids[nr_day] randomly und speichere das ergebnis in separatem vektor
    std::vector<int> shuffled_initial_agent_ids = initial_agent_ids[nr_day].at(state_ABM_idx);
    std::shuffle(shuffled_initial_agent_ids.begin(), shuffled_initial_agent_ids.end(), gen_local);
    for (int i = 0; i < number_of_exposed_agents; i++) {
      agentID_healthStatus_map.at(state_ABM_idx)[shuffled_initial_agent_ids.at(i)] = 1;
    }
    for (int i = number_of_exposed_agents; i < nbr_agents.at(state_ABM_idx); i++) {
      agentID_healthStatus_map.at(state_ABM_idx)[shuffled_initial_agent_ids.at(i)] = 0; //health status: almost everybody is susceptible
    }
  }
  
  double area;
  std::vector<std::vector<double>> corners;
  corners.resize(3, std::vector<double>(2, 0.0));

  std::vector<std::vector<int>> nbr_symptomatic_agents(nr_ABM_states, std::vector<int>(maxSteps+1,0));
  if (verbosity > 0) { 
     
    // compute initial number of symptomatic people -> new symptomatic individuals
    for (int state_ABM_idx = 0; state_ABM_idx < nr_ABM_states; state_ABM_idx++) { // iterate over all ABM domains
      for (const auto& agent_ID_healthState: agentID_healthStatus_map.at(state_ABM_idx)) {
        int healthState = agent_ID_healthState.second;
        if (healthState == 3) {
          nbr_symptomatic_agents.at(state_ABM_idx)[0]++;
        }
      }
    }

    // save ODE solution
    for (int state_ODE_idx=0; state_ODE_idx<nr_ODE_states; state_ODE_idx++) { // iterate over all ODE domains
      F_ODE.at(state_ODE_idx)[0] = x_ODE.at(state_ODE_idx)[3]; // total number of symptomatic people in ODE domain
    }

    // save PDE solution
    for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) { // iterate over all PDE domains
      std::ostringstream fn;
      fn << "../../work/output/output_data_optim_" + std::to_string(correction_step) + "/graph_abm_pde_ode/covid-slab-uv-" + stateLabels_PDE.at(state_PDE_idx) + "-" + std::to_string(run) + "-";
      fn.width(3);
      fn.fill('0');
      fn.setf(std::ios_base::right,std::ios_base::adjustfield);
      // vtkFileNo++;
      fn << init_time_step*dt_inv;
      fn.flush(); 
      #pragma omp critical(vtk_write)
      {
        writeVTK(x_PDE.at(state_PDE_idx),fn.str(),IoOptions().setPrecision(10),data_F.at(state_PDE_idx),data_DF.at(state_PDE_idx),variables_symptomatic,0,ones.at(state_PDE_idx),verbosity);
      }
      if (verbosity > 0) {
        std::cout << variableSet.at(state_PDE_idx).degreesOfFreedom() << " values written to " << fn.str() << std::endl;

      // approximate integral of symptomatic over each domain 
        for (auto const& cell: elements(x_PDE.at(state_PDE_idx).descriptions.gridView)) {
          for (int corner=0; corner<cell.subEntities(2); corner++) {
            auto xglob = cell.geometry().corner(corner);
            corners[corner][0] = xglob[0];
            corners[corner][1] = xglob[1];
          }
          area = 0.5 * std::abs(corners[0][0] * (corners[1][1]-corners[2][1]) + corners[1][0] * (corners[2][1]-corners[0][1]) + corners[2][0] * (corners[0][1]-corners[1][1]) );

          for (int corner=0; corner<cell.subEntities(2); corner++) {
            int index = x_PDE.at(state_PDE_idx).descriptions.gridView.indexSet().subIndex(cell, corner, 2); //global index
            F_PDE.at(state_PDE_idx)[0] += (area/3) * data_F.at(state_PDE_idx)[index];
          } // end for corner
        } // end for cell
      } // verbosity
    } // end for PDE domains
  }

  int filtered_event_data_size;

  std::vector<std::vector<std::vector<std::vector<double>>>> PDEcompartmentJumpingIntoABMPersons;
  for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) { // iterate over all PDE domains
    PDEcompartmentJumpingIntoABMPersons.push_back(std::vector<std::vector<std::vector<double>>>(nr_ABM_states,std::vector<std::vector<double>>(nbr_points_PDE.at(state_PDE_idx), std::vector<double>(nbr_compartments,0.0))));
  }
  std::vector<std::vector<std::vector<double>>> ODEcompartmentJumpingIntoABMPersons(nr_ODE_states, std::vector<std::vector<double>>(nr_ABM_states, std::vector<double>(nbr_compartments,0.0)));
  std::vector<std::vector<std::vector<std::vector<double>>>> ODEcompartmentJumpingIntoPDEPersons(nr_ODE_states, std::vector<std::vector<std::vector<double>>>(nr_PDE_states)); 
  for (int ode = 0; ode < nr_ODE_states; ode++) {
    for (int pde = 0; pde < nr_PDE_states; pde++) {
        ODEcompartmentJumpingIntoPDEPersons[ode][pde] =
        std::vector<std::vector<double>>(nbr_points_PDE[pde], std::vector<double>(nbr_compartments, 0.0));
    }
  }
  
  std::vector<std::vector<std::vector<double>>> ODEcompartmentJumpingIntoODEPersons(nr_ODE_states, std::vector<std::vector<double>>(nr_ODE_states, std::vector<double>(nbr_compartments,0.0)));
  std::vector<std::vector<std::vector<std::vector<double>>>> PDEcompartmentJumpingIntoODEPersons(nr_PDE_states, std::vector<std::vector<std::vector<double>>>(nr_ODE_states)); 
  for (int pde = 0; pde < nr_PDE_states; pde++) {
    for (int ode = 0; ode < nr_ODE_states; ode++) {
      PDEcompartmentJumpingIntoODEPersons[pde][ode] =
        std::vector<std::vector<double>>(nbr_points_PDE[pde], std::vector<double>(nbr_compartments, 0.0));
    }
  }

  std::vector<std::vector<std::vector<std::vector<double>>>> PDEcompartmentJumpingIntoPDEPersons(nr_PDE_states, std::vector<std::vector<std::vector<double>>>(nr_PDE_states)); 
  for (int pde = 0; pde < nr_PDE_states; pde++) {
    for (int pde2 = 0; pde2 < nr_PDE_states; pde2++) {
      PDEcompartmentJumpingIntoPDEPersons[pde][pde2] =
        std::vector<std::vector<double>>(nbr_points_PDE[pde], std::vector<double>(nbr_compartments, 0.0));
    }
  }
  int sign;

  std::vector<std::vector<int>> healthStatus_agents_in_PDE, healthStatus_agents_in_ODE;
  const auto& jump_matrix = abm.jump_matrix; 
  const auto& germany_jump_matrix_agent_IDs = abm.germany_jump_matrix_agent_IDs;
  std::vector<std::vector<std::vector<std::vector<int>>>> germany_jump_matrix_agent_IDs_with_reductions;
  germany_jump_matrix_agent_IDs_with_reductions.resize(dt_inv, std::vector<std::vector<std::vector<int>>>(nbr_domains, std::vector<std::vector<int>>(nbr_domains))); 

  std::vector<std::vector<std::vector<int>>> cnt_jumping_persons_per_compartment(nbr_domains, std::vector<std::vector<int>>(nr_ABM_states,std::vector<int>(nbr_compartments)));
  const auto& activityChangePercentage = abm.activityChangePercentage;
  int max_jumps = 0; 
  int nbr_unique_agents = 0;
  for (int steps=0; !done && steps<maxSteps; steps++) {
    std::cout << "next step " << steps+(init_time_step*dt_inv) << std::endl;

    int jump_timestep = (steps%dt_inv);
    std::cout << "jump_timestep " << jump_timestep << std::endl;
    // erwünschte ausgabe: 0 bis 47
    // steps == 0 -> jump_timestep = 0
    // steps == 47 -> jump_timestep = 47
    // steps == 48 -> jump_timestep = 0

    if (eq[0].time()>T-1.1*dt) {
      dt = T-eq[0].time();
      done = true;
    }

    // update event data for new agents
    std::tie(filtered_event_data_size, agentID_location_commuting_start, agentID_location_commuting_end) = 
      abm.filter_event_data(agent_IDs_in_ABM, filtered_event_data, gen_local, 
                        steps, nr_day, verbosity); 

    // flatten agentID_healthStatus_map for ABM step
    agentID_healthStatus_map_flatten.clear();
    for (int state_ABM_idx = 0; state_ABM_idx < nr_ABM_states; state_ABM_idx++) {
      for (const auto& [agent_id, healthState]: agentID_healthStatus_map.at(state_ABM_idx)) {
        agentID_healthStatus_map_flatten[agent_id] = healthState;
      }
    }

    // do one ABM step
    std::tie(agentID_healthStatus_map, agent_IDs_in_ABM, healthStatus_agents_in_PDE, healthStatus_agents_in_ODE) = abm.step(std::move(filtered_event_data), filtered_event_data_size, with_coupling_states, agentID_location_commuting_start, agentID_location_commuting_end, agentID_healthStatus_map_flatten, steps, param_ABMs, nr_day, population_scale, gen_local, verbosity, run, correction_step); 

    agent_IDs_in_ABM_set.clear();
    agent_IDs_in_ABM_set.reserve(agent_IDs_in_ABM.size());
    agent_IDs_in_ABM_set.insert(agent_IDs_in_ABM.begin(), agent_IDs_in_ABM.end());

    // count symptomatic agents in ABM:
    for (int state_ABM_idx = 0; state_ABM_idx < nr_ABM_states; state_ABM_idx++) {
      for (const auto& agent_ID_healthState: agentID_healthStatus_map.at(state_ABM_idx)) {
        int healthState = agent_ID_healthState.second;
        if (healthState == 3) {
          nbr_symptomatic_agents.at(state_ABM_idx)[steps+1]++;
        }
      }
    }

    // do one PDE step
    std::vector<typename VariableSet::VariableSet> dx_PDE;
    for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) { // iterate over all PDE domains
      dx_PDE.push_back(typename VariableSet::VariableSet(x_PDE.at(state_PDE_idx)));

      try {
        #pragma omp critical
        {
          dx_PDE.at(state_PDE_idx) = limex_PDE.at(state_PDE_idx).step(x_PDE.at(state_PDE_idx),dt,extrapolOrder,tolX_PDE.at(state_PDE_idx)); // no mesh adaptation! tolX_PDE[] empty
        }
      }
      catch (const Kaskade::DirectSolverException& e) {
        std::cerr << "Fehler beim Lösen der Gleichung vom Bundesland " + stateLabels_PDE.at(state_PDE_idx) + ": " << e.what() << std::endl;
        return std::make_tuple(agentID_healthStatus_map, agent_IDs_in_ABM); 
      }

      // step ahead
      x_PDE.at(state_PDE_idx) += dx_PDE.at(state_PDE_idx);
      eq.at(state_PDE_idx).setTime(eq.at(state_PDE_idx).time()+dt);


      // compute how many people are in total in the PDE domain
      std::vector<std::vector<double>> compartmentJumpingPersons_test(nbr_points_PDE.at(state_PDE_idx),std::vector<double>(nbr_compartments,0.0));
      compartmentJumpingPersons_test = getNumberOfJumps<0>(x_PDE.at(state_PDE_idx), compartmentJumpingPersons_test, state_PDE_idx); 
      compartmentJumpingPersons_test = getNumberOfJumps<1>(x_PDE.at(state_PDE_idx), compartmentJumpingPersons_test, state_PDE_idx);
      compartmentJumpingPersons_test = getNumberOfJumps<2>(x_PDE.at(state_PDE_idx), compartmentJumpingPersons_test, state_PDE_idx);
      compartmentJumpingPersons_test = getNumberOfJumps<3>(x_PDE.at(state_PDE_idx), compartmentJumpingPersons_test, state_PDE_idx);
      compartmentJumpingPersons_test = getNumberOfJumps<4>(x_PDE.at(state_PDE_idx), compartmentJumpingPersons_test, state_PDE_idx);
      compartmentJumpingPersons_test = getNumberOfJumps<5>(x_PDE.at(state_PDE_idx), compartmentJumpingPersons_test, state_PDE_idx);
      compartmentJumpingPersons_test = getNumberOfJumps<6>(x_PDE.at(state_PDE_idx), compartmentJumpingPersons_test, state_PDE_idx);
      compartmentJumpingPersons_test = getNumberOfJumps<7>(x_PDE.at(state_PDE_idx), compartmentJumpingPersons_test, state_PDE_idx);

      double total_nbr_PDE_state = 0.0;
      for (int j=0; j<compartmentJumpingPersons_test[0].size(); j++) { //compartments
        for (int i=0; i<compartmentJumpingPersons_test.size(); i++) { // grid nodes
          total_nbr_PDE_state += compartmentJumpingPersons_test[i][j];
        }
      }
      n_PDE.at(state_PDE_idx)[steps+1] = total_nbr_PDE_state;
    }

    double act_change_notAtHome = activityChangePercentage[int(t_0+steps*dt+0.1)][1]/100.0;
    double factor_beta = (1 + act_change_notAtHome) * M_PI;
    if ((1 + act_change_notAtHome) > 1) {
      factor_beta = M_PI;
    }
    for (int state_ODE_idx=0; state_ODE_idx<nr_ODE_states; state_ODE_idx++) {
      double beta = factor_beta * params_ODEs[9+state_ODE_idx] / areas_ODE.at(state_ODE_idx); 
      auto& x = x_ODE[state_ODE_idx];

      double S_ode_before = x[0];
      double E_ode_before = x[1];
      double I_ode_before = x[2];
      double SY_ode_before = x[3];
      double H_ode_before = x[4];
      double C_ode_before = x[5];
      double HC_ode_before = x[6];
      double S_times_I_plus_SY_ode_before = S_ode_before * (I_ode_before + SY_ode_before);
      x[0] =  S_ode_before + dt * (- beta * S_times_I_plus_SY_ode_before);
      x[1] =  E_ode_before + dt * (+ beta * S_times_I_plus_SY_ode_before - params_ODEs[0]*E_ode_before); 
      x[2] =  I_ode_before + dt * (+ params_ODEs[0]*E_ode_before  - (params_ODEs[5]+params_ODEs[1])*I_ode_before);
      x[3] = SY_ode_before + dt * (+ params_ODEs[1]*I_ode_before  - (params_ODEs[6]+params_ODEs[2])*SY_ode_before);
      x[4] =  H_ode_before + dt * (+ params_ODEs[2]*SY_ode_before - (params_ODEs[7]+params_ODEs[3])*H_ode_before);
      x[5] =  C_ode_before + dt * (+ params_ODEs[3]*H_ode_before  -                 params_ODEs[4] *C_ode_before);
      x[6] = HC_ode_before + dt * (+ params_ODEs[4]*C_ode_before -  params_ODEs[8]                 *HC_ode_before);
      x[7] = x[7] + dt * (+ (params_ODEs[5]*I_ode_before + params_ODEs[6]*SY_ode_before + params_ODEs[7]*H_ode_before + params_ODEs[8]*HC_ode_before));

      for (int i=0; i<x.size(); i++) {
        if (x[i] < 0) {
          std::cout << "negative solution for i=" << i << ", state_ODE_idx=" << state_ODE_idx << ", after solution value " << x[i] << std::endl; 
          std::exit(EXIT_FAILURE);
        }
      }

      F_ODE.at(state_ODE_idx)[steps+1] = x[3]; // total number of symptomatic people in ODE domain
      n_ODE.at(state_ODE_idx)[steps+1] = x[0] + x[1] + x[2] + x[3] + x[4] + x[5] + x[6] + x[7]; // total number of people in ODE domain
    }

    double total_population_in_germany = 0.0;
    int cnt_agents_in_PDE_state = 0;
    for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) { // iterate over all PDE domains
      for (int compartment=0; compartment<nbr_compartments; compartment++) {
        cnt_agents_in_PDE_state += healthStatus_agents_in_PDE.at(state_PDE_idx).at(compartment); 
      }
    }
    total_population_in_germany += cnt_agents_in_PDE_state;
    if (verbosity > 0) {
      std::cout << "agents jumping out of ABM into PDE " << total_population_in_germany << std::endl;
    }

    int cnt_agents_in_ODE_state = 0;
    for (int state_ODE_idx=0; state_ODE_idx<nr_ODE_states; state_ODE_idx++) { // iterate over all ODE domains
      for (int compartment=0; compartment<nbr_compartments; compartment++) {
        cnt_agents_in_ODE_state += healthStatus_agents_in_ODE.at(state_ODE_idx).at(compartment);
      }
    }
    total_population_in_germany += cnt_agents_in_ODE_state;
    if (verbosity > 0) {
      std::cout << "agents jumping out of ABM into ODE " << cnt_agents_in_ODE_state << std::endl;
    }

    total_nbr_agents = agent_IDs_in_ABM.size(); 
    total_population_in_germany += total_nbr_agents;

    double total_nbr_ODE_individuals = 0.0;
    for (int state_ODE_idx=0; state_ODE_idx<nr_ODE_states; state_ODE_idx++) {
      total_nbr_ODE_individuals += n_ODE.at(state_ODE_idx)[steps+1];
    }
    total_population_in_germany += total_nbr_ODE_individuals;
    if (verbosity > 0) {
      std::cout << "total population in ODE states " << (total_nbr_ODE_individuals) << std::endl;
    }

    double total_nbr_PDE_individuals = 0.0;
    for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) {
      total_nbr_PDE_individuals += n_PDE.at(state_PDE_idx)[steps+1];
    }
    total_population_in_germany += total_nbr_PDE_individuals;
    if (verbosity > 0) {
      std::cout << "total population in PDE states " << (total_nbr_PDE_individuals) << std::endl;
      std::cout << "total population in ABM states should be " << (nbr_unique_agents - cnt_agents_in_ODE_state - cnt_agents_in_PDE_state) << std::endl;
      std::cout << "total population in ABM states is        " << (total_nbr_agents) << std::endl;
      std::cout << "total population in germany " << (total_population_in_germany) << std::endl;
    }
    if (int(total_population_in_germany+0.1) != nbr_trajectories) {
      std::cout << "Total population number changed! It should be " << nbr_trajectories << " but currently is " << int(total_population_in_germany+0.1) << std::endl;
      std::exit(EXIT_FAILURE);
    }
    //////////////////////////////////////////
    // add agents to PDE or ODE states
    // PDE
    std::vector<std::vector<std::vector<double>>> compartmentAddingPersons;
    compartmentAddingPersons.resize(nr_PDE_states);

    sign = 1; // -1: from PDE (to ABM), 1: (from ABM) to PDE
    for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) { // iterate over all PDE domains
      if (!with_coupling_states[pde_idx].at(state_PDE_idx)) {
        continue;
      }

      compartmentAddingPersons.at(state_PDE_idx) = std::vector<std::vector<double>>(nbr_points_PDE.at(state_PDE_idx), std::vector<double>(nbr_compartments, 0.0));
      for (int i=0; i<nbr_compartments; i++) {
        compartmentAddingPersons.at(state_PDE_idx)[0][i] = healthStatus_agents_in_PDE.at(state_PDE_idx)[i];
      }

      processCompartments(x_PDE.at(state_PDE_idx), compartmentAddingPersons.at(state_PDE_idx), sign, state_PDE_idx);
    }

    //ODE
    std::vector<std::vector<double>> ODEcompartmentAddingPersons(nr_ODE_states, std::vector<double>(nbr_compartments,0.0));
    sign = 1; // -1: from ODE (to ABM), 1: (from ABM) to ODE
    for (int state_ODE_idx=0; state_ODE_idx<nr_ODE_states; state_ODE_idx++) { // iterate over all ODE domains
      if (!with_coupling_states[ode_idx].at(state_ODE_idx)) { 
        continue;
      }

      for (int i=0; i<nbr_compartments; i++) {
        ODEcompartmentAddingPersons.at(state_ODE_idx)[i] = healthStatus_agents_in_ODE.at(state_ODE_idx)[i];
      }
      processODECompartments(x_ODE.at(state_ODE_idx), ODEcompartmentAddingPersons.at(state_ODE_idx), sign);
    }

    // empty ODEcompartmentJumpingIntoABMPersons and PDEcompartmentJumpingIntoABMPersons completely
    for (int ode = 0; ode < nr_ODE_states; ode++) {
      for (int abm = 0; abm < nr_ABM_states; abm++) {
        std::fill(ODEcompartmentJumpingIntoABMPersons[ode][abm].begin(), ODEcompartmentJumpingIntoABMPersons[ode][abm].end(), 0.0);
      }
    }
    for (int pde = 0; pde < nr_PDE_states; pde++) {
      for (int abm = 0; abm < nr_ABM_states; abm++) {
        for (auto& row : PDEcompartmentJumpingIntoABMPersons[pde][abm]) {
          std::fill(row.begin(), row.end(), 0.0);
        }
      }
    }

    int total_jumps_from_PDE_to_ABM = 0;
    for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) { // iterate over all PDE domains
      if (!with_coupling_states[pde_idx].at(state_PDE_idx)) {
        continue;
      }  
      // PDE -> ABM
      int state_idx_from = stateIndices_PDE.at(state_PDE_idx);
      for (int state_ABM_idx=0; state_ABM_idx<nr_ABM_states; state_ABM_idx++) { 
        if (!with_coupling_states[abm_idx].at(state_ABM_idx)) {
          continue;
        }

        int state_idx_into = stateIndices_ABM.at(state_ABM_idx);
        std::vector<int> jump_vector_agent_IDs = germany_jump_matrix_agent_IDs[nr_day][jump_timestep][state_idx_from][state_idx_into];
        for (int person_idx = 0; person_idx < jump_matrix[nr_day][jump_timestep][state_idx_from][state_idx_into]; person_idx++) { // if duplicate, then we ignore the agent ID because this agent never left the ABM in the first place due to activity reductions
          int agent_ID = jump_vector_agent_IDs[person_idx];
          
          if (agent_IDs_in_ABM_set.find(agent_ID) == agent_IDs_in_ABM_set.end()) {
            agent_IDs_in_ABM_set.insert(agent_ID);
            germany_jump_matrix_agent_IDs_with_reductions[jump_timestep][state_idx_from][state_idx_into].push_back(agent_ID);
          }
        }
        //////////////
        max_jumps = germany_jump_matrix_agent_IDs_with_reductions[jump_timestep][state_idx_from][state_idx_into].size();
        if (max_jumps > 0) {
          total_jumps_from_PDE_to_ABM += max_jumps;
          sign = -1; // -1: jump from PDE (to ABM): subtracting, 1: jump (from ABM) to PDE: adding
          bool targetIsNonABM = false;
          processCompartments(x_PDE.at(state_PDE_idx), PDEcompartmentJumpingIntoABMPersons.at(state_PDE_idx).at(state_ABM_idx), sign, state_PDE_idx, max_jumps, targetIsNonABM);
        }
      }

      // process jumps out of PDE into PDE
      for (int state_PDE_into_idx=0; state_PDE_into_idx<nr_PDE_states; state_PDE_into_idx++) { 
        if (state_PDE_idx == state_PDE_into_idx) { // can not jump from and into the same state
          continue;
        }
        if (!with_coupling_states[pde_idx].at(state_PDE_into_idx)) {
          continue;
        }  

        max_jumps = jump_matrix[nr_day][jump_timestep][stateIndices_PDE.at(state_PDE_idx)][stateIndices_PDE[state_PDE_into_idx]]; 
        if (max_jumps > 0) {
          sign = -1; // -1: jump from PDE (to PDE): subtracting
          processCompartments(x_PDE.at(state_PDE_idx), PDEcompartmentJumpingIntoPDEPersons.at(state_PDE_idx)[state_PDE_into_idx], sign, state_PDE_idx, max_jumps);
          //set entries to zero: PDEcompartmentJumpingIntoPDEPersons[state_PDE_into_idx].at(state_PDE_idx)
          for (auto& row : PDEcompartmentJumpingIntoPDEPersons[state_PDE_into_idx].at(state_PDE_idx)) {
            std::fill(row.begin(), row.end(), 0.0);
          }
          // iterate over PDEcompartmentJumpingIntoPDEPersons.at(state_PDE_idx)[state_PDE_into_idx] gridnodes and compartments and add to first grid node of PDEcompartmentJumpingIntoPDEPersons[state_PDE_into_idx].at(state_PDE_idx) (other grid!)
          for (int idx=0; idx<PDEcompartmentJumpingIntoPDEPersons.at(state_PDE_idx)[state_PDE_into_idx].size(); idx++) { //grid nodes
            for (int compartment=0; compartment<nbr_compartments; compartment++) {
              PDEcompartmentJumpingIntoPDEPersons[state_PDE_into_idx].at(state_PDE_idx)[0][compartment] += PDEcompartmentJumpingIntoPDEPersons.at(state_PDE_idx)[state_PDE_into_idx][idx][compartment];
            }
          }
          sign = 1; // 1: jump (from PDE) to PDE: adding
          processCompartments(x_PDE[state_PDE_into_idx], PDEcompartmentJumpingIntoPDEPersons[state_PDE_into_idx].at(state_PDE_idx), sign, state_PDE_into_idx);
        }
      }

      // process jumps out of PDE into ODE
      for (int state_ODE_idx=0; state_ODE_idx<nr_ODE_states; state_ODE_idx++) { 
        if (!with_coupling_states[ode_idx].at(state_ODE_idx)) {
          continue;
        }
        max_jumps = jump_matrix[nr_day][jump_timestep][stateIndices_PDE.at(state_PDE_idx)][stateIndices_ODE.at(state_ODE_idx)]; 
        
        if (max_jumps > 0) {
          sign = -1; // -1: jump from PDE (to ODE): subtracting
          processCompartments(x_PDE.at(state_PDE_idx), PDEcompartmentJumpingIntoODEPersons.at(state_PDE_idx).at(state_ODE_idx), sign, state_PDE_idx, max_jumps);
          
          // iterate over PDEcompartmentJumpingIntoODEPersons.at(state_PDE_idx).at(state_ODE_idx) gridnodes and compartments and add to first grid node
          for (int compartment=0; compartment<nbr_compartments; compartment++) {
            for (int idx=1; idx<PDEcompartmentJumpingIntoODEPersons.at(state_PDE_idx).at(state_ODE_idx).size(); idx++) { //grid nodes
              PDEcompartmentJumpingIntoODEPersons.at(state_PDE_idx).at(state_ODE_idx)[0][compartment] += PDEcompartmentJumpingIntoODEPersons.at(state_PDE_idx).at(state_ODE_idx)[idx][compartment];
            }
          }

          sign = 1; // 1: jump (from PDE) to ODE: adding
          processODECompartments(x_ODE.at(state_ODE_idx), PDEcompartmentJumpingIntoODEPersons.at(state_PDE_idx).at(state_ODE_idx)[0], sign);
        }
      }
      
      if ( verbosity>0 ) {
        std::cout.flush();
        std::cout << "t= " << eq.at(state_PDE_idx).time()+dt << " dt = " << dt << '\n';

        // save solution
        std::ostringstream fm;
        fm << "../../work/output/output_data_optim_" + std::to_string(correction_step) + "/graph_abm_pde_ode/covid-slab-uv-" + stateLabels_PDE.at(state_PDE_idx) + "-" + std::to_string(run) + "-";
        fm.width(3);
        fm.fill('0');
        fm.setf(std::ios_base::right,std::ios_base::adjustfield);
        fm << (steps+1+(init_time_step*dt_inv));
        fm.flush(); 
        #pragma omp critical(vtk_write)
        {
          writeVTK(x_PDE.at(state_PDE_idx),fm.str(),IoOptions().setPrecision(10),data_F.at(state_PDE_idx),data_DF.at(state_PDE_idx),variables_symptomatic,steps+1,ones.at(state_PDE_idx),verbosity); //density
          if (run > 0 || steps % 8 != 0) { // do not delete every 8. file
            fm << ".vtu"; 
            if (std::remove(fm.str().c_str()) != 0) { // delete file
              std::perror("Fehler beim Löschen der Datei");
            } 
          }
        }

        // approximate integral of symptomatic over each domain
        for (auto const& cell: elements(x_PDE.at(state_PDE_idx).descriptions.gridView)) {
          for (int corner=0; corner<cell.subEntities(2); corner++) {
            auto xglob = cell.geometry().corner(corner);
            corners[corner][0] = xglob[0];
            corners[corner][1] = xglob[1];
          }
          area = 0.5 * std::abs(corners[0][0] * (corners[1][1]-corners[2][1]) + corners[1][0] * (corners[2][1]-corners[0][1]) + corners[2][0] * (corners[0][1]-corners[1][1]) );

          for (int corner=0; corner<cell.subEntities(2); corner++) {
            int index = x_PDE.at(state_PDE_idx).descriptions.gridView.indexSet().subIndex(cell, corner, 2); //global index
            F_PDE.at(state_PDE_idx)[steps+1] += (area/3) * data_F.at(state_PDE_idx)[index+nbr_points_PDE.at(state_PDE_idx)*(steps+1)];
          } // end for corner
        } // end for cell
        std::cout << "symptomatic in " << stateLabels_PDE.at(state_PDE_idx) << ": " << (F_PDE.at(state_PDE_idx)[steps+1]) << std::endl;
      } //verbosity
    } // end for state_PDE_idx
    if (verbosity > 0) {
      std::cout << "total_jumps_from_PDE_to_ABM " << total_jumps_from_PDE_to_ABM << std::endl;
    }

    int total_jumps_from_ODE_to_ABM = 0;
    // compute number of jumps and perform jumps outgoing from ODE states (into ABM, PDE, ODE states)
    for (int state_ODE_idx=0; state_ODE_idx<nr_ODE_states; state_ODE_idx++) { // iterate over all ODE domains
      if (!with_coupling_states[ode_idx].at(state_ODE_idx)) {
        continue;
      }
      // process jumps out of ODE (-> ABM)
      int state_idx_from = stateIndices_ODE.at(state_ODE_idx);
      for (int state_ABM_idx=0; state_ABM_idx<nr_ABM_states; state_ABM_idx++) { 
        if (!with_coupling_states[abm_idx].at(state_ABM_idx)) {
          continue;
        }

        int state_idx_into = stateIndices_ABM.at(state_ABM_idx);
        std::vector<int> jump_vector_agent_IDs = germany_jump_matrix_agent_IDs[nr_day][jump_timestep][state_idx_from][state_idx_into];
        for (int person_idx = 0; person_idx < jump_vector_agent_IDs.size(); person_idx++) { // if duplicate, then we ignore the agent ID because this agent never left the ABM in the first place due to activity reductions
          int agent_ID = jump_vector_agent_IDs[person_idx];

          if (agent_IDs_in_ABM_set.find(agent_ID) == agent_IDs_in_ABM_set.end()) {
            agent_IDs_in_ABM_set.insert(agent_ID);
            germany_jump_matrix_agent_IDs_with_reductions[jump_timestep][state_idx_from][state_idx_into].push_back(agent_ID);
          }
        }
        //////////////
        max_jumps = germany_jump_matrix_agent_IDs_with_reductions[jump_timestep][state_idx_from][state_idx_into].size();

        if (max_jumps > 0) {          
          total_jumps_from_ODE_to_ABM += max_jumps;
          sign = -1; // -1: jump from ODE (to ABM): subtracting, 1: jump (from ABM) to ODE: adding
          bool targetIsNonABM = false;
          processODECompartments(x_ODE.at(state_ODE_idx), ODEcompartmentJumpingIntoABMPersons.at(state_ODE_idx).at(state_ABM_idx), sign, max_jumps, targetIsNonABM);
        }
      }

      // process jumps out of ODE into ODE
      for (int state_ODE_into_idx=0; state_ODE_into_idx<nr_ODE_states; state_ODE_into_idx++) { 
        if (state_ODE_idx == state_ODE_into_idx) { // can not jump from and into the same state
          continue;
        }
        if (!with_coupling_states[ode_idx][state_ODE_into_idx]) {
          continue;
        }  
        max_jumps = jump_matrix[nr_day][jump_timestep][stateIndices_ODE.at(state_ODE_idx)][stateIndices_ODE[state_ODE_into_idx]]; 
        
        if (max_jumps > 0) {
          sign = -1; // -1: jump from ODE (to ODE): subtracting
          processODECompartments(x_ODE.at(state_ODE_idx), ODEcompartmentJumpingIntoODEPersons.at(state_ODE_idx)[state_ODE_into_idx], sign, max_jumps);
          sign = 1; // 1: jump (from ODE) to ODE: adding
          processODECompartments(x_ODE[state_ODE_into_idx], ODEcompartmentJumpingIntoODEPersons.at(state_ODE_idx)[state_ODE_into_idx], sign);
        }
      }

      // process jumps out of ODE into PDE
      for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) { 
        if ((!with_coupling_states[pde_idx].at(state_PDE_idx))) {
          continue;
        }
        max_jumps = jump_matrix[nr_day][jump_timestep][stateIndices_ODE.at(state_ODE_idx)][stateIndices_PDE.at(state_PDE_idx)]; 
        
        if (max_jumps > 0) {
          sign = -1; // -1: jump from ODE (to PDE): subtracting
          processODECompartments(x_ODE.at(state_ODE_idx), ODEcompartmentJumpingIntoPDEPersons.at(state_ODE_idx).at(state_PDE_idx)[0], sign, max_jumps);

          sign = 1; // 1: jump (from ODE) to PDE: adding
          processCompartments(x_PDE.at(state_PDE_idx), ODEcompartmentJumpingIntoPDEPersons.at(state_ODE_idx).at(state_PDE_idx), sign, state_PDE_idx);
        }
      }
    } // end for state_ODE_idx
    nbr_unique_agents = total_nbr_agents + total_jumps_from_ODE_to_ABM + total_jumps_from_PDE_to_ABM;
    
    if (verbosity > 0) {
      std::cout << "total_jumps_from_ODE_to_ABM " << total_jumps_from_ODE_to_ABM << std::endl;
      std::cout << "unique agents in map should be " << nbr_unique_agents << std::endl;
    }

    /////////////////////////////
    // add new agents that came from PDE and ODE domain
    for (auto& domain : cnt_jumping_persons_per_compartment) {
        for (auto& ABMstate : domain) {
            std::fill(ABMstate.begin(), ABMstate.end(), 0);
        }
    }
    // count new agents that came from ODE domain
    for (int state_ODE_idx=0; state_ODE_idx<nr_ODE_states; state_ODE_idx++) { // iterate over all ODE domains
      int state_idx = stateIndices_ODE.at(state_ODE_idx);
      for (int state_ABM_idx=0; state_ABM_idx<nr_ABM_states; state_ABM_idx++) { // iterate over all ABM domains
        for (int compartment=0; compartment<nbr_compartments; compartment++) { 
          cnt_jumping_persons_per_compartment[state_idx].at(state_ABM_idx)[compartment] += int(ODEcompartmentJumpingIntoABMPersons.at(state_ODE_idx).at(state_ABM_idx)[compartment] + 0.1); // + 0.1 to avoid rounding error when converting to int
        }
      }
    }

    // count new agents that came from PDE domain
    for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) { // iterate over all PDE domains
      int state_idx = stateIndices_PDE.at(state_PDE_idx);
      for (int state_ABM_idx=0; state_ABM_idx<nr_ABM_states; state_ABM_idx++) { // iterate over all PDE domains
        for (int compartment=0; compartment<nbr_compartments; compartment++) {
          double sum_over_grid = 0.0;
          for (int i=0; i<PDEcompartmentJumpingIntoABMPersons.at(state_PDE_idx).at(state_ABM_idx).size(); i++) { // grid nodes
            sum_over_grid += PDEcompartmentJumpingIntoABMPersons.at(state_PDE_idx).at(state_ABM_idx)[i][compartment];
          }
          cnt_jumping_persons_per_compartment[state_idx].at(state_ABM_idx)[compartment] += int(sum_over_grid + 0.1);  // + 0.1 to avoid rounding error when converting to int
        }
      }
    }

    // set health status of new agents accordingly to how many new people are in each compartment
    for (int state_ABM_idx=0; state_ABM_idx<nr_ABM_states; state_ABM_idx++) { 
      int state_idx_into = stateIndices_ABM.at(state_ABM_idx);
      for (int state_idx=0; state_idx<nbr_domains; state_idx++) { // iterate over all domains
        int person_idx = 0;
        // shuffle the agent IDs before using them! if not shuffled, the lower agent IDs will always be susceptible (if there are any susceptible)
        std::shuffle(germany_jump_matrix_agent_IDs_with_reductions[jump_timestep][state_idx][state_idx_into].begin(), germany_jump_matrix_agent_IDs_with_reductions[jump_timestep][state_idx][state_idx_into].end(), gen_local);

        for (int compartment = 0; compartment < nbr_compartments; compartment++) { 
          for (int j = 0; j < cnt_jumping_persons_per_compartment[state_idx].at(state_ABM_idx)[compartment]; j++) { 
            if (person_idx >= germany_jump_matrix_agent_IDs_with_reductions[jump_timestep][state_idx][state_idx_into].size()) {
              std::cout << "germany_jump_matrix_agent_IDs_with_reductions[jump_timestep][state_idx][state_idx_into].size() " << germany_jump_matrix_agent_IDs_with_reductions[jump_timestep][state_idx][state_idx_into].size() << std::endl;
              std::cout << "cnt_jumping_persons_per_compartment[state_idx].at(state_ABM_idx)[compartment] " << cnt_jumping_persons_per_compartment[state_idx].at(state_ABM_idx)[compartment] << std::endl;
              std::cout << "jump_matrix[nr_day][jump_timestep][state_idx][state_idx_into] " << jump_matrix[nr_day][jump_timestep][state_idx][state_idx_into] << std::endl;
              std::cout << "jump_timestep " << (jump_timestep) << ", state_idx " << state_idx << ", state_idx_into " << state_idx_into << ", compartment " << compartment << std::endl;
              std::exit(EXIT_FAILURE);
            }

            int agent_ID = germany_jump_matrix_agent_IDs_with_reductions[jump_timestep][state_idx][state_idx_into][person_idx];
              agentID_healthStatus_map.at(state_ABM_idx)[agent_ID] = compartment; // add new entries to map
              agent_IDs_in_ABM.push_back(agent_ID);
            person_idx++; // if duplicate, then we ignore the agent ID because this agent never left the ABM in the first place due to activity reductions
          }
        }
      } // end for domains
    } // end for ABM domains
    std::sort(agent_IDs_in_ABM.begin(), agent_IDs_in_ABM.end()); // sort by id
    // end jumps
    ///////////////////////////////////////// 

    // begin end of day / week / two week code (last block before next time step iteration)
    if ((steps+1) % dt_inv == 0) { // new day 
      std::cout << "day " << (int((steps+1) * dt)) << std::endl;

      // alle susceptible agent IDs sammeln, shuffeln, die ersten 0.001% zu state nr 1 ändern
      for (int state_ABM_idx = 0; state_ABM_idx < nr_ABM_states; state_ABM_idx++) {
        if (!with_coupling_states[abm_idx].at(state_ABM_idx)) {
          continue; // no disease import
        }
        std::vector<int> susceptible_agent_ids;
        for (const auto& agent_ID_healthState: agentID_healthStatus_map.at(state_ABM_idx)) {
          int healthState = agent_ID_healthState.second;
          if (healthState == 0) {
            susceptible_agent_ids.push_back(agent_ID_healthState.first);
          }
        }

        double rate_S_E = 0.00001;
        int n_change = static_cast<int>(rate_S_E * susceptible_agent_ids.size()); 
        if (verbosity > 0) {
          std::cout << "for ABM state " << stateLabels_ABM.at(state_ABM_idx) << " we change " << n_change << " susceptible(s) to exposed at rate " << rate_S_E << std::endl;
        }

        if (n_change > 0) {
          std::shuffle(susceptible_agent_ids.begin(), susceptible_agent_ids.end(), gen_local);
        }

        for (int i = 0; i < n_change; i++) {
          int agent_id = susceptible_agent_ids[i];
          agentID_healthStatus_map.at(state_ABM_idx)[agent_id] = 1;
        }
      }
    
      // because jump_timestep will be reset at the beginning of the day, germany_jump_matrix_agent_IDs_with_reductions needs to be emptied again
      germany_jump_matrix_agent_IDs_with_reductions = std::vector<std::vector<std::vector<std::vector<int>>>>(dt_inv, std::vector<std::vector<std::vector<int>>>(nbr_domains, std::vector<std::vector<int>>(nbr_domains)));    

      week_day++;
      nr_day = get_nr_day(week_day);
      std::cout << "nr_day 0/1/2 " << nr_day << std::endl;
    }
  } // end for steps
  std::cout << '\n';
  limex_PDE[0].reportTime(std::cout);

  if (verbosity<=0) {
    return std::make_tuple(agentID_healthStatus_map, agent_IDs_in_ABM); 
  }

  int days = maxSteps / dt_inv;
  std::vector<double> target(days, 0.0);
  std::vector<std::vector<double>> target_ABM(nr_ABM_states,std::vector<double>(days));
  std::vector<std::vector<double>> target_PDE(nr_PDE_states,std::vector<double>(days));
  std::vector<std::vector<double>> target_ODE(nr_ODE_states,std::vector<double>(days));
  std::vector<double> result_days(days, 0.0);
  std::vector<std::vector<double>> result_days_ABM(nr_ABM_states,std::vector<double>(days, 0.0));
  std::vector<std::vector<double>> result_days_PDE(nr_PDE_states,std::vector<double>(days, 0.0));
  std::vector<std::vector<double>> result_days_ODE(nr_ODE_states,std::vector<double>(days, 0.0));

  double norm = 0.0;
  for (int day = 0; day < days; ++day) {
    for (int state_ODE_idx=0; state_ODE_idx<nr_ODE_states; state_ODE_idx++) { // iterate over all ODE domains
      for (int hour = 0; hour < dt_inv; ++hour) {
        result_days_ODE.at(state_ODE_idx)[day] += F_ODE.at(state_ODE_idx)[day*dt_inv+hour];
      }
      result_days_ODE.at(state_ODE_idx)[day] /= dt_inv;
      std::cout << "numerical solution, only F_ODE.at(state_ODE_idx) " << result_days_ODE.at(state_ODE_idx)[day] << std::endl;

      result_days_ODE.at(state_ODE_idx)[day] *= population_scale;
      result_days[day] += result_days_ODE.at(state_ODE_idx)[day];

      target_ODE.at(state_ODE_idx)[day] = exposed[run-initial_run][day+init_time_step][ode_idx].at(state_ODE_idx);
      target_ODE.at(state_ODE_idx)[day] *= population_scale;
      target[day] += target_ODE.at(state_ODE_idx)[day];
    } // end for ODE domains

    std::cout << "numerical solution " << result_days[day] << ", target " << target[day] << ", domain " << pde_idx << std::endl;
    norm += std::pow(result_days[day]-target[day],2);
  } // end for day
  std::cout << "norm " << std::sqrt(norm) << std::endl;

  for (int day = 0; day < days; ++day) {
    for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) { // iterate over all PDE domains
      for (int hour = 0; hour < dt_inv; ++hour) {
        result_days_PDE.at(state_PDE_idx)[day] += F_PDE.at(state_PDE_idx)[day*dt_inv+hour];
      }
      result_days_PDE.at(state_PDE_idx)[day] /= dt_inv;
      std::cout << "numerical solution, only F_PDE.at(state_PDE_idx) " << result_days_PDE.at(state_PDE_idx)[day] << std::endl;

      result_days_PDE.at(state_PDE_idx)[day] *= population_scale;
      result_days[day] += result_days_PDE.at(state_PDE_idx)[day];

      target_PDE.at(state_PDE_idx)[day] = exposed[run-initial_run][day+init_time_step][pde_idx].at(state_PDE_idx);
      target_PDE.at(state_PDE_idx)[day] *= population_scale;
      target[day] += target_PDE.at(state_PDE_idx)[day];
    } // end for PDE domains

    std::cout << "numerical solution " << result_days[day] << ", target " << target[day] << ", domain " << pde_idx << std::endl;
    norm += std::pow(result_days[day]-target[day],2);
  } // end for day
  std::cout << "norm " << std::sqrt(norm) << std::endl;

  for (int day = 0; day < days; ++day) {
    for (int state_ABM_idx=0; state_ABM_idx<nr_ABM_states; state_ABM_idx++) { // iterate over all ABM domains
      for (int hour = 0; hour < dt_inv; ++hour) {
        result_days_ABM.at(state_ABM_idx)[day] += nbr_symptomatic_agents.at(state_ABM_idx)[day*dt_inv+hour];
      }
      result_days_ABM.at(state_ABM_idx)[day] /= dt_inv;
      std::cout << "numerical solution, only nbr_symptomatic_agents " << result_days_ABM.at(state_ABM_idx)[day] << std::endl;

      result_days_ABM.at(state_ABM_idx)[day] *= population_scale; // compute for 100% of population
      result_days[day] += result_days_ABM.at(state_ABM_idx)[day];

      target_ABM.at(state_ABM_idx)[day] = exposed[run-initial_run][day+init_time_step][abm_idx].at(state_ABM_idx);
      target_ABM.at(state_ABM_idx)[day] *= population_scale;
      target[day] += target_ABM.at(state_ABM_idx)[day];
    } // end for ABM domains

    std::cout << "numerical solution " << result_days[day] << ", target " << target[day] << ", domain " << abm_idx << std::endl;
    norm += std::pow(result_days[day]-target[day],2);
  } // end for day
  std::cout << "norm " << std::sqrt(norm) << std::endl;

  std::cout << "run " << run << std::endl;

  // save result in file
  std::string filename_result = "../../work/output/output_data_optim_" + std::to_string(correction_step) + "/result_0_" + std::to_string(run) + addOn + ".txt";
  std::ofstream file_result;
  file_result.open(filename_result);
  for (int day = 0; day < result_days.size(); day++) {
    file_result << result_days[day] << std::endl;
  }
  file_result.close();


  for (int state_ODE_idx=0; state_ODE_idx<nr_ODE_states; state_ODE_idx++) { // iterate over all ODE domains
    std::string filename_result_ODE = "../../work/output/output_data_optim_" + std::to_string(correction_step) + "/result_ODE_0_" + std::to_string(run) + "_" + stateLabels_ODE.at(state_ODE_idx) + addOn + ".txt";
    std::ofstream file_result_ODE;
    file_result_ODE.open(filename_result_ODE);
    for (int day = 0; day < result_days_ODE.at(state_ODE_idx).size(); day++) {
      file_result_ODE << result_days_ODE.at(state_ODE_idx)[day] << std::endl;
    }
    file_result_ODE.close();

    if (run == 0) {
      std::string filename_target_ODE = "../../work/output/output_data_optim_" + std::to_string(correction_step) + "/target_ODE_0_" + stateLabels_ODE.at(state_ODE_idx) + ".txt";
      std::ofstream file_target_ODE;
      file_target_ODE.open(filename_target_ODE);
      for (int day = 0; day < target_ODE.at(state_ODE_idx).size(); day++) {
        file_target_ODE << target_ODE.at(state_ODE_idx)[day] << std::endl;
      }
      file_target_ODE.close();
    }
  }

  for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) { // iterate over all PDE domains
    std::string filename_result_PDE = "../../work/output/output_data_optim_" + std::to_string(correction_step) + "/result_PDE_0_" + std::to_string(run) + "_" + stateLabels_PDE.at(state_PDE_idx) + addOn + ".txt";
    std::ofstream file_result_PDE;
    file_result_PDE.open(filename_result_PDE); //, std::ios_base::app);
    for (int day = 0; day < result_days_PDE.at(state_PDE_idx).size(); day++) {
      file_result_PDE << result_days_PDE.at(state_PDE_idx)[day] << std::endl;
    }
    file_result_PDE.close();

    if (run == 0) {
      std::string filename_target_PDE = "../../work/output/output_data_optim_" + std::to_string(correction_step) + "/target_PDE_0_" + stateLabels_PDE.at(state_PDE_idx) + ".txt";
      std::ofstream file_target_PDE;
      file_target_PDE.open(filename_target_PDE); //, std::ios_base::app);
      for (int day = 0; day < target_PDE.at(state_PDE_idx).size(); day++) {
        file_target_PDE << target_PDE.at(state_PDE_idx)[day] << std::endl;
      }
      file_target_PDE.close();
    }
  }

  for (int state_ABM_idx=0; state_ABM_idx<nr_ABM_states; state_ABM_idx++) { // iterate over all ABM domains
    std::string filename_result_ABM = "../../work/output/output_data_optim_" + std::to_string(correction_step) + "/result_ABM_0_" + std::to_string(run) + "_" + stateLabels_ABM.at(state_ABM_idx) + addOn + ".txt";
    std::ofstream file_result_ABM;
    file_result_ABM.open(filename_result_ABM);
    for (int day = 0; day < result_days_ABM.at(state_ABM_idx).size(); day++) {
      file_result_ABM << result_days_ABM.at(state_ABM_idx)[day] << std::endl;
    }
    file_result_ABM.close();

    if (run == 0) {
      std::string filename_target_ABM = "../../work/output/output_data_optim_" + std::to_string(correction_step) + "/target_ABM_0_" + stateLabels_ABM.at(state_ABM_idx) + ".txt";
      std::ofstream file_target_ABM;
      file_target_ABM.open(filename_target_ABM);
      for (int day = 0; day < target_ABM.at(state_ABM_idx).size(); day++) {
        file_target_ABM << target_ABM.at(state_ABM_idx)[day] << std::endl;
      }
      file_target_ABM.close();
    }
  }

  // save error/difference in file
  std::string filename_error = "../../work/output/output_data_optim_" + std::to_string(correction_step) + "/error_0_" + std::to_string(run) + addOn + ".txt";
  std::ofstream file_error;
  file_error.open(filename_error);
  for (int day = 0; day < result_days.size(); day++) {
    file_error << (result_days[day] - target[day]) << std::endl;
  }
  file_error.close();

  if (run == 0) {
    plot_result(result_days, target, correction_step, addOn);
    for (int state_ABM_idx=0; state_ABM_idx<nr_ABM_states; state_ABM_idx++) { // iterate over all ABM domains
      plot_result(result_days_ABM.at(state_ABM_idx), target_ABM.at(state_ABM_idx), correction_step, "_ABM_" + stateLabels_ABM.at(state_ABM_idx) + addOn);
    }
    for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) { // iterate over all PDE domains
      plot_result(result_days_PDE.at(state_PDE_idx), target_PDE.at(state_PDE_idx), correction_step, "_PDE_" + stateLabels_PDE.at(state_PDE_idx) + addOn);
    }
    for (int state_ODE_idx=0; state_ODE_idx<nr_ODE_states; state_ODE_idx++) { // iterate over all ODE domains
      plot_result(result_days_ODE.at(state_ODE_idx), target_ODE.at(state_ODE_idx), correction_step, "_ODE_" + stateLabels_ODE.at(state_ODE_idx) + addOn);
    }
  }
  if (!done)
  std::cout << "*** maxSteps reached ***\n";
  return std::make_tuple(agentID_healthStatus_map, agent_IDs_in_ABM);
}
#endif