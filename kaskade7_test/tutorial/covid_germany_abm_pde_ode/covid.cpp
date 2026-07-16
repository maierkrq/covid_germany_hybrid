#define FUSION_MAX_VECTOR_SIZE 20

#include <iostream>

#include <iomanip> // std::scientific, setprecision
#include <fstream> //write txt file

#include <thread>
#include <chrono>

#include <random> //default_random_engine

#include <cmath>

#include <boost/timer/timer.hpp>
#include <boost/geometry.hpp>

#include "dune/grid/config.h"
#include "dune/grid/uggrid.hh"

#include "utilities/enums.hh"
#include "fem/assemble.hh"
#include "fem/gridmanager.hh"
#include "fem/lagrangespace.hh"
#include "fem/functional_aux.hh"
#include "fem/coarsening.hh"
#include "fem/norms.hh"

#include "linalg/direct.hh"
//#include "linalg/cg.hh"               // alternative: use Kaskade7 CG
//#include "linalg/iluprecond.hh"       // with ILUK preconditioner

#include "timestepping/limexWithoutJens.hh"


#include "utilities/kaskopt.hh"
#include "utilities/gridGeneration.hh"

#include "fem/boundaryInterpolation.hh" //BoundaryEdge
#include "fem/deforminggridmanager.hh" // DeformingGridManagerDetail::Deformation

#include <omp.h>

using namespace Kaskade; 

#include "covid-integrate.hh"
#include "covid.hh"

#include "io/readPoly.hh"
#include "io/readPoly.cpp"


template<typename T>
struct Initial_1Value {
  using Scalar = double;
  static int const components = 1;
  using ValueType = Dune::FieldVector<Scalar,components>;

  Initial_1Value(int c, const T& xx, std::vector<double> VV, int sstateIdx, int rrun): component_int(c), x(xx), normalized_inv_V(VV), stateIdx(sstateIdx), run(rrun) {} 
  template <class Cell> int order(Cell const&) const { return std::numeric_limits<int>::max(); }
  template <class Cell>
  ValueType value(Cell const& cell,
                  Dune::FieldVector<typename Cell::Geometry::ctype,Cell::Geometry::coorddimension> const& localCoordinate) const 
  {
    // get index of global coordinates 
    double error_min = std::numeric_limits<int>::max();
    double error;
    int index_min;
    Dune::FieldVector<typename Cell::Geometry::ctype,Cell::Geometry::coorddimension> x_glob = cell.geometry().global(localCoordinate);
    for (int idx=0; idx<cell.subEntities(2); idx++) {
      error = (cell.geometry().corner(idx)-x_glob).two_norm(); 
      if (error < error_min) {
        error_min = error;
        index_min = idx;
      }
    }
    int glob_idx = x.descriptions.gridView.indexSet().subIndex(cell,index_min,2); //global index
    
    double s = susceptibles[run][0][pde_idx][stateIdx] * normalized_inv_V[glob_idx];
    double e = exposed[run][0][pde_idx][stateIdx] * normalized_inv_V[glob_idx];
    double r = recovered[run][0][pde_idx][stateIdx] * normalized_inv_V[glob_idx];

    double i = 0.0;
    double sy = 0.0;
    double h = 0.0;
    double c = 0.0;
    double hc = 0.0;

    if (component_int==0) return s; 
    else if (component_int==1) return e; 
    else if (component_int==2) return i;
    else if (component_int==3) return sy;
    else if (component_int==4) return h;
    else if (component_int==5) return c;
    else if (component_int==6) return hc;
    else if (component_int==7) return r;
    else
      assert("wrong index!\n"==0);
    return 0;
  }

  private:
  int component_int, stateIdx, run;
  const T& x;
  std::vector<double> normalized_inv_V;
};

void delete_non_relevant_events(const std::vector<int>& no_coupling_abm_indices, ABModel& abm) {
  const std::unordered_set<int> no_coupling_set(no_coupling_abm_indices.begin(), no_coupling_abm_indices.end());

  for (int nr_day=0; nr_day<3; nr_day++) {
    auto& dayData = abm.eventData.at(nr_day);
    const auto& facilityLabels_allCategories_nr_day = abm.facilityLabels_allCategories.at(nr_day);

    std::size_t write_idx = 0;
    int previous_agent_id = -1;
    int initial_federal_state = -1;
    for (std::size_t read_idx = 0; read_idx < dayData.size(); read_idx++) {
      const auto& row = dayData[read_idx];

      int agent_id = row[2];
      int facility_id = row[3];
      int category_id = row[4];
      int federal_state = facilityLabels_allCategories_nr_day[category_id][facility_id];

      bool keep = true;
      if (agent_id != previous_agent_id) {
        initial_federal_state = federal_state;
        previous_agent_id = agent_id;
      } 
      else if ((federal_state != initial_federal_state) && ((no_coupling_set.find(federal_state) != no_coupling_set.end()) || (no_coupling_set.find(initial_federal_state) != no_coupling_set.end()))) {
          keep = false;
      }

      if (keep) {
        if (write_idx != read_idx) {
          dayData[write_idx] = std::move(dayData[read_idx]);
        }
        write_idx++;
      }
    }

    dayData.resize(write_idx);
  }
}

int main(int argc, char *argv[]) {
  int num_threads = 1; //10; //
  std::cout << "number of threads " << num_threads << std::endl;
  omp_set_num_threads(num_threads); 

  int initial_run = 0;
  int maxRun = 1; //2; //10; //

  int dt_inv = 48;
  
  bool with_restriction = true; //false;

  int week_day = 0; // Monday, 2nd March, 2020 
  //6; // Sunday 
  int nr_day = get_nr_day(week_day);
  initial_nr_day = nr_day;
  initial_week_day = week_day;

  std::string no_coupling_state = ""; //"Berlin"; //"Nordrhein-Westfalen"; //
  std::vector<std::vector<bool>> with_coupling_states(nbr_model_types);
  for (int model_type=0; model_type<nbr_model_types; model_type++) {
    with_coupling_states[model_type].resize(stateLabels_all_models[model_type].size());
    for (int state_idx=0; state_idx<stateLabels_all_models[model_type].size(); state_idx++) {
      if (stateLabels_all_models[model_type][state_idx] == no_coupling_state) {
        with_coupling_states[model_type][state_idx] = false;
      }
      else {
        with_coupling_states[model_type][state_idx] = true;
      }
      std::cout << "model_type " << model_type << ", state_idx " << state_idx << ": " << with_coupling_states[model_type][state_idx] << std::endl;
    }
  }

  std::cout << "set local_index of domain" << std::endl;
  for (int model_type=0; model_type<nbr_model_types; model_type++) {
    for (int state_idx=0; state_idx<stateIndices_all_models[model_type].size(); state_idx++) {
      int domain = stateIndices_all_models[model_type][state_idx];
      std::cout << "domain " << domain << ", model_type " << model_type << std::endl;
      model_type_of_domain[domain] = model_type;
      local_index[domain] = state_idx;
    }
  }

  // 0: coarsest, 1: coarse, 2: fine, 3: finer(, 4: finest)
  std::vector<int> coarsity_PDEs(nr_PDE_states, 0); 

  std::vector<std::string> coarsity_str_PDEs;
  for (auto const& coarsity : coarsity_PDEs) {
    if (coarsity == 0) {
      coarsity_str_PDEs.push_back("_coarse"); //_coarsest
    }
    else if (coarsity == 1) {
      coarsity_str_PDEs.push_back(""); //_coarse
    }
    else if (coarsity == 2) {
      coarsity_str_PDEs.push_back("_fine");
    }
    else if (coarsity == 3) {
      coarsity_str_PDEs.push_back("_finer");
    }
    else if (coarsity == 4) {
      coarsity_str_PDEs.push_back("_finest");
    }
    else {
      std::cout << "Invalid coarsity level" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  double population_scale = 4.0;

  std::vector<ABModel> all_abms;
  all_abms.reserve(maxRun - initial_run);
  for (int run = initial_run; run < maxRun; run++) {
      all_abms.emplace_back(population_scale, with_restriction, dt_inv);
      all_abms.back().readEventData();
      all_abms.back().readFacilityCoordinates_allCategories();
      all_abms.back().readFacilityCoordinates_home();
      all_abms.back().set_initial_agent_ids();
      all_abms.back().set_max_number_of_agents_per_facility();
      all_abms.back().setLeavingNumberOfPeople();
      all_abms.back().readActivityChangeData();
  }

  population[pde_idx].resize(nr_PDE_states);
  for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) {
    population[pde_idx].at(state_PDE_idx) = population_pde_t_0.at(state_PDE_idx);
    std::cout << std::scientific << std::setprecision(8) << "population " << stateLabels_PDE.at(state_PDE_idx) << " " << population[pde_idx].at(state_PDE_idx) << std::endl;
  }

  population[abm_idx].resize(nr_ABM_states);
  for (int state_ABM_idx=0; state_ABM_idx<nr_ABM_states; state_ABM_idx++) {
    population[abm_idx].at(state_ABM_idx) = all_abms[0].nbr_agents[nr_day].at(state_ABM_idx); // all_abms[0] because they are the same in terms of number of agents
    std::cout << std::scientific << std::setprecision(8) << "population " << stateLabels_ABM.at(state_ABM_idx) << " " << population[abm_idx].at(state_ABM_idx) << std::endl;
  }

  population[ode_idx].resize(nr_ODE_states);
  for (int state_ODE_idx=0; state_ODE_idx<nr_ODE_states; state_ODE_idx++) {
    population[ode_idx].at(state_ODE_idx) = population_ode_t_0.at(state_ODE_idx);
    std::cout << std::scientific << std::setprecision(8) << "population " << stateLabels_ODE.at(state_ODE_idx) << " " << population[ode_idx].at(state_ODE_idx) << std::endl;
  }

  // check if we consider ABM federal states with no transitions into or out of the respective area
  bool no_coupling_abm = false;
  std::vector<int> no_coupling_abm_indices;
  for (int state_ABM_idx=0; state_ABM_idx<nr_ABM_states; state_ABM_idx++) { // iterate over all ABM domains
      if (!with_coupling_states[abm_idx].at(state_ABM_idx)) {
          no_coupling_abm = true;
          no_coupling_abm_indices.push_back(stateIndices_ABM.at(state_ABM_idx)); // 0,6, ..
          std::cout << "ABM index excluded from coupling " << stateIndices_ABM.at(state_ABM_idx) << std::endl;
      }
  }

  if (no_coupling_abm) {
    for (int run = initial_run; run < maxRun; run++) {
      std::cout << "dayData size before " << all_abms[run].eventData.at(nr_day).size() << std::endl;
      delete_non_relevant_events(no_coupling_abm_indices, all_abms[run]);
      std::cout << "dayData size after  " << all_abms[run].eventData.at(nr_day).size() << std::endl;
    }
    std::cout << "done deleting non relevant events" << std::endl;
  }

  std::vector<std::vector<double>> V_PDE;
  std::vector<std::vector<std::vector<double>>> grad_V_PDE;
  for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) { // iterate over all PDE domains
    std::string filename_V = "../../work/input/global/" + stateLabels_PDE.at(state_PDE_idx) + "_V" + coarsity_str_PDEs.at(state_PDE_idx) + ".bin";
    V_PDE.push_back(readLandscape(filename_V));

    std::string filename_grad_V = "../../work/input/global/" + stateLabels_PDE.at(state_PDE_idx) + "_grad_V" + coarsity_str_PDEs.at(state_PDE_idx) + ".bin";
    grad_V_PDE.push_back(readGradient(filename_grad_V));
  }

  std::vector<std::vector<std::vector<double>>> diseaseStatus_all;
  for (int state_idx=0; state_idx<nbr_domains; state_idx++) {
    std::string filename = stateLabels[state_idx];
    std::vector<std::vector<double>> diseaseStatus = readTargetData(filename);
    diseaseStatus_all.push_back(diseaseStatus);
  }
  std::cout << "done reading target data" << std::endl;
  
  // set values: susceptibles, exposed, recovered
  int nbrDays = diseaseStatus_all[0].size();
  int shift_days = 4; // shift from 27th February 2020 to 2nd March 2020
  nbrDays -= shift_days;
  susceptibles.resize(maxRun-initial_run, std::vector<std::vector<std::vector<double>>>(nbrDays, std::vector<std::vector<double>>(nbr_model_types)));
  exposed.resize(maxRun-initial_run, std::vector<std::vector<std::vector<double>>>(nbrDays, std::vector<std::vector<double>>(nbr_model_types)));
  recovered.resize(maxRun-initial_run, std::vector<std::vector<std::vector<double>>>(nbrDays, std::vector<std::vector<double>>(nbr_model_types)));
  for (int run=0; run<maxRun-initial_run; run++) {
    for (int idx_i=0; idx_i<nbrDays; idx_i++) { 
      for (int model_type=0; model_type<nbr_model_types; model_type++) {
        susceptibles[run][idx_i][model_type].resize(stateLabels_all_models[model_type].size());
        exposed[run][idx_i][model_type].resize(stateLabels_all_models[model_type].size());
        recovered[run][idx_i][model_type].resize(stateLabels_all_models[model_type].size());
        for (int state_idx=0; state_idx<stateLabels_all_models[model_type].size(); state_idx++) {
          exposed[run][idx_i][model_type][state_idx] = diseaseStatus_all[stateIndices_all_models[model_type][state_idx]][shift_days+idx_i][0]/population_scale; // 25% population: we model only a quarter of the population but data is given for the whole population
          recovered[run][idx_i][model_type][state_idx] = diseaseStatus_all[stateIndices_all_models[model_type][state_idx]][shift_days+idx_i][1]/population_scale; 
          susceptibles[run][idx_i][model_type][state_idx] = population[model_type][state_idx] - exposed[run][idx_i][model_type][state_idx] - recovered[run][idx_i][model_type][state_idx];
        }
      }
    }
  }

  for (int model_type=0; model_type<nbr_model_types; model_type++) {
    for (int state_idx=0; state_idx<stateLabels_all_models[model_type].size(); state_idx++) {
      std::cout << "exposed[0][0][" << stateLabels_all_models[model_type][state_idx] << "]";
      std::cout << exposed[0][0][model_type][state_idx] << std::endl;
      std::cout << "recovered[0][0][" << stateLabels_all_models[model_type][state_idx] << "]" << recovered[0][0][model_type][state_idx] << std::endl;
      std::cout << "susceptibles[0][0][" << stateLabels_all_models[model_type][state_idx] << "]" << susceptibles[0][0][model_type][state_idx] << "\n" << std::endl;
    }
  }
  //////////
  using namespace boost::fusion;

  int verbosityOpt = 1; 
  bool dump = false; 
  std::unique_ptr<boost::property_tree::ptree> pt = getKaskadeOptions(argc, argv, verbosityOpt, dump);

  std::cout << "Start program" << std::endl;

  boost::timer::cpu_timer totalTimer;

  //0: less output, no file saving; 1: more output, file savings without agents position and health states, 2: more output, all file savings
  int verbosity = 0; //1; //2; //

  int const dim = 2;
  DirectType directType = DirectType::MUMPS; 
  int order = 1;
  int extrapolOrder = 0;
  double xi = 1.0; // no scaling  

  // final times sectioned
  int nbr_time_intervals = 1;
  std::vector<int> T_vec(nbr_time_intervals);
  T_vec[0] = 14;

  double correction_term_ABM = 5.0e-2; 
  double correction_term_PDE = 5.0e+1;
  double correction_term_ODE = 5.0e+1; 

  // baseline scenario and Berlin as ABM: correction terms 
  std::vector<int> correction_term_apply_ABM(nr_ABM_states,0); 
  std::vector<int> correction_term_apply_PDE(nr_PDE_states,0);
  std::vector<int> correction_term_apply_ODE(nr_ODE_states,100);

  // // scenario Berlin as ODE: correction terms
  // std::vector<int> correction_term_apply_ABM(nr_ABM_states,0); 
  // std::vector<int> correction_term_apply_PDE(nr_PDE_states,0);
  // std::vector<int> correction_term_apply_ODE = {100,100,100,100,100,0}; 

  if ((int)correction_term_apply_ABM.size() != nr_ABM_states || (int)correction_term_apply_PDE.size() != nr_PDE_states || (int)correction_term_apply_ODE.size() != nr_ODE_states) {
    std::cerr << "Error: size of correction_term_apply != number of states for that model type";
    std::abort();
  }

  int correction_step = 0;
  
  for (int optimization_run=0; optimization_run<1; optimization_run++) { //100
    std::string sub_folder = "../../work/output/output_data_optim_" + std::to_string(correction_step) + "/graph_abm_pde_ode";
    std::filesystem::create_directories(sub_folder);

    // print all correction terms:
    std::cout << "correction terms ABM:" << std::endl;
    for (int state_ABM_idx = 0; state_ABM_idx < nr_ABM_states; state_ABM_idx++) {
      std::cout << correction_term_apply_ABM.at(state_ABM_idx) << ", ";
    }
    std::cout << std::endl;

    std::cout << "correction terms PDE:" << std::endl;
    for (int state_PDE_idx = 0; state_PDE_idx < nr_PDE_states; state_PDE_idx++) {
      std::cout << correction_term_apply_PDE.at(state_PDE_idx) << ", ";
    }
    std::cout << std::endl;

    std::cout << "correction terms ODE:" << std::endl;
    for (int state_ODE_idx = 0; state_ODE_idx < nr_ODE_states; state_ODE_idx++) {
      std::cout << correction_term_apply_ODE.at(state_ODE_idx) << ", ";
    }
    std::cout << std::endl;

    std::vector<double> beta_e_vec_pde(nbr_time_intervals); 
    beta_e_vec_pde[0] = 4.5e+2;

    std::vector<std::vector<double>> beta_e_vec_pdes(nr_PDE_states, std::vector<double>(nbr_time_intervals));
    for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) {
      for (int time_interval=0; time_interval<nbr_time_intervals; time_interval++) {
        beta_e_vec_pdes.at(state_PDE_idx)[time_interval] = beta_e_vec_pde[time_interval] + correction_term_apply_PDE.at(state_PDE_idx) * correction_term_PDE;
      }
    }

    std::vector<std::vector<double>> beta_e_vec_odes(nr_ODE_states, std::vector<double>(nbr_time_intervals));
    for (int state_ODE_idx=0; state_ODE_idx<nr_ODE_states; state_ODE_idx++) {
      for (int time_interval=0; time_interval<nbr_time_intervals; time_interval++) {
        beta_e_vec_odes.at(state_ODE_idx)[time_interval] = beta_e_vec_pde[time_interval] + correction_term_apply_ODE.at(state_ODE_idx) * correction_term_ODE;
      }
    }

    double betaFactor = 0.4;
    std::vector<double> param_abms(nr_ABM_states);
    for (int state_ABM_idx=0; state_ABM_idx<nr_ABM_states; state_ABM_idx++) {
      param_abms.at(state_ABM_idx) = 1.7e-5 * (betaFactor + correction_term_apply_ABM.at(state_ABM_idx) * correction_term_ABM);
    }

    double sigma = 1.0/3.5, gamma = 1.0/2, eta = 1.0/4, kappa = 1.0, eta_c = 1.0/21, phi_i = 1.0/4, phi_sy = 1.0/8, phi_h = 1.0/14, phi_hc = 1.0/7;

    double landscape_scale = 1e-14;
    double D_scale = 1.0e+8;

    std::vector<double> D(nbr_time_intervals);
    D[0] = D_scale;

    // scale V and grad_V
    for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) {
      for (int i=0; i<V_PDE.at(state_PDE_idx).size(); i++) { // grid nodes
        V_PDE.at(state_PDE_idx)[i] *= landscape_scale;
        grad_V_PDE.at(state_PDE_idx)[i][0] *= landscape_scale;
        grad_V_PDE.at(state_PDE_idx)[i][1] *= landscape_scale;
      }
    }

    // create V with only positive entries
    std::vector<std::vector<double>> V_pos_PDE;
    for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) {
      std::vector<double> V_pos(V_PDE.at(state_PDE_idx).size(),0.0);
      double min_V = *std::min_element(V_PDE.at(state_PDE_idx).begin(), V_PDE.at(state_PDE_idx).end());
      for (int idx=0; idx<V_PDE.at(state_PDE_idx).size(); idx++) {
        V_pos[idx] = V_PDE.at(state_PDE_idx)[idx] + min_V + 0.1; // positive entries
      }
      V_pos_PDE.push_back(V_pos);
    }

    using Grid = Dune::UGGrid<dim>; 
    std::vector<std::unique_ptr<GridManager<Grid>>> gridManagers_PDE;
    for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) {
      std::string gridFile = "domain/domain_" + stateLabels_PDE.at(state_PDE_idx) + coarsity_str_PDEs.at(state_PDE_idx) + ".1"; 
      auto gm = std::make_unique<GridManager<Grid>>(readPolyData(gridFile, nullptr, 500, 10)); // read given boundary vertices and edges to create grid
      gm->enforceConcurrentReads(true);

      std::cout << "\nGrid " << stateLabels_PDE.at(state_PDE_idx) << ": "
                << gm->grid().size(0) << " triangles,\n"
                << "      " << gm->grid().size(1) << " edges,\n"
                << "      " << gm->grid().size(dim) << " points\n\n";

      gridManagers_PDE.push_back(std::move(gm));
    }

    // construct involved spaces
    using H1Space = FEFunctionSpace<ContinuousLagrangeMapper<double,Grid::LeafGridView>>;
    using Spaces = boost::fusion::vector<H1Space const*>; 
    std::vector<H1Space> pdeSpaces;
    std::vector<Spaces> spaces_PDE;
    pdeSpaces.reserve(nr_PDE_states);
    spaces_PDE.reserve(nr_PDE_states);
    for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) {
      pdeSpaces.push_back(H1Space(*gridManagers_PDE.at(state_PDE_idx),gridManagers_PDE.at(state_PDE_idx)->grid().leafGridView(),order)); 
      spaces_PDE.push_back(Spaces(&pdeSpaces.at(state_PDE_idx)));
    }

    using VariableDescriptions = boost::fusion::vector<
                                                      Variable<SpaceIndex<0>,Components<1>,VariableId<0>>,
                                                      Variable<SpaceIndex<0>,Components<1>,VariableId<1>>, 
                                                      Variable<SpaceIndex<0>,Components<1>,VariableId<2>>,
                                                      Variable<SpaceIndex<0>,Components<1>,VariableId<3>>,
                                                      Variable<SpaceIndex<0>,Components<1>,VariableId<4>>, 
                                                      Variable<SpaceIndex<0>,Components<1>,VariableId<5>>,
                                                      Variable<SpaceIndex<0>,Components<1>,VariableId<6>>,
                                                      Variable<SpaceIndex<0>,Components<1>,VariableId<7>>
                                                      >;
    std::string varNames[nbr_compartments] = { "s", "e", "i", "sy", "h", "c", "hc", "r" };

    using VarSetDesc = VariableSetDescription<Spaces,VariableDescriptions>;
    std::vector<VarSetDesc> varSetDesc_PDE;
    std::vector<VarSetDesc::VariableSet> x_PDE;
    varSetDesc_PDE.reserve(nr_PDE_states);
    x_PDE.reserve(nr_PDE_states);
    for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) {
      varSetDesc_PDE.push_back(VarSetDesc(spaces_PDE.at(state_PDE_idx),varNames));
    }
    for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) {
      x_PDE.push_back(VarSetDesc::VariableSet(varSetDesc_PDE.at(state_PDE_idx)));
    }

    using Equation = AlievPanfilovEquation<double,VarSetDesc>;

    std::vector<int> nbr_points_PDE, nbr_cells_PDE;
    for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) {
      Kaskade::VTKWriterDetail::VTKGridInfo<typename VarSetDesc::VariableSet::Descriptions::GridView> gridInfo_PDE(x_PDE.at(state_PDE_idx).descriptions.gridView,IoOptions().setPrecision(10).setDataMode(IoOptions::conforming));
      nbr_points_PDE.push_back(gridInfo_PDE.npoints);
      nbr_cells_PDE.push_back(gridInfo_PDE.ncells);
      std::cout << "nbr_points_PDE " << stateLabels_PDE.at(state_PDE_idx) << " " << nbr_points_PDE.at(state_PDE_idx) << std::endl;
      std::cout << "nbr_cells_PDE " << stateLabels_PDE.at(state_PDE_idx) << " " << nbr_cells_PDE.at(state_PDE_idx) << "\n" << std::endl;
    }

    for (int i = 0; i < nr_PDE_states; ++i) {
      if ((int)V_PDE[i].size() != nbr_points_PDE[i]) {
        std::cerr << "Mismatch V: " << stateLabels_PDE[i]
                  << " V=" << V_PDE[i].size()
                  << " npoints=" << nbr_points_PDE[i] << "\n";
        std::abort();
      }
      if ((int)grad_V_PDE[i].size() != nbr_points_PDE[i]) {
        std::cerr << "Mismatch grad_V: " << stateLabels_PDE[i]
                  << " gradV=" << grad_V_PDE[i].size()
                  << " npoints=" << nbr_points_PDE[i] << "\n";
        std::abort();
      }
    }


    std::vector<std::vector<double>> corners;
    corners.resize(3, std::vector<double>(2, 0.0));
    double area;

    // compute density equivalent of one individual (minDensityForTransition_PDE)
    std::vector<std::vector<double>> area_of_triangles_touching_grid_point(nr_PDE_states);
    for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) {
      int nbr_grid_points = nbr_points_PDE.at(state_PDE_idx);
      area_of_triangles_touching_grid_point.at(state_PDE_idx).resize(nbr_grid_points,0.0);
      for (auto const& cell: elements(x_PDE.at(state_PDE_idx).descriptions.gridView)) {
        for (int corner=0; corner<cell.subEntities(2); corner++) {
          auto xglob = cell.geometry().corner(corner);
          corners[corner][0] = xglob[0];
          corners[corner][1] = xglob[1];
        }
        area = 0.5 * std::abs(corners[0][0] * (corners[1][1]-corners[2][1]) + corners[1][0] * (corners[2][1]-corners[0][1]) + corners[2][0] * (corners[0][1]-corners[1][1]) );
        for (int corner=0; corner<cell.subEntities(2); corner++) {
          int glob_idx = (x_PDE.at(state_PDE_idx).descriptions.gridView).indexSet().subIndex(cell,corner,2);
          area_of_triangles_touching_grid_point.at(state_PDE_idx)[glob_idx] += area;
        }
      }
      minDensityForTransition_PDE.at(state_PDE_idx).resize(nbr_grid_points);
      for (int i=0; i<nbr_grid_points; i++) {
        double inverse_area = 1.0/area_of_triangles_touching_grid_point.at(state_PDE_idx)[i];
        minDensityForTransition_PDE.at(state_PDE_idx)[i] = 3.0 * inverse_area;
      }
    }

    // compute distribution for initial values: normalized inverse of positive landscape
    double integral_of_inv_V;
    for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) {
      integral_of_inv_V = 0.0;
      for (auto const& cell: elements(x_PDE.at(state_PDE_idx).descriptions.gridView)) {
        for (int corner=0; corner<cell.subEntities(2); corner++) {
          auto xglob = cell.geometry().corner(corner);
          corners[corner][0] = xglob[0];
          corners[corner][1] = xglob[1];
        }
        area = 0.5 * std::abs(corners[0][0] * (corners[1][1]-corners[2][1]) + corners[1][0] * (corners[2][1]-corners[0][1]) + corners[2][0] * (corners[0][1]-corners[1][1]) );
        for (int corner=0; corner<cell.subEntities(2); corner++) {
          int index = x_PDE.at(state_PDE_idx).descriptions.gridView.indexSet().subIndex(cell, corner, 2); //global index
          integral_of_inv_V += (area/3.0) * 1.0/(V_pos_PDE.at(state_PDE_idx)[index]);
        }
      }
      std::cout << "integral_of_inv_V " << integral_of_inv_V << std::endl;

      // normalize
      normalized_inv_V_PDE.at(state_PDE_idx).resize(V_pos_PDE.at(state_PDE_idx).size());
      for (int idx=0; idx<V_pos_PDE.at(state_PDE_idx).size(); idx++) {
        normalized_inv_V_PDE.at(state_PDE_idx)[idx] = 1.0/(V_pos_PDE.at(state_PDE_idx)[idx] * integral_of_inv_V); 
      }
    }

    std::cout << "start runs" << std::endl;
    std::vector<std::vector<std::vector<double>>> x_ODE_all_runs(maxRun-initial_run, std::vector<std::vector<double>>(nr_ODE_states, std::vector<double>(nbr_compartments,0.0)));
    std::vector<std::vector<double>> S_ode_t_0(maxRun-initial_run);
    std::vector<std::vector<double>> E_ode_t_0(maxRun-initial_run);
    std::vector<std::vector<double>> R_ode_t_0(maxRun-initial_run);
    for (int run=0; run<maxRun-initial_run; run++) {
      for (int state_idx=0; state_idx<nr_ODE_states; state_idx++) {
        E_ode_t_0[run].push_back(exposed[run][0][ode_idx][state_idx]);
        R_ode_t_0[run].push_back(recovered[run][0][ode_idx][state_idx]);
        S_ode_t_0[run].push_back(susceptibles[run][0][ode_idx][state_idx]);
      }
    }

    std::vector<std::vector<VarSetDesc::VariableSet>> x_PDE_all_runs;
    for (int run=0; run<maxRun-initial_run; run++) { 
      x_PDE_all_runs.push_back(x_PDE);
    }

    std::vector<std::vector<std::unordered_map<int,int>>> agentID_healthStatus_all_runs(maxRun - initial_run);
    std::vector<std::vector<int>> agent_IDs_in_ABM_all_runs(maxRun - initial_run); 

    std::string filename_duration = "../../work/output/output_data_optim_" + std::to_string(correction_step) + "/duration.txt";
    {
        std::ofstream file_duration(filename_duration, std::ios::trunc);
    }

    double t_0, T;
    T = T_vec[0]; 
    t_0 = 0;
    
    int maxSteps = getParameter(pt, "maxSteps", (T-t_0)*dt_inv); 
    std::cout << "maxSteps = " << maxSteps << std::endl; 

    #pragma omp parallel for schedule(static,1) proc_bind(spread)
    for (int run=initial_run; run<maxRun; run++) {
      auto start_time = std::chrono::high_resolution_clock::now();
      std::cout << "########################## run: " << run << std::endl;      
      int week_day_local = initial_week_day;

      unsigned int run_seed = 12345 + run; // unique seed per run
      std::ranlux24_base gen_local(run_seed);

      std::vector<Equation> equations;
      for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) {
        Equation Eq(sigma, gamma, eta, kappa, eta_c, phi_i, phi_sy, phi_h, phi_hc, beta_e_vec_pdes.at(state_PDE_idx), all_abms[run-initial_run].activityChangePercentage, D, landscape_scale, grad_V_PDE.at(state_PDE_idx), T_vec, xi, varSetDesc_PDE.at(state_PDE_idx));
        equations.push_back(Eq); 
      }

      std::vector<VarSetDesc::VariableSet> x_initial_PDE;
      for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) {
        x_initial_PDE.push_back(VarSetDesc::VariableSet(varSetDesc_PDE.at(state_PDE_idx)));
      }

      std::cout << "ABM parameter choices: " << std::endl;
      for (int state_ABM_idx=0; state_ABM_idx<nr_ABM_states; state_ABM_idx++) {
        std::cout << std::scientific << std::setprecision(10) << param_abms.at(state_ABM_idx) << std::endl;
      }

      std::cout << "PDE parameter choices: " << std::endl;
      for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) {
        std::cout << std::scientific << std::setprecision(10) << beta_e_vec_pdes.at(state_PDE_idx)[0] << std::endl;
      }

      std::cout << "ODE parameter choices: " << std::endl;
      for (int state_ODE_idx=0; state_ODE_idx<nr_ODE_states; state_ODE_idx++) {
        std::cout << std::scientific << std::setprecision(10) << beta_e_vec_odes.at(state_ODE_idx)[0] << std::endl;
      }

      for (int state_ODE_idx=0; state_ODE_idx<nr_ODE_states; state_ODE_idx++) { // iterate over all ODE domains
        x_ODE_all_runs[run-initial_run].at(state_ODE_idx) = {
          S_ode_t_0[run-initial_run].at(state_ODE_idx), 
          E_ode_t_0[run-initial_run].at(state_ODE_idx), 
          0.0,
          0.0,
          0.0,
          0.0,
          0.0,
          R_ode_t_0[run-initial_run].at(state_ODE_idx)
        };
      }

      for (int state_PDE_idx=0; state_PDE_idx<nr_PDE_states; state_PDE_idx++) {
        equations.at(state_PDE_idx).setTime(t_0); 
        // set initial values
        equations.at(state_PDE_idx).scaleInitialValue<0>(Initial_1Value(0, x_PDE_all_runs[run-initial_run].at(state_PDE_idx), normalized_inv_V_PDE.at(state_PDE_idx), state_PDE_idx, run-initial_run), x_PDE_all_runs[run-initial_run].at(state_PDE_idx));
        equations.at(state_PDE_idx).scaleInitialValue<1>(Initial_1Value(1, x_PDE_all_runs[run-initial_run].at(state_PDE_idx), normalized_inv_V_PDE.at(state_PDE_idx), state_PDE_idx, run-initial_run), x_PDE_all_runs[run-initial_run].at(state_PDE_idx));
        equations.at(state_PDE_idx).scaleInitialValue<2>(Initial_1Value(2, x_PDE_all_runs[run-initial_run].at(state_PDE_idx), normalized_inv_V_PDE.at(state_PDE_idx), state_PDE_idx, run-initial_run), x_PDE_all_runs[run-initial_run].at(state_PDE_idx));
        equations.at(state_PDE_idx).scaleInitialValue<3>(Initial_1Value(3, x_PDE_all_runs[run-initial_run].at(state_PDE_idx), normalized_inv_V_PDE.at(state_PDE_idx), state_PDE_idx, run-initial_run), x_PDE_all_runs[run-initial_run].at(state_PDE_idx));
        equations.at(state_PDE_idx).scaleInitialValue<4>(Initial_1Value(4, x_PDE_all_runs[run-initial_run].at(state_PDE_idx), normalized_inv_V_PDE.at(state_PDE_idx), state_PDE_idx, run-initial_run), x_PDE_all_runs[run-initial_run].at(state_PDE_idx)); 
        equations.at(state_PDE_idx).scaleInitialValue<5>(Initial_1Value(5, x_PDE_all_runs[run-initial_run].at(state_PDE_idx), normalized_inv_V_PDE.at(state_PDE_idx), state_PDE_idx, run-initial_run), x_PDE_all_runs[run-initial_run].at(state_PDE_idx));
        equations.at(state_PDE_idx).scaleInitialValue<6>(Initial_1Value(6, x_PDE_all_runs[run-initial_run].at(state_PDE_idx), normalized_inv_V_PDE.at(state_PDE_idx), state_PDE_idx, run-initial_run), x_PDE_all_runs[run-initial_run].at(state_PDE_idx));
        equations.at(state_PDE_idx).scaleInitialValue<7>(Initial_1Value(7, x_PDE_all_runs[run-initial_run].at(state_PDE_idx), normalized_inv_V_PDE.at(state_PDE_idx), state_PDE_idx, run-initial_run), x_PDE_all_runs[run-initial_run].at(state_PDE_idx));
      }

      std::vector<double> params_ODEs = { sigma, gamma, eta, kappa, eta_c, phi_i, phi_sy, phi_h, phi_hc };
      for (int i=0; i<nr_ODE_states; i++) {
        params_ODEs.push_back(beta_e_vec_odes.at(i)[0]);
      }
      
      std::tie(agentID_healthStatus_all_runs[run-initial_run], agent_IDs_in_ABM_all_runs[run-initial_run]) = 
      integrate( 
        gridManagers_PDE, equations, varSetDesc_PDE, 
        dt_inv, population_scale, gen_local, correction_step,
        week_day_local, T, maxSteps, extrapolOrder,
        x_PDE_all_runs[run-initial_run], x_ODE_all_runs[run-initial_run], all_abms[run-initial_run], param_abms, params_ODEs,
        run, initial_run, with_coupling_states,
        directType, verbosity
      );

      auto end_time = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
      std::cout << "total duration: " << duration.count() << " milliseconds" << std::endl;
      std::cout << duration.count()/1000.0 << " seconds" << std::endl;
      std::cout << duration.count()/(1000.0*60) << " minutes" << std::endl;

      #pragma omp critical
      {
        // save time to file in milliseconds, one line for each run
        std::ofstream file_duration(filename_duration, std::ios::app);
        file_duration << duration.count() << '\n';
      }
    } // end for run
    if (verbosity == 0) {
      return 0;
    }
    //////////////////////////////////////////////////////
    int days = maxSteps / dt_inv;
    double avg_error_germany = 0.0;

    std::vector<std::vector<double>> avg_result_days_ABM(nr_ABM_states,std::vector<double>(days, 0.0));
    for (int run=initial_run; run<maxRun; run++) {
      for (int state_ABM_idx = 0; state_ABM_idx < nr_ABM_states; state_ABM_idx++) {
        std::string filename_result_ABM = "../../work/output/output_data_optim_" + std::to_string(correction_step) + "/result_ABM_0_" + std::to_string(run) + "_" + stateLabels_ABM.at(state_ABM_idx) + addOn + ".txt";
        std::ifstream file(filename_result_ABM);
        if (!file) {
          throw std::runtime_error("Missing file: " + filename_result_ABM);
        }
        std::vector<double> values;
        double val;
        while (file >> val) {
          values.push_back(val);
        }
        if ((int)values.size() < days) {
          throw std::runtime_error("Too few values in " + filename_result_ABM +
                                  " got " + std::to_string(values.size()) +
                                  " need " + std::to_string(days));
        }
        for (int day = 0; day < days; day++) {
          avg_result_days_ABM.at(state_ABM_idx)[day] += values[day];
        }
      }
    }
    // divide to get mean
    for (int state_ABM_idx = 0; state_ABM_idx < nr_ABM_states; state_ABM_idx++) {
      for (int day = 0; day < days; day++) {
        avg_result_days_ABM.at(state_ABM_idx)[day] /= (maxRun-initial_run);
      }
    }
    // read target data 
    std::vector<std::vector<double>> target_ABM_loaded(nr_ABM_states);
    for (int state_ABM_idx = 0; state_ABM_idx < nr_ABM_states; state_ABM_idx++) {
      std::string filename_target_ABM = "../../work/output/output_data_optim_" + std::to_string(correction_step) + "/target_ABM_0_" + stateLabels_ABM.at(state_ABM_idx) + ".txt";

      std::ifstream file(filename_target_ABM);
      double val;
      while (file >> val) {
        target_ABM_loaded.at(state_ABM_idx).push_back(val);
      }
    }

    // compute error
    std::vector<double> error_ABM(nr_ABM_states, 0.0);
    for (int state_ABM_idx = 0; state_ABM_idx < nr_ABM_states; state_ABM_idx++) {
      for (int day = 0; day < days; day++) {
        double diff = (avg_result_days_ABM.at(state_ABM_idx)[day] - target_ABM_loaded.at(state_ABM_idx)[day]);
        error_ABM.at(state_ABM_idx) += diff;
        avg_error_germany += std::abs(diff);
      }
      if (error_ABM.at(state_ABM_idx) > 0.0 && ((betaFactor + (correction_term_apply_ABM.at(state_ABM_idx)-1) * correction_term_ABM) > 0.0)) {
        correction_term_apply_ABM.at(state_ABM_idx) -= 1;
      }
      else if (error_ABM.at(state_ABM_idx) < 0.0) {
        correction_term_apply_ABM.at(state_ABM_idx) += 1;
      }
    }

    std::vector<std::vector<std::vector<double>>> result_days_PDE(maxRun-initial_run, std::vector<std::vector<double>>(nr_PDE_states)); 
    std::vector<std::vector<double>> avg_result_days_PDE(nr_PDE_states,std::vector<double>(days, 0.0));
    for (int run=initial_run; run<maxRun; run++) {
      for (int state_PDE_idx = 0; state_PDE_idx < nr_PDE_states; state_PDE_idx++) {
        std::string filename_result_PDE = "../../work/output/output_data_optim_" + std::to_string(correction_step) + "/result_PDE_0_" + std::to_string(run) + "_" + stateLabels_PDE.at(state_PDE_idx) + addOn + ".txt";
        std::ifstream file(filename_result_PDE);
        if (!file) {
          throw std::runtime_error("Missing file: " + filename_result_PDE);
        }
        std::vector<double> values; 
        double val;
        while (file >> val) {
          values.push_back(val);
        }
        if ((int)values.size() < days) {
          throw std::runtime_error("Too few values in " + filename_result_PDE +
                                  " got " + std::to_string(values.size()) +
                                  " need " + std::to_string(days));
        }
        for (int day = 0; day < days; day++) {
          avg_result_days_PDE.at(state_PDE_idx)[day] += values[day];
        }
      }
    }
    // divide to get mean
    for (int state_PDE_idx = 0; state_PDE_idx < nr_PDE_states; state_PDE_idx++) {
      for (int day = 0; day < days; day++) {
        avg_result_days_PDE.at(state_PDE_idx)[day] /= (maxRun-initial_run);
      }
    }
    // read target data 
    std::vector<std::vector<double>> target_PDE_loaded(nr_PDE_states);
    for (int state_PDE_idx = 0; state_PDE_idx < nr_PDE_states; state_PDE_idx++) {
      std::string filename_target_PDE = "../../work/output/output_data_optim_" + std::to_string(correction_step) + "/target_PDE_0_" + stateLabels_PDE.at(state_PDE_idx) + ".txt";

      std::ifstream file(filename_target_PDE);
      double val;
      while (file >> val) {
        target_PDE_loaded.at(state_PDE_idx).push_back(val);
      }
    }

    // compute error
    std::vector<double> error_PDE(nr_PDE_states, 0.0);
    for (int state_PDE_idx = 0; state_PDE_idx < nr_PDE_states; state_PDE_idx++) {
      for (int day = 0; day < days; day++) {
        double diff = (avg_result_days_PDE.at(state_PDE_idx)[day] - target_PDE_loaded.at(state_PDE_idx)[day]);
        error_PDE.at(state_PDE_idx) += diff;
        avg_error_germany += std::abs(diff);
      }
      if (error_PDE.at(state_PDE_idx) > 0.0 && (beta_e_vec_pde[0] + (correction_term_apply_PDE.at(state_PDE_idx)-1) * correction_term_PDE) > 0.0) {
        correction_term_apply_PDE.at(state_PDE_idx) -= 1;
      }
      else if (error_PDE.at(state_PDE_idx) < 0.0) {
        correction_term_apply_PDE.at(state_PDE_idx) += 1;
      }
    }

    std::vector<std::vector<std::vector<double>>> result_days_ODE(maxRun-initial_run, std::vector<std::vector<double>>(nr_ODE_states)); 
    std::vector<std::vector<double>> avg_result_days_ODE(nr_ODE_states,std::vector<double>(days, 0.0));
    for (int run=initial_run; run<maxRun; run++) {
      for (int state_ODE_idx = 0; state_ODE_idx < nr_ODE_states; state_ODE_idx++) {
        std::string filename_result_ODE = "../../work/output/output_data_optim_" + std::to_string(correction_step) + "/result_ODE_0_" + std::to_string(run) + "_" + stateLabels_ODE.at(state_ODE_idx) + addOn + ".txt";
        std::ifstream file(filename_result_ODE);

        if (!file) {
          throw std::runtime_error("Missing file: " + filename_result_ODE);
        }
        std::vector<double> values;
        double val;
        while (file >> val) {
          values.push_back(val);
        }
        if ((int)values.size() < days) {
          throw std::runtime_error("Too few values in " + filename_result_ODE +
                                  " got " + std::to_string(values.size()) +
                                  " need " + std::to_string(days));
        }
        for (int day = 0; day < days; day++) {
          avg_result_days_ODE.at(state_ODE_idx)[day] += values[day];
        }
      }
    }
    // divide to get mean
    for (int state_ODE_idx = 0; state_ODE_idx < nr_ODE_states; state_ODE_idx++) {
      for (int day = 0; day < days; day++) {
        avg_result_days_ODE.at(state_ODE_idx)[day] /= (maxRun-initial_run);
      }
    }
    // read target data 
    std::vector<std::vector<double>> target_ODE_loaded(nr_ODE_states);
    for (int state_ODE_idx = 0; state_ODE_idx < nr_ODE_states; state_ODE_idx++) {
      std::string filename_target_ODE = "../../work/output/output_data_optim_" + std::to_string(correction_step) + "/target_ODE_0_" + stateLabels_ODE.at(state_ODE_idx) + ".txt";

      std::ifstream file(filename_target_ODE);
      double val;
      while (file >> val) {
        target_ODE_loaded.at(state_ODE_idx).push_back(val);
      }
    }

    // compute error
    std::vector<double> error_ODE(nr_ODE_states, 0.0);
    for (int state_ODE_idx = 0; state_ODE_idx < nr_ODE_states; state_ODE_idx++) {
      for (int day = 0; day < days; day++) {
        double diff = (avg_result_days_ODE.at(state_ODE_idx)[day] - target_ODE_loaded.at(state_ODE_idx)[day]);
        error_ODE.at(state_ODE_idx) += diff;
        avg_error_germany += std::abs(diff);
      }
      if (error_ODE.at(state_ODE_idx) > 0.0 && (beta_e_vec_pde[0] + (correction_term_apply_ODE.at(state_ODE_idx)-1) * correction_term_ODE > 0.0)) {
        correction_term_apply_ODE.at(state_ODE_idx) -= 1;
      }
      else if (error_ODE.at(state_ODE_idx) < 0.0) {
        correction_term_apply_ODE.at(state_ODE_idx) += 1;
      }
    }
    avg_error_germany /= (nbr_domains*days);
    std::cout << "average error in Germany for correction step " << std::to_string(correction_step) << ": " << avg_error_germany << std::endl;

    correction_step++;
  } //end for optimization_run

  std::cout << "End program" << std::endl;
}