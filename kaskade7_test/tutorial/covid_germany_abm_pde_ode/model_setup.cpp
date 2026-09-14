#include <iostream>
#include <fstream>
#include <tuple>
#include <math.h>
#include <string> 
#include <vector> 
#include <limits>

#include <dune/common/fvector.hh>

//only important if non-ABM state becomes ABM state or vice versa:
// empty "" string means default scenario, "_Berlin_ABM" means Berlin is ABM state 
// string is needed for jump matrix in ABM.cpp
const std::string addOn = ""; //"_Berlin_ABM"; //

int initial_nr_day;
int initial_week_day;
const int nbr_day_types = 3;

const int nbr_compartments = 8;

const int nbr_model_types = 3; // PDE, ABM and ODE
const int abm_idx = 0; 
const int pde_idx = 1; 
const int ode_idx = 2; 

// all federal states in correct order: 
const std::vector<std::string> stateLabels = {
    "Schleswig-Holstein", "Hamburg", "Niedersachsen", "Bremen", "Nordrhein-Westfalen",
    "Hessen", "Rheinland-Pfalz", "Baden-Wuerttemberg", "Bayern", "Saarland",
    "Berlin", "Brandenburg", "Mecklenburg-Vorpommern", "Sachsen",
    "Sachsen-Anhalt", "Thueringen"
}; // data is given for Berlin and then Brandenburg

std::vector<int> stateIndices_ABM, stateIndices_PDE, stateIndices_ODE; 
std::vector<std::vector<int>> stateIndices_all_models = {stateIndices_ABM, stateIndices_PDE, stateIndices_ODE}; 

int nr_ABM_states, nr_PDE_states, nr_ODE_states; 

std::vector<std::string> stateLabels_ABM, stateLabels_PDE, stateLabels_ODE; 
std::vector<std::vector<std::string>> stateLabels_all_models = {stateLabels_ABM, stateLabels_PDE, stateLabels_ODE}; 

static int nbr_domains; 
static int const nbr_federal_states = 16;

std::vector<int> model_type_of_domain(nbr_federal_states, -1);
std::vector<int> local_index(nbr_federal_states, -1);

const std::vector<double> areas_ODEs = {
    1.5804e+10, 7.551e+8, 4.7710e+10, 4.199e+8, 3.4113e+10,
    2.1116e+10, 1.9858e+10, 3.5748e+10, 7.0542e+10, 2.572e+9,
    8.911e+8, 2.9654e+10, 2.3295e+10, 1.8450e+10,
    2.0555e+10, 1.6202e+10
};
std::vector<double> areas_ODE;

std::vector<std::vector<double>> minDensityForTransition_PDE; 

std::vector<std::vector<std::vector<std::vector<double>>>> susceptibles, exposed, recovered;
std::vector<std::vector<double>> population(nbr_model_types); 
std::vector<double> population_pde_t_0, population_ode_t_0;
int total_initial_population;
std::vector<std::vector<double>> normalized_inv_V_PDE; 

std::vector<double> readLandscape(std::string filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Fehler beim Öffnen der Datei " << filename << std::endl;
        return std::vector<double>();
    }

    // get filesize
    file.seekg(0, std::ios::end);
    std::streampos fileSize = file.tellg();
    file.seekg(0, std::ios::beg);

    std::size_t numDoubles = fileSize / sizeof(double);
    std::vector<double> data(numDoubles);
    file.read(reinterpret_cast<char*>(data.data()), fileSize);
    file.close();

    return data;
}

std::vector<std::vector<double>> readGradient(std::string filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Fehler beim Öffnen der Datei " << filename << std::endl;
        return std::vector<std::vector<double>>();
    }

    // get filesize
    file.seekg(0, std::ios::end);
    std::streampos fileSize = file.tellg();
    file.seekg(0, std::ios::beg);

    std::size_t numPoints = fileSize / (2 * sizeof(double));
    std::cout << "grad_V numPoints=" << numPoints << std::endl;
    std::vector<std::vector<double>> data(numPoints, std::vector<double>(2));
    for (std::size_t i = 0; i < numPoints; i++) {
        file.read(reinterpret_cast<char*>(&data[i][0]), sizeof(double));
        file.read(reinterpret_cast<char*>(&data[i][1]), sizeof(double));
    }
    file.close();

    return data;
}

std::vector<std::vector<double>> readTargetData(std::string bundesland){
    std::string filename = "../../work/input_data/global/cases_" + bundesland + ".bin";
    std::ifstream infile(filename, std::ios::binary);

    if (!infile) {
        std::cerr << "Datei " << filename << " konnte nicht geöffnet werden" << std::endl;
        return std::vector<std::vector<double>>();
    }
    // get filesize
    infile.seekg(0, std::ios::end);
    std::streampos fileSize = infile.tellg();
    infile.seekg(0, std::ios::beg);

    std::size_t numEventPoints = fileSize / (2 * sizeof(double));
    std::vector<std::vector<double>> data(numEventPoints, std::vector<double>(2)); 

    for (int i = 0; i < numEventPoints; i++) {
        double value1, value2;
        infile.read(reinterpret_cast<char*>(&value1), sizeof(value1));
        infile.read(reinterpret_cast<char*>(&value2), sizeof(value2));

        if (!infile) {
            std::cerr << "Fehler beim Lesen." << std::endl;
            break;
        }
        data[i][0] = value1;
        data[i][1] = value2;
    }
    infile.close();
    return data;
}

void plot_result(std::vector<double> result, std::vector<double> target, int correction_step, std::string fileAddOn="") {
    int numTimeSteps = (int) result.size();

    std::string colors[] = {"orange", "black"}; 

    FILE* gnuplotPipe = popen("gnuplot", "w");

    fprintf(gnuplotPipe, "set xlabel 'time (days)'\n");
    fprintf(gnuplotPipe, "set ylabel 'number of infectious'\n");
    fprintf(gnuplotPipe, "set logscale y\n"); // Y-Achse auf logarithmisch setzen

    std::string dataFileName = "data_all_provinces.dat";
    std::ofstream dataFile(dataFileName);

    for (int timeStep = 0; timeStep < numTimeSteps; timeStep++) {
        dataFile << timeStep;
        dataFile << " " << result[timeStep] << " " << target[timeStep];
        dataFile << std::endl;
    }
    dataFile.close();

    fprintf(gnuplotPipe, "set terminal png\n");
    fprintf(gnuplotPipe, "set output '../../work/output/output_data_optim_%d/Symptomatic_hybrid_interval_0_%s.png'\n", correction_step, fileAddOn.c_str());
    fprintf(gnuplotPipe, "set title ''\n");
    fprintf(gnuplotPipe, "plot ");

    fprintf(gnuplotPipe, "'%s' using 1:(column(%d*2)) with lines title sprintf('numerical, symptomatic', %d-1) lc rgb '%s', ",
            dataFileName.c_str(), 1, 1, colors[0].c_str()); 
    fprintf(gnuplotPipe, "'%s' using 1:(column(%d*2+1)) with linespoints title sprintf('target', %d-1) lc rgb '%s', ",
            dataFileName.c_str(), 1, 1, colors[1].c_str());

    pclose(gnuplotPipe);
}
