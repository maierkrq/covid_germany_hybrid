#include <iostream>
#include <fstream>
#include <tuple>
#include <math.h>
#include <string> 
#include <vector> 
#include <limits>

#include <dune/common/fvector.hh>

//only important if ABM states change:
const std::string addOn = ""; //"_Berlin_ABM"; //

int initial_nr_day;
int initial_week_day;

const int nbr_compartments = 8;

const int nbr_model_types = 3; // PDE, ABM and ODE
const int pde_idx = 0;
const int abm_idx = 1; 
const int ode_idx = 2; 

const std::vector<int> stateIndices_PDE = {1,3,9,10}; // Berlin is number 10, can be moved to stateIndices_ABM or stateIndices_ODE
const std::vector<int> stateIndices_ABM = {0,6,11,12,13,14,15}; 
const std::vector<int> stateIndices_ODE = {2,4,5,7,8}; 
const std::vector<std::vector<int>> stateIndices_all_models = {stateIndices_PDE, stateIndices_ABM, stateIndices_ODE}; 

const int nr_PDE_states = stateIndices_PDE.size();
const int nr_ABM_states = stateIndices_ABM.size();
const int nr_ODE_states = stateIndices_ODE.size();

const std::vector<std::string> stateLabels = {"Schleswig-Holstein", "Hamburg", "Niedersachsen", "Bremen", "Nordrhein-Westfalen",
    "Hessen", "Rheinland-Pfalz", "Baden-Wuerttemberg", "Bayern", "Saarland",
    "Berlin", "Brandenburg", "Mecklenburg-Vorpommern", "Sachsen",
    "Sachsen-Anhalt", "Thueringen"
}; // data is given for Berlin and then Brandenburg


const std::vector<std::string> stateLabels_PDE = {"Hamburg","Bremen","Saarland","Berlin"}; // Berlin can be moved to stateLabels_ABM or stateLabels_ODE
const std::vector<std::string> stateLabels_ABM = {"Schleswig-Holstein","Rheinland-Pfalz","Brandenburg","Mecklenburg-Vorpommern","Sachsen","Sachsen-Anhalt","Thueringen"};
const std::vector<std::string> stateLabels_ODE = {"Niedersachsen","Nordrhein-Westfalen","Hessen","Baden-Wuerttemberg","Bayern"};
const std::vector<std::vector<std::string>> stateLabels_all_models = {stateLabels_PDE, stateLabels_ABM, stateLabels_ODE}; 

static int const nbr_domains = stateLabels.size();

std::vector<int> model_type_of_domain(nbr_domains, -1);
std::vector<int> local_index(nbr_domains, -1);

const std::vector<double> areas_ODE = {4.7710e+10, 3.4113e+10, 2.1116e+10, 3.5748e+10, 7.0542e+10 };
// const std::vector<double> areas_ODE = {4.7710e+10, 3.4113e+10, 2.1116e+10, 3.5748e+10, 7.0542e+10, 8.911e+8 }; // if Berlin is in ODE, then we need to use this areas_ODE

std::vector<std::vector<double>> minDensityForTransition_PDE(nr_PDE_states);

std::vector<std::vector<std::vector<std::vector<double>>>> susceptibles;
std::vector<std::vector<std::vector<std::vector<double>>>> exposed;
std::vector<std::vector<std::vector<std::vector<double>>>> recovered;
std::vector<std::vector<double>> population(nbr_model_types); 
std::vector<double> population_pde_t_0(nr_PDE_states); 
std::vector<double> population_ode_t_0(nr_ODE_states);
std::vector<std::vector<double>> normalized_inv_V_PDE(nr_PDE_states);

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
    for (std::size_t i = 0; i < numPoints; ++i) {
        file.read(reinterpret_cast<char*>(&data[i][0]), sizeof(double));
        file.read(reinterpret_cast<char*>(&data[i][1]), sizeof(double));
    }
    file.close();

    return data;
}

std::vector<std::vector<double>> readTargetData(std::string bundesland){
    std::string filename = "../../work/input/global/cases_" + bundesland + ".bin";
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

    for (int i = 0; i < numEventPoints; ++i) {
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
    fprintf(gnuplotPipe, "set output 'output_data_optim_%d/Symptomatic_hybrid_interval_0_%s.png'\n", correction_step, fileAddOn.c_str());
    fprintf(gnuplotPipe, "set title ''\n");
    fprintf(gnuplotPipe, "plot ");

    fprintf(gnuplotPipe, "'%s' using 1:(column(%d*2)) with lines title sprintf('numerical, symptomatic', %d-1) lc rgb '%s', ",
            dataFileName.c_str(), 1, 1, colors[0].c_str()); 
    fprintf(gnuplotPipe, "'%s' using 1:(column(%d*2+1)) with linespoints title sprintf('target', %d-1) lc rgb '%s', ",
            dataFileName.c_str(), 1, 1, colors[1].c_str());

    pclose(gnuplotPipe);
}