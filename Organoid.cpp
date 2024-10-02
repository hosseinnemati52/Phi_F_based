#include <iostream>
#include <iomanip>
#include <math.h>
#include <tgmath.h> 
#include <unistd.h>
// #include <array>
#include <vector>
// #include <cstdlib>
// #include <ctime>
// #include <utility>
// #include <tuple>
// #include <cmath>
// #include <map>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <cstdlib>
#include <stdio.h>
// #include <cstring>
#include <string>
#include <random>
#include <chrono>
#include <sys/stat.h>
// #include <filesystem>
#include <string.h>
#include <algorithm>
#include <map>
#include <set>
#include <utility> // for std::pair
#include <variant>


///////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////#DEFINITIONS///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
#define PI              (3.14159265359)
#define CYCLING_STATE   (1)
#define G1_ARR_STATE    (-1)
#define G0_STATE        (-2)
#define DIFF_STATE      (-3)
#define APOP_STATE      (-4)
///////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////#DEFINITIONS///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
////////////////////////////// USING /////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace std::chrono;
//////////////////////////////////////////////////////////////////////////
////////////////////////////// USING /////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
////////////////////////////// STRUCTS////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
struct Config {
    
    //############### GENERAL ##############
    int N_UpperLim; // upper limit for defining vectors
    int NTypes; // number of cell types.
    std::vector<double> typeR0; // R for each type @ phi = 0
    std::vector<double> typeR2PI; // R for each type @ phi = 2*PI
    std::vector<std::vector<double>> typeTypeEpsilon; // This is a matrix that shows synchronization strength of the cell phase for each pair of cell types which are neighbors;
    //############### GENERAL ##############

    //############### FORCE AND FRICTION ##############
    std::vector<double> typeGamma;
    std::vector<std::vector<double>> typeTypeGammaCC;
    std::vector<std::vector<double>> typeTypeF_rep_max;
    std::vector<std::vector<double>> typeTypeF_abs_max;
    double R_eq_coef;
    double R_cut_coef_force;
    //############### FORCE AND FRICTION ##############
    
    //############### SELF-PROPULSION ##############
    std::vector<double> typeFm;
    std::vector<double> typeDr;
    //############### SELF-PROPULSION ##############
    
    // //############### Phi Noise ##############
    // double PhiNoiseSigma;
    // //############### Phi Noise ##############

    //############### G1 phase border / (2*PI) ##############
    double G1Border;
    //############### G1 phase border / (2*PI) ##############
    
    // //############### THRESHOLD OMEGAS ##############
    // double omegaThG1_arr;
    // double omegaThG0;
    // double omegaThDiff;
    // double omegaThApop;
    // //############### THRESHOLD OMEGAS ##############

    //############### THRESHOLD FITNESSES ##############
    std::vector<double> typeFit0;
    double Fit_Th_G1_arr;
    double Fit_Th_G0;
    double Fit_Th_Diff;
    double Fit_Th_Apop;
    //############### THRESHOLD FITNESSES ##############
    
    //############### TIMING & SAMPLING ##############
    double maxTime;
    double dt;
    double dt_sample;
    int samplesPerWrite;
    int writePerZip;
    double printingTimeInterval;
    //############### TIMING & SAMPLING ##############

    // //############### FITNESS TO OMEGA MAPPING ##############
    // std::vector<double> typeOmega0;
    // std::vector<double> typeOmegaLim;
    // std::vector<double> typeFit0;
    // std::vector<double> typeFitLim;
    // //############### FITNESS TO OMEGA MAPPING ##############

    //############### GAME ##############
    double R_cut_coef_game;
    double GameNoiseSigma;
    double tau;
    std::vector<std::vector<double>> typeTypePayOff_mat_real_C;
    std::vector<std::vector<double>> typeTypePayOff_mat_real_F1;
    std::vector<std::vector<double>> typeTypePayOff_mat_real_F2;
    std::vector<std::vector<double>> typeTypePayOff_mat_imag_C;
    std::vector<std::vector<double>> typeTypePayOff_mat_imag_F1;
    std::vector<std::vector<double>> typeTypePayOff_mat_imag_F2;
    //############### GAME ##############

    //############### INITIALIZAION ##############
    std::string initConfig;
    //############### INITIALIZAION ##############

};

enum class ObjectType {
    VECTOR,
    MATRIX
};
//////////////////////////////////////////////////////////////////////////
////////////////////////////// STRUCTS////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
////////////////////////////// PROTOTYPES ////////////////////////////////
//////////////////////////////////////////////////////////////////////////
void simulationDataReader(int* NSitesPtr, double* LxPtr, double* LyPtr, double* AlphaPtr, double* KbPtr, double* TemPtr,  int* NumCellsPtr, \
                          double* AvgCellAreaPtr, double* LambdaPtr, long* maxMCStepsPtr, int* samplesPerWritePtr, \
                          int* printingTimeIntervalPtr, int* numLinkedListPtr, string* initConfigPtr);

std::vector<double> parseVector(const std::string& str);

std::vector<std::vector<double>> parseMatrix(const std::string& str);

void trim(std::string& str);

void readConfigFile(const std::string& filename, Config& config);

void initializer(const int N_UpperLim, int* NCellsPtr, vector<int>& NCellsPerType,
                 const vector<double> typeR0, const vector<double> typeR2PI, 
                 vector<int>& cellType, vector<double>& cellX, vector<double>& cellY, 
                 vector<double>& cellVx, vector<double>& cellVy, 
                 vector<double>& cellPhi, vector<int>& cellState, vector<double>& cellTheta, vector<double>& cellR, vector<double>& cellArea,
                 vector<vector<double>>& cellFitness, const vector<double>& typeFit0);

void initial_read(const int N_UpperLim, int* NCellsPtr, vector<int>& NCellsPerType,
                 const vector<double> typeR0, const vector<double> typeR2PI, 
                 vector<int>& cellType, vector<double>& cellX, vector<double>& cellY, 
                 vector<double>& cellVx, vector<double>& cellVy, 
                 vector<double>& cellPhi, vector<int>& cellState, vector<double>& cellTheta, vector<double>& cellR, vector<double>& cellArea,
                 vector<vector<double>>& cellFitness);

void writeIntVectorToFile(const std::vector<int>& vec, int NCells, const std::string& filename);

void writeIntMatrixToFile(const std::vector<std::vector<int>>& mat, int NCells, int NCols, const std::string& filename);

void writeDoubleVectorToFile(const std::vector<double>& vec, int NCells, const std::string& filename);

void writeDoubleMatrixToFile(const std::vector<std::vector<double>>& mat, int NCells, int NCols, const std::string& filename);

void readIntVectorFromFile(const std::string& filename, std::vector<int>& data);

void readDoubleVectorFromFile(const std::string& filename, std::vector<double>& data);

void readIntMatrixFromFile(const std::string& filename, std::vector<std::vector<int>>& data);

void readDoubleMatrixFromFile(const std::string& filename, std::vector<std::vector<double>>& data);

//////////////////////////////////////////////////////////////////////////
////////////////////////////// PROTOTYPES ////////////////////////////////
//////////////////////////////////////////////////////////////////////////

int main()
{

    /////////////////// FOLDERS NAMES ////////////////////
    std::string dataFolderName = "data";
    std::string initFolderName = "init";
    std::string mainResumeFolderName = "main_resume";
    std::string backupResumeFolderName = "backup_resume";
    std::string loadFolderName;
    /////////////////// FOLDERS NAMES ////////////////////

    /////////////////// MAKING SUB DIRECTORIES /////////////////
    // This block is for windows:
    // mkdir(dataFolderName.c_str()); //making data folder
    // mkdir(initFolderName.c_str()); //making init folder
    // mkdir(mainResumeFolderName.c_str()); //making main_resume folder
    // mkdir(backupResumeFolderName.c_str()); //making backup_resume folder

    // This block is for Linux:
    // mkdir(dataFolderName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); //making data folder
    // mkdir(initFolderName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); //making init folder
    mkdir(mainResumeFolderName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); //making main_resume folder
    mkdir(backupResumeFolderName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); //making backup_resume folder
    /////////////////// MAKING SUB DIRECTORIES /////////////////


    //////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////
    ///////////// Parameters reading /////////////////////////////////////////////////////////////
    Config config;
    // readConfigFile("params_test.csv", config);
    readConfigFile("params.csv", config);

    ///////////////////////////////// VARIABLES DEFINITION /////////////////////////////
    //############### GENERAL ##############
    int N_UpperLim; // upper limit for defining vectors
    int NTypes; // number of cell types.
    vector<double> typeR0(NTypes); // The radius of each type at \phi = 0.0; 
    vector<double> typeR2PI(NTypes); // The radius of each type at \phi = 2*PI; 
    vector<vector<double>> typeTypeEpsilon(NTypes, vector<double>(NTypes)); // This is a matrix that shows synchronization strength of the cell phase for each pair of cell types which are neighbors;
    //############### GENERAL ##############

    //############### FORCE AND FRICTION ##############
    vector<double> typeGamma(NTypes); // velocity = (1/gamma) * total_force;
    vector<vector<double>> typeTypeGammaCC(NTypes, vector<double>(NTypes)); // This is a matrix that shows friction coefficient between each pair of cell types;
    vector<vector<double>> typeTypeF_rep_max(NTypes, vector<double>(NTypes)); // This is a matrix that shows F_rep_max between each pair of cell types;
    vector<vector<double>> typeTypeF_abs_max(NTypes, vector<double>(NTypes)); // This is a matrix that shows F_abs_max between each pair of cell types;
    double R_eq_coef; // R_eq_coef = (R_eq)/ (R[cell_1]+R[cell_2]);
    double R_cut_coef_force; //  R_cut_coef_force = (R_cut_off_force)/ (R[cell_1]+R[cell_2]); 
    //############### FORCE AND FRICTION ##############
    
    //############### SELF-PROPULSION ##############
    vector<double> typeFm(NTypes); // self-propulsion force amplitude of each cell type;
    vector<double> typeDr(NTypes); // Diffusion coef of self-propulsion direction of each cell type;
    //############### SELF-PROPULSION ##############

    //############### G1 phase border / (2*PI) ##############
    double G1Border; // The fraction of the whole period at which G1 ends (G1Border = (G1-end phase)/(2*PI) )
    //############### G1 phase border / (2*PI) ##############
    
    //############### THRESHOLD FITNESSES ##############
    vector<double> typeFit0(NTypes); // F_0^WT and F_0^C, which are exactly omega_0 values.
    double Fit_Th_G1_arr; // the fitness where WT cells go into G1-arrest. (only happens when they are already in G1).
    double Fit_Th_G0; // the fitness where WT cells go from G1-arrest into G0.
    double Fit_Th_Diff; // the fitness where WT cells go from G0 into differentiated state.
    double Fit_Th_Apop; // the fitness where WT cells go from differentiated state into apoptosis.
    //############### THRESHOLD FITNESSES ##############
    
    //############### TIMING & SAMPLING ##############
    double maxTime; // total simulaion time
    double dt; // simulation and integration time step
    double dt_sample; // sampling intervals
    int samplesPerWrite; // how many samples per each data writing operation
    int writePerZip; // how many written set of data for each zipping operation
    double printingTimeInterval; // for terminal show 
    //############### TIMING & SAMPLING ##############

    //############### GAME ##############
    double R_cut_coef_game;   //R_cut_coef_game = (R_cut_game)/ (R[cell_1]+R[cell_2]); 
    double GameNoiseSigma; // std of the noise term in payoff matrix entries
    double tau; // the characteristic memory time for cells fitnesses
    vector<vector<double>> typeTypePayOff_mat_real_C(NTypes, vector<double>(NTypes)); // constant term
    vector<vector<double>> typeTypePayOff_mat_real_F1(NTypes, vector<double>(NTypes)); // coefficient of fitness of player no.1 (the one that GAINS the value in the payoff matrix)
    vector<vector<double>> typeTypePayOff_mat_real_F2(NTypes, vector<double>(NTypes)); // coefficient of fitness of player no.2 (the one that LOSES the value in the payoff matrix)
    vector<vector<double>> typeTypePayOff_mat_imag_C(NTypes, vector<double>(NTypes)); // imaginary values like above
    vector<vector<double>> typeTypePayOff_mat_imag_F1(NTypes, vector<double>(NTypes));
    vector<vector<double>> typeTypePayOff_mat_imag_F2(NTypes, vector<double>(NTypes));
    //############### GAME ##############

    //############### INITIALIZAION ##############
    std::string initConfig;
    //############### INITIALIZAION ##############
    ///////////////////////////////// VARIABLES DEFINITION /////////////////////////////


    ///////////////////////////////// VARIABLES VALUE ASSIGHNMENT /////////////////////
    //############### GENERAL ##############
    N_UpperLim = config.N_UpperLim; // upper limit for defining vectors
    NTypes = config.NTypes; // number of cell types.
    typeR0 = config.typeR0; // The radius of each type at \phi = 0.0; 
    typeR2PI = config.typeR2PI; // The radius of each type at \phi = 2*PI; 
    typeTypeEpsilon = config.typeTypeEpsilon; // This is a matrix that shows synchronization strength of the cell phase for each pair of cell types which are neighbors;
    //############### GENERAL ##############

    //############### FORCE AND FRICTION ##############
    typeGamma = config.typeGamma; // velocity = (1/gamma) * total_force;
    typeTypeGammaCC = config.typeTypeGammaCC; // This is a matrix that shows friction coefficient between each pair of cell types;
    typeTypeF_rep_max = config.typeTypeF_rep_max; // This is a matrix that shows F_rep_max between each pair of cell types;
    typeTypeF_abs_max = config.typeTypeF_abs_max; // This is a matrix that shows F_abs_max between each pair of cell types;
    R_eq_coef = config.R_eq_coef; // R_eq_coef = (R_eq)/ (R[cell_1]+R[cell_2]);
    R_cut_coef_force = config.R_cut_coef_force; //  R_cut_coef_force = (R_cut_off_force)/ (R[cell_1]+R[cell_2]); 
    //############### FORCE AND FRICTION ##############
    
    //############### SELF-PROPULSION ##############
    typeFm = config.typeFm; // self-propulsion force amplitude of each cell type;
    typeDr = config.typeDr; // Diffusion coef of self-propulsion direction of each cell type;
    //############### SELF-PROPULSION ##############

    //############### G1 phase border / (2*PI) ##############
    G1Border = config.G1Border; // The fraction of the whole period at which G1 ends (G1Border = (G1-end phase)/(2*PI) )
    //############### G1 phase border / (2*PI) ##############
    
    //############### THRESHOLD FITNESSES ##############
    typeFit0 = config.typeFit0; // F_0^WT and F_0^C, which are exactly omega_0 values.
    Fit_Th_G1_arr = config.Fit_Th_G1_arr; // the fitness where WT cells go into G1-arrest. (only happens when they are already in G1).
    Fit_Th_G0 = config.Fit_Th_G0; // the fitness where WT cells go from G1-arrest into G0.
    Fit_Th_Diff = config.Fit_Th_Diff; // the fitness where WT cells go from G0 into differentiated state.
    Fit_Th_Apop = config.Fit_Th_Apop; // the fitness where WT cells go from differentiated state into apoptosis.
    //############### THRESHOLD FITNESSES ##############
    
    //############### TIMING & SAMPLING ##############
    maxTime = config.maxTime; // total simulaion time
    dt = config.dt; // simulation and integration time step
    dt_sample = config.dt_sample; // sampling intervals
    samplesPerWrite = config.samplesPerWrite; // how many samples per each data writing operation
    writePerZip = config.writePerZip; // how many written set of data for each zipping operation
    printingTimeInterval = config.printingTimeInterval; // for terminal show 
    //############### TIMING & SAMPLING ##############

    //############### GAME ##############
    R_cut_coef_game = config.R_cut_coef_game;   //R_cut_coef_game = (R_cut_game)/ (R[cell_1]+R[cell_2]); 
    GameNoiseSigma = config.GameNoiseSigma; // std of the noise term in payoff matrix entries
    tau = config.tau;
    typeTypePayOff_mat_real_C = config.typeTypePayOff_mat_real_C; // constant term
    typeTypePayOff_mat_real_F1 = config.typeTypePayOff_mat_real_F1; // coefficient of fitness of player no.1 (the one that GAINS the value in the payoff matrix)
    typeTypePayOff_mat_real_F2 = config.typeTypePayOff_mat_real_F2; // coefficient of fitness of player no.2 (the one that LOSES the value in the payoff matrix)
    typeTypePayOff_mat_imag_C = config.typeTypePayOff_mat_imag_C; // imaginary values like above
    typeTypePayOff_mat_imag_F1 = config.typeTypePayOff_mat_imag_F1;
    typeTypePayOff_mat_imag_F2 = config.typeTypePayOff_mat_imag_F2;
    //############### GAME ##############

    //############### INITIALIZAION ##############
    initConfig =  config.initConfig;
    //############### INITIALIZAION ##############
    ///////////////////////////////// VARIABLES VALUE ASSIGHNMENT /////////////////////


    vector<double> typeA_min(NTypes); // minimum area of each cell type
    vector<double> typeA_max(NTypes); // maximum area of each cell type
    for (int type_c = 0; type_c < NTypes; type_c++)
    {
        typeA_min[type_c] = PI * typeR0[type_c] * typeR0[type_c];
        typeA_max[type_c] = PI * typeR2PI[type_c] * typeR2PI[type_c];
    }
    ///////////// Parameters reading /////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////

    


    int NCells;     // total number of cells (changes during the simulation)
    vector<int> NCellsPerType(NTypes);
    vector<int> cellType(N_UpperLim);
    
    vector<double> cellX(N_UpperLim);
    vector<double> cellY(N_UpperLim);
    vector<double> cellVx(N_UpperLim);
    vector<double> cellVy(N_UpperLim);
    vector<double> cellFx(N_UpperLim); // force x
    vector<double> cellFy(N_UpperLim); // force y
    vector<double> cellTheta(N_UpperLim); // self-propulsion direction of cells, \in [0, 2*pi]

    vector<double> cellPhi(N_UpperLim); // phase in cell cycle, \in [0, 2*pi]
    vector<int> cellState(N_UpperLim); // Cell state : {cycling:CYCLING_STATE, G1_arr:G1_ARR_STATE, G0:G0_STATE, differentiated:DIFF_STATE, apop:APOP_STATE, does not exist: 0}
    // vector<double> cellOmega(N_UpperLim); // cell fitness. Changes by the game. equals to d phi / dt for cycling cells
    
    vector<double> cellR(N_UpperLim); // The radius of each cell. It may change time to time.
    vector<double> cellArea(N_UpperLim); // The area of each cell. It may change time to time.
    vector<vector<double>> cellFitness(N_UpperLim, vector<double>(2)); // This stores the complex fitness of cells (Real and Imaginary parts).

    ///////////// BUNCHES FOR WRITING SAMPLES /////////////////
    vector<double> tBunch;
    vector<vector<int>> cellTypeBunch;
    vector<vector<double>> cellXBunch;
    vector<vector<double>> cellYBunch;
    vector<vector<double>> cellVxBunch;
    vector<vector<double>> cellVyBunch;
    vector<vector<double>> cellPhiBunch;
    vector<vector<int>> cellStateBunch;
    vector<vector<double>> cellRBunch;
    vector<vector<vector<double>>> cellFitnessBunch;
    int writeCounter = 0;
    ///////////// BUNCHES FOR WRITING SAMPLES /////////////////

    ////////////////////// DEFINING RANDOM GENERATOR ///////////////
    std::mt19937 mt_rand;
    unsigned long MT_MAX = mt_rand.max();
    unsigned long MT_MIN = mt_rand.min();
    unsigned long randState;
    ////////////////////// DEFINING RANDOM GENERATOR ///////////////

    /////////////////// RANDOM GENERATOR SEEDING /////////////////
    random_device rd; // random seed creation
    int mt_rand_seed = rd();
    // mt_rand_seed = 3758017832;

    //Seeding
    mt_rand.seed(mt_rand_seed);
    
    ///////////////// INITIALIZATION ////////////////////
    if ( initConfig == "read" )
    {
    initial_read(N_UpperLim, &NCells, NCellsPerType,
                typeR0, typeR2PI, 
                cellType, cellX, cellY, 
                cellVx, cellVy, 
                cellPhi, cellState, cellTheta, cellR, cellArea,
                cellFitness);
    } else if ( initConfig == "gen" )
    {
    initializer(N_UpperLim, &NCells, NCellsPerType,
                 typeR0, typeR2PI, 
                 cellType, cellX, cellY, 
                 cellVx, cellVy, 
                 cellPhi, cellState, cellTheta, cellR, cellArea,
                 cellFitness, typeFit0);
    }

    ///////////////// INITIALIZATION ////////////////////

    // double dt = 0.01;
    double t = dt;
    double t_eps = 1e-8;
    double tLastSampling = 0.0;
    int bunchInd = 1;
    double Fx, Fy, F; // these are forces
    int cellC_1, cellC_2;
    double distance, distance2; // distance2 means distance^2
    double delta_x, delta_y; // that of 2 minus that of 1
    double R_cut_force,R_cut_game, R_eq;
    int cellType_1, cellType_2; // type of cell 1 and type of cell 2
    double deltaOmega, KuramotoTerm;
    double FaddTermX, FaddTermY;
    int newBornCells, newBornInd;
    double A_min, A_max;
    double rand_divison_angle;
    double daughterX_1, daughterY_1, daughterX_2, daughterY_2;
    double max_interaction_r2, max_interaction_r; // maximum interaction distance, and its square.
    max_interaction_r =  (max(R_cut_coef_force, R_cut_coef_game)) * 2.0 * (*std::max_element(typeR2PI.begin(), typeR2PI.end()));
    max_interaction_r2 = max_interaction_r * max_interaction_r;
    double gain_cell_1_real, gain_cell_1_imag;
    double gain_noise_real, gain_noise_imag;
    int int_rand_1, int_rand_2;
    double uniform_rand_1, uniform_rand_2, gauss_rand_1, gauss_rand_2; // for Box-Muller transform


    
    
    // Update Auxiliary properties vectors
    vector<vector<double>> cellFitnessUpdated(N_UpperLim, vector<double>(2)); // This stores the updated values of fitness of cells (Real and Imaginary parts).
    vector<double> cellPhiUpdated(N_UpperLim); // Updated phases in cell cycle, \in [0, 2*pi]
    vector<int> cellStateUpdated(N_UpperLim); // Cell state Updated values: {cycling:CYCLING_STATE, G1_arr:G1_ARR_STATE, G0:G0_STATE, differentiated:DIFF_STATE, apop:APOP_STATE, does not exist: 0}
    for (cellC_1 = 0; cellC_1 < N_UpperLim; cellC_1++) // initializing: setting everything to zero
        {
            cellFitnessUpdated[cellC_1][0] = 0.0;
            cellFitnessUpdated[cellC_1][1] = 0.0;

            cellPhiUpdated[cellC_1] = 0.0;

            cellStateUpdated[cellC_1] = 0;
        }
    // Update Auxiliary properties vectors


    /////// SIMULATION LOOP /////////////
    while (t < maxTime)
    {
        // setting forces to zero
        for (cellC_1 = 0; cellC_1 < NCells; cellC_1++)
        {
            cellFx[cellC_1]= 0.0;
            cellFy[cellC_1]= 0.0;
        }
        // setting forces to zero
        
        // calculating Fx , Fy, and Updated fitnesses, without changing X, Y, Vx, Vy, and fitnesses
        for (cellC_1 = 0; cellC_1 < NCells; cellC_1++) // loop on cellC_1
        {
            cellType_1 = cellType[cellC_1];

            // cellOmega[cellC_1] = typeOmega0[cellType_1];
            // cellOmega[cellC_1] = typeFit0[cellType_1];
            KuramotoTerm = 0.0;

            for (cellC_2 = cellC_1 + 1 ; cellC_2 < NCells; cellC_2++) // loop on cellC_2, for interactions (force and game)
            {
                
                delta_x = cellX[cellC_2] - cellX[cellC_1];
                delta_y = cellY[cellC_2] - cellY[cellC_1];
                distance2 = delta_x * delta_x + delta_y * delta_y;

                if (distance2 > max_interaction_r2) // if the distance is larger than the maximum interaction distance, do nothing and go to the next cellC_2. No mutual interaction!
                {
                    continue;
                }
                else // if (distance2 <= max_interaction_r2), they MAY interact. 
                {
                    distance = pow(distance2, 0.5);

                    R_cut_force = R_cut_coef_force * (cellR[cellC_2] + cellR[cellC_1]);
                    R_eq =   R_eq_coef * (cellR[cellC_2] + cellR[cellC_1]);

                    cellType_2 = cellType[cellC_2];

                    /////////// Forces ////////////////
                    if (distance < R_cut_force ) // They do interact forcewise
                    {
                        
                        if (distance < R_eq )
                        {
                            F = typeTypeF_rep_max[cellType_1][cellType_2] * (distance - R_eq) / R_eq;
                        } else
                        {
                            F = typeTypeF_abs_max[cellType_1][cellType_2] * (distance - R_eq) / (R_cut_force - R_eq);
                        }
                        
                        FaddTermX = ( F * (delta_x / distance)  + typeTypeGammaCC[cellType_1][cellType_2] * (cellVx[cellC_2] - cellVx[cellC_1]) );
                        FaddTermY = ( F * (delta_y / distance)  + typeTypeGammaCC[cellType_1][cellType_2] * (cellVy[cellC_2] - cellVy[cellC_1]) );

                        cellFx[cellC_1] += FaddTermX;
                        cellFy[cellC_1] += FaddTermY;

                        cellFx[cellC_2] -= FaddTermX;
                        cellFy[cellC_2] -= FaddTermY;
                    }
                    /////////// Forces ////////////////

                    /////////// Game-interaction ////////////////
                    R_cut_game = R_cut_coef_game * (cellR[cellC_2] + cellR[cellC_1]);
                    if (distance < R_cut_game ) // They do interact gamewise
                    {   

                        int_rand_1 = mt_rand();
                        while(int_rand_1 == MT_MIN || int_rand_1 == MT_MAX){int_rand_1 = mt_rand();}
                        int_rand_2 = mt_rand();
                        while(int_rand_2 == MT_MIN || int_rand_2 == MT_MAX){int_rand_2 = mt_rand();}

                        uniform_rand_1 = (((long double)(int_rand_1)-MT_MIN)/((long double)MT_MAX-MT_MIN));
                        uniform_rand_2 = (((long double)(int_rand_2)-MT_MIN)/((long double)MT_MAX-MT_MIN));

                        gauss_rand_1 =  pow( (-2.0 * log(uniform_rand_1)) , 0.5) * cos(2.0 * PI * uniform_rand_2); // Box-Muller transform
                        // gauss_rand_2 = 

                        gain_noise_real = (GameNoiseSigma / dt) * gauss_rand_1;
                        // gain_noise_imag = 

                        gain_cell_1_real = gain_noise_real + \
                                           2;
                                           INJAAAAAAAAAAAAAAAAAAAAAAAAAA
                        // gain_cell_1_imag = 


                        
                        // cellFitness[cellC_2][0] += (typeTypePayOff_mat_real[cellType_2][cellType_1] * dt) ;
                        // cellFitness[cellC_2][1] += (typeTypePayOff_mat_imag[cellType_2][cellType_1] * dt) ;
                        // cellFitness[cellC_1][0] += (typeTypePayOff_mat_real[cellType_1][cellType_2] * dt) ;
                        // cellFitness[cellC_1][1] += (typeTypePayOff_mat_imag[cellType_1][cellType_2] * dt) ;
                        /////////// Game (Fitness update) ////////////////
                        KuramotoTerm += (typeTypeEpsilon[cellType_1][cellType_2] * sin (cellPhi[cellC_2] - cellPhi[cellC_1]));
                    }
                    /////////// Game-interaction ////////////////


                }

                

                

                
                

                
                
                

                
            }

            

            /////////// delta Omega ////////////////
            if (cellType_1==0)
            {
                if (cellFitness[cellC_1][0]>=0.0)
                {
                    deltaOmega = 0.0;
                } else if (cellFitness[cellC_1][0] < 0.0)
                {
                    // deltaOmega = (typeOmegaLim[cellType_1]- typeOmega0[cellType_1]) * (1.0 - exp(cellFitness[cellC_1][0]));
                }


            } else if (cellType_1==1)
            {
                if (cellFitness[cellC_1][0] <= 0.0)
                {
                    deltaOmega = 0.0;
                } else if (cellFitness[cellC_1][0] > 0.0)
                {
                    // deltaOmega = (typeOmegaLim[cellType_1]- typeOmega0[cellType_1]) * (1.0 - exp(-cellFitness[cellC_1][0]));
                }
            }

            cellOmega[cellC_1] += deltaOmega;
            cellOmega[cellC_1] += KuramotoTerm;
            /// !!! Updating cellPhi must be the final step of loop, because it may change the NCells !!! ///
            
            /////////// delta Omega ////////////////
            
        }

        newBornCells = 0;
        newBornInd = NCells;
        for (cellC_1 = 0; cellC_1 < NCells; cellC_1++)
        {
            cellType_1 = cellType[cellC_1];

            cellVx[cellC_1] = cellFx[cellC_1] / typeGamma[cellType_1];
            cellVy[cellC_1] = cellFy[cellC_1] / typeGamma[cellType_1];

            cellX[cellC_1] += cellVx[cellC_1] * dt;
            cellY[cellC_1] += cellVy[cellC_1] * dt;

            cellPhi[cellC_1] += cellOmega[cellC_1] * dt;
            
            A_min = typeA_min[cellType_1];
            A_max = typeA_max[cellType_1];

            if (cellPhi[cellC_1] < 2.0 *PI)
            {
                
                cellArea[cellC_1] = A_min + (A_max - A_min) * 0.5 * (1 - cos(cellPhi[cellC_1]/2.0));
                cellR[cellC_1] = pow(cellArea[cellC_1] / PI, 0.5);

            } else // CELL DIVISION
            {
                rand_divison_angle = (((long double)(mt_rand())-MT_MIN)/((long double)MT_MAX-MT_MIN)) * (2.0 * PI);

                daughterX_1 = cellX[cellC_1] + (cellR[cellC_1] / 1.4142) * cos(rand_divison_angle);
                daughterY_1 = cellY[cellC_1] + (cellR[cellC_1] / 1.4142) * sin(rand_divison_angle);
                daughterX_2 = cellX[cellC_1] - (cellR[cellC_1] / 1.4142) * cos(rand_divison_angle);
                daughterY_2 = cellY[cellC_1] - (cellR[cellC_1] / 1.4142) * sin(rand_divison_angle);
                
                // new Born cell
                cellType[newBornInd] = cellType_1;

                cellX[newBornInd] = daughterX_2;
                cellY[newBornInd] = daughterY_2;

                cellVx[newBornInd] = cellVx[cellC_1] / 2.0;
                cellVy[newBornInd] = cellVy[cellC_1] / 2.0;

                cellArea[newBornInd] = A_min;
                cellR[newBornInd] = typeR0[cellType_1];

                cellFitness[newBornInd][0] = 0.0;
                cellFitness[newBornInd][1] = 0.0;
                // new Born cell

                // original cell
                cellX[cellC_1] = daughterX_1;
                cellY[cellC_1] = daughterY_1;

                cellVx[cellC_1] = cellVx[cellC_1] / 2.0;
                cellVy[cellC_1] = cellVy[cellC_1] / 2.0;

                cellArea[cellC_1] = A_min;
                cellR[cellC_1] = typeR0[cellType_1];

                cellPhi[cellC_1] = 0.0 ;
                // original cell

                newBornInd++;
                newBornCells++;
            }
        }

        NCells += newBornCells;

        
        ///// Sampling operation //////
        if ((t-tLastSampling) > dt_sample - t_eps)
        {

            //  take sample and add it to the current bunch
            tBunch.push_back(t);
            cellTypeBunch.push_back(cellType);
            cellXBunch.push_back(cellX);
            cellYBunch.push_back(cellY);
            cellVxBunch.push_back(cellVx);
            cellVyBunch.push_back(cellVy);
            cellPhiBunch.push_back(cellPhi);
            cellRBunch.push_back(cellR);
            cellFitnessBunch.push_back(cellFitness);

            ////// writing the bunch /////////
            if (cellTypeBunch.size() == samplesPerWrite)
            {
                
                // write the bunch
                

                tBunch.clear();
                cellTypeBunch.clear();
                cellXBunch.clear();
                cellYBunch.clear();
                cellVxBunch.clear();
                cellVyBunch.clear();
                cellPhiBunch.clear();
                cellRBunch.clear();
                cellFitnessBunch.clear();

                bunchInd++;
                
                ///////////// NUMBERED DATA ZIPPING ///////////////////
                if ((writeCounter+1)%writePerZip ==0)
                {
                    // string argv1 = 'data_'+to_string((int)(writeCounter/writePerZip))+'.zip';
                    std::string argv1 = "data_" + std::to_string(static_cast<int>(writeCounter / writePerZip)+1) + ".zip";
                    std::string argv2 = "data";
                    std::string command = "python3 dataZipperNum.py " + argv1 + " " + argv2;
                    system(command.c_str());
                }
                ///////////// NUMBERED DATA ZIPPING ///////////////////



                ///////////////////////// LS SAVING ///////////////////////
                for(int resumeFolderCounter=1; resumeFolderCounter<=2; resumeFolderCounter++)
                {
                if (resumeFolderCounter==1)
                {
                    loadFolderName = mainResumeFolderName;
                }
                else if (resumeFolderCounter==2)
                {
                    loadFolderName = backupResumeFolderName;
                }

                // writing the final time
                ofstream tLSCheck;
                tLSCheck.open(loadFolderName+"/"+"tLSCheck.csv");
                tLSCheck << t;
                tLSCheck.close();

                writeIntVectorToFile(cellType, NCells, loadFolderName+"/Type_LS.txt");
                writeDoubleVectorToFile(cellX, NCells, loadFolderName+"/X_LS.txt");
                writeDoubleVectorToFile(cellY, NCells, loadFolderName+"/Y_LS.txt");
                writeDoubleVectorToFile(cellVx, NCells, loadFolderName+"/Vx_LS.txt");
                writeDoubleVectorToFile(cellVy, NCells, loadFolderName+"/Vy_LS.txt");
                writeDoubleVectorToFile(cellPhi, NCells, loadFolderName+"/Phi_LS.txt");
                writeDoubleVectorToFile(cellR, NCells, loadFolderName+"/R_LS.txt");
                writeDoubleMatrixToFile(cellFitness, NCells, 2,  loadFolderName+"/Fitness_LS.txt");

                // writing the final state of random generator
                std::ofstream randStateLS(loadFolderName+"/"+"randStateLS.csv");
                randStateLS << mt_rand;
                randStateLS.close();

                // writing round counter
                ofstream roundCounterLS;
                roundCounterLS.open(loadFolderName+"/"+"roundCounterLS.csv");
                roundCounterLS << writeCounter;
                roundCounterLS.close();

                // writing the final time
                ofstream tLS;
                tLS.open(loadFolderName+"/"+"tLS.csv");
                tLS << t;
                tLS.close();

                // cout << "\033[2J\033[1;1H";
                }
                ///////////////////////// LS SAVING ///////////////////////


                writeCounter++;
            }
            ////// writing the bunch /////////

            // writeIntVectorToFile(cellType, NCells, "data/Type_"+ to_string(ind) + ".txt");
            // writeDoubleVectorToFile(cellX, NCells, "data/X_"+ to_string(ind) + ".txt");
            // writeDoubleVectorToFile(cellY, NCells, "data/Y_"+ to_string(ind) + ".txt");
            // writeDoubleVectorToFile(cellVx, NCells, "data/Vx_"+ to_string(ind) + ".txt");
            // writeDoubleVectorToFile(cellVy, NCells, "data/Vy_"+ to_string(ind) + ".txt");
            // writeDoubleVectorToFile(cellPhi, NCells, "data/Phi_"+ to_string(ind) + ".txt");
            // writeDoubleVectorToFile(cellR, NCells, "data/R_"+ to_string(ind) + ".txt");
            // // writeDoubleVectorToFile(cellTheta, NCells, "init/cellTheta_init.txt");
            // writeDoubleMatrixToFile(cellFitness, NCells, 2,  "data/Fitness_"+ to_string(ind) + ".txt");
            tLastSampling = t;
        }
        ///// Sampling operation //////


        t += dt;
    }
    /////// SIMULATION LOOP /////////////
    






    return 0;
}


//////////////////////////////////////////////////////////////////////////
////////////////////////////// FUNCTIONS /////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void simulationDataReader(int* NSitesPtr, double* LxPtr, double* LyPtr, double* AlphaPtr, double* KbPtr, double* TemPtr,  int* NumCellsPtr, \
                          double* AvgCellAreaPtr, double* LambdaPtr, long* maxMCStepsPtr, int* samplesPerWritePtr, \
                          int* printingTimeIntervalPtr, int* numLinkedListPtr, string* initConfigPtr)
{   
    // int L_read;
    // int NumCells_read;
    // int samplesPerWrite_read;

    fstream newfile;
    newfile.open("simulationData_vec.csv",ios::in); //open a file to perform read operation using file object
    if (newfile.is_open()){ //checking whether the file is open
        string tp;
        while(getline(newfile, tp)){ //read data from file object and put it into string.
            // cout << tp << "\n"; //print the data of the string
            // if (strstr((const char)tp, "NSites = "))
            const char* tpChar = tp.c_str();
            if (strstr(tpChar, "NSites = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *NSitesPtr = std::stoi(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "Lx = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *LxPtr = std::stod(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "Ly = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *LyPtr = std::stod(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "Alpha = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *AlphaPtr = std::stod(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "Kb = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *KbPtr = std::stod(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "Tem = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *TemPtr = std::stod(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "NumCells = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *NumCellsPtr = std::stoi(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "AvgCellArea = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *AvgCellAreaPtr = std::stod(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "Lambda = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *LambdaPtr = std::stod(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "SweepLength = "))
            {
                // Do nothing
                continue;
            }
            if (strstr(tpChar, "maxMCSteps = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *maxMCStepsPtr = std::stoi(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "samplesPerWrite = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *samplesPerWritePtr = std::stoi(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "printingTimeInterval = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *printingTimeIntervalPtr = std::stoi(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "numLinkedList = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *numLinkedListPtr = std::stoi(num_str);  // convert substring to int
                continue;
            }
            if (strstr(tpChar, "initConfig = "))
            {
                std::size_t pos = tp.find('=');  // find position of '='
                std::string num_str = tp.substr(pos+2);  // extract substring starting from position after '='
                *initConfigPtr = num_str;
                continue;
            }
        }
    }
    newfile.close(); //close the file object.

    // *L_read_ptr = L_read;
    // *NumCells_read_ptr = NumCells_read;
    // *samplesPerWrite_read_ptr = samplesPerWrite_read;
}

std::vector<double> parseVector(const std::string& str) {
    std::vector<double> result;
    std::stringstream ss(str.substr(1, str.size() - 2)); // Remove the square brackets
    std::string item;
    while (std::getline(ss, item, ',')) {
        result.push_back(std::stod(item));
    }
    return result;
}

std::vector<std::vector<double>> parseMatrix(const std::string& str) {
    std::vector<std::vector<double>> result;
    std::stringstream ss(str.substr(1, str.size() - 2)); // Remove the square brackets
    std::string row;
    while (std::getline(ss, row, ';')) {
        result.push_back(parseVector(row));
    }
    return result;
}

void trim(std::string& str) {
    str.erase(str.begin(), std::find_if(str.begin(), str.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
    str.erase(std::find_if(str.rbegin(), str.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), str.end());
}

void readConfigFile(const std::string& filename, Config& config) {
    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string key;
        std::getline(ss, key, '=');
        std::string value;
        std::getline(ss, value);

        trim(key);
        trim(value);

        if (key == "N_UpperLim") {
            config.N_UpperLim = std::stoi(value);
        } else if (key == "NTypes") {
            config.NTypes = std::stoi(value);
        } else if (key == "typeR0") {
            config.typeR0 = parseVector(value);
        } else if (key == "typeR2PI") {
            config.typeR2PI = parseVector(value);
        } else if (key == "typeGamma") {
            config.typeGamma = parseVector(value);
        } else if (key == "typeTypeGammaCC") {
            config.typeTypeGammaCC = parseMatrix(value);
        } else if (key == "typeTypeEpsilon") {
            config.typeTypeEpsilon = parseMatrix(value);
        } else if (key == "typeFm") {
            config.typeFm = parseVector(value);
        } else if (key == "typeDr") {
            config.typeDr = parseVector(value);
        } else if (key == "R_eq_coef") {
            config.R_eq_coef = std::stod(value);
        } else if (key == "R_cut_coef_force") {
            config.R_cut_coef_force = std::stod(value);
        } else if (key == "R_cut_coef_game") {
            config.R_cut_coef_game = std::stod(value);
        } else if (key == "typeTypeF_rep_max") {
            config.typeTypeF_rep_max = parseMatrix(value);
        } else if (key == "typeTypeF_abs_max") {
            config.typeTypeF_abs_max = parseMatrix(value);
        } else if (key == "typeTypePayOff_mat_real_C") {
            config.typeTypePayOff_mat_real_C = parseMatrix(value);
        } else if (key == "typeTypePayOff_mat_real_F1") {
            config.typeTypePayOff_mat_real_F1 = parseMatrix(value);
        } else if (key == "typeTypePayOff_mat_real_F2") {
            config.typeTypePayOff_mat_real_F2 = parseMatrix(value);
        } else if (key == "typeTypePayOff_mat_imag_C") {
            config.typeTypePayOff_mat_imag_C = parseMatrix(value);
        } else if (key == "typeTypePayOff_mat_imag_F1") {
            config.typeTypePayOff_mat_imag_F1 = parseMatrix(value);
        } else if (key == "typeTypePayOff_mat_imag_F2") {
            config.typeTypePayOff_mat_imag_F2 = parseMatrix(value);
        } else if (key == "typeOmega0") {
            // config.typeOmega0 = parseVector(value);
            {}
        } else if (key == "typeOmegaLim") {
            // config.typeOmegaLim = parseVector(value);
            {}
        } else if (key == "typeFit0") {
            config.typeFit0 = parseVector(value);
        } else if (key == "typeFitLim") {
            // config.typeFitLim = parseVector(value);
            {}
        } else if (key == "maxTime") {
            config.maxTime = std::stod(value);
        } else if (key == "dt") {
            config.dt = std::stod(value);
        } else if (key == "dt_sample") {
            config.dt_sample = std::stod(value);
        } else if (key == "samplesPerWrite") {
            config.samplesPerWrite = std::stoi(value);
        } else if (key == "writePerZip") {
            config.writePerZip = std::stoi(value);
        } else if (key == "printingTimeInterval") {
            config.printingTimeInterval = std::stod(value);
        } else if (key == "GameNoiseSigma") {
            config.GameNoiseSigma = std::stod(value);
        } else if (key == "PhiNoiseSigma") {
            // config.PhiNoiseSigma = std::stod(value);
            {}
        } else if (key == "G1Border") {
            config.G1Border = std::stod(value);
        } else if (key == "omegaThG1_arr") {
            // config.omegaThG1_arr = std::stod(value);
            {}
        } else if (key == "omegaThG0") {
            // config.omegaThG0 = std::stod(value);
            {}
        } else if (key == "omegaThDiff") {
            // config.omegaThDiff = std::stod(value);
            {}
        } else if (key == "omegaThApop") {
            // config.omegaThApop = std::stod(value);
            {}
        } else if (key == "Fit_Th_G1_arr") {
            config.Fit_Th_G1_arr = std::stod(value);
        } else if (key == "Fit_Th_G0") {
            config.Fit_Th_G0 = std::stod(value);
        } else if (key == "Fit_Th_Diff") {
            config.Fit_Th_Diff = std::stod(value);
        } else if (key == "Fit_Th_Apop") {
            config.Fit_Th_Apop = std::stod(value);
        } else if (key == "initConfig") {
            config.initConfig = value;
        } else if (key == "tau") {
            config.tau = std::stod(value);
        }
    }
}

void initializer(const int N_UpperLim, int* NCellsPtr, vector<int>& NCellsPerType,
                 const vector<double> typeR0, const vector<double> typeR2PI, 
                 vector<int>& cellType, vector<double>& cellX, vector<double>& cellY, 
                 vector<double>& cellVx, vector<double>& cellVy, 
                 vector<double>& cellPhi, vector<int>& cellState, vector<double>& cellTheta, vector<double>& cellR, vector<double>& cellArea,
                 vector<vector<double>>& cellFitness, const vector<double>& typeFit0)
{
    NCellsPerType[0] = 500;
    NCellsPerType[1] = 50;

    int NCells = NCellsPerType[0] + NCellsPerType[1];
    
    *NCellsPtr = NCells;

    std::srand(static_cast<unsigned int>(std::time(0)));
    
    int random_int;
    float random_float;

    vector<double> A_tot_sq(2);

    int cellC = 0;

    int typeInd = 0;
    double A_min = PI * typeR0[typeInd] * typeR0[typeInd];
    double A_max = PI * typeR2PI[typeInd] * typeR2PI[typeInd];
    while (cellC < NCellsPerType[0])
    {   
        cellType[cellC] = typeInd;

        random_float = static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX);
        cellPhi[cellC] = random_float * (2*PI);

        cellState[cellC] = CYCLING_STATE; // all the cells start from cycling state

        // cellArea[cellC] = A_min + (A_max - A_min) * 0.5 * (1 - cos(cellPhi[cellC]/2.0)); // cosine area independency to phi
        cellArea[cellC] = A_min + (A_max - A_min) * cellPhi[cellC] / (2 * PI); // linear area independency to phi
        cellR[cellC] = pow(cellArea[cellC] / PI, 0.5);
        A_tot_sq[typeInd] += 4 * cellR[cellC] * cellR[cellC];

        random_float = static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX);
        cellTheta[cellC] = random_float * (2*PI);

        cellVx[cellC] = 0.0;
        cellVy[cellC] = 0.0;

        cellFitness[cellC][0] = typeFit0[typeInd];
        cellFitness[cellC][1] = 0.0;

        cellC++;
    }

    typeInd = 1;
    A_min = PI * typeR0[typeInd] * typeR0[typeInd];
    A_max = PI * typeR2PI[typeInd] * typeR2PI[typeInd];
    while (cellC < NCellsPerType[0] + NCellsPerType[1])
    {   
        cellType[cellC] = typeInd;

        random_float = static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX);
        cellPhi[cellC] = random_float * (2*PI);

        cellState[cellC] = CYCLING_STATE; // all the cells start from cycling state

        // cellArea[cellC] = A_min + (A_max - A_min) * 0.5 * (1 - cos(cellPhi[cellC]/2.0));  // cosine area independency to phi
        cellArea[cellC] = A_min + (A_max - A_min) * cellPhi[cellC] / (2 * PI); // linear area independency to phi
        cellR[cellC] = pow(cellArea[cellC] / PI, 0.5);
        A_tot_sq[typeInd] += 4 * cellR[cellC] * cellR[cellC];

        random_float = static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX);
        cellTheta[cellC] = random_float * (2*PI);

        cellVx[cellC] = 0.0;
        cellVy[cellC] = 0.0;

        cellFitness[cellC][0] = typeFit0[typeInd];
        cellFitness[cellC][1] = 0.0;

        cellC++;
    }

    while (cellC < N_UpperLim)
    {
        cellPhi[cellC] = 0;
        cellState[cellC] = 0;

        cellArea[cellC] = 0;
        cellR[cellC] = 0;
        cellTheta[cellC] = 0;

        cellX[cellC] = 0;
        cellY[cellC] = 0;

        cellVx[cellC] = 0;
        cellVy[cellC] = 0;

        cellFitness[cellC][0] = 0;
        cellFitness[cellC][1] = 0;

        cellC++;
    }

    ////// initialization of X and Y //////
    double Lx , Ly;
    double R_tot;
    R_tot = pow( (A_tot_sq[0] + A_tot_sq[1]) / PI, 0.5);
    Lx = 2.0 * R_tot;
    Ly = 2.0 * R_tot;

    // finding h (border of WT and Cancer cells)
    double h = R_tot;
    double dh = 0.001 * R_tot;
    double a = 0.0;
    while (a < A_tot_sq[1])
    {
        a = R_tot * R_tot * acos(h/R_tot) - h * pow(R_tot * R_tot - h * h , 0.5);
        h = h - dh;
    }
    // finding h (border of WT and Cancer cells)

    cellC = 0;
    while (cellC < NCellsPerType[0] + NCellsPerType[1])
    {   
        typeInd = cellType[cellC];
        
        int repeat_cond, out_cond, not_part_cond;
        double x, y;

        do
        {
            random_float = static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX);
            x = -R_tot + (2 * R_tot) * random_float;

            random_float = static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX);
            y = -R_tot + (2 * R_tot) * random_float;

            out_cond =  ( (x*x + y*y) > (R_tot*R_tot) );
            if (typeInd == 0)
            {
                not_part_cond = (y > h);
            }
            else if (typeInd == 1)
            {
                not_part_cond = (y <= h);
            }

            repeat_cond = out_cond || not_part_cond;

        } while (repeat_cond);
        
        
        cellX[cellC] = x;
        cellY[cellC] = y;

        cellC++;
    }
    ////// initialization of X and Y //////

    //// This block is for windows:
    // mkdir(ppDataFolderName.c_str()); //making data folder
    //// This block is for Linux:
    mkdir("init", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); //making backup_resume folder

    writeIntVectorToFile(cellType, NCells, "init/cellType_init.txt");
    writeDoubleVectorToFile(cellX, NCells, "init/cellX_init.txt");
    writeDoubleVectorToFile(cellY, NCells, "init/cellY_init.txt");
    writeDoubleVectorToFile(cellVx, NCells, "init/cellVx_init.txt");
    writeDoubleVectorToFile(cellVy, NCells, "init/cellVy_init.txt");
    writeDoubleVectorToFile(cellPhi, NCells, "init/cellPhi_init.txt");
    writeIntVectorToFile(cellState, NCells, "init/cellState_init.txt");
    writeDoubleVectorToFile(cellR, NCells, "init/cellR_init.txt");
    writeDoubleVectorToFile(cellTheta, NCells, "init/cellTheta_init.txt");
    writeDoubleMatrixToFile(cellFitness, NCells, 2,  "init/cellFitness_init.txt");
}

void initial_read(const int N_UpperLim, int* NCellsPtr, vector<int>& NCellsPerType,
                 const vector<double> typeR0, const vector<double> typeR2PI, 
                 vector<int>& cellType, vector<double>& cellX, vector<double>& cellY, 
                 vector<double>& cellVx, vector<double>& cellVy, 
                 vector<double>& cellPhi, vector<int>& cellState, vector<double>& cellTheta, vector<double>& cellR, vector<double>& cellArea,
                 vector<vector<double>>& cellFitness)
{
    
    vector<int> cellType_read;
    vector<double> cellX_read;
    vector<double> cellY_read;
    vector<double> cellVx_read;
    vector<double> cellVy_read;
    vector<double> cellPhi_read;
    vector<int> cellState_read;
    vector<double> cellTheta_read;
    vector<double> cellR_read;
    vector<vector<double>> cellFitness_read;

    readIntVectorFromFile("init/Type_init.txt", cellType_read);
    readDoubleVectorFromFile("init/X_init.txt", cellX_read);
    readDoubleVectorFromFile("init/Y_init.txt", cellY_read);
    readDoubleVectorFromFile("init/R_init.txt", cellR_read);
    readDoubleVectorFromFile("init/Vx_init.txt", cellVx_read);
    readDoubleVectorFromFile("init/Vy_init.txt", cellVy_read);
    readDoubleVectorFromFile("init/Phi_init.txt", cellPhi_read);
    readIntVectorFromFile("init/State_init.txt", cellState_read);
    readDoubleVectorFromFile("init/Theta_init.txt", cellTheta_read);
    readDoubleMatrixFromFile("init/Fitness_init.txt", cellFitness_read);

    int NCells = cellType_read.size();
    *NCellsPtr = NCells;

    NCellsPerType[0] = 0;
    NCellsPerType[1] = 0;

    for (int cellC = 0; cellC < NCells; cellC++)
    {
        cellType[cellC] = cellType_read[cellC];
        cellX[cellC] = cellX_read[cellC];
        cellY[cellC] = cellY_read[cellC];
        cellR[cellC] = cellR_read[cellC];
        cellVx[cellC] = cellVx_read[cellC];
        cellVy[cellC] = cellVy_read[cellC];
        cellPhi[cellC] = cellPhi_read[cellC];
        cellState[cellC] = cellState_read[cellC];
        // cellTheta[cellC] = cellTheta_read[cellC];
        cellTheta[cellC] = 0.0;

        cellFitness[cellC][0]  = cellFitness_read[cellC][0];
        cellFitness[cellC][1]  = cellFitness_read[cellC][1];

        NCellsPerType[cellType[cellC]]++;

        cellArea[cellC] = (PI * cellR[cellC] * cellR[cellC]);
    }

}

void writeIntVectorToFile(const std::vector<int>& vec, int NCells, const std::string& filename) {
    std::ofstream outFile(filename);

    if (!outFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    int vecSize = vec.size();
    NCells = std::min(NCells, vecSize);

    for (int i = 0; i < NCells; ++i) {
        outFile << vec[i] << std::endl;
    }

    outFile.close();
    std::cout << "Data written to file: " << filename << std::endl;
}

void writeIntMatrixToFile(const std::vector<std::vector<int>>& mat, int NCells, int NCols, const std::string& filename) {
    std::ofstream outFile(filename);

    if (!outFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    int matRows = mat.size();
    int matCols = mat[0].size();
    NCells = std::min(NCells, matRows);
    NCols = std::min(NCols, matCols);

    for (int i = 0; i < NCells; ++i) {
        for (int j = 0; j < NCols; ++j) {
            outFile << mat[i][j];
            if (j < NCols - 1) {
                outFile << ", ";
            }
        }
        outFile << std::endl;
    }

    outFile.close();
    std::cout << "Data written to file: " << filename << std::endl;
}

void writeDoubleVectorToFile(const std::vector<double>& vec, int NCells, const std::string& filename) {
    std::ofstream outFile(filename);

    if (!outFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    outFile << std::fixed << std::setprecision(8);

    int vecSize = vec.size();
    NCells = std::min(NCells, vecSize);

    for (int i = 0; i < NCells; ++i) {
        outFile << vec[i] << std::endl;
    }

    outFile.close();
    std::cout << "Data written to file: " << filename << std::endl;
}

void writeDoubleMatrixToFile(const std::vector<std::vector<double>>& mat, int NCells, int NCols, const std::string& filename) {
    std::ofstream outFile(filename);

    if (!outFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    outFile << std::fixed << std::setprecision(8);

    int matRows = mat.size();
    int matCols = mat[0].size();
    NCells = std::min(NCells, matRows);
    NCols = std::min(NCols, matCols);

    for (int i = 0; i < NCells; ++i) {
        for (int j = 0; j < NCols; ++j) {
            outFile << mat[i][j];
            if (j < NCols - 1) {
                outFile << ", ";
            }
        }
        outFile << std::endl;
    }

    outFile.close();
    std::cout << "Data written to file: " << filename << std::endl;
}

void readIntVectorFromFile(const std::string& filename, std::vector<int>& data) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line;
    int value;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        while (ss >> value) {
            data.push_back(value);
            if (ss.peek() == ',') ss.ignore();
        }
    }

    file.close();
}

void readDoubleVectorFromFile(const std::string& filename, std::vector<double>& data) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line;
    double value;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        while (ss >> value) {
            data.push_back(value);
            if (ss.peek() == ',') ss.ignore();
        }
    }

    file.close();
}

void readIntMatrixFromFile(const std::string& filename, std::vector<std::vector<int>>& data) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::vector<int> row;
        int value;
        char comma;

        while (ss >> value) {
            row.push_back(value);
            if (ss.peek() == ',') ss.ignore();
        }

        data.push_back(row);
    }

    file.close();
}

void readDoubleMatrixFromFile(const std::string& filename, std::vector<std::vector<double>>& data) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::vector<double> row;
        double value;
        char comma;

        while (ss >> value) {
            row.push_back(value);
            if (ss.peek() == ',') ss.ignore();
        }

        data.push_back(row);
    }

    file.close();
}
//////////////////////////////////////////////////////////////////////////
////////////////////////////// FUNCTIONS /////////////////////////////////
//////////////////////////////////////////////////////////////////////////