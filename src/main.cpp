#include "swap.h"

// Default run parameters
int N = 5;
double Size = pow(N/density, 1/3.);
double T = 0.04; 
int tau = 100000;
int tw = 1;
int cycles = 1;
int steps = tw*(cycles-1)+tau;
int linPoints = 50;
int logPoints = 50;
double p_flip = 0.2;

// Setting arrays
configuration cfg;
// int *mol_index = nullptr;
// double *X = nullptr, *Y = nullptr, *Z = nullptr, *S = nullptr, *Sref = nullptr, 
//        *X0 = nullptr, *Y0 = nullptr, *Z0 = nullptr;
// double *Xfull = nullptr, *Yfull = nullptr, *Zfull = nullptr, 
//        *Xref = nullptr, *Yref = nullptr, *Zref = nullptr;
// std::vector < std::vector <double>> Xtw, Ytw, Ztw;
// std::vector < std::vector<int> > NL, NN, BN;
// std::vector < std::vector < std::vector <int>>> NN_tw, RL;
std::vector < std::string > allObs;

std::string input;
std::string outdir;

//-----------------------------------------------------------------------------
//  main.cpp
int main(int argc, const char * argv[]) {
    
    // Random number generator
    srand(31);

    // Define the command-line options
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("input", po::value<std::string>(&input)->required(), "set input file")
        ("outdir", po::value<std::string>(&outdir)->required(), "set out directory")
        ("N", po::value<int>(&N)->default_value(N), "set system size")
        ("T", po::value<double>(&T)->default_value(T), "set temperature")
        ("tau", po::value<int>(&tau)->default_value(tau), "set single-run time")
        ("tw", po::value<int>(&tw)->default_value(tw), "set waiting time")
        ("cycles", po::value<int>(&cycles)->default_value(cycles), "set number of cycles")
        ("lin", po::value<int>(&linPoints)->default_value(linPoints), "set number of lin-spaced snapshots")
        ("log", po::value<int>(&logPoints)->default_value(logPoints), "set number of log-spaced snapshots")
        ("p_flip", po::value<double>(&p_flip)->default_value(p_flip), "set flip-attempt probability")
        ("MSD", "Flag to compute MSD")
        ("Cb", "Flag to compute Cb")
        ("Fs", "Flag to compute Fs")
        ("U", "Flag to compute U");
    // std::string input = motherdir + argv[1];
    // std::string outdir = motherdir + argv[2] + "results/";

    // Parse the command-line arguments
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
    } catch (const po::error &ex) {
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    // Handle the help option
    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 0;
    }

    // Parsing the observables in order of appearance
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--MSD") {
            allObs.push_back("MSD");
        } else if (arg == "--Cb") {
            allObs.push_back("Cb");
        } else if (arg == "--Fs") {
            allObs.push_back("Fs");
        } else if (arg == "--U") {
            allObs.push_back("U");
        }
    }

    // Resizing arrays
    Size = pow(N/density, 1./3.);
    steps = tw*(cycles-1)+tau;
    // mol_index = new int[N];
    // X = new double[N]; Y = new double[N]; Z = new double[N]; 
    // S = new double[N]; Sref = new double[N]; 
    // X0 = new double[N]; Y0 = new double[N]; Z0 = new double[N];
    // Xfull = new double[N]; Yfull = new double[N]; Zfull = new double[N]; 
    // Xref = new double[N]; Yref = new double[N]; Zref = new double[N];

    // creating outdir if not existing
    fs::path out_path = outdir;
    if(!fs::is_directory(out_path)){
        
        fs::create_directory(outdir);
    }
    
     // Writing params.txt file
    std::ofstream params;
    params.open(outdir + "params.txt");
    params << "rootdir" << " " << "N" << " " << "T" << " " 
           << "tau" << " " << "tw" << " " << "cycles" << " " 
           << "logPoints" << " " << "linPoints" << " " << "p_flip" << std::endl;
    params << outdir << " " << N << " " << T << " " << tau << " " << tw << " "
           << cycles << " " << logPoints << " " << linPoints << " " << p_flip << std::endl;
    params.close();

    // Read init config
    cfg = ReadTrimCFG(input);
    cfg.GetBonds(); cfg.UpdateNL();
    // Do simulation with timer
    double t0 = time(NULL); // Timer
    MC(outdir, logPoints, linPoints); 
    std::cout << "Time taken: " << (time(NULL) - t0) << "s" << std::endl; 
    std::cout << "Done" << std::endl;

    // Freeing allocated memory
    // delete[] mol_index;
    // delete[] X; delete[] Y; delete[] Z; delete[] S; delete[] Sref; 
    // delete[] X0; delete[] Y0; delete[] Z0;
    // delete[] Xfull; delete[] Yfull; delete[] Zfull;
    // delete[] Xref; delete[] Yref; delete[] Zref;

    return 0;
}