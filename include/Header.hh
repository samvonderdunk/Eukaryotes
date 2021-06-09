#ifndef GeneralHeader
#define GeneralHeader

using namespace std;
//using namespace __gnu_cxx;

#include <stdio.h>
#include <stdlib.h>
#include <ext/numeric>
#include <cmath>
#include <list>
#include <vector>
#include <set>
#include <iostream>
#include <iomanip>
#include <map>
#include <iterator>
#include <algorithm>
#include <boost/next_prior.hpp>
#include <cstring>
#include <sstream>
#include <fstream>
#include "dSFMT.h"	// Kirsten used it for some functions, so lets just stick with it
//#include "Xmgrace.hh"	// Grace pipeline added by Brem
#include <sys/time.h>
#include <ctime>
#include <typeinfo>
#include <string>

// #include "/home/sam/Programmes/nvwa-1.1/nvwa/debug_new.h"

#define toDigit(c) (c-'0')  // Converts char to digit
#define REGULATOR 0
#define BSITE 1
#define HOUSE 2

//Bead variables
const int sequence_length =	20;
const int signalp_length =	1;	//As long as I am not using it, make it small.

//Grid size and host size, setting the dimensions of the model.
const int NR = 50;
const int NC = 50;	//Gradient is over columns.

//Main settings
const bool relative_replication = false;	//Remove selection against genome size by scaling replication length with genome length.
const int rel_repl_full = 65;	//If we're doing relative replication, how many nutrients are considered to be needed for replication of the entire genome.
const bool gene_replication = false;	//Only genes take time to replicate.

//Genome parameters
const int nr_household_genes =			50;
const double leakage_to_host =			0.01;
const double leakage_to_symbiont =	0.01;

//Regulatory parameters
const double k_zero =		0.0000001;
const double epsilon =	1.00;

//Runtime and output parameters
const int TimeZero =						0;
const int default_SimTime =			1000000;
const int TimeTerminalOutput =	1;
const int TimeSaveGrid =				100;
const int TimePruneFossils =		100;
const int TimeOutputFossils =		100;
const int TimeSaveBackup =			100;

//Population parameters
const double death_rate_host = 				0.001;
const double death_rate_symbiont = 		0.001;
const double nutrient_abundance = 		100.;
const double max_organelle_density = 	100.;	//Only used with the relative nutrient function; k, number of organelles at which point nutrients will be completely depleted.

/* MUTATION PARAMETERS */
const int WeightRange = 3;  //Weights range from -WeightRange to +WeightRange.

const double regulator_threshold_mu = 		0.0005;
const double regulator_activity_mu = 			0.0005;
const double regulator_sequence_mu = 			0.0001;

const double bsite_sequence_mu = 					0.0001;
const double bsite_activity_mu = 					0.0005;

const double regulator_duplication_mu = 	0.0005;
const double regulator_deletion_mu = 			0.0005;
const double regulator_innovation_mu = 		0.0005;
const double regulator_shuffle_mu = 			0.001;

const double bsite_duplication_mu = 			0.0005;
const double bsite_deletion_mu = 					0.0005;
const double bsite_innovation_mu = 				0.005;
const double bsite_shuffle_mu = 					0.001;

const double house_duplication_mu = 			0.0001;
const double house_deletion_mu = 					0.0001;
const double house_innovation_mu = 				0.0;
const double house_shuffle_mu = 					0.001;

//Transfer mutations.
const double regulator_transfer_mu_HtoS = 0.0001;
const double regulator_transfer_mu_StoH = 0.0001;
const double bsite_transfer_mu_HtoS = 		0.0001;
const double bsite_transfer_mu_StoH = 		0.0001;
const double house_transfer_mu_HtoS = 		0.0001;
const double house_transfer_mu_StoH = 		0.0001;

//Variables defined in World.cc
extern int Time;
extern int initial_seed;
extern unsigned long long seed_draws;
extern string folder;
extern string genome_initialisation;
extern string expression_initialisation;
extern string backup_reboot;
extern string anctrace_reboot;
extern int SimTime;
extern bool mutations_on;

extern dsfmt_t dsfmt;
inline double uniform()
{
  seed_draws ++;
  return dsfmt_genrand_close_open(&dsfmt);
}

inline int uniform_shuffle (int i)
{
	return (int)(RAND_MAX*uniform()) % i;
}

const string genome_file="";
const string expression_file="";
const string backup_file="";
const string anctrace_file="";

//The current definition of the stages.
const bool StageTargets[4][5] = {
  true, false, false, true, true,       // 1 0 0 1 1    G1
  false, false, true, false, true,      // 0 0 1 0 1    S
  false, true, false, false, false,     // 0 1 0 0 0    G2
  true, false, false, false, false,     // 1 0 0 0 0    M
                                        // 1-CtrA 2-GcrA 3-DnaA 4-CcrM 5-SciP
};

#endif
