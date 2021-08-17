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

//Input files.
const string init_genome_file="";
const string init_expression_file="";
const string init_backup_file="";
const string init_anctrace_file="";
const string init_lineage_file="";

//Bead variables
const int sequence_length =	20;
const int signalp_length =	1;	//As long as I am not using it, make it small.

//Grid size and host size, setting the dimensions of the model.
const int NR = 50;
const int NC = 50;	//Gradient is over columns.

//Main settings
const int relative_replication = -1;	// -1 to turn off. If we're doing relative replication (i.e. not -1), how many nutrients are considered to be needed for replication of the entire genome. This removes selection against genome size by scaling replication length with genome length.
const bool gene_replication = false;	//Only genes take time to replicate.
const int moran_symbionts = 5;	// -1 for no Moran process, positive number defines constant number of symbionts inside each host.

//Genome parameters
const int nr_household_genes =			50;
const double leakage_to_host =			0.00;
const double leakage_to_symbiont =	0.00;

//Regulatory parameters
const double k_zero =		0.0000001;
const double epsilon =	1.00;

//Runtime and output parameters
const int default_TimeZero =						00000000;
const int default_SimTime =							10000000;
const int default_TimeTerminalOutput =	100;
const int default_TimeSaveGrid =				100;
const int default_TimePruneFossils =		1000;
const int default_TimeOutputFossils =		10000;
const int default_TimeSaveBackup =			10000;

//Population parameters
const double death_rate_host =								0.001;
const double death_rate_symbiont =						0.001;
const double default_nutrient_abundance =			30.;
const int default_nutrient_competition =			2;
//Options for nutrient_competition:
// 1, classic nutrient function (e.g. Paramecium tetraurelia):		n_ij = ( n_tot - (x_nei-x_i) ) / x_i
// 2, first smooth nutrient function (e.g. Paramecium caudatum):	n_ij = n_tot / x_nei
// 3, distribute by cell, by organelle (e.g. Volvox carteri IV):	n_ij = (n_tot / c_nei) / x_i
// where n_ij is nutrients at site i for organelle j,
//       n_tot is total nutrient_abundance (see par above), i.e. influx per site,
//       x_nei is total number of organelles in the neighbourhood,
//       x_i is total number of organelles at site i,
//       c_nei is number of cells (or hosts) in the neighbourhood.

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
const double regulator_transfer_mu_HtoS = 0.0000;
const double regulator_transfer_mu_StoH = 0.0000;
const double bsite_transfer_mu_HtoS = 		0.0000;
const double bsite_transfer_mu_StoH = 		0.0000;
const double house_transfer_mu_HtoS = 		0.0000;
const double house_transfer_mu_StoH = 		0.0000;

//Variables defined in World.cc
extern int Time;
extern int seed;
extern unsigned long long seed_draws;
extern string folder;

extern string genome_file;
extern string expression_file;
extern string backup_file;
extern string anctrace_file;
extern string lineage_file;

extern int TimeZero;
extern int SimTime;
extern int TimeTerminalOutput;
extern int TimeSaveGrid;
extern int TimePruneFossils;
extern int TimeOutputFossils;
extern int TimeSaveBackup;

extern bool follow_single_individual;
extern bool follow_with_fixed_symbionts;
extern bool trace_lineage;
extern bool log_lineage;
extern bool mutations_on;
extern bool well_mixing;

extern double nutrient_abundance;
extern int nutrient_competition;

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

//The current definition of the stages.
const bool StageTargets[4][5] = {
  true, false, false, true, true,       // 1 0 0 1 1    G1
  false, false, true, false, true,      // 0 0 1 0 1    S
  false, true, false, false, false,     // 0 1 0 0 0    G2
  true, false, false, false, false,     // 1 0 0 0 0    M
                                        // 1-CtrA 2-GcrA 3-DnaA 4-CcrM 5-SciP
};

#endif
