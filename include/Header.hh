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

#define HOST 0
#define SYMBIONT 1

#define HOUSE 0
#define BSITE 1
#define REGULATOR 2
#define EFFECTOR 3

const int max_input_files=20;

//Bead variables
const int signalp_length =		1;	//As long as I am not using it, make it small.
const int regulator_length =	20;
const int effector_length =		10;

//Grid size and host size, setting the dimensions of the model.
const int NR_max = 100;
const int NC_max = 100;	//Gradient is over columns.

//Main settings
const int relative_replication = -1;	// -1 to turn off. If we're doing relative replication (i.e. not -1), how many nutrients are considered to be needed for replication of the entire genome. This removes selection against genome size by scaling replication length with a given genome length.
const bool gene_replication = false;	//Only genes take time to replicate.
const bool moran_symbionts = false;	// If true, hosts evolve with constant symbiont numbers (i.e. Moran process), at the level that they are initialised.
const bool safe_symbiont_distribution = false;	//If true, it means that daughter cell can end up with no less than 1 symbiont (if initial division says 0, then copy one of the symbionts from the other daughter). Symbionts can still be lost by their own fault (basal death and failed division). Only makes sense without Moran process.
const int symbiont_overgrowth = -1;	//The maximal number of symbiont spots. If value is set to -1, there is no max. and no overgrowth (new symbionts are always appended). If value is positive, there is overgrowth (symbiont division always carries a probability to overgrow a colleague, and this probability equals 1 when nr_symbionts==symbiont_overgrowth).
const int host_growth = 0;
//Options for host growth / cell-cycle fitness criterion.
// 0, as in Prokaryotes: hosts overgrow one another, dividing as soon as they reach M.
// 1, hosts wait for empty sites, but waiting is free (expression unimportant once they reach M).
// 2, hosts wait for empty sites, but need to actively maintain M expression (i.e. making the whole cell-cycle more complex again). This option may be implemented later...
const bool perfect_transport = false;	//Genes with a 0 in their signalp, are always moved to the host; genes with a 1 in their signalp always get targetted to the symbionts. Moving here means that they do not stick around in the compartment where they are created.
const bool nutshare_evolve = false;	//Host can evolve how much nutrients it claims from the environment, passing on the remaining fraction to its symbionts (equally divided among these). Each host has an identical claim on environmental nutrients, i.e. independent of how many symbionts it has. This corresponds to nutrient competition 4.
const bool mutation_epochs = false;	//Turn up mutation rate during specified time periods (see parameters below).
const int seq_hdist = 0;	//Maximal hamming distance to still be called this particular gene type. When this value is set to 0, we have the same case as in Prokaryotes.
const int eff_hdist = 0;	//Maximal hamming distance for defining effectors.
const bool empty_division_killing = true;	//Hosts first divide, potentially killing a neighbour, before realising that the new cell does not get symbionts and dies right away.

//Genome parameters
const int nr_household_genes =			50;
const double leakage_to_host =			0.00;
const double leakage_to_symbiont =	0.00;

//Regulation parameters
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
// 0, constant nutrient level, unaffected by cells.
// 1, classic nutrient function (e.g. Paramecium tetraurelia):		n_ij = ( n_tot - (x_nei-x_i) ) / x_i
// 2, first smooth nutrient function (e.g. Paramecium caudatum):	n_ij = n_tot / x_nei
// 3, distribute by cell, by organelle (e.g. Volvox carteri IV):	n_ij = (n_tot / c_nei) / x_i
// 4, distribute by cell, let host claim its share and then divide among symbionts:	n_iH = (n_tot / c_nei)*claim_H and n_iS = (n_tot / c_nei)*(1 - claim_H) / (x_i - 1)
// 5, distribute by cell, let host give a fraction nutrient_claim to each of its symbionts, itself taking whatever's left: n_iH = max(0, 1 - (x_i-1)*claim_H) * (n_tot / c_nei) and n_iS = min(1, (x_i-1)*claim_H)/(x_i-1) * (n_tot / c_nei)
// where n_ij is nutrients at site i for organelle j,
//       n_tot is total nutrient_abundance (see par above), i.e. influx per site,
//       x_nei is total number of organelles in the neighbourhood,
//       x_i is total number of organelles at site i,
//       c_nei is number of cells (or hosts) in the neighbourhood,
//       claim_H is the fraction of nutrients claimed by the host over its symbionts,
//       n_iH is the number of nutrients ending up in the host.
//       n_iS is the number of nutrients ending up in each of the symbionts.
const int default_strain_competition = 1;
//Options for strain competition (i.e. initial distribution):
// 1, each strain in its own sector (vertical stripes).
// 2, all strains mixed (i.e. each site has equal probability to be any of the strains).
// 3, field divided into blocks (both horizontal and vertical division, assumed to start with a square number of strains).

//INVASION parameters.
const int equilibration_time = 1000;
const int add_finish_time = 1000;

/* MUTATION PARAMETERS */

// Used by nutshare_evolve option.
const double init_nutrient_claim =				0.1;	//If we start with 4 symbionts and 1 hosts, that means they initially share fairly.
const double nutrient_claim_mu =					0.001;
const double nutrient_claim_mu_delta =		0.05;

const double high_mu_factor =							3;
const int high_mu_interval =							40000;
const int high_mu_period =								10000;

const int WeightRange = 3;  //Weights range from -WeightRange to +WeightRange.

//Mutation rates are specified like this:
// mu[HOST][DUPLICATION][BSITE]
// mut[HOST][EFFECTOR] specifies the transfer mutation rate FROM HOST (to symbiont) for effectors.

#define DUPLICATION 0
#define DELETION 1
#define SHUFFLE 2
#define INVENTION 3
#define THRESHOLD 4
#define SIGNALP 5
#define SEQUENCE 6
#define ACTIVITY 7

extern double mu[2][8][4];
extern double muT[2][4];

//Effector definitions. Previously used to make hard-coded regulatory types.
const bool effector_types[5][effector_length] =
{
	true, false, true, false, true, false, true, false, true, false,
	false, false, true, true, false, false, true, true, false, false,
	false, false, false, true, true, true, false, false, false, true,
	true, true, true, true, false, false, false, false, true, true,
	false, false, false, false, false, true, true, true, true, true
};

//The current definition of the stages.
const bool StageTargets[4][5] = {
  true, false, false, true, true,       // 1 0 0 1 1    G1
  false, false, true, false, true,      // 0 0 1 0 1    S
  false, true, false, false, false,     // 0 1 0 0 0    G2
  true, false, false, false, false,     // 1 0 0 0 0    M
                                        // 1-CtrA 2-GcrA 3-DnaA 4-CcrM 5-SciP
};

//Variables defined in World.cc
extern int Time;
extern int seed;
extern string folder;
extern int NR;
extern int NC;
extern int NCfull;

extern string genome_files[max_input_files];
extern string expression_files[max_input_files];
extern string definition_files[max_input_files];
extern string mutation_file;
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

extern bool invasion_experiment;	//If set, do an invasion experiment, and stop simulation when the first cell hits the right-most column.
extern int invasion_complete;
extern bool follow_single_individual;
extern bool follow_with_fixed_symbionts;
extern bool trace_lineage;
extern bool log_lineage;
extern bool mutations_on;
extern bool well_mixing;

extern int init_stage;

extern double nutrient_abundance;
extern int nutrient_competition;
extern int strain_competition;

extern dsfmt_t dsfmt;
inline double uniform() { return dsfmt_genrand_close_open(&dsfmt); }
inline int uniform_shuffle (int i) { return (int)(RAND_MAX*uniform()) % i; }


#endif
