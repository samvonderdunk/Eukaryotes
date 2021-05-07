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

//Bead variables
const int sequence_length = 20;

//Grid size and host size, setting the dimensions of the model.
const int NR=50;
const int NC=550;	//Gradient is over columns.
const int HS=100;	//Host size, i.e. maximal number of symbionts inside a cell beside the host.

//Genome parameters
const int WeightRange = 3;  //Weights range from -WeightRange to +WeightRange.
const int nr_household_genes = 50;
const double leakage_to_host = 0.01;
const double leakage_to_symbiont = 0.01;

//Mutation parameters
const double gene_threshold_mu = 0.0005;
const double gene_activity_mu = 0.0005;
const double gene_binding_domain_mu = 0.0001;

const double tfbs_binding_site_mu = 0.0001;
const double tfbs_activity_mu = 0.0005;

const double gene_duplication_mu = 0.0005;
const double gene_deletion_mu = 0.0005;
const double gene_innovation_mu = 0.0005; //I set this 10x lower than other mutation rates on purpose.
const double gene_shuffle_mu = 0.001;

const double tfbs_duplication_mu = 0.0005;
const double tfbs_deletion_mu = 0.0005;
const double tfbs_innovation_mu = 0.005;
const double tfbs_shuffle_mu = 0.001;

const double house_duplication_mu = 0.0001;
const double house_deletion_mu = 0.0001;
const double house_shuffle_mu = 0.001;

//Regulatory parameters
const double k_zero = 0.0000001;
const double epsilon = 1.00;

//Runtime and output parameters
const int TimeZero=0;
const int default_SimTime=1000000;

//Population parameters
const double diffusion_rate = 0.;  // >1: multiple Margolus steps per time step, <1: probability of single Margolus step each time step.
const double death_rate = 0.001;
const double repl_rate = 1.0;
const int repl_scale = 3;  //Neighbourhood for replication, standard=3x3.

//Variables defined in World.cc
extern int Time;
extern int initial_seed;
extern unsigned long long seed_draws;
extern string folder;
extern string genome_init;
extern string genestate_init;
extern string backup_reboot;
extern string anctrace_reboot;
extern int SimTime;

extern dsfmt_t dsfmt;
inline double uniform()
{
  seed_draws ++;
  return dsfmt_genrand_close_open(&dsfmt);
}

const string genome_file="";
const string genestate_file="";
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
