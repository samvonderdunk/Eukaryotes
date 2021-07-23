#include "Header.hh"
#include "Population.hh"
#include "Cell.hh"
#include "Organelle.hh"
#include "Genome.hh"
#include "Regulator.hh"
#include "Bsite.hh"
#include "House.hh"
#include "Bead.hh"
#include "dSFMT.h"

dsfmt_t dsfmt;
int Time;
int seed = time(0);	//Unless changed at command line.
unsigned long long seed_draws = 0;	//May be obsolete in the future, but still needed for old backups (before 23/07/2021).

string folder = "/linuxhome/tmp/sam/Eukaryotes/";
string genome_file = init_genome_file;
string expression_file = init_expression_file;
string backup_file = init_backup_file;
string anctrace_file = init_anctrace_file;
string lineage_file = init_lineage_file;

int TimeZero = default_TimeZero;
int SimTime = default_SimTime;
int TimeTerminalOutput = default_TimeTerminalOutput;
int TimeSaveGrid = default_TimeSaveGrid;
int TimePruneFossils = default_TimePruneFossils;
int TimeOutputFossils = default_TimeOutputFossils;
int TimeSaveBackup = default_TimeSaveBackup;

bool mutations_on = true;
bool follow_single_individual = false;
bool follow_with_fixed_symbionts = false;
bool trace_lineage = false;
bool log_lineage = false;

double nutrient_abundance = default_nutrient_abundance;
int nutrient_competition = default_nutrient_competition;

void Setup(int argc, char** argv);

int main(int argc, char** argv) {

	/* ############## Setup ############## */
	printf("\n\033[93m### Setup ###\033[0m\n");
	Population* P;
	for(int q=0; q<argc; q++)	printf("%s ", argv[q]);
	Setup(argc, argv);
	dsfmt_init_gen_rand(&dsfmt, seed);	//Used to seed uniform().
	srand(seed);												//Used to seed random_shuffle(...).
	printf("\b\nSetup completed...\n\n");

	if (follow_single_individual)
	{
		printf("\033[93m### Start ###\033[0m\n");
		P = new Population();
		P->FollowSingleCell();
	}
	else
	{
		/* ############## Initialisation ############## */
		printf("\033[93m### Initialisation ###\033[0m\n");
		P = new Population();
		if(backup_file != "")	P->ContinuePopulationFromBackup();
		else	P->InitialisePopulation();
		if(lineage_file != "")	P->ReadLineageFile();
		printf("Initialisation completed...\n\n");

		/* ############## Simulation ############## */

		printf("\033[93m### Simulation ###\033[0m\n");
		for(Time=TimeZero; Time<SimTime+1; Time++){	//We do one extra step, because output is generated at the beginning of a step, such that time=0 is the field as it is initialised.
			P->UpdatePopulation();		//Main next-state function, updating the population.
		}
		printf("Simulation completed...\n\n");
	}

	/* ############## End ############## */

	printf("\033[93m### End ###\033[0m\n");
	delete P;
	P = NULL;
	printf("Eukaryotes completed...\n\n");

}

void PrintUsage(bool full)
{
	printf("\n\033[93m### Eukaryotes --- usage ###\033[0m\nArguments:\n   -s [seed]\t\t\tSet seed for random number generator (e.g. 211)\n   -p [project title]\t\tDefines folder for local storage\n   -g [genomes file]\t\tSee World.cc for format (e.g. CellX.g)\n   -e [expressions file]\tSee World.cc for format (e.g. CellX_G1.g)\n   -b [backup file]\t\tStart from backup (e.g. /path/backup00090000.txt)\n   -a [ancestor file]\t\tContinue ancestor trace (e.g. /path/anctrace00090000.txt)\n   -l [lineage file]\t\tLineage record (e.g. Host_MRCA_t1960k.out, see Programmes -L1 and -L2)\n   -t0 [start time]\t\tSet starting time (default: 0)\n   -tN [end time]\t\tSet simulation time (default: 10M)\n   -tT [term time]\t\tSet interval for terminal output (default: 100)\n   -tS [snap time]\t\tSet interval for saving snapshot (default: 100)\n   -tP [prune time]\t\tSet time interval for fossil pruning (default: 1000)\n   -tF [fossil time]\t\tSet time interval for saving fossil record (default: 10k)\n   -tB [backup time]\t\tSet interval for saving backup (default: 10k)\n   -nA [nutrient abundance]\tSet nutrient influx per site\n   -nC [nutrient competition]\tChoose type of nutrient competition (1: subtract/divide, 2: divide all, 3: divide/divide)\n\nFlags:\n   --nomut\t\t\tNo mutations\n   --help\t\t\tPrint full usage info\n\nProgrammes:\n   -S1\t\t\t\tFollow single cell with growing symbionts\n   -S2\t\t\t\tFollow single with fixed symbiont numbers\n   -L1\t\t\t\tTrace complete lineage (every organelle is considered a mutant)\n   -L2\t\t\t\tLog complete lineage (obtained using -L1)\n");
	if (full)
	{
		printf("\n\033[93m### Eukaryotes --- formats ###\033[0m\n\n<genomes file>\t\tHost genome on first line, each next line a symbiont.\n   (G2:0:-3:1:10010100010101010001).(H).(0:01010101000110010100).(...\n   (G4:1:-1:2:01010101100101000001).(H).(2:10100011001010001010).(...\n   (G4:1:-1:2:01010101100101000001).(H).(2:10100011001010001010).(...\n\n<expressions file>\tHost expression on first line, each next line a symbiont (matching with genomes file)\n   {10101110}\n   {...\n   {...\n");
	}
	exit(1);
}

void Setup(int argc, char** argv) {

	string ReadOut, command;
	bool project_name_found = false;

	for(int i=1;i<argc;i++)	//Loop through input arguments.
	{
		ReadOut = (char*) argv[i];	//There does not seem to be a quicker way to compare the input arguments with a string.

		/* ######### */
		/* ARGUMENTS */
		/* ######### */

		if(ReadOut=="-s" && (i+1)!=argc)
		{
			seed = atoi(argv[i+1]);
			i++;
			continue;
		}

		else if(ReadOut=="-p" && (i+1)!=argc)
		{
			folder += argv[i+1];
			project_name_found = true;
			i++;
			continue;
		}

		else if(ReadOut=="-g" && (i+1)!=argc)
		{
			genome_file = argv[i+1];
			i++;
			continue;
		}

		else if(ReadOut=="-e" && (i+1)!=argc)
		{
			expression_file = argv[i+1];
			i++;
			continue;
		}

		else if(ReadOut=="-b" && (i+1)!=argc)
		{
			backup_file = argv[i+1];
			i++;
			continue;
		}

		else if(ReadOut=="-a" && (i+1)!=argc)
		{
			anctrace_file = argv[i+1];
			i++;
			continue;
		}

		else if(ReadOut=="-l" && (i+1)!=argc)
		{
			lineage_file = argv[i+1];
			i++;
			continue;
		}

		else if(ReadOut=="-t0" && (i+1)!=argc)
		{
			TimeZero = atoi(argv[i+1]);
			i++;
			continue;
		}

		else if(ReadOut=="-tN" && (i+1)!=argc)
		{
			SimTime = atoi(argv[i+1]);
			i++;
			continue;
		}

		else if(ReadOut=="-tT" && (i+1)!=argc)
		{
			TimeTerminalOutput = atoi(argv[i+1]);
			i++;
			continue;
		}

		else if(ReadOut=="-tS" && (i+1)!=argc)
		{
			TimeSaveGrid = atoi(argv[i+1]);
			i++;
			continue;
		}

		else if(ReadOut=="-tP" && (i+1)!=argc)
		{
			TimePruneFossils = atoi(argv[i+1]);
			i++;
			continue;
		}

		else if(ReadOut=="-tF" && (i+1)!=argc)
		{
			TimeOutputFossils = atoi(argv[i+1]);
			i++;
			continue;
		}

		else if(ReadOut=="-tB" && (i+1)!=argc)
		{
			TimeSaveBackup = atoi(argv[i+1]);
			i++;
			continue;
		}

		else if(ReadOut=="-nA" && (i+1)!=argc)
		{
			nutrient_abundance = atof(argv[i+1]);
			i++;
			continue;
		}

		else if(ReadOut=="-nC" && (i+1)!=argc)
		{
			nutrient_competition = atoi(argv[i+1]);
			i++;
			continue;
		}

		/* ##### */
		/* FLAGS */
		/* ##### */

		else if(ReadOut=="--nomut")
		{
			mutations_on = false;
		}

		else if(ReadOut=="--help")
		{
			PrintUsage(true);
		}

		/* ########## */
		/* PROGRAMMES */
		/* ########## */

		else if(ReadOut=="-S1")
		{
			follow_single_individual = true;
			mutations_on = false;
			follow_with_fixed_symbionts = false;
		}
		else if(ReadOut=="-S2")
		{
			follow_single_individual = true;
			mutations_on = false;
			follow_with_fixed_symbionts = true;
		}
		else if(ReadOut=="-L1")
		{
			trace_lineage = true;
		}
		else if(ReadOut=="-L2")
		{
			log_lineage = true;
		}

		else
		{
			PrintUsage(false);
		}
	}

	if (!project_name_found)	folder += "Project_Name";	//I did not manage to give the date as an extension to the folder.

	printf("\n\nFolder:\t\t\t%s\nGenomes:\t\t%s\nExpression:\t\t%s\nBackup:\t\t\t%s\nAnctrace:\t\t%s\nLineage:\t\t%s\nStart time:\t\t%d\nEnd time:\t\t%d\nt-Terminal:\t\t%d\nt-Snap:\t\t\t%d\nt-Prune:\t\t%d\nt-Ancestry:\t\t%d\nt-Backup:\t\t%d\nNutrient abundance:\t%f\nNutrient comp.:\t\t%d\n\nMutations:\t\t%s\nFollow var. host:\t%s\nFollow fixed host:\t%s\nTrace lineage:\t\t%s\nLog lineage:\t\t%s\n", folder.c_str(), genome_file.c_str(), expression_file.c_str(), backup_file.c_str(), anctrace_file.c_str(), lineage_file.c_str(), TimeZero, SimTime, TimeTerminalOutput, TimeSaveGrid, TimePruneFossils, TimeOutputFossils, TimeSaveBackup, nutrient_abundance, nutrient_competition, mutations_on?"Yes":"No", follow_single_individual?"Yes":"No", follow_with_fixed_symbionts?"Yes":"No", trace_lineage?"Yes":"No", log_lineage?"Yes":"No");


	//Set up all directories for data.
	command = "mkdir -p " + folder;
	system(command.c_str());
	command = "mkdir -p " + folder + "/snapsamples";
	system(command.c_str());
	command = "mkdir -p " + folder + "/backups";
	system(command.c_str());
	command = "mkdir -p " + folder + "/ancestors";
	system(command.c_str());

}
