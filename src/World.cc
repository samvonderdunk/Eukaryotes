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
int initial_seed = time(0);
unsigned long long seed_draws = 0;
string folder = "/linuxhome/tmp/sam/Eukaryotes/";
string genome_initialisation = genome_file;
string expression_initialisation = expression_file;
string backup_reboot = backup_file;
string anctrace_reboot = anctrace_file;
int SimTime = default_SimTime;
bool mutations_on = true;
bool follow_single_individual = false;

void Setup(int argc, char** argv);

int main(int argc, char** argv) {

	/* ############## Setup ############## */
	printf("\n\033[93m### Setup ###\033[0m\n");
	Population* P;
	Setup(argc, argv);
	printf("Function call: ");
	for(int q=0; q<argc; q++)	printf("%s ", argv[q]);
	dsfmt_init_gen_rand(&dsfmt, initial_seed);	//Used to seed uniform().
	srand(initial_seed);	//Used to seed random_shuffle(...).
	printf("\b\nSetup completed...\n\n");

	if (follow_single_individual)
	{
		printf("\033[93m### Start ###\033[0m\n");
		// P = new Population();
		// P->FollowSingleIndividual();
	}
	else
	{
		/* ############## Initialisation ############## */
		printf("\033[93m### Initialisation ###\033[0m\n");
		P = new Population();
		if(backup_reboot != "")	P->ContinuePopulationFromBackup();
		else	P->InitialisePopulation();
		printf("Initialisation completed...\n\n");

		/* ############## Simulation ############## */

		printf("\033[93m### Simulation ###\033[0m\n");
		for(Time=TimeZero; Time<SimTime+1; Time++){	//We do one extra step, because output is generated at the beginning of a step, such that time=0 is the field as it is initialised.
			P->UpdatePopulation();		//Main next-state function, updating the population.
		}
		//Make sure that you save all possible things in the last timestep, if you did not already choose your parameters such.
		// Time--;
		printf("Simulation completed...\n\n");
	}

	/* ############## End ############## */

	printf("\033[93m### End ###\033[0m\n");
	delete P;
	P = NULL;
	printf("Eukaryotes completed...\n\n");

}



void Setup(int argc, char** argv) {

	string ReadOut, command;
	bool project_name_found = false;
	bool initial_seed_set = false;

	for(int i=1;i<argc;i++)	//Loop through input arguments.
	{
		ReadOut = (char*) argv[i];	//There does not seem to be a quicker way to compare the input arguments with a string.

		/* ############## */
		/* INPUT SETTINGS */
		/* ############## */

		if(ReadOut=="-s" && (i+1)!=argc)
		{
			initial_seed = atoi(argv[i+1]);
			initial_seed_set = true;
			printf("Seed = %i\n", initial_seed);
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

		/* Format of genome input file:
			HOST::(G2:3:-2:11101010101... etc...)
			SYMB::(G4:1:1:0001010100101... etc..)
			SYMB::(G1:0:-1:000100101010... etc..)
			etc..
		*/

		else if(ReadOut=="-g" && (i+1)!=argc)
		{
			genome_initialisation = argv[i+1];
			printf("Genome input: %s\n", genome_initialisation.c_str());
			i++;
			continue;
		}

		/* Format of expression input file:
			HOST::[1,0,1,0,0,0,1,1,0]
			SYMB::[...]
			SYMB::[...]
			etc..
		*/

		else if(ReadOut=="-e" && (i+1)!=argc)
		{
			expression_initialisation = argv[i+1];
			printf("Expression input: %s\n", expression_initialisation.c_str());
			i++;
			continue;
		}

		else if(ReadOut=="-b" && (i+1)!=argc)
		{
			backup_reboot = argv[i+1];
			printf("Backup-file input: %s\n", backup_reboot.c_str());
			i++;
			continue;
		}

		else if(ReadOut=="-a" && (i+1)!=argc)
		{
			anctrace_reboot = argv[i+1];
			printf("Anctrace input: %s\n", anctrace_reboot.c_str());
			i++;
			continue;
		}

		else if(ReadOut=="-t" && (i+1)!=argc)
		{
			SimTime = atoi(argv[i+1]);
			printf("Simulation time: %d\n", SimTime);
			i++;
			continue;
		}

		/* ########### */
		/* RUN OPTIONS */
		/* ########### */

		else if(ReadOut=="-nomut")
		{
			mutations_on = false;
			printf("Simulating without mutations.\n");
		}

		//Follow a single immortal individual for many time steps.
		else if(ReadOut=="-S")
		{
			follow_single_individual = true;
			printf("Following a single, immortal individual through time.\n");
		}

		else	//Print usage/help.
		{
			printf("\nWARNING: OLD HELP MESSAGE\033[93m### Prokaryotes --- usage ###\033[0m\nArgument options:\n   -p [project title]\t\tDefines folder for local storage\n   -s [seed]\t\t\tSet seed for random number generator (e.g. 211)\n   -i [genome initialisation file]\t\te.g. CellX.g\n   -g [expression initialisation file]\te.g. CellX_Expr.g\n   -e [env]\t\t\tInitial environment (e.g. -3)\n   -b [backup file]\t\tStart from backup (e.g. /path/backup00090000.txt)\n   -a [ancestor file]\t\tContinue ancestor trace (e.g. /path/anctrace00090000.txt)\n   -nomut\t\t\tNo mutations\n   -t [max. time]\t\tSet simulation time (e.g. 100)\n Programmes:\n   -M [nr_mutants]\tGenerate mutants\n   -MS\t\t\tScan neutral mutational path\n   -A [nr_states]\tSimulate state-space transitions [nr. of initial states = max(nr_states, total nr. unique states)]\n   -S\t\t\tFollow single immortal individual/lineage through time [simulating until Time==SimTime]\n");
			exit(1);
		}
	}

	if (!initial_seed_set)	printf("Seed = %li\n", time(0));
	if (!project_name_found)	folder += "Project_Name";	//I did not manage to give the date as an extension to the folder.

	//Set up all directories for data.
	command = "mkdir -p " + folder;
	system(command.c_str());
	printf("Folder = %s\n", folder.c_str());
	command = "mkdir -p " + folder + "/snapsamples";
	system(command.c_str());
	command = "mkdir -p " + folder + "/backups";
	system(command.c_str());
	command = "mkdir -p " + folder + "/ancestors";
	system(command.c_str());

}
