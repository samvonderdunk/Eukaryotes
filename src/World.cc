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
	printf("\n\033[93m### Eukaryotes --- usage ###\033[0m\nArguments:\n   -s [seed]\t\t\tSet seed for random number generator (e.g. 211)\n   -p [project title]\t\tDefines folder for local storage\n   -g [genomes file]\t\tSee World.cc for format (e.g. CellX.g)\n   -e [expressions file]\tSee World.cc for format (e.g. CellX_G1.g)\n   -b [backup file]\t\tStart from backup (e.g. /path/backup00090000.txt)\n   -a [ancestor file]\t\tContinue ancestor trace (e.g. /path/anctrace00090000.txt)\n   -t [max. time]\t\tSet simulation time (e.g. 100)\n\nFlags:\n   --nomut\t\t\tNo mutations\n   --help\t\t\tPrint full usage info\n\nProgrammes:\n   -S\t\t\t\tFollow single immortal cell/lineage through time (simulating until Time==SimTime)\n");
	if (full)
	{
		printf("\n\033[93m### Eukaryotes --- formats ###\033[0m\n\n<genomes file>\t\tHost genome on first line, each next line a symbiont.\n   (G2:0:-3:1:10010100010101010001).(H).(0:01010101000110010100).(...\n   (G4:1:-1:2:01010101100101000001).(H).(2:10100011001010001010).(...\n   (G4:1:-1:2:01010101100101000001).(H).(2:10100011001010001010).(...\n\n<expressions file>\tHost expression on first line, each next line a symbiont (matching with genomes file)\n   {10101110}\n   {...\n   {...\n");
	}
	exit(1);
}

void Setup(int argc, char** argv) {

	string ReadOut, command;
	bool project_name_found = false;
	bool initial_seed_set = false;

	for(int i=1;i<argc;i++)	//Loop through input arguments.
	{
		ReadOut = (char*) argv[i];	//There does not seem to be a quicker way to compare the input arguments with a string.

		/* ######### */
		/* ARGUMENTS */
		/* ######### */

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

		else if(ReadOut=="-g" && (i+1)!=argc)
		{
			genome_initialisation = argv[i+1];
			printf("Genome input: %s\n", genome_initialisation.c_str());
			i++;
			continue;
		}

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

		/* ##### */
		/* FLAGS */
		/* ##### */

		else if(ReadOut=="--nomut")
		{
			mutations_on = false;
			printf("Simulating without mutations.\n");
		}

		else if(ReadOut=="--help")
		{
			PrintUsage(true);
		}

		/* ########## */
		/* PROGRAMMES */
		/* ########## */

		else if(ReadOut=="-S")	//Follow a single immortal individual for many time steps.
		{
			follow_single_individual = true;
			printf("Following a single, immortal individual through time.\n");
		}

		else
		{
			PrintUsage(false);
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
