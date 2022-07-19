#include "Header.hh"
#include "Population.hh"
#include "Cell.hh"
#include "Organelle.hh"
#include "Genome.hh"
#include "Gene.hh"
#include "Bsite.hh"
#include "House.hh"
#include "Bead.hh"
#include "dSFMT.h"

dsfmt_t dsfmt;
int Time;
int seed = time(0);	//Unless changed at command line.

string folder = "/linuxhome/tmp/sam/Eukaryotes/";
string genome_files[max_input_files];	//Apparently strings don't need to be explicitly set to "".
string expression_files[max_input_files];
string definition_files[max_input_files];
string mutation_file;
string backup_file;
string anctrace_file;

int TimeZero = default_TimeZero;
int SimTime = default_SimTime;
int TimeTerminalOutput = default_TimeTerminalOutput;
int TimeSaveGrid = default_TimeSaveGrid;
int TimePruneFossils = default_TimePruneFossils;
int TimeOutputFossils = default_TimeOutputFossils;
int TimeSaveBackup = default_TimeSaveBackup;
int NR = NR_max;
int NC = NC_max;
int NCfull = NC_max;

bool mutations_on = true;
bool well_mixing = false;

int init_stage = 0;

double nutrient_condition[nr_sectors] = {0.};
int nutrient_competition = default_nutrient_competition;
int strain_competition = default_strain_competition;

double mu[8][4] = {0.};
double muWGD = 0.;

void Setup(int argc, char** argv);
void SetMutationRates();

int main(int argc, char** argv) {

	int q;
	/* ############## Setup ############## */
	printf("\n\033[93m### Setup ###\033[0m\n");
	Population* P;
	for(q=0; q<argc; q++)	printf("%s ", argv[q]);
	Setup(argc, argv);
	SetMutationRates();
	dsfmt_init_gen_rand(&dsfmt, seed);	//Used to seed uniform().
	srand(seed);												//Used to seed random_shuffle(...).
	printf("\b\nSetup completed...\n\n");

	/* ############## Initialisation ############## */
	printf("\033[93m### Initialisation ###\033[0m\n");
	P = new Population();
	if(backup_file != "")	P->ContinuePopulationFromBackup();
	else	P->InitialisePopulation();
	printf("Initialisation completed...\n\n");

	/* ############## Simulation ############## */

	printf("\033[93m### Simulation ###\033[0m\n");
	for(Time=TimeZero; Time<SimTime+1; Time++){	//We do one extra step, because output is generated at the beginning of a step, such that time=0 is the field as it is initialised.
		P->UpdatePopulation();		//Main next-state function, updating the population.
	}
	printf("Simulation completed...\n\n");

	/* ############## End ############## */

	printf("\033[93m### End ###\033[0m\n");
	delete P;
	P = NULL;
	printf("Eukaryotes completed...\n\n");

}

void PrintUsage(bool full)
{
	printf("\n\033[93m### TFEvolution --- usage ###\033[0m\nArguments:\n   -s [seed]\t\t\tSet seed for random number generator (e.g. 211)\n   -p [project title]\t\tDefines folder for local storage\n   -g [genomes file]\t\tSee World.cc or --help for format\n   -e [expressions file]\tSee World.cc or --help for format\n   -d [definitions file]\tSee World.cc or --help for format\n   -m [mutations file]\tSee World.cc or --help for format\n   -b [backup file]\t\tStart from backup (e.g. /path/backup00090000.txt)\n   -a [ancestor file]\t\tContinue ancestor trace (e.g. /path/anctrace00090000.txt)\n   -t0 [start time]\t\tSet starting time (default: 0)\n   -tN [end time]\t\tSet simulation time (default: 10M)\n   -tT [term time]\t\tSet interval for terminal output (default: 100)\n   -tS [snap time]\t\tSet interval for saving snapshot (default: 100)\n   -tP [prune time]\t\tSet time interval for fossil pruning (default: 1000)\n   -tF [fossil time]\t\tSet time interval for saving fossil record (default: 10k)\n   -tB [backup time]\t\tSet interval for saving backup (default: 10k)\n   -nA [nutrient conditions 80,70,60]\tSet nutrient influx per site\n   -nC [nutrient competition]\tChoose type of nutrient competition (1: subtract/divide, 2: divide all, 3: divide/divide)\n   -sC [strain competition]\tChoose distribution of strains (1: vertical stripes, 2: mixed)\n   -r [number of rows]\n   -c [number of columns]\n   -st [initial stage]\n\nFlags:\n   --nomut\t\t\tNo mutations\n   --mixed\t\t\tWell-mixing\n   --help\t\t\tPrint full usage info\n\nProgrammes:\n   [empty]\n");
	if (full)
	{
		printf("\n\033[93m### TFEvolution --- formats ###\033[0m\n\n<genomes file>\t\tHost genome on first line, each next line a symbiont.\n   (R2:0:-3:1:10010100010101010001).(H).(0:01010101000110010100).(...\n   (R4:1:-1:2:01010101100101000001).(H).(2:10100011001010001010).(...\n   (R4:1:-1:2:01010101100101000001).(H).(2:10100011001010001010).(...\n\n<expressions file>\tHost expression on first line, each next line a symbiont (matching with genomes file)\n   {10101110}\n   {...\n   {...\n\n<definitions file>\tHost type definitions on first line, each next line a symbiont\n   (R1:1:10010100010101010001);(R2:-5:01010100101010111000);(...\n   (R1:2:10100011001010001010);(...\n   (R1:2:10100011001010001010);(...\n<mutations file> Per line a different type of mutation; then 2 sets of 4 columns for the 4 bead types\n   #DUPLICATION 0.0001 0.0003 0.0001 0.0001\n   #DELETION 0.0001 0.0003 0.0001 0.0001\n");
	}
	exit(1);
}


void PrintLog()
{
	int i;

	printf("\n\nFolder:\t\t\t%s\nGenomes:\t\t", folder.c_str());
	for (i=0; i<max_input_files; i++)
	{
		if (genome_files[i] != "")
		{
			if (i > 0)	printf(", ");
			printf("%s", genome_files[i].c_str());
		}
	}
	printf("\nExpression:\t\t");
	for (i=0; i<max_input_files; i++)
	{
		if (expression_files[i] != "")
		{
			if (i > 0)	printf(", ");
			printf("%s", expression_files[i].c_str());
		}
	}
	printf("\nDefinition:\t\t");
	for (i=0; i<max_input_files; i++)
	{
		if (definition_files[i] != "")
		{
			if (i > 0)	printf(", ");
			printf("%s", definition_files[i].c_str());
		}
	}
	printf("\nMutations:\t\t%s\nBackup:\t\t\t%s\nAnctrace:\t\t%s\nStart time:\t\t%d\nEnd time:\t\t%d\nt-Terminal:\t\t%d\nt-Snap:\t\t\t%d\nt-Prune:\t\t%d\nt-Ancestry:\t\t%d\nt-Backup:\t\t%d\n", mutation_file.c_str(), backup_file.c_str(), anctrace_file.c_str(), TimeZero, SimTime, TimeTerminalOutput, TimeSaveGrid, TimePruneFossils, TimeOutputFossils, TimeSaveBackup);
	printf("Nutrient conditions:\t");
	for (i=0; i<nr_sectors; i++)
	{
		printf("%f", nutrient_condition[i]);
		if (i != nr_sectors-1)	printf(", ");
	}
	printf("\nNutrient comp.:\t\t%d\nStrain comp.:\t\t%d\nNR:\t\t\t%d\nNC:\t\t\t%d\n\nMutations:\t\t%s\nMixing:\t\t\t%s\n", nutrient_competition, strain_competition, NR, NC, mutations_on?"Yes":"No", well_mixing?"Yes":"No");
}


void Setup(int argc, char** argv) {

	string ReadOut, ReadOutN, command;
	bool project_name_found = false;
	int i, j, i_condition=0;

	for (i=0; i<nr_sectors; i++)	nutrient_condition[i]=default_conditions[i];

	for(i=1; i<argc; i++)	//Loop through input arguments.
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

		else if(ReadOut=="-g")
		{
			j=1;
			while (i+j != argc)
			{
				ReadOutN = (char*) argv[i+j];
				if (ReadOutN.substr(0,1)=="-" || j == max_input_files+1){
					i += j-1;	//The main for-loop will also add 1 to i.
					j=1;
					break;
				}
				else	//Haven't gotten to the next different command-line argument (or max. number of input files), so keep interpreting as another genome file.
				{
					genome_files[j-1] = argv[i+j];
					j++;
				}
			}
			if (i+j == argc)
			{
				i += j-1;
				break;
			}
		}

		else if(ReadOut=="-e")
		{
			j=1;
			while (i+j != argc)
			{
				ReadOutN = (char*) argv[i+j];
				if (ReadOutN.substr(0,1)=="-" || j == max_input_files+1){
					i += j-1;
					j=1;	//Otherwise the statement below might think that i+j==argc, even though that is not really true.
					break;
				}
				else
				{
					expression_files[j-1] = argv[i+j];
					j++;
				}
			}
			if (i+j == argc)
			{
				i += j-1;
				break;
			}
		}

		else if(ReadOut=="-d")
		{
			j=1;
			while (i+j != argc)
			{
				ReadOutN = (char*) argv[i+j];
				if (ReadOutN.substr(0,1)=="-" || j == max_input_files+1){
					i += j-1;
					j=1;	//Otherwise the statement below might think that i+j==argc, even though that is not really true.
					break;
				}
				else
				{
					definition_files[j-1] = argv[i+j];
					j++;
				}
			}
			if (i+j == argc)
			{
				i += j-1;
				break;
			}
		}

		else if(ReadOut=="-m" && (i+1)!=argc)
		{
			mutation_file = argv[i+1];
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
			nutrient_condition[i_condition] = atof(argv[i+1]);
			i_condition++;
			i++;
			continue;
		}

		else if(ReadOut=="-nC" && (i+1)!=argc)
		{
			nutrient_competition = atoi(argv[i+1]);
			i++;
			continue;
		}

		else if(ReadOut=="-sC" && (i+1)!=argc)
		{
			strain_competition = atoi(argv[i+1]);
			i++;
			continue;
		}

		else if(ReadOut=="-r" && (i+1)!=argc)
		{
			NR = atoi(argv[i+1]);
			i++;
			continue;
		}

		else if(ReadOut=="-c" && (i+1)!=argc)
		{
			NC = atoi(argv[i+1]);
			i++;
			continue;
		}

		else if(ReadOut=="-st" && (i+1)!=argc)
		{
			init_stage = atoi(argv[i+1]);
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

		else if(ReadOut=="--mixed")
		{
			well_mixing = true;
		}

		else if(ReadOut=="--help")
		{
			PrintUsage(true);
		}

		else
		{
			PrintUsage(false);
		}
	}

	if (!project_name_found)	folder += "Project_Name";	//I did not manage to give the date as an extension to the folder.

	PrintLog();

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

void SetMutationRates()
{
	string line;
	int count_lines=0, success;
	char* data;
	char buffer[20];
	ifstream infile(mutation_file.c_str());

	if (!infile.is_open())
	{
		printf("Mutation-file could not be opened.\n");
		exit(1);
	}

	printf("Reading mutation rates from file: %s\n", mutation_file.c_str());
	while(getline(infile,line))
	{
		data = (char*)line.c_str();
		if (count_lines < 8)
		{
			success = sscanf(data, "#%s\t%lf %lf %lf %lf", buffer, &mu[count_lines][HOUSE], &mu[count_lines][BSITE], &mu[count_lines][REGULATOR], &mu[count_lines][EFFECTOR]);
		}

		if(count_lines < 8 && success != 5)
		{
			cerr << "Error: mutation file potentially corrupt.\n" << endl;
			exit(1);
		}

		count_lines++;
	}
}
