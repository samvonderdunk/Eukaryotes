### Readme for TFEvol branch of Eukaryotes ###

This branch is a simpler version / stripped down from the Eukaryotes code, which means that it includes the speed improvements implemented in Eukaryotes with respect to Prokaryotes, while no longer really encoding eukaryotes.

# RNA & Beyond: Evolution of TF specificity in prokaryotes and eukaryotes.
The idea is to change the binding mechanism: regulatory genes have binding motifs that specify binding probabilities to binding sites (i.e. sequences). Potentially the target expression and fitness functions of cells will also be altered in this project.

# Installation.
Fork / clone the whole git repository and make by running "make" on the command line in the main directory. You will compile the binary called Eukaryotes and can run this using several input files (see below).

# Preparation.
1) You will need input genomes, expressions, definitions (may become obselete), and mutation rates all as separate files. An example for each is included in the repository in the folder example_input. This example is an evolved cell that executes a prokaryotic cell-cycle very efficiently. It should do well under a wide range of nutrient conditions and changes in genome length (as long as the regulatory network is not meddled with).
2) You will also need to change the base folder where output is directed to. Do this by changing the "string folder = " in World.cc at line 16 to a local folder on your pc. Always make the whole project again after making such changes.

# Execution.
Run from the command-line as follows (runs 100 timesteps on a 50x50 grid, showing the population size at each timestep):
# ./Eukaryotes -s 11 -p test -g example_input/genome_8.g -e example_input/expression_8.g -d example_input/definition_8.g -m example_input/mutations.g -r 50 -c 50 -tT 1 -tN 100
To see what other run options there are, just give a "faulty" argument, e.g.:
# ./Eukaryotes -X
