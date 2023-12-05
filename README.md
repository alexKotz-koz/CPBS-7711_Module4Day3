# README

## WARNING: This project will take days to run with no test data, modify and run at your own discresion


### Installation Requirements:

### Project Info/ System Requirements:
- Python version 3.11.
- Linux or Unix based Operating System required to run the project using the instructions below.
- This project was built using Python 3.11 on macOS Sonoma v14.0.
- This project utilizes the concurrent module for threading from pythons standard library. 
- Development machine: 2.8 GHz Quad-Core Intel 1.7 Processor with 16 GB memory. 
- Libraries for visualizaiton:
    - Matplotlib (v3.8.1), NetworkX (v3.2.1), netgraph (v4.12.11)

### Ouput files in M4D3 sub-directories:

#### created_at_runtime:

**averageGeneScores.json**: An object containing every gene in STRING 1.txt and Input.gmt.txt with its associated average gene score and locusId

**faNetwork.txt**: A flat network file containing FA-FA gene connections

**finalFaPopulation.json**: An object containing the final optimized fa subnetwork population with its associated average density

**generations.json**: A series of objects containing information and data regarding each generation of the genetic algorithm optimization routine

**random_fa_subnetworks.json**: The initial population of random fa subnetowrks

**random_non_fa_subnetworks.json**: A series of non fa subnetwork populations

**scoring_geneScores.txt**: A raw file containing a list of all gene scores calculated in main.py


#### static:

**Input.gmt.txt**: A list of the 12 known genes associated with Fanconi Anemia with there cooresponding loci.

**STRING 1.txt**: A network of connected genes and their edge measurement.

**faNetwork.txt**: A subnetwork of connected fa genes, derived from Module 1 Day 3 project. 

#### supplemental_documents_for_submission:

**Day3_Output.gmt**: A list of average gene scores, calculated from the final optimized fa population, each gene is listed in its associated locus

# Attention: Format: gene<space>gene...<density score (becuase the pvalues will all be 0, not giving any statistical significance)><space><pval>
**Day3_Output_Network{DD}_pval<p-value>.txt**: A series of files of the top ten scoring subnetworks and thier associated p-values. Naming schema: DD from 01 to 10 and the corresponding calculated p-value for the population.

**generation_statistics.txt**: A file with statistics of each generation of the FA genetic algorithm routine.



### Setup and Configuration:

- Download and extract zip folder containing the two source text files and main.py, into a known directory on host machine.

<hr>

## Run

1. Open an instance of the terminal.

**Note**: If using a conda environment, please activate the conda environment prior to running the main.py file.

2. Navigate to the directory in which the project was extracted to.

    Example: 
        
        cd ~\Desktop\M4D3\

3. Run the main.py file using python interpreter.

    **Note**: This project requires python version 3.11 or higher, to run the main.py file, please assure you are using the correct interpreter. 

    Example:

        python3 main.py
_***The project will take ~ 5 hours to execute without adjustments to main.py(), as specified below.***_

