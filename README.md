
# Welcome
![SNIPR Logo](SNIPR.png)

SNIPR is a software to design single nucleotide specific programmable riboregulator. 

See published paper: (Precise and Programmable Detection of Mutations Using Ultraspecific Riboregulators)[https://www.cell.com/cell/fulltext/S0092-8674(20)30155-0]

# Requirements

The SNIPR software makes use of the multi-objective sequence design features and energy calculation capabilities of NUPACK. Users must download the (NUPACK source code)[http://www.nupack.org/downloads] and compile it on their computers. Please follow the NUPACK manual to install NUPACK. 

After installation, it is necessary to add NUPACK to the system’s environmental variables to ensure that it can be located by the SNIPR design software. 

From a terminal window, open or create a .bash_profile file in your home directory with the command:

`open ~/.bash_profile`

This file will be run whenever a bash shell is opened and can be used to set environment variables. In the resulting text editor environment, add the NUPACK directory as the following example.

`export NUPACKHOME=/home/Documents/nupack3.2.0`


`export PATH=$PATH":${NUPACKHOME}/bin"`

[For nano, save your changes by typing ctrl + o and hit return to save. Exit nano by typing ctrl + x.]

After setting the environment variables, Matlab will need to be opened from a terminal window to use NUPACK based functions. A sample command to open Matlab is below and will vary depending on the Matlab install directory:

`/Applications/MATLAB_R2017b.app/bin/matlab`


# SNIPR sequence input

In the input folder of the SNIPR install folder, fill the target_input.csv file with the design name, the wild-type target sequence, the mutated target sequence, and the first 29 bases of the output gene. For the target sequence input, the mutation position should be at 15 bases away from the 5’ end of the transcript. Sequences should be specified in the 5’ to 3’ direction.

Full details describing parameter input requirements can be found in the comments in the source code. 

# Running the SNIPR software
Before running the main code SNIPR.m，there are several parameters that can be specified by the user in the main function. The initial version of SNIPR.m is set to default parameters that we have found work well for device generation.

The size of library that the code will generate for each designs can be specified. For each target, there will be 9 conformation designs with different combinations of forward and reverse toehold length. The total number of designs generated by the code will be product of the library size and 9. We recommend setting the library size to at least 10 to generate a sufficient number of designs for later screening. 

`library_num = 10;`

If your computer has multiple cores, the code run time can be reduced by using the parallel implementation of the algorithm by setting:

`IS_PARALLEL = 1;`

Set IS_PARALLEL to 0 to use the single-core implementation.

The parameter select_num is used to specify how many designs you want to try experimentally. SNIPR design is challenging since it targets tiny differences in sequence. We recommend trying at least 5 designs to increase the likelihood of getting a successful sensor. 

`select_num = 9;`

Lastly, it is necessary to specify which of the two potential targets – the wild-type (WT) or the mutant (SNP) – to detect. To design the SNIPR to detect the mutant target:

`SNP_TARGET = 1;`
 
For detecting the wild-type, set SNP_TARGET to 0.

After finishing setting the parameters, you can run the main function SNIPR.m to generate the requested designs. The design process may take some time to finish, but the Matlab command window will keep you apprised of the progress of the algorithm.
