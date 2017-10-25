# Novo_Stitch
Novo&Stitch is an assembly reconciliation tool which takes advantage of optical maps to accurately carry out assembly reconciliation. One or more optical maps are used to obtain coordinates for the contigs, which are then stitched based on their alignments. The presence of the optical map dramatically reduces the complexityof the problem and the possibility of a misjoin. Combinatorial optimization models and technologies such as graph model, dynamic programming, weighted vertex cover model on hypergraph, greedy strategy, linear programming are used for solving some subproblems like data reduction, error correction and post-processing. 
Extensive experimental results demonstrate that Novo&Stitch can significantly improve the contiguity of de novo genome assemblies without introducing misassembles or reducing completeness.
Until now, Novo&Stitch can only run on Unix/Linux systems.  

# DEPENDENCY
1. python
The majoy part of Novo&Stitch is written in Python, so Python has to be installed. 
Python2.7(or above) is suggested.  

2. perl
In Novo$Stitch, we use one perl script "fa2cmap_multi.pl" of a scaffolding tool Irys-scaffolding.
So perl has to be installed.

3. g++
To speed up the running of blastn, OpenMP is used for parallelization. That part of program is written in C/C++ language and executable programs are called by python codes. So g++ has to be installed for compiling those C/C++ codes
Make sure your g++ supports compiling OpenMP programs. 

4. GLPK
GLPK is a tool for solving Linear Programming (LP) model. Since the false alignments removal module of Novo&Stitch uses a standalone tool of GLPK called "glpsol" to solve LP model, GLPK has to be installed. 
GLPK can be found from https://www.gnu.org/software/glpk/#TOCdownloading.

5. RefAligner
RefAligner is a tool developed by BioNano company to align contigs to optical maps. It is called by Novo&Stitch, so it has to be installed. 
RefAligner can be found from http://www.bnxinstall.com/RefalignerAssembler/Linux/ 

6. blastn
blastn is called by Novo&Stitch to get the sequence alignments between overlapped contigs, so it has to been installed. 

7. fa2cmap_multi.pl file
In Novo$Stitch, we use one perl script "fa2cmap_multi.pl" of a scaffolding tool Irys-scaffolding. 
The users need to download this script from https://github.com/i5K-KINBRE-script-share/Irys-scaffolding/blob/e8e8f177dce2bf59421bd00c517ab7dc683e25d4/KSU_bioinfo_lab/assemble_XeonPhi/third-party/fa2cmap_multi.pl
and then put it in our ./novo_stitch/tools directory.

# INSTALLATION
In Novo$Stitch, only the part written by C/C++ language needs to be compiled. 
Plaase follow the steps as follow:
$cd ./novo_stitch
$g++ ./scripts/Script_mapping_contigs2contigs.cc -o ./tools/Script_mapping_contigs2contigs -fopenmp
$g++ ./scripts/Script_mapping_contigs2stitched_contigs.cc -o ./tools/Script_mapping_contigs2stitched_contigs -fopenmp

And as we said in DEPENDENCY, don't forget to download fa2cmap_multi.pl script and put it in ./novo_stitch/tools directory. 

