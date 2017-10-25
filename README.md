

# Novo&Stitch Software




***********
# DESCRIPTION
***********

Novo&Stitch is an assembly reconciliation tool which takes advantage of optical maps to accurately carry out assembly reconciliation. One or more optical maps are used to obtain coordinates for the contigs, which are then stitched based on their alignments. The presence of the optical map dramatically reduces the complexityof the problem and the possibility of a misjoin. Combinatorial optimization models and technologies such as graph model, dynamic programming, weighted vertex cover model on hypergraph, greedy strategy, linear programming are used for solving some subproblems like data reduction, error correction and post-processing. 
Extensive experimental results demonstrate that Novo&Stitch can significantly improve the contiguity of de novo genome assemblies without introducing misassembles or reducing completeness.
Until now, Novo&Stitch can only run on Unix/Linux systems.  


**********
# DEPENDENCY
**********

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

7. fa2cmap_multi.pl
In Novo$Stitch, we use one perl script "fa2cmap_multi.pl" of a scaffolding tool Irys-scaffolding. 
The users need to download this script from https://github.com/i5K-KINBRE-script-share/Irys-scaffolding/blob/e8e8f177dce2bf59421bd00c517ab7dc683e25d4/KSU_bioinfo_lab/assemble_XeonPhi/third-party/fa2cmap_multi.pl
and then put it in our ./Novo_Stitch/tools directory.   




************
# INSTALLATION
************

In Novo$Stitch, only the part written by C/C++ language needs to be compiled. 
Plaase follow the steps as follow:
$cd ./Novo_Stitch
$g++ ./scripts/Script_mapping_contigs2contigs.cc -o ./tools/Script_mapping_contigs2contigs -fopenmp
$g++ ./scripts/Script_mapping_contigs2stitched_contigs.cc -o ./tools/Script_mapping_contigs2stitched_contigs -fopenmp

And as we said in DEPENDENCY, don't forget to download fa2cmap_multi.pl script and put it in ./Novo_Stitch/tools directory. 




**********
# PARAMETERS
**********

There are two kinds of parameters, functional parameters and performance related parameters. 
Functional parameters are the ones that users may have to set by themselves. For performance related parameters, most users may not be able to understand the meaning of them, unless they are very familiar with every step of Novo&Stitch algorithm. Therefore, we offer two sets of suggested parameters for users in "performance related parameters" section. 

1. functional parameters

-f: specify the parameter file which lists all the input paramerters.
Novo$Stitch offers two ways to input parameters for users as follow:
(1) first way
e.g. 
$cd ./Novo_Stitch
$python ./scripts/main.py -x /home/stelo/BIONANO_in_progress/tools/RefAligner -i /home/weihua/cowpea/fastas_cowpea_eight.txt -o /home/weihua/cowpea/eight_stitch_loose_false_and_contained_BssSI -t BssSI -m /home/stelo/BIONANO_in_progress/vu_bsss1_102.cmap -p 32 -a 0 -b 0.2 -c 5000 -d 0.5 -e 0.8

(2) second way
e.g.
$cd ./Novo_Stitch
$python ./scripts/main.py -f ~/phytophthora/parameters.txt
Then in parameters.txt, the parameters are listed line by line as follows:
-x /home/stelo/BIONANO_in_progress/tools/RefAligner
-i /home/weihua/cowpea/fastas_cowpea_eight.txt
-o /home/weihua/cowpea/eight_stitch_loose_false_and_contained_BssSI
-t BssSI
-m /home/stelo/BIONANO_in_progress/vu_bsss1_102.cmap
-p 32
-a 0
-b 0.2
-c 5000
-d 0.5
-e 0.8


-i: specify the fasta list file which lists the address of input fasta files line by line.
e.g. fasta_list.txt 
/home/stelo/phytophthora/seven_2/seven_2.fasta
/home/weihua/canu.fasta
./falcon.fasta
This parameter is required to be specified by users. 

-o: specify the output dirctory which will contains all of the output files and intermediate files generated. This parameter is required to be specified by users. 

-m: specify the optical map file. This parameter is required to be specified by users. 

-t: specify the restriction enzyme (e.g. BspQI). This parameter is required to be specified by users. 

-p: specify the number of threads. The default value is "32".

-x: specify the address of Refaligner executable program on your machine. e.g. /home/weihua/tools/RefAligner
If you can run RefAligner in any directory of your machine by command like "RefAligner [parameters]", then you can put "RefAligner" as this parameter after "-x". Anyway, this parameter should be the first term of the command your input when you run RefAligner on your machine. The default value is "RefAligner". 

-g: specify the address of glpsol executable prgram on your machine. e.g. /home/weihua/tools/glpsol
You should be able to find it after installing GLPK
Same as Refaligner, you may also be able to use just "glpsol". The default value is "glpsol".

-y: specify the address of blastn executable program on your machine. e.g. /home/weihua/tools/blastn
Same as Refaligner, you may also be able to use just "blastn". The default value is "blastn".


2. performance related parameters

-a: specify threshold_1 which is used when stitching overlapped contigs. When stitching, only the alignments whose length is larger than threshold_1 are considered to use. The default value is "3000".

-b: specify threshold_2 which is used when stitching overlapped contigs. threshold_2 is the threshold of ratio between difference of two contigs' coordinates on optical map and smaller contigs's length. When stitching, only when the ratio is smaller than threshold_2, the alignment is considered to use.  
The default value is "0.1".

-c: specify threshold_3 which is used in post-stitch processing. In post-processing, only the alignments whose length is larger than threshold_3 are considered to use. The default value is "10000".

-d: specify threshold_4 which is used in post-stitch processing. 
threshold_4 is a threshold of ratio between length of two alignments' overlap and the length of smaller alignments. In post-processing, if the ratio is larger than threshold_4, we don't consider that alignment. The default value is "0.5".

-e: specify threshold_5 which is used in post-stitch processing. 
threshold_5 is a threshold for the ratio between the total length of parts of contig mapped to stitched contig and the total length of this contig. In post-stitching, if the ratio is smaller than threshold_5, we will cancel this stitch. The default value is "0.9".

-h: specify threshold_6 which is used in pre-stitch processing. Only the alignments whose confidence from Refaligner is larger than threshold_6 are used. The default value is "25".

-r: specify threshold_7 which is used in false alignments removal. When building hypergraph, if the diffence between the two distance between two contigs' middle point coordinates on two optical map molecules is larger than threshold_7 of larger contig's length, one hyperedge is built on the four alignments. The default value is "0.2".

As we said before, here we offer two sets of performance related parameters for users who cannot understand the meaning of them as follow:
strict (default): -a 3000 -b 0.1 -c 10000 -d 0.5 -e 0.9 -h 25 -r 0.2
loose: -a 0 -b 0.2 -c 5000 -d 0.5 -e 0.8 -h 25 -r 0.2
The strict set is more conservative than loose set, which means it stitches less but has less chance to make mistakes than loose set. 




************
# OUTPUT FILES
************

There are two kinds of output files, stitched contigs file and log files. Stitched contigs file contains the final stitched contigs, while the log files contain the information which shows the whole process of running. All of the output files are stored in output directory specified by users. There are some other intermediate files inside, which could be ignored by users.
 
1. stitched contigs file

final_stitched.fasta: It's a fasta file which contains the stitched contigs.

2. log files and some intermediate files

(1) in main output directory 

parameters.log: shows the values of all parameters used.

input.log: shows all the input fasta files

(2) in each iteration_* directory

alms_0_initial.log: all alignments 

alms_1_removed_lowconf.log: alignments after removing low quilty alignments

alms_2_removed_false_alms.log: alignments after removing false alignments

alms_3_removed_contained_locally.log: alignments after removing contained alignments locally

alms_4_unified.log: alignments after coordinates unification

alms_5_removed_contained_globally.log: alignments after removing contained alignments globally

alms_6_mtp.log: alignments of mtp contigs

ugraph_1.log: first association graph

ugraph_2.log: second association graph

forest_1.log: msf of first association graph

forest_2.log: msf of second association graph

shifts.log: shifts information for coordinates unification

alm_graph_edges.txt: edges of alignment graph

alm_graph_vertices.txt: vertices of alignment graph

cover.log: optimal vertex cover got for false alignments removal

mtp_nodes.log: mtp contigs

stitch_info.log: stitch information

confidences.txt: the confidences of the contigs mapped to stitched contigs

cancel_id_list.txt: the information of cancelled stitch

input_contigs.fasta: input contigs for this iteration

input_contigs_mtp.fasta: mtp contigs for this iteration

adjusted_contigs_total.fasta: ouput stitched contigs for this iteration



