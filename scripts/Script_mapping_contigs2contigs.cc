#include <stdio.h>
#include <vector>
#include <string>
#include <omp.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
using namespace std;

int main(int argc, char* argv[])
{
	string BLASTN(argv[1]);
	string GROUP_FILE(argv[2]);
	string CONTIGS_REPO(argv[3]);

	vector<string> ids;	
	ifstream fi;
	fi.open(GROUP_FILE.c_str());
	string line;
	while(getline(fi, line))
	{
		if (line=="") continue;
		ids.push_back(line);
		 
	}
	fi.close();

	#pragma omp parallel
	{
		int num_threads = omp_get_num_threads();
		int id = omp_get_thread_num();

		for (int j=0; j<ids.size(); j++)
		{
			if (j%num_threads != id) continue;
			printf("thread%d,j=%d\n",id, j);
			string line = ids[j];

			string command1 = "formatdb -i " + CONTIGS_REPO + "/molecule_" + line + ".fasta -p F -o F";
			const char * c_1 = command1.c_str();
			string command2 = BLASTN + " -query " + CONTIGS_REPO + "/molecule_" + line + ".fasta -db " + CONTIGS_REPO + "/molecule_" + line + ".fasta -out " + CONTIGS_REPO + "/self_molecule_" + line + ".blast -outfmt 6";
			const char * c_2 = command2.c_str();
			system(c_1);
			system(c_2);
		
		}
	}
}



