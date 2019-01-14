#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <typeinfo>
#include <random>
#include <algorithm>

#include "monteCarlo.h"

using namespace std;

string build_fasta_file(vector<tuple<string, int> > final_path,
                        map<string, map<string, vector<vector<float> > > > groupedCR,
                        vector<string> keysCR,
                        map<string, map<string, vector<vector<float> > > > groupedRR,
                        vector<string> keysRR,
                        map<string, string> fastaReads,
                        map<string, string> fastaContigs) {
    bool firstContig = true;
    string fastaString = ">final_path\n";
    int lengthOfFinalPath;
    lengthOfFinalPath = final_path.size();
    tuple<string, int> first;
    tuple<string, int> second;
    for (int i = 0; i < lengthOfFinalPath - 1; i++) {
        first = final_path[i];
        string firstInPath = get<0>(first);
        second = final_path[i + 1];
        string secondInPath = get<0>(second);
        int secondInPathIndex = get<1>(second);
        int firstInPathIndex = get<1>(first);
        if ((firstInPath[0] == 'C' or firstInPath[0] == 'c') and (firstInPath[1] == 't') and (firstInPath[2] == 'g')) {
            vector<float> readCont = groupedCR[firstInPath][secondInPath][secondInPathIndex];
            string fastaOfRead = fastaReads[secondInPath];
            int OH2 = readCont[3];
            int EL1 = readCont[4];
            if (firstContig) {
                string fastaOfContig = fastaContigs[firstInPath];
                fastaString.append(fastaOfContig);
                firstContig = false;
            }
            if (OH2 != 0) {
                fastaString = fastaString.substr(0, fastaString.size() - OH2);
            }
            fastaString.append(fastaOfRead.substr(fastaOfRead.size() - EL1, fastaOfRead.size()));

        } else if ((secondInPath[0] == 'C' or secondInPath[0] == 'c') and (secondInPath[1] == 't') and (secondInPath[2] == 'g')) {
            vector<float> readCont = groupedCR[secondInPath][firstInPath][firstInPathIndex];
            string fastaOfRead = fastaReads[firstInPath];
            int OH1 = readCont[2];
            int EL2 = readCont[5];
            if (OH1 != 0) {
                fastaString = fastaString.substr(0, fastaString.size() - OH1);
            }
            fastaString.append(fastaOfRead.substr(fastaOfRead.size() - EL2, fastaOfRead.size()));
        } else {
            vector<float> readCont = groupedRR[firstInPath][secondInPath][secondInPathIndex];
            string fastaOfRead = fastaReads[secondInPath];
            int OH2 = readCont[3];
            int EL1 = readCont[4];
            if (OH2 != 0) {
                fastaString = fastaString.substr(0, fastaString.size() - OH2);
            }
            fastaString.append(fastaOfRead.substr(fastaOfRead.size() - EL1, fastaOfRead.size()));
        }
    }
    return fastaString;

}

void saveFasta(string fasta, string path) {
      ofstream myfile;
      myfile.open (path);
      myfile << fasta;
      myfile.close();
}
