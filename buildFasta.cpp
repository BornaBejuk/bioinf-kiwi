#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <typeinfo>
#include <random>
#include <algorithm>

#include "monteCarlo.h"

using namespace std;

// builds fasta file, currently not in use
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

// builds fasta string by given scaffold order
string buildFastaString(vector<vector<tuple<string, int> > > finalOrder, map<string, map<string, vector<vector<float> > > > groupedCR,
                        map<string, map<string, vector<vector<float> > > > groupedRR, map<string, string> fastaReads, map<string, string> fastaContigs, vector<string> keysCR) {

    string fastaString;
    string currentTarget;
    string currentQuery;
    int currentTargetIndex;
    int currentQueryIndex;

    int flagFirstEver = 1;
    int flagMiddle = 1;
    string middle;
    int flagMiddleAfter = 0;


    int leftIndex;
    int rightIndex;

    string substringToInsert;

    vector<string> contigsAppended;

    // string begin;
    // string end;
    // float ctgBeginLen;
    // float ctgEndLen;
    // float pathLen = 0.0;
    // string currentTargetRead;
    // string currentQueryRead;
    for( auto pair : finalOrder) {
        currentTarget = get<0>(pair[0]);
        currentTargetIndex = get<1>(pair[1]);
        currentQuery = get<0>(pair[1]);
        currentQueryIndex = get<1>(pair[1]);
        // cout << "TARGET: " << currentTarget << " Query: " << currentQuery << endl;

        if ((currentTarget[0] == 'C' or currentTarget[0] == 'c') and (currentTarget[1] == 't') and (currentTarget[2] == 'g')) {
            if( flagFirstEver == 1 ){
                contigsAppended.push_back(currentTarget);
                flagFirstEver = 0;

                leftIndex = groupedCR[currentTarget][currentQuery][currentQueryIndex][3]; // OH2
                rightIndex = fastaContigs[currentTarget].size(); // end

                // cout << "LEn  ctg " << fastaContigs[currentTarget].size() << endl;
                substringToInsert = fastaContigs[currentTarget].substr(leftIndex, rightIndex);
                fastaString = substringToInsert + fastaString;

                leftIndex = 0;
                rightIndex = groupedCR[currentTarget][currentQuery][currentQueryIndex][4]; // EL1
                // cout << "LEn  read " << fastaReads[currentQuery].size() << endl;

                substringToInsert = fastaReads[currentQuery].substr(leftIndex, rightIndex);
                fastaString = substringToInsert + fastaString;

            } else if ( flagMiddle == 1) {
                contigsAppended.push_back(currentTarget);
                flagMiddle = 0;
                flagMiddleAfter = 1;
                middle = currentTarget;

                leftIndex = 0;
                rightIndex = groupedCR[currentTarget][currentQuery][currentQueryIndex][5]; // EL2

                substringToInsert = fastaContigs[currentTarget].substr(leftIndex, rightIndex);
                fastaString = substringToInsert + fastaString;

                // this is removing OH1
                // leftIndex = 0;
                rightIndex = groupedCR[currentTarget][currentQuery][currentQueryIndex][4]; // OH1

                // substringToInsert = fastaReads[currentQueryRead].substr(leftIndex, rightIndex);
                fastaString = fastaString.substr(rightIndex);

            } else if( flagMiddleAfter == 1) {
                flagMiddleAfter = 0;
                flagMiddle = 1;
                if( currentTarget != middle ) {
                    // in this case we dont have the same contig so we cant continue building but start again with new contig
                    contigsAppended.push_back(currentTarget);

                    leftIndex = groupedCR[currentTarget][currentQuery][currentQueryIndex][3]; // OH2
                    rightIndex = fastaContigs[currentTarget].size(); // end

                    substringToInsert = fastaContigs[currentTarget].substr(leftIndex, rightIndex);
                    fastaString = substringToInsert + fastaString;

                    leftIndex = 0;
                    rightIndex = groupedCR[currentTarget][currentQuery][currentQueryIndex][4]; // EL1

                    substringToInsert = fastaReads[currentQuery].substr(leftIndex, rightIndex);
                    fastaString = substringToInsert + fastaString;
                } else {
                    leftIndex = 0;
                    rightIndex = groupedCR[currentTarget][currentQuery][currentQueryIndex][3]; // OH2

                    // remove OH2
                    fastaString = fastaString.substr(rightIndex);

                    leftIndex = 0;
                    rightIndex = groupedCR[currentTarget][currentQuery][currentQueryIndex][4]; // EL1

                    substringToInsert = fastaReads[currentQuery].substr(leftIndex, rightIndex);
                    fastaString = substringToInsert + fastaString;
                }
            }
        } else {
            // read read
            leftIndex = 0;
            rightIndex = groupedRR[currentTarget][currentQuery][currentQueryIndex][3]; // OH2

            // substringToInsert = fastaReads[currentTarget].substr(leftIndex, rightIndex);
            // remove left OH2
            // TODO check if this is goood, else fastaString.erase(rightIndex);
            fastaString = fastaString.substr(rightIndex);

            leftIndex = 0;
            rightIndex = groupedRR[currentTarget][currentQuery][currentQueryIndex][4]; // EL1

            substringToInsert = fastaReads[currentQuery].substr(leftIndex, rightIndex);
            fastaString = substringToInsert + fastaString;

        }
    }

    // cout << "FIRst part over" << endl;

    if( contigsAppended.size() < keysCR.size()) {
        for( auto ctg : keysCR) {
            if( !(std::find(contigsAppended.begin(), contigsAppended.end(), ctg) != contigsAppended.end())){
                fastaString = fastaString + fastaContigs[ctg];
                cout << "This contig is not already appended: " << ctg << endl;
            }
        }
    }

    fastaString = ">final_path\n" + fastaString;
    return fastaString;
}

// saves fasta string in file
void saveFasta(string fasta, string path) {
      ofstream myfile;
      myfile.open (path);
      myfile << fasta;
      myfile.close();
}
