#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <typeinfo>
#include <random>
#include <algorithm>

#include "monteCarlo.hpp"

using namespace std;

// author: Karlo Brajdic

string switchStrand(string s){
    string out = "";
    for( auto c : s){
        if( c == 'C') {
            out += ("G");
        }
        if( c == 'G') {
            out += ("C");
        }
        if( c == 'A') {
            out += ("T");
        }
        if( c == 'T') {
            out += ("A");
        }
    }
    reverse(out.begin(), out.end());
    return out;
}

// author: Karlo Brajdic
// builds fasta string by given scaffold order
string buildFastaString(vector<vector<tuple<string, int> > > finalOrder, map<string, map<string, vector<vector<float> > > > &groupedCR,
                        map<string, map<string, vector<vector<float> > > > &groupedRR, map<string, string> &fastaReads, map<string, string> &fastaContigs, vector<string> keysCR) {

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

    int ctgCounter = 0;
    for( auto pair : finalOrder) {
        currentTarget = get<0>(pair[0]);
        currentTargetIndex = get<1>(pair[1]);
        currentQuery = get<0>(pair[1]);
        currentQueryIndex = get<1>(pair[1]);

        if ((currentTarget[0] == 'C' or currentTarget[0] == 'c') and (currentTarget[1] == 't') and (currentTarget[2] == 'g')) {
            if( flagFirstEver == 1 ){
                contigsAppended.push_back(currentTarget);
                flagFirstEver = 0;
                leftIndex = groupedCR[currentTarget][currentQuery][currentQueryIndex][3]; // OH2
                rightIndex = fastaContigs[currentTarget].size(); // end

                substringToInsert = fastaContigs[currentTarget].substr(leftIndex, rightIndex);
                if( groupedCR[currentTarget][currentQuery][currentQueryIndex][9] == 1.0 ){
                    substringToInsert = switchStrand(substringToInsert);
                }
                fastaString = substringToInsert + fastaString;

                leftIndex = 0;
                rightIndex = groupedCR[currentTarget][currentQuery][currentQueryIndex][4]; // EL1

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
                if( groupedCR[currentTarget][currentQuery][currentQueryIndex][9] == 1.0 ){
                    substringToInsert = switchStrand(substringToInsert);
                }
                fastaString = substringToInsert + fastaString;

                // this is removing OH1
                // leftIndex = 0;
                rightIndex = groupedCR[currentTarget][currentQuery][currentQueryIndex][4]; // OH1

                fastaString = fastaString.substr(rightIndex);

            } else if( flagMiddleAfter == 1) {
                flagMiddleAfter = 0;
                flagMiddle = 1;
                if( currentTarget != middle ) {
                    fastaString = "\n>ctg" + to_string(ctgCounter) + "\n" + fastaString;
                    ctgCounter += 1;
                    // in this case we dont have the same contig so we cant continue building but start again with new contig
                    contigsAppended.push_back(currentTarget);

                    leftIndex = groupedCR[currentTarget][currentQuery][currentQueryIndex][3]; // OH2
                    rightIndex = fastaContigs[currentTarget].size(); // end

                    substringToInsert = fastaContigs[currentTarget].substr(leftIndex, rightIndex);
                    if( groupedCR[currentTarget][currentQuery][currentQueryIndex][9] == 1.0 ){
                        substringToInsert = switchStrand(substringToInsert);
                    }
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
                    if( groupedCR[currentTarget][currentQuery][currentQueryIndex][9] == 1.0 ){
                        substringToInsert = switchStrand(substringToInsert);
                    }
                    fastaString = substringToInsert + fastaString;
                }
            }
        } else {
            // read read
            leftIndex = 0;
            rightIndex = groupedRR[currentTarget][currentQuery][currentQueryIndex][3]; // OH2

            // remove left OH2
            // TODO check if this is goood, else fastaString.erase(rightIndex);
            fastaString = fastaString.substr(rightIndex);

            leftIndex = 0;
            rightIndex = groupedRR[currentTarget][currentQuery][currentQueryIndex][4]; // EL1

            substringToInsert = fastaReads[currentQuery].substr(leftIndex, rightIndex);
            if( groupedRR[currentTarget][currentQuery][currentQueryIndex][9] == 1.0 ){
                substringToInsert = switchStrand(substringToInsert);
            }
            fastaString = substringToInsert + fastaString;

        }
    }

    if( contigsAppended.size() < keysCR.size()) {
        for( auto ctg : keysCR) {
            if( !(std::find(contigsAppended.begin(), contigsAppended.end(), ctg) != contigsAppended.end())){
                fastaString = "\n>ctg" + to_string(ctgCounter) + "\n" + fastaString;
                ctgCounter += 1;
                fastaString = fastaString + fastaContigs[ctg];
            }
        }
    }
    fastaString = ">ctg" + to_string(ctgCounter) + "\n" + fastaString;
    return fastaString;
}

// author: Karlo Brajdic
// saves fasta string in file
void saveFasta(string fasta, string path) {
      ofstream myfile;
      myfile.open (path);
      myfile << fasta;
      myfile.close();
}
