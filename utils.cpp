#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>
#include <chrono>


#include "utils.hpp"

using namespace std;
typedef std::chrono::high_resolution_clock Clock;

// author: Karlo Brajidc

// loads overlap data files
void loadData(string path, vector<string> &queryNames, vector<int> &queryLens, vector<float> &queryStarts,
            vector<float> &queryEnds, vector<string> &targetNames, vector<int> &targetLens,
            vector<float> &targetStarts, vector<float> &targetEnds, vector<float> &resMatches,
            vector<float> &blockLens, vector<float> &SI, float SImin, vector<float> &extensionSides, vector<float> &strands) {

    string line;

    ifstream inputFile(path);
    float sisum = 0.0;
    float counter = 0.0;
    if (inputFile.good()) {
        while (getline(inputFile, line)) {
            stringstream linestream(line);
            string qName;
            int qLen;
            float qStart;
            float qEnd;
            string strand;
            string tName;
            int tLen;
            float tStart;
            float tEnd;
            float resMatch;
            float bLen;
            int mapQuality;

            getline(linestream, qName, '\t');

            linestream >> qLen >> qStart >> qEnd >> strand >> tName >> tLen >> tStart >> tEnd >> resMatch >> bLen;

            if( qName == tName) {
                continue;
            }

            float si = resMatch / bLen;
            sisum += si;
            counter += 1;
            if( si > SImin) {
                if( extendRight(qEnd, qLen, tEnd, tLen)) {
                    extensionSides.push_back(1);
                } else if( extendLeft(qStart, tStart)) {
                    extensionSides.push_back(0);
                } else {
                    continue;
                }
                queryNames.push_back(qName);
                queryLens.push_back(qLen);
                queryStarts.push_back(qStart);
                queryEnds.push_back(qEnd);
                if( strand == "+"){
                    strands.push_back(0.0);
                } else {
                    strands.push_back(1.0);
                }
                targetLens.push_back(tLen);
                targetNames.push_back(tName);
                targetStarts.push_back(tStart);
                targetEnds.push_back(tEnd);
                resMatches.push_back(resMatch);
                blockLens.push_back(bLen);
                SI.push_back(si);
            }
        }
    }
    cout << "AVG SI: " << sisum / counter << endl;
}

// determines if target is extended to right
bool extendRight(float queryEnd, int queryLen, float targetEnd, int targetLen) {
    return (targetLen - targetEnd - 1) < (queryLen - queryEnd - 1);
}

// determines if target is extended to left
bool extendLeft(float queryStart, float targetStart) {
    return targetStart < queryStart;
}

// filters reads that are fully contained inside contig
void filterContained(vector<string> &queryNames, vector<int> &queryLens, vector<float> &queryStarts,
            vector<float> &queryEnds, vector<string> &targetNames, vector<int> &targetLens,
            vector<float> &targetStarts, vector<float> &targetEnds, vector<float> &resMatches,
            vector<float> &blockLens, vector<float> &extensionSides) {

    for( int i = 0; i < queryNames.size(); i++) {
        if( extendRight(queryEnds[i], queryLens[i], targetEnds[i], targetLens[i])) {
            extensionSides.push_back(1);
        } else if( extendLeft(queryStarts[i], targetStarts[i])) {
            extensionSides.push_back(0);
        } else {
            queryNames.erase(queryNames.begin() + i);
            queryLens.erase(queryLens.begin() + i);
            queryStarts.erase(queryStarts.begin() + i);
            queryEnds.erase(queryEnds.begin() + i);
            targetLens.erase(targetLens.begin() + i);
            targetNames.erase(targetNames.begin() + i);
            targetStarts.erase(targetStarts.begin() + i);
            targetEnds.erase(targetEnds.begin() + i);
            resMatches.erase(resMatches.begin() + i);
            blockLens.erase(blockLens.begin() + i);

            --i;
        }
    }
}

// calculates sequence index
void calculateSI(vector<float> &SI, vector<float> &resMatches, vector<float> &blockLens) {

    for( int i = 0; i < resMatches.size(); i++) {
        SI.push_back(resMatches[i] / blockLens[i]);
    }
}

// filters by SImin, currently not in use
void filterBySI(float SImin, vector<string> &queryNames, vector<int> &queryLens, vector<float> &queryStarts,
            vector<float> &queryEnds, vector<string> &targetNames, vector<int> &targetLens,
            vector<float> &targetStarts, vector<float> &targetEnds, vector<float> &resMatches,
            vector<float> &blockLens, vector<float> &extensionSides, vector<float> SI) {

    for( int i = 0; i < queryStarts.size(); i++) {
        if( SI[i] < SImin ) {
            queryNames.erase(queryNames.begin() + i);
            queryLens.erase(queryLens.begin() + i);
            queryStarts.erase(queryStarts.begin() + i);
            queryEnds.erase(queryEnds.begin() + i);
            targetLens.erase(targetLens.begin() + i);
            targetNames.erase(targetNames.begin() + i);
            targetStarts.erase(targetStarts.begin() + i);
            targetEnds.erase(targetEnds.begin() + i);
            // resMatches.erase(resMatches.begin() + i);
            // blockLens.erase(blockLens.begin() + i);
            extensionSides.erase(extensionSides.begin() + i);
            SI.erase(SI.begin() + i);
            --i;
        }
    }
}

// calculates overlap length
void calculateOL(vector<float> &OL1, vector<float> &OL2, vector<float> &queryStarts, vector<float> &queryEnds, vector<float> &targetStarts, vector<float> &targetEnds) {

    for( int i = 0; i < queryStarts.size(); i++) {
        // TODO add - 1 to each because both are included, end and start
        OL1.push_back(queryEnds[i] - queryStarts[i]-1);
        OL2.push_back(targetEnds[i] - targetStarts[i]-1);
    }
}

// calculates overhang length
void calculateOH(vector<float> &OH1, vector<float> &OH2, vector<int> &queryLens, vector<float> &queryStarts, vector<float> &queryEnds,
    vector<int> &targetLens, vector<float> &targetStarts, vector<float> &targetEnds, vector<float> &extensionSides) {

    for( int i = 0; i < targetLens.size(); i++) {
        if( extensionSides[i] == 1) {
            OH1.push_back(queryStarts[i]);
            OH2.push_back(targetLens[i] - targetEnds[i] - 1);
        } else {
            OH1.push_back(queryLens[i] - queryEnds[i] - 1);
            OH2.push_back(targetStarts[i]);
        }
    }
}

// calculates extension length
void calculateEL(vector<float> &EL1, vector<float> &EL2, vector<int> &queryLens, vector<float> &queryStarts, vector<float> &queryEnds,
    vector<int> &targetLens, vector<float> &targetStarts, vector<float> &targetEnds, vector<float> &extensionSides) {

    for( int i = 0; i < targetLens.size(); i++) {
        if( extensionSides[i] == 1) {
            EL1.push_back(queryLens[i] - queryEnds[i] - 1);
            EL2.push_back(targetStarts[i]);
        } else {
            EL1.push_back(queryStarts[i]);
            EL2.push_back(targetLens[i] - targetEnds[i] - 1);
        }
    }
}

// calculates overlap score
void calculateOS(vector<float> &OS, vector<float> &OL1, vector<float> &OL2, vector<float> &SI) {

    for( int i = 0; i < OL1.size(); i++) {
        OS.push_back(0.5 * (OL1[i] + OL2[i]) * SI[i]);
    }
}

// calculates extension score
void calculateES(vector<float> &ES1, vector<float> &ES2, vector<float> &OS, vector<float> &EL1, vector<float> &EL2, vector<float> &OH1, vector<float> &OH2) {

    for( int i = 0; i < EL1.size(); i++) {
        ES1.push_back(OS[i] + 0.5 * EL2[i] - 0.5 * (OH1[i] + OH2[i]));
        ES2.push_back(OS[i] + 0.5 * EL1[i] - 0.5 * (OH1[i] + OH2[i]));
    }
}

// calculates all scores
void calculateScores(vector<float> &OS, vector<float> &ES, vector<float> &SI, vector<float> &EL,
                    vector<float> &OH, vector<float> &resMatches, vector<float> &blockLens) {

}

// loads fasta files
map<string, string> loadFasta(string path) {

    map<string, string> fasta;
    string line;
    ifstream inputFile(path);
    string name;
    if (inputFile.good()) {
        while (getline(inputFile, line)) {
            if( line.at(0) == '>' or line.at(0) == '@') {
                stringstream linestream(line);
                getline(linestream, name);
                if( name[1] == 'C' or name[1] == 'c') {
                    name = name.substr(0,5);
                }
                name = name.substr(1,name.size());
            } else if( line.at(0) == 'A' or line.at(0) == 'C' or line.at(0) == 'T' or line.at(0) == 'G'){
                fasta[name] = line;
            }
        }
    }
    return fasta;
}

void timeIt(std::chrono::_V2::system_clock::time_point& start, std::string action){
  std::chrono::duration<double> elapsed = Clock::now() - start;
  cout << "\033[92m" << action << " done in: " << elapsed.count() << " seconds.\033[0m\n\n";
  start = Clock::now();
}
