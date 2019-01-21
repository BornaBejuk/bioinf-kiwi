#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>
#include <string>
#include<stdio.h>
#include<stdlib.h>

#include "utils.hpp"
#include "monteCarlo.hpp"
#include "selectPaths.hpp"
#include "scaffolding.hpp"
#include "buildFasta.hpp"
#include "dfsSearch.hpp"

typedef std::chrono::high_resolution_clock Clock;

using namespace std;

// authors: Karlo Brajdic

int main(int argc, char **argv) {
	    string pathCR;
    string pathRR;
    string pathFastaCtgs;
    string pathFastaReads;
    string pathFastaOut;
    float SImin;
    int maxDepth;
    int nTimes;
    if (argc == 2) {
        char* filename = argv[1];
        cout << filename << endl;
        std::ifstream file(filename, std::ifstream::in);

        if (!file) {
            cout << "parameter filename error" << endl;
        }
        string value;
        while (1) {
            file >> value;
            if (value == "pathCR"){
                file >> pathCR;
            } else if (value == "pathRR"){
                file >> pathRR;
            } else if (value == "pathFastaCtgs"){
                file >> pathFastaCtgs;
            } else if (value == "pathFastaReads"){
                file >> pathFastaReads;
            } else if (value == "pathFastaOut"){
                file >> pathFastaOut;
            } else if (value == "SImin"){
                file >> SImin;
            } else if (value == "maxDepth"){
                file >> maxDepth;
            } else if (value == "nTimes"){
                file >> nTimes;
            } else if (value == "end") {
                break;
            }

        }
    } else if (argc == 17){
        string value;
        for (int i=1; i <= 15; i+=2) {
            value = argv[i];
            if (value == "pathCR"){
                pathCR = argv[i+1];
            } else if (value == "pathRR"){
                pathRR = argv[i+1];
            } else if (value == "pathFastaCtgs"){
                pathFastaCtgs = argv[i+1];
            } else if (value == "pathFastaReads"){
                pathFastaReads = argv[i+1];
            } else if (value == "pathFastaOut"){
                pathFastaOut = argv[i+1];
            } else if (value == "SImin"){
                SImin = atof(argv[i+1]);
            } else if (value == "maxDepth"){
                maxDepth = atoi(argv[i+1]);
            } else if (value == "nTimes"){
                nTimes = atoi(argv[i+1]);
            }
        }
    } else {
    pathCR =  "data/EColi-synthetic/overlaps-c-r.paf";
    pathRR = "data/EColi-synthetic/overlaps-r-r.paf";
    pathFastaCtgs = "data/EColi-synthetic/ecoli_test_contigs.fasta";
    pathFastaReads = "data/EColi-synthetic/ecoli_test_reads.fasta";
    pathFastaOut = "data/EColi-synthetic/final.fasta";
    SImin = 0.9;
    maxDepth = 40;
    nTimes = 50;
	}

    //string pathCR = "data/EColi-synthetic/overlaps-c-r.paf";
    //string pathRR = "data/EColi-synthetic/overlaps-r-r.paf";
    //string pathFastaCtgs = "data/EColi-synthetic/ecoli_test_contigs.fasta";
    //string pathFastaReads = "data/EColi-synthetic/ecoli_test_reads.fasta";
    //string pathFastaOut = "data/EColi-synthetic/final.fasta";

    // string pathCR = "data/CJejuni-real/overlaps-c-r.paf";
    // string pathRR = "data/CJejuni-real/overlaps-r-r.paf";
    // string pathFastaCtgs = "data/CJejuni-real/CJejuni-contigs.fasta";
    // string pathFastaReads = "data/CJejuni-real/CJejuni-reads.fastq";
    // string pathFastaOut = "data/CJejuni-real/final.fasta";

    // string pathCR = "data/BGrahamii-real/overlaps-c-r.paf";
    // string pathRR = "data/BGrahamii-real/overlaps-r-r.paf";
    // string pathFastaCtgs = "data/BGrahamii-real/BGrahamii-contigs.fasta";
    // string pathFastaReads = "data/BGrahamii-real/BGrahamii-reads.fastq";
    // string pathFastaOut = "data/BGrahamii-real/final.fasta";

    vector<string> queryNames;
    vector<int> queryLens;
    vector<float> queryStarts;
    vector<float> queryEnds;

    vector<string> targetNames;
    vector<int> targetLens;
    vector<float> targetStarts;
    vector<float> targetEnds;

    vector<float> resMatches;
    vector<float> blockLens;

    vector<float> extensionSides; // 1 is right, 0 is left
    vector<float> SI;
    float SImin = 0.9;
    vector<float> OL1;
    vector<float> OL2;
    vector<float> OH1;
    vector<float> OH2;
    vector<float> EL1;
    vector<float> EL2;
    vector<float> OS;
    vector<float> ES1;
    vector<float> ES2;

    loadData(pathCR, queryNames, queryLens, queryStarts, queryEnds, targetNames, targetLens, targetStarts, targetEnds, resMatches, blockLens, SI, SImin, extensionSides);
    calculateOL(OL1, OL2, queryStarts, queryEnds, targetStarts, targetEnds);
    calculateOH(OH1, OH2, queryLens, queryStarts, queryEnds, targetLens, targetStarts, targetEnds, extensionSides);
    calculateEL(EL1, EL2, queryLens, queryStarts, queryEnds, targetLens, targetStarts, targetEnds, extensionSides);
    calculateOS(OS, OL1, OL2, SI);
    calculateES(ES1, ES2, OS, EL1, EL2, OH1, OH2);

    map<string, map<string, vector<vector<float> > > > groupedCR;
    for( int i = 0; i < queryNames.size(); i++) {
        vector<float> tmp = {extensionSides[i], ES2[i], OH1[i], OH2[i], EL1[i], EL2[i], OL2[i], OS[i], SI[i]};
        groupedCR[targetNames[i]][queryNames[i]].push_back(tmp);
    }
    vector<string> keysCR;
    for( auto key : groupedCR) {
        keysCR.push_back(key.first);
    }

    cout << "CR loaded" << endl;

    queryNames.clear();
    queryLens.clear();
    queryStarts.clear();
    queryEnds.clear();

    targetNames.clear();
    targetLens.clear();
    targetStarts.clear();
    targetEnds.clear();

    resMatches.clear();
    blockLens.clear();

    extensionSides.clear();
    SI.clear();
    OL1.clear();
    OL2.clear();
    OH1.clear();
    OH2.clear();
    EL1.clear();
    EL2.clear();
    OS.clear();
    ES1.clear();
    ES2.clear();

    loadData(pathRR, queryNames, queryLens, queryStarts, queryEnds, targetNames, targetLens, targetStarts, targetEnds, resMatches, blockLens, SI, SImin, extensionSides);
    calculateOL(OL1, OL2, queryStarts, queryEnds, targetStarts, targetEnds);
    calculateOH(OH1, OH2, queryLens, queryStarts, queryEnds, targetLens, targetStarts, targetEnds, extensionSides);
    calculateEL(EL1, EL2, queryLens, queryStarts, queryEnds, targetLens, targetStarts, targetEnds, extensionSides);
    calculateOS(OS, OL1, OL2, SI);
    calculateES(ES1, ES2, OS, EL1, EL2, OH1, OH2);

    map<string, map<string, vector<vector<float> > > > groupedRR;
    for( int i = 0; i < queryNames.size(); i++) {
        vector<float> tmp = {extensionSides[i], ES2[i], OH1[i], OH2[i], EL1[i], EL2[i], OS[i], SI[i]};
        groupedRR[targetNames[i]][queryNames[i]].push_back(tmp);
    }

    vector<string> keysRR;
    for( auto key : groupedRR) {
        keysRR.push_back(key.first);
    }

    cout << "RR loaded" << endl;

    map<float, vector<vector<tuple<string, int> > > > paths;
    map<float, vector<vector<tuple<string, int> > > > pathsTmp;

    //int maxDepth = 40;

    //int nTimes = 50;
    paths = monteCarloWrapper(keysCR, groupedCR, keysRR, groupedRR, maxDepth, nTimes);

    // int branchingFactor = 5;
    // int measureIndex = 6; // overlap score
    // pathsTmp = dfsApproach(keysCR, groupedCR, keysRR, groupedRR, maxDepth, branchingFactor, measureIndex);
    // for( auto side : paths) {
    //     for( auto path : pathsTmp[side.first]){
    //         side.second.push_back(path);
    //     }
    // }
    //
    // measureIndex = 1; // extension score
    // pathsTmp = dfsApproach(keysCR, groupedCR, keysRR, groupedRR, maxDepth, branchingFactor, measureIndex);
    // for( auto side : paths) {
    //     for( auto path : pathsTmp[side.first]){
    //         side.second.push_back(path);
    //     }
    // }


    // TODO concatenate them or change mapPaths to work with both paths
    map<tuple<string, string>, vector<vector<tuple<string, int> > > >  pathsMapLeft;
    map<tuple<string, string>, vector<vector<tuple<string, int> > > >  pathsMapRight;
    pathsMapLeft = mapPaths(0.0, paths);
    pathsMapRight = mapPaths(1.0, paths);

    cout << "Paths left: " << endl;
    for( auto key : pathsMapLeft){
        cout << get<0>(key.first) << " " << get<1>(key.first) << " paths:" << key.second.size() << endl;
    }
    cout << "Paths right: " << endl;
    for( auto key : pathsMapRight){
        cout << get<0>(key.first) << " " << get<1>(key.first) << " paths: " << key.second.size() << endl;
    }

    map<tuple<string, string>, vector<tuple<vector<tuple<string, int> >, float> > > pathLengthsMap;
    pathLengthsMap = calculatePathLengths(pathsMapLeft, groupedCR, groupedRR);

    // for( auto key : pathLengthsMap) {
    //     for( auto tapl : key.second) {
    //         float length = get<1>(tapl);
    //         cout << get<0>(key.first) << "->" << get<1>(key.first) << ":" << (int) length << endl;
    //     }
    // }



    vector<tuple<string, string> > scaffoldContigs;
    scaffoldContigs = getScaffoldContigs2(keysCR.size(), pathsMapLeft, pathsMapRight);

    // scaffoldContigs.clear();
    // scaffoldContigs.push_back(make_tuple("Ctg0", "Ctg3"));
    // scaffoldContigs.push_back(make_tuple("Ctg3", "Ctg2"));
    // scaffoldContigs.push_back(make_tuple("Ctg4", "Ctg5"));
    // scaffoldContigs.push_back(make_tuple("Ctg5", "Ctg1"));

    // for( auto tapl : scaffoldContigs) {
    //     cout << "da" << get<0>(tapl) << " " << get<1>(tapl) << endl;
    // }
    map<tuple<string, string>, vector<tuple<string, int> > > chosenPaths;
    chosenPaths = mapConsensusPath(dividePathsIntoGroups(pathLengthsMap, 10));
    // for( auto key : scaffoldContigs) {
    //     cout << get<0>(key) << " " << get<1>(key) << endl;
    //     chosenPaths[key] = pathsMapLeft[key][0];
    // }
    // for( auto tapl : chosenPaths) {
    //     cout << "ne" << get<0>(tapl.first) << " " << get<1>(tapl.first) << endl;
    // }
    vector<vector<tuple<string, int> > > finalOrder;
    finalOrder = buildFinalScaffoldOrder(chosenPaths, scaffoldContigs);

    string currentTarget;
    string currentQuery;
    int currentTargetIndex;
    int currentQueryIndex;
    cout << "Scaffold order:" << endl;
    for( auto pair : finalOrder) {
        currentTarget = get<0>(pair[0]);
        currentTargetIndex = get<1>(pair[1]);
        currentQuery = get<0>(pair[1]);
        currentQueryIndex = get<1>(pair[1]);

        cout << currentTarget << " " << currentTargetIndex << " " << currentQuery << " " << currentQueryIndex << endl;
    }
    cout << "Scaffold order end." << endl;

    map<string, string> fastaReads;
    fastaReads = loadFasta(pathFastaReads);

    map<string, string> fastaContigs;
    fastaContigs = loadFasta(pathFastaCtgs);

    string fastaString;
    fastaString = buildFastaString(finalOrder, groupedCR, groupedRR, fastaReads, fastaContigs, keysCR);

    saveFasta(fastaString, pathFastaOut);

    return 0;
}
