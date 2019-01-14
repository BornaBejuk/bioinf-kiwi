#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>

#include "utils.h"
#include "monteCarlo.h"
#include "selectPaths.h"
#include "scaffolding.h"

using namespace std;

int main() {

    string pathCR = "data/EColi-synthetic/overlaps-c-r.paf";
    string pathRR = "data/EColi-synthetic/overlaps-r-r.paf";
    string pathFastaCtgs = "data/EColi-synthetic/ecoli_test_contigs.fasta";
    string pathFastaReads = "data/EColi-synthetic/ecoli_test_reads.fasta";

    // string pathCR = "data/CJejuni-real/overlaps-c-r.paf";
    // string pathRR = "data/CJejuni-real/overlaps-r-r.paf";
    // string pathFastaCtgs = "data/CJejuni-real/CJejuni-contigs.fasta";
    // string pathFastaReads = "data/CJejuni-real/CJejuni-contigs.fasta";

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
    float SImin = 0.7;
    vector<float> OL1;
    vector<float> OL2;
    vector<float> OH1;
    vector<float> OH2;
    vector<float> EL1;
    vector<float> EL2;
    vector<float> OS;
    vector<float> ES1;
    vector<float> ES2;

    loadData(pathCR, queryNames, queryLens, queryStarts, queryEnds, targetNames, targetLens, targetStarts, targetEnds, resMatches, blockLens, SI, SImin);
    filterContained(queryNames, queryLens, queryStarts, queryEnds, targetNames, targetLens, targetStarts, targetEnds, resMatches, blockLens, extensionSides);

    calculateSI(SI, resMatches, blockLens);
    filterBySI(SImin, queryNames, queryLens, queryStarts, queryEnds, targetNames, targetLens, targetStarts, targetEnds, resMatches, blockLens, extensionSides, SI);
    calculateOL(OL1, OL2, queryStarts, queryEnds, targetStarts, targetEnds);
    calculateOH(OH1, OH2, queryLens, queryStarts, queryEnds, targetLens, targetStarts, targetEnds, extensionSides);
    calculateEL(EL1, EL2, queryLens, queryStarts, queryEnds, targetLens, targetStarts, targetEnds, extensionSides);
    calculateOS(OS, OL1, OL2, SI);
    calculateES(ES1, ES2, OS, EL1, EL2, OH1, OH2);

    map<string, map<string, vector<vector<float> > > > groupedCR;
    for( int i = 0; i < queryNames.size(); i++) {
        vector<float> tmp = {extensionSides[i], ES2[i], OH1[i], OH2[i], EL1[i], EL2[i], OL2[i]};
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

    loadData(pathRR, queryNames, queryLens, queryStarts, queryEnds, targetNames, targetLens, targetStarts, targetEnds, resMatches, blockLens, SI, SImin);
    // cout << queryNames.size() << endl;
    filterContained(queryNames, queryLens, queryStarts, queryEnds, targetNames, targetLens, targetStarts, targetEnds, resMatches, blockLens, extensionSides);
    // cout << queryNames.size() << endl;
    // calculateSI(SI, resMatches, blockLens);
    // cout << "SI loaded" << endl;
    // filterBySI(0.3, queryNames, queryLens, queryStarts, queryEnds, targetNames, targetLens, targetStarts, targetEnds, resMatches, blockLens, extensionSides, SI);

    calculateOL(OL1, OL2, queryStarts, queryEnds, targetStarts, targetEnds);
    calculateOH(OH1, OH2, queryLens, queryStarts, queryEnds, targetLens, targetStarts, targetEnds, extensionSides);
    calculateEL(EL1, EL2, queryLens, queryStarts, queryEnds, targetLens, targetStarts, targetEnds, extensionSides);
    calculateOS(OS, OL1, OL2, SI);
    calculateES(ES1, ES2, OS, EL1, EL2, OH1, OH2);

    map<string, map<string, vector<vector<float> > > > groupedRR;
    for( int i = 0; i < queryNames.size(); i++) {
        vector<float> tmp = {extensionSides[i], ES2[i], OH1[i], OH2[i], EL1[i], EL2[i]};
        groupedRR[targetNames[i]][queryNames[i]].push_back(tmp);
    }

    vector<string> keysRR;
    for( auto key : groupedRR) {
        keysRR.push_back(key.first);
    }

    cout << "RR loaded" << endl;
    // cout << groupedRR["m161108_211237_00127_c101051402550000001823235612291637_s1_p0/100336/0_19619"]["m161103_175158_00127_c101051712550000001823235612291635_s1_p0/140004/0_2390"][0][0] << endl;
    // for( auto x : groupedCR){
    //     cout << x.first << " contains:" << endl;
    //     for( auto y : x.second){
    //         cout << y.first << ':' << y.second[0][1] << '\n';
    //         // cout << y.second << '\n';
    //         break;
    //     }
    //     break;
    // //     cout << x.first << ' '  << x.second << ' ' << '\n';
    // }
    //
    // for( auto x : groupedRR){
    //     cout << x.first << " contains:" << endl;
    //     for( auto y : x.second){
    //         cout << y.first << ':' << y.second[0][1] << '\n';
    //         // cout << y.second << '\n';
    //         break;
    //     }
    //     break;
    // //     cout << x.first << ' '  << x.second << ' ' << '\n';
    // }

    map<float, vector<vector<tuple<string, int> > > > paths;
    int maxDepth = 50;
    int nTimes = 50;
    // keysCR.clear();
    // keysCR.push_back("ctg1");
    // monteCarlo("ctg2", 0.0, keysCR, groupedCR, keysRR, groupedRR, maxDepth, nTimes);

    paths = monteCarloWrapper(keysCR, groupedCR, keysRR, groupedRR, maxDepth, nTimes);
    // string read;
    // int number;
    // for( auto side : paths){
    //     cout << side.first << " contains:" << endl;
    //     for( auto p : side.second){
    //         for( auto element : p) {
    //             tie(read, number) = element;
    //             cout << read << " " << number << '\n';
    //         }
    //         cout << endl;
    //     }
    // }


    // TODO concatenate them or change mapPaths to work with both paths
    map<tuple<string, string>, vector<vector<tuple<string, int> > > >  pathsMapLeft;
    map<tuple<string, string>, vector<vector<tuple<string, int> > > >  pathsMapRight;
    pathsMapLeft = mapPaths(0.0, paths);
    pathsMapRight = mapPaths(1.0, paths);

    // paths are now in one map, regardless of extension side
    map<tuple<string, string>, vector<tuple<vector<tuple<string, int> >, float> > > pathLengthsMap;
    pathLengthsMap = calculatePathLengths(pathsMapLeft, groupedCR, groupedRR);

    // for( auto key : pathLengthsMap) {
    //     for( auto tapl : key.second) {
    //         float length = get<1>(tapl);
    //         cout << get<0>(key.first) << "->" << get<1>(key.first) << ":" << (int) length << endl;
    //     }
    // }

    // string read;
    // int number;
    for( auto key : pathsMapLeft){
        cout << get<0>(key.first) << " " << get<1>(key.first) << " paths:" << endl;
        // for( auto p : key.second){
        //     for( auto element : p) {
        //         tie(read, number) = element;
        //         cout << read << " " << number << '\n';
        //     }
        //     cout << endl;
        // }
    }


    vector<tuple<string, string> > scaffoldContigs;
    scaffoldContigs = getScaffoldContigs(keysCR.size(), pathsMapLeft);

    map<tuple<string, string>, vector<tuple<string, int> > > chosenPaths;
    for( auto key : scaffoldContigs) {
        chosenPaths[key] = pathsMapLeft[key][0];
    }

    for( auto tapl : scaffoldContigs) {
        cout << "da" << get<0>(tapl) << " " << get<1>(tapl) << endl;
    }
    for( auto tapl : chosenPaths) {
        cout << "ne" << get<0>(tapl.first) << " " << get<1>(tapl.first) << endl;
    }

    vector<vector<tuple<string, int> > > finalOrder;
    finalOrder = buildFinalScaffoldOrder(chosenPaths);

    map<string, string> fastaReads;
    // fastaReads = loadFasta(pathFastaReads);

    // allPaths = monteCarlo(keysCR[2], 0, keysCR, groupedCR, keysRR, groupedRR, maxDepth);
    // cout << proba[targetNames[0]][queryNames[0]][0] << '\n';
    // cout << queryLens[0] << '\n';
    // cout << targetLens[0] << '\n';
    // cout << extensionSides[0] << '\n';
    // cout << OL1[0] << '\n';
    // cout << OH1[0] << ' ' << OH2[0] << '\n';
    // cout << EL1[0] << ' ' << EL2[0] << '\n';
    // cout << SI[0] << '\n';

    return 0;
}
