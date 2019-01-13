#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>

#include "utils.h"
#include "monteCarlo.h"
#include "selectPaths.h"

using namespace std;

int main() {

    string pathCR = "data/EColi-synthetic/overlaps-c-r.paf";
    string pathRR = "data/EColi-synthetic/overlaps-r-r.paf";

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

    filterContained(queryNames, queryLens, queryStarts, queryEnds, targetNames, targetLens, targetStarts, targetEnds, resMatches, blockLens, extensionSides);
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
    int maxDepth = 30;
    int nTimes = 1;
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

    map<tuple<string, string>, vector<vector<tuple<string, int> > > >  pathsMap;
    pathsMap = mapPaths(0.0, paths);

    map<tuple<string, string>, vector<tuple<vector<tuple<string, int> >, float> > > pathLengthsMap;
    pathLengthsMap = calculatePathLengths(pathsMap, groupedCR, groupedRR);

    // for( auto key : pathLengthsMap) {
    //     for( auto tapl : key.second) {
    //         float length = get<1>(tapl);
    //         cout << get<0>(key.first) << "->" << get<1>(key.first) << ":" << (int) length << endl;
    //     }
    // }

    // string read;
    // int number;
    // for( auto key : pathsMap){
    //     cout << get<0>(key.first) << " " << get<1>(key.first) << " paths:" << endl;
    //     for( auto p : key.second){
    //         for( auto element : p) {
    //             tie(read, number) = element;
    //             cout << read << " " << number << '\n';
    //         }
    //         cout << endl;
    //     }
    // }

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
