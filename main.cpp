#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>

#include "utils.h"

using namespace std;

int main() {
    // cout << "Hellou" << endl;
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

    vector<string> extensionSides;
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

    loadData(queryNames, queryLens, queryStarts, queryEnds, targetNames, targetLens, targetStarts, targetEnds, resMatches, blockLens);
    filterContained(queryNames, queryLens, queryStarts, queryEnds, targetNames, targetLens, targetStarts, targetEnds, resMatches, blockLens, extensionSides);

    calculateSI(SI, resMatches, blockLens);
    filterBySI(SImin, queryNames, queryLens, queryStarts, queryEnds, targetNames, targetLens, targetStarts, targetEnds, resMatches, blockLens, extensionSides, SI);
    calculateOL(OL1, OL2, queryStarts, queryEnds, targetStarts, targetEnds);
    calculateOH(OH1, OH2, queryLens, queryStarts, queryEnds, targetLens, targetStarts, targetEnds, extensionSides);
    calculateEL(EL1, EL2, queryLens, queryStarts, queryEnds, targetLens, targetStarts, targetEnds, extensionSides);
    calculateOS(OS, OL1, OL2, SI);
    calculateES(ES1, ES2, OS, EL1, EL2, OH1, OH2);

    map<string, map<string, vector<float> > > proba;
    proba[targetNames[0]][queryNames[0]].push_back(ES1[0]);
    for( auto x : proba){
        cout << x.first << " contains:" << endl;
        for( auto y : x.second){
            cout << y.first << ':' << y.second[0] << '\n';
            // cout << y.second << '\n';
        }
    //     cout << x.first << ' '  << x.second << ' ' << '\n';
    }
    cout << proba[targetNames[0]][queryNames[0]][0] << '\n';
    // cout << queryLens[0] << '\n';
    // cout << targetLens[0] << '\n';
    // cout << extensionSides[0] << '\n';
    // cout << OL1[0] << '\n';
    // cout << OH1[0] << ' ' << OH2[0] << '\n';
    // cout << EL1[0] << ' ' << EL2[0] << '\n';
    // cout << SI[0] << '\n';

    return 0;
}
