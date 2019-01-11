#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
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
    cout << queryEnds.size() << '\n';
    filterContained(queryNames, queryLens, queryStarts, queryEnds, targetNames, targetLens, targetStarts, targetEnds, resMatches, blockLens, extensionSides);
    cout << queryEnds.size() << '\n';
    calculateSI(SI, resMatches, blockLens);
    filterBySI(SImin, queryNames, queryLens, queryStarts, queryEnds, targetNames, targetLens, targetStarts, targetEnds, resMatches, blockLens, extensionSides, SI);
    calculateOL(OL1, OL2, queryStarts, queryEnds, targetStarts, targetEnds);
    calculateOH(OH1, OH2, queryLens, queryStarts, queryEnds, targetLens, targetStarts, targetEnds, extensionSides);
    calculateOH(EL1, EL2, queryLens, queryStarts, queryEnds, targetLens, targetStarts, targetEnds, extensionSides);
    calculateOS(OS, OL1, OL2, SI);
    calculateES(ES1, ES2, OS, EL1, EL2, OH1, OH2);
    // cout << queryLens[0] << '\n';
    // cout << targetLens[0] << '\n';
    // cout << extensionSides[0] << '\n';
    // cout << OL1[0] << '\n';
    // cout << OH1[0] << '\n';
    // cout << EL1[0] << '\n';
    // cout << SI[0] << '\n';

    return 0;
}
