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


    loadData(queryNames, queryLens, queryStarts, queryEnds, targetNames, targetLens, targetStarts, targetEnds, resMatches, blockLens);

    cout << queryNames[0] << '\n';
    return 0;
}
