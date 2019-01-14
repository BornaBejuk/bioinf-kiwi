#include <vector>

using namespace std;

void loadData(string path, vector<string> &queryNames, vector<int> &queryLens, vector<float> &queryStarts,
            vector<float> &queryEnds, vector<string> &targetNames, vector<int> &targetLens,
            vector<float> &targetStarts, vector<float> &targetEnds, vector<float> &resMatches,
            vector<float> &blockLens, vector<float> &SI, float SImin);

void filterContained(vector<string> &queryNames, vector<int> &queryLens, vector<float> &queryStarts,
            vector<float> &queryEnds, vector<string> &targetNames, vector<int> &targetLens,
            vector<float> &targetStarts, vector<float> &targetEnds, vector<float> &resMatches,
            vector<float> &blockLens, vector<float> &extensionSides);

void calculateScores(vector<float> &OS, vector<float> &ES, vector<float> &SI, vector<float> &EL,
                    vector<float> &OH, vector<float> &resMatches, vector<float> &blockLens);

void calculateSI(vector<float> &SI, vector<float> &resMatches, vector<float> &blockLens);

void filterBySI(float SImin, vector<string> &queryNames, vector<int> &queryLens, vector<float> &queryStarts,
            vector<float> &queryEnds, vector<string> &targetNames, vector<int> &targetLens,
            vector<float> &targetStarts, vector<float> &targetEnds, vector<float> &resMatches,
            vector<float> &blockLens, vector<float> &extensionSides, vector<float> SI);

void calculateOL(vector<float> &OL1, vector<float> &OL2, vector<float> &queryStarts, vector<float> &queryEnds, vector<float> &targetStarts, vector<float> &targetEnds);

void calculateOH(vector<float> &OH1, vector<float> &OH2, vector<int> &queryLens, vector<float> &queryStarts, vector<float> &queryEnds,
    vector<int> &targetLens, vector<float> &targetStarts, vector<float> &targetEnds, vector<float> &extensionSides);

void calculateEL(vector<float> &EL1, vector<float> &EL2, vector<int> &queryLens, vector<float> &queryStarts, vector<float> &queryEnds,
        vector<int> &targetLens, vector<float> &targetStarts, vector<float> &targetEnds, vector<float> &extensionSides);

void calculateOS(vector<float> &OS, vector<float> &OL1, vector<float> &OL2, vector<float> &SI);

void calculateES(vector<float> &ES1, vector<float> &ES2, vector<float> &OS, vector<float> &EL1, vector<float> &EL2, vector<float> &OH1, vector<float> &OH2);

map<string, string> loadFasta(string path);
