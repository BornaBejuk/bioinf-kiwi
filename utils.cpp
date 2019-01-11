#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

using namespace std;

void loadData(vector<string> &queryNames, vector<int> &queryLens, vector<float> &queryStarts,
            vector<float> &queryEnds, vector<string> &targetNames, vector<int> &targetLens,
            vector<float> &targetStarts, vector<float> &targetEnds, vector<float> &resMatches,
            vector<float> &blockLens) {

    string line;

    ifstream inputFile("data/EColi-synthetic/overlaps-c-r.paf");

    if (inputFile.good()) {
        int current_number = 0;
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

            queryNames.push_back(qName);
            queryLens.push_back(qLen);
            queryStarts.push_back(qStart);
            queryEnds.push_back(qEnd);
            targetLens.push_back(tLen);
            targetNames.push_back(tName);
            targetStarts.push_back(tStart);
            targetEnds.push_back(tEnd);
            resMatches.push_back(resMatch);
            blockLens.push_back(bLen);

        }
    }
}
