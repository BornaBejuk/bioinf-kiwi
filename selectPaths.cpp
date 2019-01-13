#include <iostream>
#include <vector>
#include <map>
#include <typeinfo>

#include "selectPaths.h"

using namespace std;

map<tuple<string, string>, vector<vector<tuple<string, int> > > > mapPaths(float extensionSide, map<float, vector<vector<tuple<string, int> > > > paths) {

    map<tuple<string, string>, vector<vector<tuple<string, int> > > > pathsMap;
    string read;
    int number;
    string start;
    string end;
    tuple<string, string> key;
    for( auto path : paths[extensionSide]){
        // make key
        tie(start, number) = path.front();
        tie(end, number) = path.back();
        path.erase(path.begin());
        path.pop_back();
        key = make_tuple(start, end);
        // cout << get<0>(path[0]) << endl;
        pathsMap[key].push_back(path);
    }

    return pathsMap;
}

// map[(ctg1,ctg2)] = [([(read1, 0), (read2, 0), (read3, 0)], 12345)]
map<tuple<string, string>, vector<tuple<vector<tuple<string, int> >, float> > > calculatePathLengths(map<tuple<string, string>, vector<vector<tuple<string, int> > > > pathsMap, map<string, map<string, vector<vector<float> > > > groupedCR, map<string, map<string, vector<vector<float> > > > groupedRR) {

    map<tuple<string, string>, vector<tuple<vector<tuple<string, int> >, float> > > pathLengthsMap;

    string begin;
    string end;
    float ctgBeginLen;
    float ctgEndLen;
    float pathLen = 0.0;
    string currentTargetRead;
    int currentTargetIndex;
    string currentQueryRead;
    int currentQueryIndex;
    for( auto key : pathsMap) {
        begin = get<0>(key.first);
        end = get<1>(key.first);
        // cout << begin << " " << end << endl;

        for( auto path : key.second) {
            pathLen = 0.0;
            currentQueryRead = get<0>(path.front());
            currentQueryIndex = get<1>(path.front());
            // len = OL2 + EL2 + EL1
            ctgBeginLen = groupedCR[begin][currentQueryRead][currentQueryIndex][6] + groupedCR[begin][currentQueryRead][currentQueryIndex][5] + groupedCR[begin][currentQueryRead][currentQueryIndex][4];

            currentQueryRead = get<0>(path.back());
            currentQueryIndex = get<1>(path.back());

            ctgEndLen = groupedCR[end][currentQueryRead][currentQueryIndex][6] + groupedCR[end][currentQueryRead][currentQueryIndex][5] + groupedCR[end][currentQueryRead][currentQueryIndex][4];

            pathLen += ctgBeginLen + ctgEndLen;

            for( int i = 0; i < (path.size() - 1); i++){
                currentTargetRead = get<0>(path[i]);
                currentQueryRead = get<0>(path[i+1]);
                currentQueryIndex = get<1>(path[i+1]);

                // - OH2 + EL1
                pathLen += groupedRR[currentTargetRead][currentQueryRead][currentQueryIndex][3] + groupedRR[currentTargetRead][currentQueryRead][currentQueryIndex][4];
            }
            // cout << (int) pathLen << endl;

            pathLengthsMap[key.first].push_back(make_tuple(path, pathLen));
        }

    }

    return pathLengthsMap;
}
