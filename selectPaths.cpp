#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <typeinfo>

#include "selectPaths.hpp"

using namespace std;

// author: Karlo Brajdic
// maps paths according to side they were extended
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

// author: Karlo Brajdic
// calculates path length for all paths
map<tuple<string, string>, vector<tuple<vector<tuple<string, int> >, float> > > calculatePathLengths(map<tuple<string, string>, vector<vector<tuple<string, int> > > > pathsMap, map<string, map<string, vector<vector<float> > > > &groupedCR, map<string, map<string, vector<vector<float> > > > &groupedRR) {

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
                pathLen += groupedRR[currentTargetRead][currentQueryRead][currentQueryIndex][4] - groupedRR[currentTargetRead][currentQueryRead][currentQueryIndex][3];
            }
            // cout << (int) pathLen << endl;

            pathLengthsMap[key.first].push_back(make_tuple(path, pathLen));
        }

    }

    return pathLengthsMap;
}


// author: Borna Bejuk
// divides paths into groups for each key in map pathLengthsMap
map<tuple<string, string>, vector<vector<tuple<vector<tuple<string, int>>, float> > > > dividePathsIntoGroups(map<tuple<string, string>, vector<tuple<vector<tuple<string, int>>, float> > > pathLengthsMap, int smallestGroupNumber) {

    map<tuple<string, string>, vector<vector<tuple<vector<tuple<string, int>>, float>>>> mapOfGroups;

    for (auto key : pathLengthsMap) {
        vector<float> lengthsOfPaths;
        for (auto path : key.second) {
            lengthsOfPaths.push_back(get<1>(path));
        }
        sort(lengthsOfPaths.begin(), lengthsOfPaths.end());
        float lowest = lengthsOfPaths[0];
        float highest = lengthsOfPaths[lengthsOfPaths.size() - 1];
        vector<vector<tuple<vector<tuple<string, int>>, float>>> groups;

        float groupWindow = (float) (highest - lowest) / lengthsOfPaths.size();
        if (key.second.size() < smallestGroupNumber) {
            vector<tuple<vector<tuple<string, int>>, float>> first_group;
            first_group = pathLengthsMap[key.first];
            mapOfGroups[key.first].push_back(first_group);
            continue;
        }
        for (int i = 0; i < lengthsOfPaths.size(); i++) {
            vector<tuple<vector<tuple<string, int>>, float>> temp;
            groups.push_back(temp);
        }
        for (auto path : key.second) {
            float length = get<1>(path);
            if (length == highest) {
                length -= 1;
            }
            length = length - lowest;
            int position = (int) (length / groupWindow);
            groups[position].push_back(path);
        }
        mapOfGroups[key.first] = groups;
    }
    return mapOfGroups;
}


// author: Borna Bejuk
// compares two tuples
bool sortbysec(const tuple<vector<tuple<string, int> >, float>& a, const tuple<vector<tuple<string, int> >, float>& b) {
    return (get<1>(a) < get<1>(b));
}


// author: Borna Bejuk
// chooses the path for each pair of contigs
map<tuple<string,string>, vector<tuple<string, int> > > mapConsensusPath(map<tuple<string, string>, vector<vector<tuple<vector<tuple<string, int>>, float> > > > mapOfGroups) {

    map<tuple<string,string>, vector<tuple<string, int>>> mapOfChosenPaths;
    for (auto key: mapOfGroups) {
        //vector<vector<tuple<vector<tuple<string, int>>, float>>> groups = mapOfGroups[key.first];
        vector<int> vectorOfGroupSizes;
        for (auto group : key.second) {
            vectorOfGroupSizes.push_back(group.size());
        }
        int max_index = distance(vectorOfGroupSizes.begin(), max_element(vectorOfGroupSizes.begin(), vectorOfGroupSizes.end()));
        vector<tuple<vector<tuple<string, int> >, float> > paths = key.second[max_index];
        int numOfPaths = paths.size();
        int index = (int) numOfPaths / 2;
        sort(paths.begin(), paths.end(), sortbysec);
        mapOfChosenPaths[key.first] = get<0>(paths[index]);
    }
    return mapOfChosenPaths;

}
