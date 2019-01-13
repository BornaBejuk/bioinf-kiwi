#include <iostream>
#include <vector>
#include <map>

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

        pathsMap[key].push_back(path);
    }

    return pathsMap;
}
