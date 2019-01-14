#include <iostream>
#include <vector>
#include <map>


using namespace std;

map<tuple<string, string>, vector<vector<tuple<string, int> > > > mapPaths(float extensionSide, map<float, vector<vector<tuple<string, int> > > > paths);

map<tuple<string, string>, vector<tuple<vector<tuple<string, int> >, float> > > calculatePathLengths(map<tuple<string, string>, vector<vector<tuple<string, int> > > > pathsMap, map<string, map<string, vector<vector<float> > > > groupedCR, map<string, map<string, vector<vector<float> > > > groupedRR);