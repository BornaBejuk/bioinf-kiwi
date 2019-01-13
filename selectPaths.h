#include <iostream>
#include <vector>
#include <map>


using namespace std;

map<tuple<string, string>, vector<vector<tuple<string, int> > > > mapPaths(float extensionSide, map<float, vector<vector<tuple<string, int> > > > paths);
