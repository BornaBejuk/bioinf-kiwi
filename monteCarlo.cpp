#include <vector>
#include <map>
#include <iostream>


#include "monteCarlo.h"

using namespace std;

vector<vector<string> > monteCarlo(string start, float side, vector<string> keysCR, map<string, map<string, vector<vector<float> > > > groupedCR,
                        vector<string> keysRR, map<string, map<string, vector<vector<float> > > > groupedRR, int maxDepth){

    vector<vector<string> > allPaths;
    int nTimes = 100;
    for( int i = 0; i < nTimes; i++){
        vector<vector<string> > paths = depthFirstSearch(start, side, keysCR, groupedCR, keysRR, groupedRR, maxDepth);
        vector<vector<string> >::iterator row;
        // vector<string>::iterator col;
        for (row = paths.begin(); row != paths.end(); row++) {
            // cout << paths[row][0] << '\n';
            // allPaths.push_back(paths[row]);
            // for (col = row->begin(); col != row->end(); col++) {
                // do stuff ...
            // }
        }
    }

    return allPaths;
}


vector<vector<string> > depthFirstSearch(string start, float side, vector<string> keysCR, map<string, map<string, vector<vector<float> > > > groupedCR,
                        vector<string> keysRR, map<string, map<string, vector<vector<float> > > > groupedRR, int maxDepth) {

    vector<vector<string> > test;
    return test;
}
