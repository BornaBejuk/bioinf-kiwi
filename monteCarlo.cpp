#include <vector>
#include <map>

using namespace std;

vector<vector<string> > monteCarlo(string start, string side, vector<string> keysCR, map<string, map<string, vector<float> > > groupedCR,
                        vector<string> keysRR, map<string, map<string, vector<float> > > groupedRR, int maxDepth){

    vector<vector<string> > allPaths;
    int nTimes = 100;
    for( int i = 0; i < nTimes; i++){
        vector<vector<string> > paths = depthFirstSearch(start, side, keysCR, groupedCR, keysRR, groupedRR, maxDepth);
        vector<vector<string> >::iterator row;
        // vector<string>::iterator col;
        for (row = paths.begin(); row != paths.end(); row++) {
            allPaths.push_back(paths[row]);
            // for (col = row->begin(); col != row->end(); col++) {
            //     // do stuff ...
            // }
        }
    }

    return allPaths;
}
