#include <vector>
#include <map>

using namespace std;

vector<vector<string> > monteCarlo(string start, float side, vector<string> keysCR, map<string, map<string, vector<vector<float> > > > groupedCR,
                        vector<string> keysRR, map<string, map<string, vector<vector<float> > > > groupedRR, int maxDepth);

vector<vector<string> > depthFirstSearch(string start, float side, vector<string> keysCR, map<string, map<string, vector<vector<float> > > > groupedCR,
                        vector<string> keysRR, map<string, map<string, vector<vector<float> > > > groupedRR, int maxDepth);
