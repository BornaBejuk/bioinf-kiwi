#include <vector>
#include <map>

using namespace std;

vector<tuple<string, int> > monteCarlo(string start, float side, vector<string> keysCR, map<string, map<string, vector<vector<float> > > > groupedCR,
                        vector<string> keysRR, map<string, map<string, vector<vector<float> > > > groupedRR, int maxDepth);

vector<tuple<string, int> > depthFirstSearch(string start, float side, vector<string> keysCR, map<string, map<string, vector<vector<float> > > > groupedCR,
                        vector<string> keysRR, map<string, map<string, vector<vector<float> > > > groupedRR, int maxDepth);

tuple<string, int> getMCReadForContig(string contig, float side, map<string, map<string, vector<vector<float> > > > groupedCR);

tuple<string, int> getMCReadForRead(string read, float side, map<string, map<string, vector<vector<float> > > > groupedCR, map<string, map<string, vector<vector<float> > > > groupedRR);
