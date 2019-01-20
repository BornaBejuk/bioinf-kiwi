#include <vector>
#include <map>

using namespace std;

map<float, vector<vector<tuple<string, int> > > > monteCarloWrapper(vector<string> keysCR, map<string, map<string, vector<vector<float> > > > &groupedCR, vector<string> keysRR,
                                                                    map<string, map<string, vector<vector<float> > > > &groupedRR, int maxDepth, int nTimes);

vector<vector<tuple<string, int> > > monteCarlo(string start, float side, vector<string> keysCR, map<string, map<string, vector<vector<float> > > > &groupedCR,
                        vector<string> keysRR, map<string, map<string, vector<vector<float> > > > &groupedRR, int maxDepth, int nTimes);

vector<tuple<string, int> > mcSearch(string start, float side, vector<string> keysCR, map<string, map<string, vector<vector<float> > > > &groupedCR,
                        vector<string> keysRR, map<string, map<string, vector<vector<float> > > > &groupedRR, int maxDepth);

tuple<string, int> getMCReadForContig(string contig, float side, map<string, map<string, vector<vector<float> > > > &groupedCR, int measureIndex);

tuple<string, int> getMCReadForRead(string read, float side, string startContig, map<string, map<string, vector<vector<float> > > > &groupedCR, map<string, map<string, vector<vector<float> > > > &groupedRR, vector<tuple<string, int> > &path, int measureIndex);
