#include <vector>
#include <map>

using namespace std;

map<float, vector<vector<tuple<string, int> > > > dfsApproach(vector<string> keysCR, map<string, map<string, vector<vector<float> > > > &groupedCR, vector<string> keysRR, map<string, map<string, vector<vector<float> > > > &groupedRR, int maxDepth, int branchingFactor, int measureIndex);

vector<vector<tuple<string, int> > > dfsWrapper(string start, float side, vector<string> keysCR, map<string, map<string, vector<vector<float> > > > &groupedCR, vector<string> keysRR, map<string, map<string, vector<vector<float> > > > &groupedRR, int maxDepth, int branchingFactor, int measureIndex);

vector<tuple<string, int> > dfsSearch(string start, float side, vector<string> keysCR, map<string, map<string, vector<vector<float> > > > &groupedCR, vector<string> keysRR, map<string, map<string, vector<vector<float> > > > &groupedRR, int maxDepth, int measureIndex, int skip, int branchingFactor);

tuple<string, int> getReadForContig(string contig, float side, map<string, map<string, vector<vector<float> > > > &groupedCR, int skip);

tuple<string, int> getReadForRead(string read, float side, string startContig, map<string, map<string, vector<vector<float> > > > &groupedCR, map<string, map<string, vector<vector<float> > > > &groupedRR, vector<tuple<string, int> > &path, int measureIndex, vector<tuple<string, int> > &skipping);
