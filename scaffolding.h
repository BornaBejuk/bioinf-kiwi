#include <map>
#include <vector>

using namespace std;

vector<tuple<string, string> > getScaffoldContigs(int ctgNumber, map<tuple<string, string>, vector<vector<tuple<string, int> > > >  pathsMap);
