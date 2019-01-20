#include <map>
#include <vector>

using namespace std;

vector<tuple<string, string> > getScaffoldContigs2(int ctgNumber, map<tuple<string, string>, vector<vector<tuple<string, int> > > > &pathsMapLeft, map<tuple<string, string>, vector<vector<tuple<string, int> > > > &pathsMapRight);

vector<tuple<string, string> > getScaffoldContigs(int ctgNumber, map<tuple<string, string>, vector<vector<tuple<string, int> > > >  pathsMap);

vector<vector<tuple<string, int> > > buildFinalScaffoldOrder(map<tuple<string, string>, vector<tuple<string, int> > > chosenPaths, vector<tuple<string, string> > scaffoldContigs);
