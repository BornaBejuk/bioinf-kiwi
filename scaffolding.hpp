#include <map>
#include <vector>

using namespace std;

bool checkVector0(vector<tuple<string, string> > vec, string value);

bool checkVector1(vector<tuple<string, string> > vec, string value);

vector<tuple<string, string> > getScaffoldContigs2(int ctgNumber, map<tuple<string, string>, vector<vector<tuple<string, int> > > > &pathsMapLeft, map<tuple<string, string>, vector<vector<tuple<string, int> > > > &pathsMapRight);

vector<vector<tuple<string, int> > > buildFinalScaffoldOrder(map<tuple<string, string>, vector<tuple<string, int> > > chosenPaths, vector<tuple<string, string> > scaffoldContigs);
