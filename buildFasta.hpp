#include <vector>
#include <map>

using namespace std;

string buildFastaString(vector<vector<tuple<string, int> > > finalOrder, map<string, map<string, vector<vector<float> > > > &groupedCR,
                        map<string, map<string, vector<vector<float> > > > &groupedRR, map<string, string> &fastaReads, map<string, string> &fastaContigs, vector<string> keysCR);


void saveFasta(string fasta, string path);
