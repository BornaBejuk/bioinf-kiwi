#include <vector>
#include <map>

using namespace std;

string build_fasta_file(vector<tuple<string, int> > final_path,
                        map<string, map<string, vector<vector<float> > > > groupedCR, 
                        vector<string> keysCR, 
                        map<string, map<string, vector<vector<float> > > > groupedRR,
                        vector<string> keysRR,
                        map<string, string> fastaReads,
                        map<string, string> fastaContigs);