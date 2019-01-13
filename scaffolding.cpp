#include <vector>
#include <map>
#include <iostream>
#include <typeinfo>
#include <random>
#include <algorithm>

#include "scaffolding.h"


// returns one path, map[ctg1,ctg2] = [(read1,0),(read2,0)...]
// returns tuples for path, e.g. [(ctg1,ctg2), (ctg2,ctg3)]
vector<tuple<string, string> > getScaffoldContigs(int ctgNumber, map<tuple<string, string>, vector<vector<tuple<string, int> > > > pathsMap) {

    vector<tuple<string, string> > contigs;
    map<tuple<string, string>, int> pathCount;

    for( auto key : pathsMap) {
        pathCount[key.first] = key.second.size();
    }

    if( pathsMap.size() < ctgNumber) {
        // return all
        for( auto key : pathsMap) {
            contigs.push_back(key.first);
        }
    } else {
        // find lowest branch and dont include it
        int lowest = 1000000;
        tuple<string,string> lowestTuple;
        for( auto key : pathCount) {
            if( key.second < lowest) {
                lowest = key.second;
                lowestTuple = key.first;
            }
        }

        for( auto key : pathsMap) {
            if( key.first != lowestTuple){
                contigs.push_back(key.first);
            }
        }
    }

    return contigs;
}
