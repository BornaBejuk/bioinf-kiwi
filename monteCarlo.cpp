#include <vector>
#include <map>
#include <iostream>
#include <typeinfo>
#include <random>
#include <algorithm>

#include "monteCarlo.h"

using namespace std;

vector<tuple<string, int> > monteCarlo(string start, float side, vector<string> keysCR, map<string, map<string, vector<vector<float> > > > groupedCR,
                        vector<string> keysRR, map<string, map<string, vector<vector<float> > > > groupedRR, int maxDepth){

    vector<tuple<string, int> > allPaths;
    int nTimes = 1000;
    for( int i = 0; i < nTimes; i++){
        cout << "i:" << i << endl;
        vector<tuple<string, int> > paths = depthFirstSearch(start, side, keysCR, groupedCR, keysRR, groupedRR, maxDepth);
        // vector<tuple<string, int> >::iterator row;
        // vector<string>::iterator col;
        // for (row = paths.begin(); row != paths.end(); row++) {
            // cout << paths[row][0] << '\n';
            // allPaths.push_back(paths[row]);
            // for (col = row->begin(); col != row->end(); col++) {
                // do stuff ...
            // }
        // }
    }

    return allPaths;
}


vector<tuple<string, int> > depthFirstSearch(string start, float side, vector<string> keysCR, map<string, map<string, vector<vector<float> > > > groupedCR,
                        vector<string> keysRR, map<string, map<string, vector<vector<float> > > > groupedRR, int maxDepth) {

    vector<tuple<string, int> > path;
    vector<string> stack;
    string read;
    int number;

    path.push_back(make_tuple(start, -1));

    tie(read, number) = getMCReadForContig(start, side, groupedCR);
    // cout << read << ' ' << number << endl;
    path.push_back(make_tuple(read, number));

    stack.push_back(read);
    while( !stack.empty()) {

        string currentTarget = stack.back();
        stack.pop_back();
        // cout << path.size() << endl;
        if( path.size() <= maxDepth) {
            if ( std::find(keysCR.begin(), keysCR.end(), currentTarget) != keysCR.end() ) {
                path.push_back(make_tuple(currentTarget, -1));
                cout << currentTarget << endl;
                return path;
            }
            else {
                tie(read, number) = getMCReadForRead(currentTarget, side, groupedCR, groupedRR);
                if( number == -1) {
                    break;
                }
                // cout << "Nasel sam:" << currentTarget << ' ' << read << number << endl;
                path.push_back(make_tuple(read, number));
                stack.push_back(read);
            }
        }
    }
    path = {make_tuple("", -2)};
    return path;
}

tuple<string, int> getMCReadForContig(string contig, float side, map<string, map<string, vector<vector<float> > > > groupedCR) {

    float sum = 0.0;

    for( auto query : groupedCR[contig]){
        for( int i = 0; i < query.second.size(); i++) {
            // pick ES
            if( side == query.second[i][0]){
                if( query.second[i][1] > 0) {
                    sum += query.second[i][1];
                }
            }
        }
    }

    float offset = 0.0;
    random_device rand_dev;
    mt19937 generator(rand_dev());
    uniform_int_distribution<int> distr(0, sum);
    int pick = distr(generator);

    for( auto query : groupedCR[contig]){
        for( int i = 0; i < query.second.size(); i++) {
            // pick ES
            if( side == query.second[i][0]){
                if( query.second[i][1] > 0) {
                    offset += query.second[i][1];
                    if( offset > pick) {
                        return make_tuple(query.first, i);
                    }
                }
            }
        }
    }

    return make_tuple("", -1);
}

tuple<string, int> getMCReadForRead(string read, float side, map<string, map<string, vector<vector<float> > > > groupedCR, map<string, map<string, vector<vector<float> > > > groupedRR) {

    // try to find contig == goal
    for( auto target : groupedCR) {
        for( auto query : target.second){
            if( query.first != read) {
                continue;
            }
            for( int i = 0; i < query.second.size(); i++) {
                // if extends from other side
                if( side != query.second[i][0]){
                    cout << "GOAL FOUND!" << target.first << endl;
                    return make_tuple(target.first, -1);
                }
            }
        }
    }


    float sum = 0.0;
    float offset = 0.0;
    random_device rand_dev;
    mt19937 generator(rand_dev());
    uniform_int_distribution<int> distr(0, sum);
    int pick = distr(generator);

    for( auto query : groupedRR[read]){
        for( int i = 0; i < query.second.size(); i++) {
            // cout << read << query.first << endl;
            // cout << side << " " << query.second[i][0] << endl;
            if( side == query.second[i][0]){
                if( query.second[i][1] > 0) {
                    offset += query.second[i][1];
                    if( offset > pick) {
                        // cout << "Read found!" << query.first << endl;
                        return make_tuple(query.first, i);
                    }
                }
            }
        }
    }

    return make_tuple("", -1);
}
