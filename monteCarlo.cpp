#include <vector>
#include <map>
#include <iostream>
#include <typeinfo>
#include <random>
#include <algorithm>

#include "monteCarlo.hpp"
#include "utils.hpp"

typedef std::chrono::high_resolution_clock Clock;

using namespace std;

// author: Karlo Brajdic

// wrapper function which calls monte-carlo for every contig
map<float, vector<vector<tuple<string, int> > > > monteCarloWrapper(vector<string> keysCR, map<string, map<string, vector<vector<float> > > > &groupedCR, vector<string> keysRR,
                                                                    map<string, map<string, vector<vector<float> > > > &groupedRR, int maxDepth, int nTimes) {

    map<float, vector<vector<tuple<string, int> > > > pathsMap;
    vector<vector<tuple<string, int> > > paths;

    float extensionSide = 0.0;

    for( auto key : keysCR) {
        paths = monteCarlo(key, extensionSide, keysCR, groupedCR, keysRR, groupedRR, maxDepth, nTimes);
        for( auto path : paths) {
            pathsMap[extensionSide].push_back(path);
        }
    }

    extensionSide = 1.0;
    for( auto key : keysCR) {
        paths = monteCarlo(key, extensionSide, keysCR, groupedCR, keysRR, groupedRR, maxDepth, nTimes);
        for( auto path : paths) {
            pathsMap[extensionSide].push_back(path);
        }
    }

    return pathsMap;
}

// function which repeats monte-carlo search nTimes for given starting contig
vector<vector<tuple<string, int> > > monteCarlo(string start, float side, vector<string> keysCR, map<string, map<string, vector<vector<float> > > > &groupedCR,
                        vector<string> keysRR, map<string, map<string, vector<vector<float> > > > &groupedRR, int maxDepth, int nTimes){

    vector<vector<tuple<string, int> > > paths;
    vector<tuple<string, int> > path;

    int flag = 0;

    for( int i = 0; i < nTimes; i++){
        flag = 0;
        if( groupedCR[start].size() == 0 ) {
            break;
        }
        path = mcSearch(start, side, keysCR, groupedCR, keysRR, groupedRR, maxDepth);
        for( auto p1 : paths) {
            // check if path already exists
            if( equal(p1.begin(), p1.end(), path.begin()) != 0) {
                flag = 1;
                break;
            }
        }
        // if path isn't a duplicate and size is bigger than 2, append it
        if( path.size() > 1 && flag == 0){
            paths.push_back(path);
        }
    }

    return paths;
}

// function which tries to find path between two contigs
vector<tuple<string, int> > mcSearch(string start, float side, vector<string> keysCR, map<string, map<string, vector<vector<float> > > > &groupedCR,
                        vector<string> keysRR, map<string, map<string, vector<vector<float> > > > &groupedRR, int maxDepth) {

    vector<tuple<string, int> > path;
    vector<string> stack;
    string read;
    int number;

    path.push_back(make_tuple(start, -1));

    tie(read, number) = getMCReadForContig(start, side, groupedCR, 1);

    path.push_back(make_tuple(read, number));
    stack.push_back(read);

    while( !stack.empty()) {
        string currentTarget = stack.back();
        stack.pop_back();
        if( path.size() <= maxDepth) {
            if ( std::find(keysCR.begin(), keysCR.end(), currentTarget) != keysCR.end() ) {
                // cout << "Current path length: " << path.size() << endl;
                return path;
            }
            else {
                while(true) {
                    tie(read, number) = getMCReadForRead(currentTarget, side, start, groupedCR, groupedRR, path, 1);
                    if( number == -2) {
                        break;
                    }
                    // check if we already visited that read
                    if ( std::find(path.begin(), path.end(), make_tuple(read ,number)) != path.end() ) {
                        continue;
                    } else {
                        break;
                    }
                }
                path.push_back(make_tuple(read, number));
                stack.push_back(read);
            }
        }
    }
    path = {make_tuple("", -2)};
    return path;
}

// function which returns read for given contig using roulette wheel probability scheme
tuple<string, int> getMCReadForContig(string contig, float side, map<string, map<string, vector<vector<float> > > > &groupedCR, int measureIndex) {

    float sum = 0.0;
    for( auto query : groupedCR[contig]){
        for( int i = 0; i < query.second.size(); i++) {
            // pick ES
            if( side == query.second[i][0]){
                if( query.second[i][measureIndex] > 0) {
                    sum += query.second[i][measureIndex];
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
                if( query.second[i][measureIndex] > 0) {
                    offset += query.second[i][measureIndex];
                    if( offset > pick) {
                        return make_tuple(query.first, i);
                    }
                }
            }
        }
    }

    return make_tuple("", -1);
}

// function which returns read for given read using roulette wheel probability scheme
tuple<string, int> getMCReadForRead(string read, float side, string startContig, map<string, map<string, vector<vector<float> > > > &groupedCR, map<string, map<string, vector<vector<float> > > > &groupedRR, vector<tuple<string, int> > &path, int measureIndex) {

    // try to find contig == goal
    for( auto target : groupedCR) {
        if( target.first == startContig) {
            continue;
        }
        for( auto query : target.second){
            if( query.first != read) {
                continue;
            }
            for( int i = 0; i < query.second.size(); i++) {
                // if extends from other side
                if( side != query.second[i][0]){
                    // remove the last read from path
                    path.pop_back();
                    // append the last read with corrected queryIndex
                    path.push_back(make_tuple(query.first, i));
                    return make_tuple(target.first, -1);
                }
            }
        }
    }

    float sum = 0.0;

    for( auto query : groupedRR[read]){
        for( int i = 0; i < query.second.size(); i++) {
            // pick ES
            if( side == query.second[i][0]){
                if( query.second[i][measureIndex] > 0) {
                    sum += query.second[i][measureIndex];
                }
            }
        }
    }

    float offset = 0.0;
    random_device rand_dev;
    mt19937 generator(rand_dev());
    uniform_int_distribution<int> distr(0, sum);
    int pick = distr(generator);

    for( auto query : groupedRR[read]){
        for( int i = 0; i < query.second.size(); i++) {
            if( side == query.second[i][0]){
                if( query.second[i][measureIndex] > 0) {
                    offset += query.second[i][measureIndex];
                    if( offset > pick) {
                        return make_tuple(query.first, i);
                    }
                }
            }
        }
    }

    return make_tuple("", -2);
}
