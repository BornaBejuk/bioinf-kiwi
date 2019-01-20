#include <vector>
#include <map>
#include <iostream>
#include <typeinfo>
#include <random>
#include <algorithm>

#include "monteCarlo.hpp"
#include "utils.hpp"
#include "dfsSearch.hpp"

typedef std::chrono::high_resolution_clock Clock;

using namespace std;

map<float, vector<vector<tuple<string, int> > > > dfsApproach(vector<string> keysCR, map<string, map<string, vector<vector<float> > > > &groupedCR, vector<string> keysRR,
                                                                    map<string, map<string, vector<vector<float> > > > &groupedRR, int maxDepth, int branchingFactor, int measureIndex) {

    map<float, vector<vector<tuple<string, int> > > > pathsMap;
    vector<vector<tuple<string, int> > > paths;

    float extensionSide = 0.0;

    for( auto key : keysCR) {
        paths = dfsWrapper(key, extensionSide, keysCR, groupedCR, keysRR, groupedRR, maxDepth, branchingFactor, measureIndex);
        for( auto path : paths) {
            pathsMap[extensionSide].push_back(path);
        }
    }

    // extensionSide = 1.0;
    // for( auto key : keysCR) {
    // paths = dfsWrapper(key, extensionSide, keysCR, groupedCR, keysRR, groupedRR, maxDepth, branchingFactor, measureIndex);
    //     for( auto path : paths) {
    //         pathsMap[extensionSide].push_back(path);
    //     }
    // }

    return pathsMap;
}


vector<vector<tuple<string, int> > > dfsWrapper(string start, float side, vector<string> keysCR, map<string, map<string, vector<vector<float> > > > &groupedCR,
                        vector<string> keysRR, map<string, map<string, vector<vector<float> > > > &groupedRR, int maxDepth, int branchingFactor, int measureIndex){

    vector<vector<tuple<string, int> > > paths;
    vector<tuple<string, int> > path;

    int flag = 0;
    if( groupedCR[start].size() == 0 ) {
        return paths;
    }
    for( int skip = 0; skip < groupedCR[start].size(); skip++){
        flag = 0;
        path = dfsSearch(start, side, keysCR, groupedCR, keysRR, groupedRR, maxDepth, measureIndex, skip, branchingFactor);
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

vector<tuple<string, int> > dfsSearch(string start, float side, vector<string> keysCR, map<string, map<string, vector<vector<float> > > > &groupedCR,
                        vector<string> keysRR, map<string, map<string, vector<vector<float> > > > &groupedRR, int maxDepth, int measureIndex, int skip, int branchingFactor) {

    vector<tuple<string, int> > path;
    vector<tuple<string, int> > alreadyLooked;
    vector<string> stack;
    string read;
    int number;
    vector<tuple<string, int> > skipping;

    path.push_back(make_tuple(start, -1));

    tie(read, number) = getReadForContig(start, side, groupedCR, skip);

    path.push_back(make_tuple(read, number));
    alreadyLooked.push_back(make_tuple(read, number));
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
                skipping.clear();
                while( skipping.size() < branchingFactor) {
                    tie(read, number) = getReadForRead(currentTarget, side, start, groupedCR, groupedRR, path, measureIndex, skipping);
                    if( number == -2 || number == -1) {
                        skipping.push_back(make_tuple(read, number));
                        break;
                    }
                    // check if we already visited that read
                    if ( std::find(alreadyLooked.begin(), alreadyLooked.end(), make_tuple(read ,number)) != alreadyLooked.end() ) {
                        continue;
                    } else {
                        skipping.push_back(make_tuple(read, number));
                    }
                }
                // if read doesnt have any
                if( skipping.size() == 0) {
                    path.pop_back();
                } else {
                    alreadyLooked.push_back(skipping.front());
                    path.push_back(skipping.front());
                    reverse(skipping.begin(), skipping.end());
                    for( auto val : skipping){
                        stack.push_back(get<0>(val));
                    }
                }
            }
        }
    }
    path = {make_tuple("", -2)};
    return path;
}

tuple<string, int> getReadForContig(string contig, float side, map<string, map<string, vector<vector<float> > > > &groupedCR, int skip) {

    tuple<string, int> maxTuple = make_tuple("", -1);
    for( auto query : groupedCR[contig]){
        for( int i = 0; i < query.second.size(); i++) {
            if( side == query.second[i][0]){
                if ( skip > 0) {
                    skip -= 1;
                    continue;
                } else {
                    return make_tuple(query.first, i);
                }
            }
        }
    }
    return maxTuple;
}

tuple<string, int> getReadForRead(string read, float side, string startContig, map<string, map<string, vector<vector<float> > > > &groupedCR, map<string, map<string, vector<vector<float> > > > &groupedRR, vector<tuple<string, int> > &path, int measureIndex, vector<tuple<string, int> > &skipping) {

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


    float maximum = 0.0;
    tuple<string, int> maxTuple = make_tuple("", -1);

    for( auto query : groupedRR[read]){
        for( int i = 0; i < query.second.size(); i++) {
            // pick measure
            if( side == query.second[i][0]){
                if ( std::find(skipping.begin(), skipping.end(), make_tuple(query.first, i)) != skipping.end() ) {
                    continue;
                }
                if( query.second[i][measureIndex] > maximum) {
                    maximum = query.second[i][measureIndex];
                    maxTuple = make_tuple(query.first, i);
                }
            }
        }
    }
    return maxTuple;
}
