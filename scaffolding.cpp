#include <vector>
#include <map>
#include <iostream>
#include <typeinfo>
#include <random>
#include <algorithm>

#include "scaffolding.h"


// returns one path, map[ctg1,ctg2] = [(read1,0),(read2,0)...]
// returns tuples for path, e.g. [(ctg1,ctg2), (ctg2,ctg3)]
// vector<tuple<string, string> > getScaffoldContigs(int ctgNumber, map<tuple<string, string>, vector<vector<tuple<string, int> > > > pathsMap) {
//
//     vector<tuple<string, string> > contigs;
//     vector<tuple<string, string> > tmp;
//     map<tuple<string, string>, int> pathCount;
//     cout << "Get scaffold contigs" << endl;
//     for( auto key : pathsMap) {
//         cout << "key1,2 in pathsmap: " << get<0>(key.first) << " " << get<1>(key.first) << endl;
//         // check if key is already in path count
//         if( std::find(tmp.begin(), tmp.end(), make_tuple(get<1>(key.first), get<0>(key.first))) != tmp.end()){
//             cout << "jabadabaduuuu" << endl;
//             pathCount[make_tuple(get<1>(key.first), get<0>(key.first))] += key.second.size();
//         } else {
//             tmp.push_back(key.first);
//             pathCount[key.first] = key.second.size();
//         }
//     }
//
//     if( tmp.size() < ctgNumber) {
//         // return all
//         for( auto key : pathCount) {
//             contigs.push_back(key.first);
//         }
//     } else {
//         // find lowest branch and dont include it
//         int lowest = 1000000;
//         tuple<string,string> lowestTuple;
//         for( auto key : pathCount) {
//             if( key.second < lowest) {
//                 lowest = key.second;
//                 lowestTuple = key.first;
//             }
//         }
//
//         for( auto key : pathCount) {
//             if( key.first != lowestTuple){
//                 contigs.push_back(key.first);
//             }
//         }
//     }
//
//     return contigs;
// }

vector<tuple<string, string> > getScaffoldContigs(int ctgNumber, map<tuple<string, string>, vector<vector<tuple<string, int> > > > pathsMap) {

    vector<tuple<string, string> > contigs;
    vector<tuple<string, string> > tmp;
    map<tuple<string, string>, int> pathCount;
    cout << "Get scaffold contigs" << endl;
    for( auto key : pathsMap) {
        cout << "key1,2 in pathsmap: " << get<0>(key.first) << " " << get<1>(key.first) << endl;
        // check if key is already in path count
        if( std::find(tmp.begin(), tmp.end(), make_tuple(get<1>(key.first), get<0>(key.first))) != tmp.end()){
            cout << "jabadabaduuuu" << endl;
            pathCount[make_tuple(get<1>(key.first), get<0>(key.first))] += key.second.size();
        } else {
            tmp.push_back(key.first);
            pathCount[key.first] = key.second.size();
        }
    }

    if( tmp.size() < ctgNumber) {
        // return all
        for( auto key : pathCount) {
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

        for( auto key : pathCount) {
            if( key.first != lowestTuple){
                contigs.push_back(key.first);
            }
        }
    }

    return contigs;
}


vector<vector<tuple<string, int> > > buildFinalScaffoldOrder(map<tuple<string, string>, vector<tuple<string, int> > > chosenPaths) {

    vector<vector<tuple<string, int> > > finalOrder;
    vector<tuple<string, int> > tmp;
    tuple<string, int> last;
    string begin;
    string end;
    tuple<string, int> target;
    tuple<string, int> query;

    for( auto key : chosenPaths) {
        begin = get<0>(key.first);
        end = get<1>(key.first);

        target = make_tuple(begin,-1);
        for( auto element : key.second) {
            tmp.clear();
            query = element;
            tmp.push_back(target);
            tmp.push_back(query);
            finalOrder.push_back(tmp);
            target = query;
        }
        tmp.clear();
        target = make_tuple(end,-1);;
        query = key.second.back();
        tmp.push_back(target);
        tmp.push_back(query);
        finalOrder.push_back(tmp);
    }

    for( auto vec : finalOrder) {
        for( auto element : vec) {
            cout << "Element:" << get<0>(element) << " " << get<1>(element) << endl;
        }
    }

    return finalOrder;
}
