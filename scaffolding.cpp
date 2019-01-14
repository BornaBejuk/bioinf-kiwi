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

// Determines which paths should be included in scaffold
vector<tuple<string, string> > getScaffoldContigs(int ctgNumber, map<tuple<string, string>, vector<vector<tuple<string, int> > > > pathsMap) {


    vector<tuple<string, string> > contigs;
    vector<tuple<string, string> > tmp;
    map<tuple<string, string>, int> pathCount;
    cout << "Get scaffold contigs" << endl;
    string name;
    int maxNumber = 0;
    int minNumber = 100;
    int number;
    int number2;
    tuple<string, string> currentKey;
    for( auto key : pathsMap) {
        name = get<0>(key.first).substr(0,get<0>(key.first).size()-1);
        // cout << "contig name: " << name << endl;

        // cout << "key1,2 in pathsmap: " << get<0>(key.first) << " " << get<1>(key.first) << endl;

        number = stoi(get<0>(key.first).substr(get<0>(key.first).size()-1, get<0>(key.first).size()));
        number2 = stoi(get<1>(key.first).substr(get<1>(key.first).size()-1, get<0>(key.first).size()));

        if( number > maxNumber ) {
            maxNumber = number;
        }
        if( number2 > maxNumber ) {
            maxNumber = number2;
        }
        if( number < minNumber ) {
            minNumber = number;
        }
        if( number2 < minNumber ) {
            minNumber = number2;
        }
        // cout << "contig number: " << number << endl;
        // break;
        tmp.push_back(key.first);
        // if( stoi(get<0>(key.first)) )
    }

    // cout << "Contig name" << name << " min,max number: " << minNumber << " " << maxNumber << endl;

    string lastCtg = name + to_string(maxNumber);
    string currentCtg;
    for( int i = maxNumber - 1; i >= minNumber; i--) {
        currentCtg = name + to_string(i);
        currentKey = make_tuple(lastCtg,currentCtg);
        // cout << "Current key: " << get<0>(currentKey) << " " << get<1>(currentKey) << endl;
        if( std::find(tmp.begin(), tmp.end(), currentKey) != tmp.end()){
            contigs.push_back(currentKey);
            cout << "current key appended" << endl;
        }
        lastCtg = currentCtg;
    }

    vector<string> queries;
    if( contigs.size() < ( maxNumber - minNumber + 1)) {
        for( auto key : pathsMap) {
            // check if key is already in path count
            if( std::find(contigs.begin(), contigs.end(), make_tuple(get<1>(key.first), get<0>(key.first))) != contigs.end() || std::find(queries.begin(), queries.end(), get<1>(key.first)) != queries.end()){
                cout << "Same key or end so dont append it" << endl;
                // pathCount[make_tuple(get<1>(key.first), get<0>(key.first))] += key.second.size();
            } else {
                queries.push_back(get<1>(key.first));
                contigs.push_back(key.first);
            }
        }
    }
    return contigs;
}

// builds pairs and orders scaffold so we can easily access maps by these keys in buidling scaffold
vector<vector<tuple<string, int> > > buildFinalScaffoldOrder(map<tuple<string, string>, vector<tuple<string, int> > > chosenPaths, vector<tuple<string, string> > scaffoldContigs) {

    vector<vector<tuple<string, int> > > finalOrder;
    vector<tuple<string, int> > tmp;
    tuple<string, int> last;
    string begin;
    string end;
    tuple<string, int> target;
    tuple<string, int> query;

    for( auto key : scaffoldContigs) {
        begin = get<0>(key);
        end = get<1>(key);

        target = make_tuple(begin,-1);
        for( auto element : chosenPaths[key]) {
            tmp.clear();
            query = element;
            tmp.push_back(target);
            tmp.push_back(query);
            finalOrder.push_back(tmp);
            target = query;
        }
        tmp.clear();
        target = make_tuple(end,-1);;
        query = chosenPaths[key].back();
        tmp.push_back(target);
        tmp.push_back(query);
        finalOrder.push_back(tmp);
    }

    // for( auto vec : finalOrder) {
    //     for( auto element : vec) {
    //         cout << "Element:" << get<0>(element) << " " << get<1>(element) << endl;
    //     }
    // }

    return finalOrder;
}
