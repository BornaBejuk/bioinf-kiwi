#include <vector>
#include <map>
#include <iostream>
#include <typeinfo>
#include <random>
#include <algorithm>

#include "scaffolding.hpp"

// author: Karlo Brajdic

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

bool checkVector0(vector<tuple<string, string> > vec, string value) {
    for( auto key : vec) {
        if( value == get<0>(key)) {
            return true;
        }
    }
}

bool checkVector1(vector<tuple<string, string> > vec, string value) {
    for( auto key : vec) {
        if( value == get<1>(key)) {
            return true;
        }
    }
}

vector<tuple<string, string> > getScaffoldContigs2(int ctgNumber, map<tuple<string, string>, vector<vector<tuple<string, int> > > > &pathsMapLeft, map<tuple<string, string>, vector<vector<tuple<string, int> > > > &pathsMapRight) {

    vector<tuple<string, string> > contigs;
    map<tuple<string, string>, int> pathNumber;
    int maximum = 0;
    tuple<string, string> maximumKey;

    // save number of paths for each pair of contigs
    for( auto key : pathsMapLeft){
        // cout << "Counting " << get<0>(key.first) << " " << get<1>(key.first) << " paths:" << key.second.size() << endl;
        pathNumber[key.first] = key.second.size();
    }

    while( ctgNumber > (contigs.size() + 1)) {
        // find the pair with highest number of paths between them
        for( auto key : pathNumber) {
            cout << "Path considered " << get<0>(key.first) << " " << get<1>(key.first) << endl;
            if( key.second > maximum) {
                maximum = key.second;
                maximumKey = key.first;
            }
        }
        if( maximum == 0){
            break;
        }
        cout << "Max key " << get<0>(maximumKey) << " " << get<1>(maximumKey) << endl;
        cout << "Maximum " << maximum << endl;
        contigs.push_back(maximumKey);
        // exit(1);
        pathNumber.erase(maximumKey);
        for (auto it = pathNumber.cbegin(); it != pathNumber.cend() /* not hoisted */; /* no increment */) {

            // cout << "spinning: " << get<0>((*it).first) << " " << get<1>((*it).first) << " " << (*it).second << endl;
            if( get<0>(maximumKey) == get<0>((*it).first) || get<1>(maximumKey) == get<1>((*it).first)) {
                pathNumber.erase((*it).first);
                it++;
            } else {
                ++it;
            }
        }
        for (auto it = pathNumber.cbegin(); it != pathNumber.cend() /* not hoisted */; /* no increment */) {

            // cout << "spinning: " << get<0>((*it).first) << " " << get<1>((*it).first) << " " << (*it).second << endl;
            if( get<0>(maximumKey) == get<0>((*it).first) || get<1>(maximumKey) == get<1>((*it).first)) {
                pathNumber.erase((*it).first);
                it++;
            } else {
                ++it;
            }
        }
        maximum = 0;
    }

    int dontAppend = 0;
    int forDelete = -1;
    string first;
    string last;
    vector<tuple<string, string> > tmpVec;
    vector<tuple<string, string> > contigsFiltered;
    vector<vector<tuple<string, string> > > contigsVec;
    for( auto key : contigs) {
        cout << "----------------------------------" << endl;
        cout << "Key considered " << get<0>(key) << " " << get<1>(key) << endl;
        // if there are no my contigs yet, make one
        if( contigsVec.size() == 0) {
            tmpVec.clear();
            tmpVec.push_back(key);
            contigsVec.push_back(tmpVec);
            cout << "yup1 " << get<0>(key) << " " << get<1>(key) << endl;
        } else {
            // first check if it belongs to some contig
            for( int i = 0; i < contigsVec.size(); i++) {
                if( checkVector0(contigsVec[i], get<1>(key))){
                    cout << "number 1 if" << endl;
                    if( checkVector1(contigsVec[i], get<0>(key))) {
                        // if both are inside, dont append it, otherwise we will get circle
                        cout << "already inside " << "yup " << get<0>(key) << " " << get<1>(key) << endl;
                        dontAppend = 1;
                        // break;
                    } else {
                        cout << "this else" << endl;
                        // push current key and look if second ctg in key already has contig
                        for( int j = 0; j < contigsVec.size(); j++) {
                            if( i == j) {
                                continue;
                            }
                            if( checkVector0(contigsVec[j], get<1>(key)) || checkVector1(contigsVec[j], get<0>(key))) {
                                cout << "this if" << endl;
                                for( auto key2 : contigsVec[j]) {
                                    cout << "found it in another so copy " << get<0>(key) << " " << get<1>(key) << endl;
                                    cout << "found it in another so copy2 " << get<0>(key2) << " " << get<1>(key2) << endl;
                                    contigsVec[i].push_back(key2);
                                }
                                forDelete = j;
                                dontAppend = 1;
                                break;
                            }
                        }
                        contigsVec[i].push_back(key);
                        dontAppend = 1;
                    }
                } else if( checkVector1(contigsVec[i], get<0>(key))){
                    cout << "2number 1 if" << endl;
                    if( checkVector0(contigsVec[i], get<1>(key))) {
                        // if both are inside, dont append it, otherwise we will get circle
                        cout << "2already inside " << "yup " << get<0>(key) << " " << get<1>(key) << endl;
                        dontAppend = 1;
                        // break;
                    } else {
                        cout << "2this else" << endl;
                        // push current key and look if second ctg in key already has contig
                        for( int j = 0; j < contigsVec.size(); j++) {
                            if( i == j) {
                                continue;
                            }
                            if( checkVector0(contigsVec[j], get<1>(key)) || checkVector1(contigsVec[j], get<0>(key))) {
                                cout << "2this if" << endl;
                                for( auto key2 : contigsVec[j]) {
                                    cout << "2found it in another so copy " << get<0>(key) << " " << get<1>(key) << endl;
                                    cout << "2found it in another so copy2 " << get<0>(key2) << " " << get<1>(key2) << endl;
                                    contigsVec[i].push_back(key2);
                                }
                                forDelete = j;
                                dontAppend = 1;
                                break;
                            }
                        }
                        contigsVec[i].push_back(key);
                        dontAppend = 1;
                    }
                }
                if( dontAppend == 1) {
                    break;
                }
            }
            if( dontAppend == 1){
                dontAppend = 0;
            } else {
                // create new contig
                tmpVec.clear();
                tmpVec.push_back(key);
                contigsVec.push_back(tmpVec);
            }
            if( forDelete > -1) {
                contigsVec.erase(contigsVec.begin() + forDelete);
                forDelete = -1;
            }
        }

        for( auto vec : contigsVec) {
            cout << "contig: " << endl;
            for( auto key : vec) {
                cout << "yup " << get<0>(key) << " " << get<1>(key) << endl;
            }
        }

    }

    bool isFirst = true;
    vector<vector<tuple<string, string> > > contigsVec2;
    vector<tuple<string, string> > tmp2;
    cout << "EL finalo" << endl;
    for( auto vec : contigsVec) {
        cout << "contig: " << endl;
        for( auto key : vec) {
            first = get<0>(key);
            // last = get<1>(key);
            for( auto key2 : vec) {
                if( key == key2) {
                    continue;
                }
                if( get<1>(key2) == first) {
                    isFirst = false;
                    break;
                }
            }
            if( isFirst) {
                break;
            }
            cout << "yup " << get<0>(key) << " " << get<1>(key) << endl;
        }
        cout << "FIRST " << first << endl;
        while( tmp2.size() < vec.size() ) {
            for(auto key : vec) {
                if( std::find(tmp2.begin(), tmp2.end(), key) != tmp2.end()) {
                    continue;
                } else {
                    if( first == get<0>(key)) {
                        tmp2.push_back(key);
                        first = get<1>(key);
                    }
                }
            }
        }
        contigsVec2.push_back(tmp2);
        tmp2.clear();
    }
    contigs.clear();
    cout << "EL finalo 2" << endl;
    for( auto vec : contigsVec2) {
        cout << "contig: " << endl;
        for( auto key : vec) {
            cout << "yup " << get<0>(key) << " " << get<1>(key) << endl;
            contigs.push_back(key);
        }
    }
    return contigs;
}

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
    cout << "Contigs " << contigs.size() << endl;
    vector<string> queries;
    if( contigs.size() < ( maxNumber - minNumber)) {
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
