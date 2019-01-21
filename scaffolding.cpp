#include <vector>
#include <map>
#include <iostream>
#include <typeinfo>
#include <random>
#include <algorithm>

#include "scaffolding.hpp"

// author: Karlo Brajdic

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

vector<tuple<string, string> > getScaffoldContigs(int ctgNumber, map<tuple<string, string>, vector<vector<tuple<string, int> > > > &pathsMapLeft, map<tuple<string, string>, vector<vector<tuple<string, int> > > > &pathsMapRight) {

    vector<tuple<string, string> > contigs;
    map<tuple<string, string>, int> pathNumber;
    int maximum = 0;
    tuple<string, string> maximumKey;

    // save number of paths for each pair of contigs
    for( auto key : pathsMapLeft){
        pathNumber[key.first] = key.second.size();
    }
    for( auto key : pathsMapRight) {
        tuple<string, string> key2 = make_tuple(get<1>(key.first), get<0>(key.first));
        if( pathNumber.find(key2) == pathNumber.end()) {
            continue;
        } else {
            pathNumber[key2] += key.second.size();
        }
    }

    while( ctgNumber > (contigs.size() + 1)) {
        // find the pair with highest number of paths between them
        for( auto key : pathNumber) {
            if( key.second > maximum) {
                maximum = key.second;
                maximumKey = key.first;
            }
        }
        if( maximum == 0){
            break;
        }
        contigs.push_back(maximumKey);
        pathNumber.erase(maximumKey);
        for (auto it = pathNumber.cbegin(); it != pathNumber.cend(); ) {
            if( get<0>(maximumKey) == get<0>((*it).first) || get<1>(maximumKey) == get<1>((*it).first)) {
                pathNumber.erase((*it).first);
                it++;
            } else {
                ++it;
            }
        }
        for (auto it = pathNumber.cbegin(); it != pathNumber.cend() ;) {
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
        // if there are no my contigs yet, make one
        if( contigsVec.size() == 0) {
            tmpVec.clear();
            tmpVec.push_back(key);
            contigsVec.push_back(tmpVec);
        } else {
            // first check if it belongs to some contig
            for( int i = 0; i < contigsVec.size(); i++) {
                if( checkVector0(contigsVec[i], get<1>(key))){
                    if( checkVector1(contigsVec[i], get<0>(key))) {
                        // if both are inside, dont append it, otherwise we will get circle
                        dontAppend = 1;
                    } else {
                        // push current key and look if second ctg in key already has contig
                        for( int j = 0; j < contigsVec.size(); j++) {
                            if( i == j) {
                                continue;
                            }
                            if( checkVector0(contigsVec[j], get<1>(key)) || checkVector1(contigsVec[j], get<0>(key))) {
                                for( auto key2 : contigsVec[j]) {
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
                    if( checkVector0(contigsVec[i], get<1>(key))) {
                        // if both are inside, dont append it, otherwise we will get circle
                        dontAppend = 1;
                    } else {
                        // push current key and look if second ctg in key already has contig
                        for( int j = 0; j < contigsVec.size(); j++) {
                            if( i == j) {
                                continue;
                            }
                            if( checkVector0(contigsVec[j], get<1>(key)) || checkVector1(contigsVec[j], get<0>(key))) {
                                for( auto key2 : contigsVec[j]) {
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
    }

    bool isFirst = true;
    vector<vector<tuple<string, string> > > contigsVec2;
    vector<tuple<string, string> > tmp2;
    for( auto vec : contigsVec) {
        for( auto key : vec) {
            isFirst = true;
            first = get<0>(key);
            for( auto key2 : vec) {
                if( key == key2) {
                    continue;
                }
                if( get<1>(key2) == first || get<0>(key2) == first) {
                    isFirst = false;
                    break;
                }
            }
            if( isFirst) {
                break;
            }
        }
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
    for( auto vec : contigsVec2) {
        for( auto key : vec) {
            contigs.push_back(key);
        }
    }
    return contigs;
}

// builds pairs and orders scaffold so we can easily access maps by these keys in buidling scaffold
vector<vector<tuple<string, int> > > buildFinalScaffoldOrder(map<tuple<string, string>, vector<tuple<string, int> > > chosenPaths, vector<tuple<string, string> > scaffoldContigs) {

    vector<vector<tuple<string, int> > > finalOrder;
    vector<tuple<string, int> > tmp;
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

    return finalOrder;
}
