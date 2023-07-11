//
// Created by qqq on 2021/9/15.
//
#include "common.h"
#include "GenOperator.h"
#include "GenerateAchrom.h"
#include "tools.h"
using namespace std;

void CalculateLevelList() {
    LevelIdOfTask.resize(comConst.NumOfTsk);
    vector<int> InDegree;   //variables for recording the number of parent tasks whose level have not been calculated;用于记录尚未计算其级别的父任务数的变量
    vector<int> stk;        //a set for recording the index of tasks whose inDegree is equal to 0;
    InDegree.assign(comConst.NumOfTsk, 0);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        InDegree[i] = Tasks[i].parents.size();
    }
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        if (InDegree[i] == 0) stk.push_back(i);
    }
    int MaxLevel = 0;
    while (!stk.empty()) {
        int v = stk[0];
        LevelIdOfTask[v] = 0;
        for (int i = 0; i < Tasks[v].parents.size(); ++i) {
            if (LevelIdOfTask[Tasks[v].parents[i]] >= LevelIdOfTask[v]) {
                LevelIdOfTask[v] = LevelIdOfTask[Tasks[v].parents[i]] + 1;
            }
        }
        if(LevelIdOfTask[v] + 1> MaxLevel) {
            MaxLevel = LevelIdOfTask[v] + 1;
            TskLstInLvl.resize(MaxLevel);
        }
        TskLstInLvl[LevelIdOfTask[v]].push_back(v);
        stk.erase(stk.begin());
        for (int i = 0; i < Tasks[v].children.size(); ++i) {
            InDegree[Tasks[v].children[i]]--;
            if (InDegree[Tasks[v].children[i]] == 0) {
                stk.push_back(Tasks[v].children[i]);
            }
        }
    }
}

void CalculateDescendants() {
    Descendants.resize(comConst.NumOfTsk);
    for(int i = TskLstInLvl.size()-2; i >= 0; --i) {
        for(int taskId : TskLstInLvl[i]) {
            for(int childId : Tasks[taskId].children) {
                Descendants[taskId].insert(childId);
                Descendants[taskId].insert(Descendants[childId].begin(),Descendants[childId].end());
            }
        }
    }
}


void CalculateAncestors() {
    Ancestors.resize(comConst.NumOfTsk);
    for(int i = 1; i < TskLstInLvl.size(); ++i) {
        for(int taskId : TskLstInLvl[i]) {
            for(int parentId : Tasks[taskId].parents) {
                Ancestors[taskId].insert(parentId);
                Ancestors[taskId].insert(Ancestors[parentId].begin(),Ancestors[parentId].end());
            }
        }
    }
}

void IndexSortByValueOnAscend(vector<int>& ind, vector<double>& value) {
    for (int i = 0; i < ind.size(); ++i) {
        ind[i] = i;
    }
    sort(ind.begin(), ind.end(), [&value](int v1, int v2) { return value[v1] < value[v2]; });
}

void IndexSortByValueOnAscend(vector<int>& ind, vector<int>& value) {
    for (int i = 0; i < ind.size(); ++i) {
        ind[i] = i;
    }
    sort(ind.begin(), ind.end(), [&value](int v1, int v2) { return value[v1] < value[v2]; });
}

void IndexSortByValueOnDescend(vector<int>& ind, vector<double>& value) {
    vector<double> result;
    for (int i = 0; i < ind.size(); ++i) {
        ind[i] = i;
    }
    sort(ind.begin(), ind.end(), [&value](int v1, int v2) { return value[v1] > value[v2]; });
}


bool SortValueByDescend(pair<int,double>& a, pair<int,double>& b) {
    return a.second > b.second + PrecisionValue;
}

bool SortByEnergyConsumption(chromosome& a, chromosome& b) {
    return a.EnergyConsumption + PrecisionValue < b.EnergyConsumption;
}

//{generate a random number in [Start,End) }
double RandomDouble(int start, int end) {
    double ret = rand() % ((end - start) * 1000) / 1000.0 + start;
    return ret;
}

//{generate a random number in [Start,End] }
double RandomDouble2(int start, int end) {
    double ret = start + rand() % (end - start) +  rand() %  1001 / 1000.0;
    return ret;
}


double CalculatePowerByLoad(double ld, int HTid) {
    if(HTid == 0) {
        if (ld <= 0.1) {
            return ld*46.7 + 9.33;
        } else if (ld <= 0.2) {
            return ld*32 + 10.8;
        } else if (ld <= 0.3) {
            return ld*27 + 11.8;
        } else if (ld <= 0.4) {
            return ld*26 + 12.1;
        } else if (ld <= 0.5) {
            return ld*30 + 10.5;
        } else if (ld <= 0.6) {
            return ld*40 + 5.5;
        } else if (ld <= 0.7) {
            return ld*51 - 1.1;
        } else if (ld <= 0.8) {
            return ld*63 - 9.5;
        } else if (ld <= 0.9) {
            return ld*70 - 15.1;
        } else if (ld <= 1.0 + PrecisionValue) {
            return ld*39 + 12.8;
        } else {
            cout << "load can not is larger than 1" << endl; exit(0);
        }
    } else {
        if (ld <= 0.1) {
            return ld*55 + 10.0;
        } else if (ld <= 0.2) {
            return ld*29 + 12.6;
        } else if (ld <= 0.3) {
            return ld*34 + 11.6;
        } else if (ld <= 0.4) {
            return ld*37 + 10.7;
        } else if (ld <= 0.5) {
            return ld*45 + 7.5;
        } else if (ld <= 0.6) {
            return ld*55 + 2.5;
        } else if (ld <= 0.7) {
            return ld*64 - 2.9;
        } else if (ld <= 0.8) {
            return ld*81 - 14.8;
        } else if (ld <= 0.9) {
            return ld*106 - 34.8;
        } else if (ld <= 1.0 + PrecisionValue) {
            return ld*111 - 39.3;
        } else {
            cout << "load can not is larger than 1" << endl; exit(0);
        }
    }
}
