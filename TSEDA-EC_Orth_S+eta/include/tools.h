//
// Created by qqq on 2021/9/15.
//

#ifndef FRAME_TOOLS_H
#include "common.h"
#define FRAME_TOOLS_H
void CalculateLevelList();
void CalculateDescendants();
void CalculateAncestors();
bool SortValueByDescend(pair<int,double>& a, pair<int,double>& b);
double CalculatePowerByLoad(double wl, int HTid);
bool SortByEnergyConsumption(chromosome& a, chromosome& b);
void IndexSortByValueOnAscend(vector<int>& ind, vector<double>& value);

#endif //FRAME_TOOLS_H
