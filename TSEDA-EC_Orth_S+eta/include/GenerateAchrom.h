//
// Created by qqq on 2021/9/15.
//

#ifndef FRAME_GENERATEACHROM_H
#include "common.h"
#define FRAME_GENERATEACHROM_H

void UpdateITL(set<double>& ITLofRscId, double& StartTime, double& EndTime);
void W_Cal_Average_S(vector<double>& w);
void C_Cal_Average(vector<vector<double>>& c);
void Calculate_Rank_b_S(vector<double>& RankList, vector<double>& ExeTime);
void IntChr(chromosome& chrom);
chromosome GnrChr_HEFT_S(vector<double> Rank_b);
chromosome GnrChr_HMEC_S(vector<double> Rank_b);
void GnrRscLstOfChr(chromosome& chrom, vector<vector<double>>& PMR);
double GnrML_Evl_EFT_S(chromosome& ch);
double GnrML_Evl_MEC_S(chromosome& ch);
double DcdEvl_S(chromosome& ch, bool IsFrw);
void SelectRsc_EFT_S(chromosome& ch, vector<set<double>>& ITL, int& TaskIndex, int& RscIndex, double& FinalStartTime, double& FinalEndTime);
double FindIdleTimeSlot(set<double>& ITLofRscId, double& ExeTime, double& ReadyTime);
double CalculateEnergy(chromosome& Chrom);
double CalculateECByDelta2(chromosome& Chrom, int& NumBer);
double IFBD_S(chromosome& ch);
void LBCA_S(chromosome& chrom);
void InitProModelOfResAlc(vector<vector<double>>& PMR);
void InitProModelOfTskSch(vector<vector<double>>& PMS, vector<int>& NumOfAncestors, vector<int>& NumOfNonDescendants);
void GnrRscLstOfChr(chromosome& chrom, vector<vector<double>>& PMR);
chromosome GnrTskLstOfChr(vector<vector<double> >& PMS, vector<double>& eta_TSO);
void UpdatePMR(vector<vector<double>>& PMR, chromosome& bstChrom);
void UpdatePMS(vector<vector<double>>& PMS, chromosome& bstChrom);

#endif //FRAME_GENERATEACHROM_H
