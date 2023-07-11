//
// Created by qqq on 2021/10/15.
//

#include "TSEDA.h"
#include "tools.h"
#include "config.h"
#include "GenerateAchrom.h"

double runTSEDA(string XmlFile, string RscAlcFile,double& SchTime, int& iteration) {
    clock_t start = clock();
    ReadFile(XmlFile, RscAlcFile);
    ConfigParameter_TSEDA();
    CalculateLevelList();
    CalculateDescendants();
    CalculateAncestors();
    vector<double> Rank_b(comConst.NumOfTsk, 0);
    vector<double> ww(comConst.NumOfTsk, 0);
//    vector<vector<double>> cc(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk, 0));
    W_Cal_Average_S(ww);
//    C_Cal_Average(cc);
    Calculate_Rank_b_S(Rank_b, ww);

    double MaxRank_b = Rank_b[0];
    for (int i = 1; i < comConst.NumOfTsk; ++i) {
        if (MaxRank_b + PrecisionValue < Rank_b[i]) {
            MaxRank_b = Rank_b[i];
        }
    }

    vector<int> NumOfAncestors(comConst.NumOfTsk);
    vector<int> NumOfDescendants(comConst.NumOfTsk);
    vector<int> NumOfNonDescendants(comConst.NumOfTsk);
//    #pragma omp parallel for
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        NumOfAncestors[i] = Ancestors[i].size();
        NumOfDescendants[i] = Descendants[i].size();
        NumOfNonDescendants[i] = comConst.NumOfTsk - NumOfDescendants[i];
    }

    vector<vector<double>> PMR(comConst.NumOfTsk, vector<double>(comConst.NumOfRsc, 0));
    InitProModelOfResAlc(PMR);                      //PMR[i][j] represents the probability that task i is assigned to resource j;
    vector<vector<double>> PMS(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk, 0));
    InitProModelOfTskSch(PMS, NumOfAncestors, NumOfNonDescendants); //PMS[i][k] represents the probability that the k-th scheduled task is task i

    vector<chromosome> Population(Parameter_TSEDA.NumOfChromPerPop);
    chromosome BstChrom = GnrChr_HMEC_S(Rank_b);
    chromosome Chrom_HEFT = GnrChr_HEFT_S(Rank_b);
    if (Chrom_HEFT.EnergyConsumption + PrecisionValue < BstChrom.EnergyConsumption)
        BstChrom = Chrom_HEFT;

    vector<double> eta_TSO(comConst.NumOfTsk);
    double RunTime = (double) (clock() - start) / CLOCKS_PER_SEC;

    while (RunTime + PrecisionValue < Parameter_TSEDA.RunTimeRatioOfStg1 * SchTime) {
        #pragma omp parallel for
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            eta_TSO[i] = pow(Rank_b[i]/MaxRank_b, (1-RunTime/SchTime) * Parameter_TSEDA.fdhi);
        }

        #pragma omp parallel for
        for (int n = 0; n < Parameter_TSEDA.NumOfChromPerPop; ++n) {
            Population[n] = GnrTskLstOfChr(PMS, eta_TSO);
            GnrML_Evl_MEC_S(Population[n]);
        }

        for (int n = 0; n < Parameter_TSEDA.NumOfChromPerPop; ++n) {
            if (Population[n].EnergyConsumption + PrecisionValue < BstChrom.EnergyConsumption)
                BstChrom = Population[n];
        }

        UpdatePMR(PMR, BstChrom);  UpdatePMS(PMS, BstChrom);
        ++iteration;
        RunTime = (double) (clock() - start) / CLOCKS_PER_SEC;
    }

    while (RunTime + PrecisionValue < SchTime){
        #pragma omp parallel for
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            eta_TSO[i] = pow(Rank_b[i]/MaxRank_b, (1-RunTime/SchTime) * Parameter_TSEDA.fdhi);
        }

        #pragma omp parallel for
        for (int n = 0; n < Parameter_TSEDA.NumOfChromPerPop; ++n) {
            Population[n] = GnrTskLstOfChr(PMS, eta_TSO);
            GnrRscLstOfChr(Population[n], PMR);
            DcdEvl_S(Population[n],true);
            CalculateEnergy(Population[n]);
        }

        sort(Population.begin(), Population.end(), SortByEnergyConsumption);
        #pragma omp parallel for
        for (int n = 0; n < Parameter_TSEDA.NumOfImproveOfPop; ++n) {
            IFBD_S(Population[n]); //IFBSS;
            LBCA_S(Population[n]); //LBCAS;
        }

        for (int n = 0; n < Parameter_TSEDA.NumOfImproveOfPop; ++n) {
            if (Population[n].EnergyConsumption + PrecisionValue < BstChrom.EnergyConsumption)
                BstChrom = Population[n];
        }

        UpdatePMR(PMR, BstChrom);  UpdatePMS(PMS, BstChrom);
        ++iteration;
        RunTime = (double) (clock() - start) / CLOCKS_PER_SEC;
    }
    SchTime = RunTime;
    return BstChrom.EnergyConsumption;
}
