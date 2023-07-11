//
// Created by yx on 22-5-7.
//

#include "ADBRKGA.h"
#include <common.h>
#include "tools.h"
#include "config.h"
#include "GenerateAchrom.h"
#include "GenOperator.h"

double runADBRKGA(string XmlFile, string RscAlcFile, double& SchTime, int& iteration){
    double CurNum = 0; //double RunTime = 0;
    clock_t start = clock();
    ReadFile(XmlFile, RscAlcFile);
    ConfigParameter_ADBRKGA();

    CalculateLevelList();
    vector<double> ww(comConst.NumOfTsk,0.0);
//    vector<vector<double>> cc(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk,0.0));
    W_Cal_Average_S(ww);
//    C_Cal_Average(cc);
    vector<double> Rank_b(comConst.NumOfTsk,0.0);
    Calculate_Rank_b_S(Rank_b, ww);
    vector<double> Rank_t(comConst.NumOfTsk,0.0);
    Calculate_Rank_t_S(Rank_t, ww);

    set<chromosome> TemSubPop;
    chromosome chrom_DIHEFT = GnrChr_DIHEFT_S(Rank_b);
    chromosome chrom_HEFT_b = GnrChr_HEFT_b_ADBRKGA_S(Rank_b);
    chromosome chrom_HEFT_t = GnrChr_HEFT_t_ADBRKGA_S(Rank_t);
    chromosome chrom_IHEFT3_b = GnrChr_IHEFT3_b_S(Rank_b);
    chromosome chrom_IHEFT3_t = GnrChr_IHEFT3_t_S(Rank_t);
//    chromosome chrom_HEMC = GnrChr_HMEC(Rank_b);
//    GnrCode_RK(chrom_HEMC);
    TemSubPop.insert(chrom_DIHEFT);
    TemSubPop.insert(chrom_HEFT_b);
    TemSubPop.insert(chrom_HEFT_t);
    TemSubPop.insert(chrom_IHEFT3_b);
    TemSubPop.insert(chrom_IHEFT3_t);
//    TemSubPop.insert(chrom_HEMC);

    int terminationNum = ceil(CDT * sqrt(ModelScale)/Parameter_ADBRKGA.NumOfChromPerPop); double TotalNum = terminationNum;
    int NumOfNoImpGen = 0, MaxNoImpGen = 0;

    while (TemSubPop.size() < Parameter_ADBRKGA.NumOfChromPerPop){
        chromosome chrom_Ran = GnrChr_Lvl_Ran();
        AdpDcd_S(chrom_Ran,CurNum,TotalNum);
        CalculateEnergy(chrom_Ran);
        TemSubPop.insert(chrom_Ran);
    }

    vector<chromosome> population;
    population.assign(TemSubPop.begin(),TemSubPop.end());
    double bestFitness = population[0].EnergyConsumption;

    vector<double> A(Parameter_ADBRKGA.NumOfChromPerPop);
    CalSlctProb_Rank(1+1.0/Parameter_ADBRKGA.NumOfChromPerPop, A ,Parameter_ADBRKGA.NumOfChromPerPop); //calculate the cumulative probabilities
    int ImprovementNumber = ceil(Parameter_ADBRKGA.ImprovementRate * Parameter_ADBRKGA.NumOfChromPerPop);
    int ImmigrationNumber = ceil(Parameter_ADBRKGA.ImmigrationRate * Parameter_ADBRKGA.NumOfChromPerPop);

    while (1) {
        ++iteration;
        CurNum = MaxNoImpGen + 1;
        vector<chromosome> NewPopulation(Parameter_ADBRKGA.NumOfChromPerPop);
        //#pragma omp parallel for
        for (int n = 0; n < Parameter_ADBRKGA.NumOfChromPerPop; ++n) {
            int ind1 = SelectChrom(A);
            int ind2 = SelectChrom(A);
            while (ind1 == ind2) {
                ind2 = SelectChrom(A);
            }
            chromosome chrom1 ;
            chromosome chrom2 ;
            if (population[ind1].MakeSpan + PrecisionValue < population[ind2].MakeSpan) {
                chrom1 = population[ind1];
                chrom2 = population[ind2];
            } else {
                chrom1 = population[ind2];
                chrom2 = population[ind1];
            }
            NewPopulation[n] = Crs_BPUC(chrom1,chrom2);
        }

        #pragma omp parallel for
        for ( int n = 0; n < Parameter_ADBRKGA.NumOfChromPerPop; ++n ) {
            AdpDcd_S(NewPopulation[n],CurNum,TotalNum);
            CalculateEnergy(NewPopulation[n]);
        }

        sort(NewPopulation.begin(), NewPopulation.end(), SortByEnergyConsumption);
        #pragma omp parallel for
        for (int n = 0; n < ImprovementNumber; ++n) {
            LBCA_IFBS(NewPopulation[n]);
//            CalculateEnergy(NewPopulation[n]);
        }

        set<chromosome> NxtPop;
        NxtPop.insert(population.begin(),population.end());
        NxtPop.insert(NewPopulation.begin(),NewPopulation.end());
        set<chromosome>::iterator iter = NxtPop.begin();
        advance(iter,(Parameter_ADBRKGA.NumOfChromPerPop - ImmigrationNumber));
        NxtPop.erase(iter,NxtPop.end());
        while (NxtPop.size() < Parameter_ADBRKGA.NumOfChromPerPop) {
            chromosome chrom_Ran = GnrChr_Lvl_Ran();
            AdpDcd_S(chrom_Ran,CurNum,TotalNum);
            CalculateEnergy(chrom_Ran);
            NxtPop.insert(chrom_Ran);
        }
        population.assign(NxtPop.begin(),NxtPop.end());
        ++NumOfNoImpGen;
        if( population[0].EnergyConsumption + PrecisionValue < bestFitness ){
            bestFitness = population[0].EnergyConsumption;
            NumOfNoImpGen = 0;
        }
        MaxNoImpGen = XY_MAX(MaxNoImpGen,NumOfNoImpGen);
        if (NumOfNoImpGen == terminationNum) {
            SchTime = (double)(clock()-start)/CLOCKS_PER_SEC;
            break;
        }
    }
    return population[0].EnergyConsumption;
}
