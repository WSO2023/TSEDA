//
// Created by qqq on 2021/9/15.
//

#include "GenerateAchrom.h"
#include <cstdlib>
#include "tools.h"
#include "common.h"

//{calculate the average execution time of tasks}
void W_Cal_Average_S(vector<double>& w) {
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        double sum = 0, AllTransferData = Tasks[i].IFileSizeSum + Tasks[i].OFileSizeSum;
        for (int RscId : Tasks[i].ElgRsc)
            sum += Tasks[i].length / Rscs[RscId].pc + AllTransferData / VALUE * 8 / Rscs[RscId].bw;
        w[i] = sum / Tasks[i].ElgRsc.size();
    }
}
//{calculate the average transfer time among tasks}
void C_Cal_Average(vector<vector<double>>& c) {
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        if(Tasks[i].parents.size() == 0){
            continue;
        }
        for (int j = 0; j < Tasks[i].parents.size(); ++j) {
            int parent = Tasks[i].parents[j];
            double sum1 = 0;
            double sum = 0;
            sum = ParChildTranFileSizeSum[parent][i] / VALUE;
            for (int k = 0; k < Tasks[i].ElgRsc.size(); ++k) {
                for (int y = 0; y < Tasks[parent].ElgRsc.size(); ++y) {
                    if (Tasks[i].ElgRsc[k] == Tasks[parent].ElgRsc[y]) {
                        continue;
                    } else {
                        sum1 += sum * 8 / XY_MIN(Rscs[Tasks[i].ElgRsc[k]].bw, Rscs[Tasks[parent].ElgRsc[y]].bw);
                    }
                }
            }
            c[parent][i] = sum1 / (double) (Tasks[i].ElgRsc.size() * Tasks[parent].ElgRsc.size());
        }
    }
}

//calculate the rank of tasks based on shared database
void Calculate_Rank_b_S(vector<double>& RankList, vector<double>& ExeTime){
    for (int i = 0; i < TskLstInLvl[TskLstInLvl.size()-1].size(); ++i) {
        int TaskId=TskLstInLvl[TskLstInLvl.size()-1][i];
        RankList[TaskId] = ExeTime[TaskId];
    }
    for(int i =TskLstInLvl.size()-2 ;i >=0 ;--i){
        for (int TaskId: TskLstInLvl[i]) {
            for (int Child: Tasks[TaskId].children) {
                if(RankList[TaskId] + PrecisionValue < RankList[Child] ){
                    RankList[TaskId] = RankList[Child];
                }
            }
            RankList[TaskId] = RankList[TaskId] + ExeTime[TaskId];
        }
    }
}

chromosome GnrChr_HEFT_S(vector<double> Rank_b) {
    vector<int> ind(comConst.NumOfTsk);
    chromosome TemChrom;
    IntChr(TemChrom);
    IndexSortByValueOnAscend(ind, Rank_b);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        TemChrom.TskSchLst[i] = ind[comConst.NumOfTsk - i - 1];
    }
    GnrML_Evl_EFT_S(TemChrom);
    CalculateEnergy(TemChrom);
    return TemChrom;
}

chromosome GnrChr_HMEC_S(vector<double> Rank_b) {
    vector<int> ind(comConst.NumOfTsk);
    chromosome TemChrom;
    IntChr(TemChrom);
    IndexSortByValueOnAscend(ind, Rank_b);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        TemChrom.TskSchLst[i] = ind[comConst.NumOfTsk - i - 1];
    }
    GnrML_Evl_MEC_S(TemChrom);
    return TemChrom;
}

double GnrML_Evl_MEC_S(chromosome& ch){
    for (int i = 0; i < comConst.NumOfTsk; ++i)
        ch.RscAlcLst[i] = -1;
    vector<set<double> > ITL;                           //the idle time-slot lists  for all resources
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0); a.insert(9999999999 * 1.0);
        ITL.push_back(a);
    }
    ch.MakeSpan = 0;
    ch.EnergyConsumption = 0;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscIndex = -1,TaskIndex = ch.TskSchLst[i];
        double FinalEC = 9999999999, FinalEndTime = 9999999999, FinalStartTime = 0;
        for (int v: Tasks[TaskIndex].ElgRsc) {
            double ReadyTime = 0, TransferDataSize = Tasks[TaskIndex].ExternalInputFileSizeSum;
            for (int ParentIndex : Tasks[TaskIndex].parents) {
                int ParentRscIndex = ch.RscAlcLst[ParentIndex];
                if(v != ParentRscIndex){
                    TransferDataSize += ParChildTranFileSizeSum[ParentIndex][TaskIndex];
                }
                if (ReadyTime + PrecisionValue < ch.EndTime[ParentIndex]){
                    ReadyTime = ch.EndTime[ParentIndex];
                }
            }
            double ExeTime = Tasks[TaskIndex].length  / Rscs[v].pc + (TransferDataSize + Tasks[TaskIndex].OFileSizeSum) / VALUE * 8 / Rscs[v].bw;
            ch.StartTime[TaskIndex] = FindIdleTimeSlot(ITL[v], ExeTime, ReadyTime);
            ch.EndTime[TaskIndex] = ch.StartTime[TaskIndex] + ExeTime;
            ch.RscAlcLst[TaskIndex] = v;
            double CurEC = CalculateECByDelta2(ch,i);
            //{find/record the min EC}
            if (CurEC + PrecisionValue < FinalEC) {
                FinalStartTime = ch.StartTime[TaskIndex];
                FinalEndTime = ch.EndTime[TaskIndex];
                FinalEC = CurEC;
                RscIndex = v;
            }
        }
        ch.EnergyConsumption = FinalEC;
        ch.StartTime[TaskIndex] = FinalStartTime;
        ch.EndTime[TaskIndex] = FinalEndTime;
        ch.RscAlcLst[TaskIndex] = RscIndex;
        ch.MakeSpan = XY_MAX(ch.MakeSpan, FinalEndTime);
        UpdateITL(ITL[RscIndex],FinalStartTime,FinalEndTime);
    }
    return ch.EnergyConsumption;
}

// I/O independent
double DcdEvl_S(chromosome& ch, bool IsFrw) {
    ch.MakeSpan = 0;
    vector<set<double> > ITL;                   //record the idle time-slot of all resources
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0);  a.insert(9999999999 * 1.0);
        ITL.push_back(a);
    }
    //startDecode
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int TaskIndex = ch.TskSchLst[i], RscIndex = ch.RscAlcLst[TaskIndex];  //obtain the task and resource (Rsc) allocated to this task
        double ReadyTime = 0;
        double TransferSize = Tasks[TaskIndex].ExternalInputFileSizeSum; //the size of input files that need to obtain from share database;
        for (int ParentTask: Tasks[TaskIndex].parents) {
            if(RscIndex != ch.RscAlcLst[ParentTask]) {
                TransferSize = TransferSize + ParChildTranFileSizeSum[ParentTask][TaskIndex];
            }
        }
        if(IsFrw) {                              //forward-loading
            for (int ParentTask: Tasks[TaskIndex].parents) {
                if (ReadyTime + PrecisionValue < ch.EndTime[ParentTask]) {
                    ReadyTime = ch.EndTime[ParentTask];
                }
            }
        } else {                                //backward-loading
            for (int ChildTask: Tasks[TaskIndex].children) {
                if (ReadyTime + PrecisionValue < ch.EndTime[ChildTask]) {
                    ReadyTime = ch.EndTime[ChildTask];
                }
            }
        }
        double ExecutionTime = (TransferSize + Tasks[TaskIndex].OFileSizeSum) / VALUE * 8 / Rscs[RscIndex].bw + Tasks[TaskIndex].length / Rscs[RscIndex].pc;
        ch.StartTime[TaskIndex] = FindIdleTimeSlot(ITL[RscIndex],ExecutionTime,ReadyTime);
        ch.EndTime[TaskIndex] = ch.StartTime[TaskIndex] + ExecutionTime;
        if (ch.MakeSpan + PrecisionValue < ch.EndTime[TaskIndex]) {
            ch.MakeSpan = ch.EndTime[TaskIndex];
        }
        //{update ITL}
        UpdateITL(ITL[RscIndex],ch.StartTime[TaskIndex],ch.EndTime[TaskIndex]);
    }
    return ch.MakeSpan;
}

//{initialize chromosome to allocate spaces}
void IntChr(chromosome& chrom) {
    chrom.TskSchLst.resize(comConst.NumOfTsk, -1);
    chrom.Code_RK.resize(comConst.NumOfTsk);
    chrom.RscAlcLst.resize(comConst.NumOfTsk, -1);
    chrom.EndTime.resize(comConst.NumOfTsk);
    chrom.StartTime.resize(comConst.NumOfTsk);
    chrom.Code_TD.resize(comConst.NumOfRsc);
    chrom.TskSchPart.resize(comConst.NumOfTsk);
    chrom.RscAlcPart.resize(comConst.NumOfTsk);
    chrom.VTskSchPart.resize(comConst.NumOfTsk,0);
    chrom.VRscAlcPart.resize(comConst.NumOfTsk,0);
}

void UpdateITL(set<double>& ITLofRscId,double& StartTime,double& EndTime){
    if(ITLofRscId.find(StartTime) != ITLofRscId.end()) {
        ITLofRscId.erase(StartTime);
    } else {
        ITLofRscId.insert(StartTime);
    }
    if(ITLofRscId.find(EndTime) != ITLofRscId.end()) {
        ITLofRscId.erase(EndTime);
    } else {
        ITLofRscId.insert(EndTime);
    }
}

double GnrML_Evl_EFT_S(chromosome& ch) {
    for (int i = 0; i < comConst.NumOfTsk; ++i)
        ch.RscAlcLst[i] = -1;
    vector<set<double> > ITL;                           //the idle time-slot lists  for all resources
    double makespan = 0;
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0); a.insert(9999999999 * 1.0);
        ITL.push_back(a);
    }
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscId = -1, TaskId = ch.TskSchLst[i];
        double FinalEndTime = 9999999999, FinalStartTime = 0;
        SelectRsc_EFT_S(ch,ITL,TaskId,RscId,FinalStartTime,FinalEndTime);  //Find the resource that can finish the task earliest
        ch.EndTime[TaskId] = FinalEndTime;
        ch.StartTime[TaskId] = FinalStartTime;
        ch.RscAlcLst[TaskId] = RscId;
        UpdateITL(ITL[RscId],FinalStartTime,FinalEndTime);              //update ITL
        makespan = XY_MAX(makespan, FinalEndTime);
    }
    ch.MakeSpan = makespan;
    return makespan;
}

void SelectRsc_EFT_S(chromosome& ch, vector<set<double>>& ITL, int& TaskIndex, int& RscIndex, double& FinalStartTime, double& FinalEndTime) {
    for (int RscIdOfCrnTsk : Tasks[TaskIndex].ElgRsc) {
        double ReadyTime = 0;
        double TransferData = Tasks[TaskIndex].ExternalInputFileSizeSum;
        for (int ParentIndex: Tasks[TaskIndex].parents) { //calculate the ready time and transfer data of the task
            int RscIdOfPrnTsk = ch.RscAlcLst[ParentIndex];
            if(RscIdOfCrnTsk != RscIdOfPrnTsk){
                TransferData = TransferData + ParChildTranFileSizeSum[ParentIndex][TaskIndex];
            }
            if (ReadyTime + PrecisionValue < ch.EndTime[ParentIndex]){
                ReadyTime = ch.EndTime[ParentIndex];
            }
        }
        double ExeTime = Tasks[TaskIndex].length / Rscs[RscIdOfCrnTsk].pc + (TransferData + Tasks[TaskIndex].OFileSizeSum) / VALUE * 8 / Rscs[RscIdOfCrnTsk].bw;
        double StartTime = FindIdleTimeSlot(ITL[RscIdOfCrnTsk],ExeTime,ReadyTime); //Find an idle time-slot as early as possible from ITL
        double EndTime = StartTime + ExeTime;
        //{find/record the earliest finish time}
        if (EndTime + PrecisionValue < FinalEndTime) {
            FinalStartTime = StartTime;
            FinalEndTime = EndTime;
            RscIndex = RscIdOfCrnTsk;
        }
    }
}

double FindIdleTimeSlot(set<double>& ITLofRscId,double& ExeTime,double& ReadyTime){
    set<double>::iterator pre  = ITLofRscId.begin();
    set<double>::iterator post = ITLofRscId.begin();
    ++post;
    while(post != ITLofRscId.end()) {
        if((*post - *pre) > ExeTime - PrecisionValue && ReadyTime - PrecisionValue < (*post)-ExeTime) {
            return  XY_MAX(*pre, ReadyTime);
        } else {
            ++pre; ++pre; ++post; ++post;
        }
    }
}

void InitProModelOfResAlc(vector<vector<double> >& PMR) {
    for(int i = 0; i < comConst.NumOfTsk; ++i) {
        for(int j : Tasks[i].ElgRsc) {
            PMR[i][j] =  1.0 / Tasks[i].ElgRsc.size();
        }
    }
}

void InitProModelOfTskSch(vector<vector<double>>& PMS, vector<int>& NumOfAncestors, vector<int>& NumOfNonDescendants) {
    vector<int> STS(comConst.NumOfTsk, 0);
    for(int i = 0; i < comConst.NumOfTsk; ++i) {
        for(int k = NumOfAncestors[i]; k < NumOfNonDescendants[i]; ++k) {
            PMS[i][k] = 1;
            ++STS[k];
        }
    }//PMS[i][k] represents the probability that the k-th scheduled task is task i

    for(int k = 0; k < comConst.NumOfTsk; ++k) {
        for(int i = 0; i < comConst.NumOfTsk; ++i) {
            PMS[i][k] = PMS[i][k] / STS[k];
        }
    }
}

chromosome GnrTskLstOfChr(vector<vector<double> >& PMS, vector<double>& eta_TSO) {
    chromosome chrom;
    IntChr(chrom);
    vector<int > upr(comConst.NumOfTsk,0);
    list<int> RTI;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i]=Tasks[i].parents.size();
        if (upr[i]==0)  RTI.push_back(i);
    }
    for (int i = 0; i < comConst.NumOfTsk; i++) {
        double sum = 0;
        for(int k : RTI){
            sum += PMS[k][i] * eta_TSO[k];
        }
        vector<double> SltProb(comConst.NumOfTsk);
        for (int k : RTI) {
            SltProb[k] = PMS[k][i] * eta_TSO[k] / sum;
        }
        double ProbSum = 0, rnd = double(rand()%100) / 100;
        for (int k : RTI) {
            ProbSum += SltProb[k];
            if (rnd + PrecisionValue < ProbSum) {
                chrom.TskSchLst[i] = k;
                break;
            }
        }
        RTI.erase(find(RTI.begin(), RTI.end(), chrom.TskSchLst[i]));
        for (int k = 0; k < Tasks[chrom.TskSchLst[i]].children.size(); ++k) {
            upr[Tasks[chrom.TskSchLst[i]].children[k]]--;
            if (upr[Tasks[chrom.TskSchLst[i]].children[k]] == 0){
                RTI.push_back(Tasks[chrom.TskSchLst[i]].children[k]);
            }
        }
    }
    return chrom;
}

void GnrRscLstOfChr(chromosome& chrom, vector<vector<double>>& PMR) {
    for (int i = 0; i < comConst.NumOfTsk; i++) {
        double rnd = double(rand()%100) / 100;
        double sum = 0;
        for(int j: Tasks[i].ElgRsc){
            sum += PMR[i][j];
            if(rnd < sum) {
                chrom.RscAlcLst[i] = j;
                break;
            }
        }
    }
}

double IFBD_S(chromosome& ch) {
    bool IsFrw = false;
    chromosome NewChrom = ch;
    chromosome OldChrom;
    do {
        OldChrom = NewChrom;
        vector<int> ind(comConst.NumOfTsk);
        IndexSortByValueOnAscend(ind, OldChrom.EndTime);
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            NewChrom.TskSchLst[comConst.NumOfTsk - 1 - i] = ind[i];
        }
        DcdEvl_S(NewChrom, IsFrw);
        CalculateEnergy(NewChrom);
        IsFrw = !IsFrw;
    } while (NewChrom.EnergyConsumption + PrecisionValue < OldChrom.EnergyConsumption);
    if (IsFrw) {
        ch = OldChrom;
    } else {
        ch = NewChrom;
    }
    return ch.EnergyConsumption;
}

//{ Load Balancing with Communication Reduction Improvement (LBCRI)}
void LBCA_S(chromosome& ch) {
    chromosome OldCh = ch;
    vector<double> Id(comConst.NumOfRsc,0);
    //{calculate the loads of resources,find out the set TSK[j] of tasks allocated to resources j; }
    vector<vector<int> > TSK(comConst.NumOfRsc);
    vector<double> TransferDataSize;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscIndex = ch.RscAlcLst[i];
        double TransferDataSize = Tasks[i].ExternalInputFileSizeSum;
        for (int Parent: Tasks[i].parents) {
            int PrnRscId = ch.RscAlcLst[Parent];
            if (RscIndex != PrnRscId) {
                TransferDataSize += ParChildTranFileSizeSum[Parent][i];
            }
        }
        Id[RscIndex] += Tasks[i].length / Rscs[RscIndex].pc + (TransferDataSize + Tasks[i].OFileSizeSum) / VALUE * 8 / Rscs[RscIndex].bw;
        TSK[RscIndex].push_back(i);
    }
    vector<int> ind(comConst.NumOfRsc);
    IndexSortByValueOnAscend(ind, Id);         //sorting according to loads
    int RscWithMinLd = ind[0];          //find out the resource (Rsc) with the lowest load;
    set<int> ST;
    if (abs(Id[RscWithMinLd]) < PrecisionValue) {
        ST.insert(Rscs[RscWithMinLd].ElgTsk.begin(), Rscs[RscWithMinLd].ElgTsk.end());
    } else {
        //traverse the tasks allocated to the resource with lowest load and add their parents and children to set ST
        for (int TaskIndex: TSK[RscWithMinLd]) {
            ST.insert(Tasks[TaskIndex].children.begin(),Tasks[TaskIndex].children.end());
            ST.insert(Tasks[TaskIndex].parents.begin(),Tasks[TaskIndex].parents.end());
        }
        //delete the tasks which have been allocated the resource with lowest load
        for (int i = 0; i < TSK[RscWithMinLd].size(); ++i) {
            ST.erase(TSK[RscWithMinLd][i]);
        }
        //delete the tasks which can not be performed by the resource with lowest load
        for (auto iter = ST.begin(); iter != ST.end();) {
            if (find(Rscs[RscWithMinLd].ElgTsk.begin(), Rscs[RscWithMinLd].ElgTsk.end(), *iter) ==
                Rscs[RscWithMinLd].ElgTsk.end())
                iter = ST.erase(iter);
            else
                ++iter;
        }
        if(ST.empty()){
            ST.insert(Rscs[RscWithMinLd].ElgTsk.begin(), Rscs[RscWithMinLd].ElgTsk.end());
        }
    }
    //Sort the tasks in ST according to the load of the resource to which the task is allocated
    vector<pair<int, double >> t;
    for (auto s:ST) {
        t.push_back(pair<int, double>(s, Id[ch.RscAlcLst[s]]));
    }
    sort(t.begin(), t.end(), SortValueByDescend);
    ch.RscAlcLst[t[0].first] = RscWithMinLd;
    DcdEvl_S(ch, true);
    CalculateEnergy(ch);
    IFBD_S(ch);
    if (OldCh.EnergyConsumption+ PrecisionValue < ch.EnergyConsumption) {
        ch = OldCh;
    }
}

double CalculateEnergy(chromosome &Chrom){
    Chrom.EnergyConsumption = 0;
    set<double>AllTime;
    AllTime.insert( Chrom.EndTime.begin(),Chrom.EndTime.end());
    AllTime.insert(Chrom.StartTime.begin(),Chrom.StartTime.end());
    vector<double>T;
    T.assign(AllTime.begin(),AllTime.end());
    for(int x = 0; x < T.size()-1; ++x) {
        if (fabs(T[x+1]-T[x]) < 1e-5) {
            continue;
        }
        vector<double> HstLd(HstSet.size(),0);
        for(int i = 0; i < comConst.NumOfTsk; ++i) {
            int taskId = Chrom.TskSchLst[i];
            int RscId = Chrom.RscAlcLst[taskId];
            if(Chrom.StartTime[taskId] - PrecisionValue < T[x] && T[x+1]  - PrecisionValue < Chrom.EndTime[taskId]) {
                int taskHT = Rscs[RscId].Hostid;
                double TemLd = Rscs[RscId].pc / HstSet[taskHT].pc;
                HstLd[taskHT] += TemLd;
            }
        }
        for(int k = 0; k < HstSet.size(); ++k) {
            Chrom.EnergyConsumption += (T[x+1] - T[x]) * CalculatePowerByLoad(HstLd[k],k);
        }
    }
    return Chrom.EnergyConsumption;
}

double CalculateECByDelta2(chromosome &ch, int &Number){
    double DeltaEC = 0;
    vector<int>STn;
    int CurTask = ch.TskSchLst[Number];
    int CurHT = Rscs[ch.RscAlcLst[CurTask]].Hostid;
    for (int i = 0; i < Number; ++i) {
        int taskId = ch.TskSchLst[i];
        int RscId = ch.RscAlcLst[taskId];
        int taskHT = Rscs[RscId].Hostid;
        if((ch.EndTime[taskId] > ch.StartTime[CurTask] + PrecisionValue && ch.EndTime[taskId] - PrecisionValue < ch.EndTime[CurTask] && taskHT == CurHT) ||
           (ch.StartTime[taskId] > ch.StartTime[CurTask] - PrecisionValue && ch.StartTime[taskId] + PrecisionValue < ch.EndTime[CurTask] && taskHT == CurHT) ||
           (ch.StartTime[taskId] + PrecisionValue < ch.StartTime[CurTask] && ch.EndTime[taskId] > ch.EndTime[CurTask] + PrecisionValue && taskHT == CurHT)){ //加精度控制-xy-已改-qmq
            STn.push_back(taskId);
        }
    }
    set<double>AllTime;
    vector<double>T;
    AllTime.insert(ch.StartTime[CurTask]);
    AllTime.insert(ch.EndTime[CurTask]);
    for (int i = 0; i < STn.size(); ++i) {
        if (ch.StartTime[STn[i]] > ch.StartTime[CurTask] + PrecisionValue ){
            AllTime.insert(ch.StartTime[STn[i]]);
        }
        if (ch.EndTime[STn[i]] + PrecisionValue < ch.EndTime[CurTask]){
            AllTime.insert(ch.EndTime[STn[i]]);
        }
    }
    T.assign(AllTime.begin(),AllTime.end());
    for(int x = 0; x < T.size()-1; ++x) {
        if (fabs(T[x+1]-T[x]) < 1e-5) {
            continue;
        }
        double HstLd = 0;
        for(int i = 0; i < STn.size(); ++i) {
            int taskId = STn[i];
            int RscId = ch.RscAlcLst[taskId];
            if(ch.StartTime[taskId] - PrecisionValue < T[x] && T[x+1]  - PrecisionValue < ch.EndTime[taskId]) {
                double TemLd = Rscs[RscId].pc / HstSet[CurHT].pc;
                HstLd += TemLd;
            }
        }
        double CurLd = HstLd + Rscs[ch.RscAlcLst[CurTask]].pc / HstSet[CurHT].pc;
        double TemEc = (T[x+1] - T[x]) * (CalculatePowerByLoad(CurLd,CurHT) - CalculatePowerByLoad(HstLd,CurHT));
        DeltaEC += TemEc;
    }
    if (ch.EndTime[CurTask] > ch.MakeSpan + PrecisionValue ){
        for (int i = 0; i < HstSet.size(); ++i) {
            DeltaEC += (ch.EndTime[CurTask] - ch.MakeSpan) * CalculatePowerByLoad(0,i);
        }
    }
    return ch.EnergyConsumption + DeltaEC;
}

void UpdatePMR(vector<vector<double>>& PMR, chromosome& bstChrom){
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        for (int j = 0; j < comConst.NumOfRsc; ++j) {
            int count = 0;
            if (bstChrom.RscAlcLst[i] == j) {
                count = 1;
            }
            PMR[i][j] = (1 - Parameter_TSEDA.theta1) * PMR[i][j] + Parameter_TSEDA.theta1 * count;
        }
    }
}

void UpdatePMS(vector<vector<double>>& PMS, chromosome& bstChrom){
    for(int i = 0; i < comConst.NumOfTsk; ++i) {
        for(int j = 0; j < comConst.NumOfTsk; ++j) {
            int count = 0;
            if(bstChrom.TskSchLst[i] == j) {
                count = 1;
            }
            PMS[j][i] = (1 - Parameter_TSEDA.theta2) * PMS[j][i] + Parameter_TSEDA.theta2 * count;
        }
    }
}

