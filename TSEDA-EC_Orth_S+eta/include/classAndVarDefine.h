//
// Created by qqq on 2021/9/15.
//

#ifndef FRAME_CLASSANDVARDEFINE_H
#include <vector>
#include <string>
#include <set>
#include <map>
#include <list>
#define FRAME_CLASSANDVARDEFINE_H
using namespace std;

class vfile {
public:
    string FileName;     //file name
    int source;          //the source of file, -1:from the shared server; i: from task i
    double size;         //the size of file
};

class Task{
public:
    double length;
    vector<int> ElgRsc;
    vector<int> parents;
    vector<int> children;
    vector<vfile> IFile;
    vector<vfile> OFile;
    double IFileSizeSum = 0.0;
    double OrignalInputFileSizeSum = 0.0;
    double ExternalInputFileSizeSum = 0.0; //the size of all input files which are not generated by their parent tasks
    double OFileSizeSum = 0.0;
};

class Resource{
public:
    int Hostid;
    vector<int> ElgTsk;
    double pc,bw;
    Resource(){

    };
    Resource(int id,double pc,double bw){
        this->Hostid = id;
        this->pc = pc;
        this->bw = bw;
    }
};

class HT {
public:
    vector<int> RscslnHt;
    double pc, bw;
    HT(vector<int> RscslnHt, double pc, double bw) {
        this->RscslnHt = RscslnHt;
        this->pc = pc;
        this->bw = bw;
    }
};

class chromosome{
public:
    vector<int> RscAlcLst;       //resources allocation, task-to-resource mapping, (Match String)
    vector<int> TskSchLst;       //task scheduling order (Scheduling String)
    vector<double> Code_RK;
    vector<list<int>> Code_TD;
    vector<double> RscAlcPart;
    vector<double> TskSchPart;
    vector<double> VTskSchPart;
    vector<double> VRscAlcPart;
    vector<double> StartTime, EndTime;
    double MakeSpan;
    double EnergyConsumption;

    bool operator<(const chromosome& otherChrom) const{
        return this->EnergyConsumption + 1e-6 < otherChrom.EnergyConsumption;
    }
};

class ComConst{
public:
    int NumOfTsk;
    int NumOfRsc;
};

class Paramet_TSEDA {
public:
    int NumOfChromPerPop;
    double theta1;
    double theta2;
    double fdhi;
    int NumOfImproveOfPop;
    double RunTimeRatioOfStg1;
};

class Orthogonal{
public:
    double PopSizeFactor ;       //the number of chromosomes in each population
    double theta1;
    double theta2;
    double fdhi;
    double ImprovementRate;
    double RunTimeRatioOfStg1;
};

extern vector<Task> Tasks;
extern vector<vector<int>> TskLstInLvl; //task list (set) in each level;
extern vector<vector<double> > ParChildTranFileSizeSum;
extern vector<set<int>> Descendants;
extern vector<set<int>> Ancestors;
extern vector<int> LevelIdOfTask;
extern vector<Resource> Rscs;
extern ComConst comConst;
extern vector<HT> HstSet;
extern double ModelScale;
extern Paramet_TSEDA Parameter_TSEDA;
extern vector<Orthogonal> orthogonal;

#endif //FRAME_CLASSANDVARDEFINE_H
