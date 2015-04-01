/*************************************************************************
  > File Name: flowAsyTwist.h
 ************************************************************************/
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "TComplex.h"
#include "TMath.h"
#include "TProfile.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include <Event.h>

using namespace std;

enum{ 
    NCENT = 8,
    nz=50,

    HALF_TRK = 7,
    HALF_CAL = 8,
    REF_FCAL = 7, //4+3

};

const double PI=acos(-1.0);

float  TrkLowerEta[] =  {-2.5, -2.0, -1.5, -1.0, -0.5, -0.5, -0.2, 0.0, 0.0, 0.2, 0.5, 1.0, 1.5, 2.0};
float  TrkUpperEta[] =  {-2.0, -1.5, -1.0, -0.5, -0.2,  0.0,  0.0, 0.2, 0.5, 0.5, 1.0, 1.5, 2.0, 2.5};
float  CalLowerEta[] =  {-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,-0.5,-0.2,0.0,0.0,0.2,0.5,1.0,1.5,2.0,2.5,3.0,   -4.9,-4.4,-4.0,-3.6,3.2,3.6,4.0,4.4,  -4.9,-4.0,3.2,4.0,  -4.9,3.2};
float  CalUpperEta[] =  {-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,-0.2, 0.0, 0.0,0.2,0.5,0.5,1.0,1.5,2.0,2.5,3.0,3.5,   -4.4,-4.0,-3.6,-3.2,3.6,4.0,4.4,4.9,  -4.0,-3.2,4.0,4.9,  -3.2,4.9};

int Trk_For[HALF_TRK] = { 7, 8, 9, 10, 11, 12, 13 }; 
int Trk_Bac[HALF_TRK] = { 6, 5, 4,  3,  2,  1,  0 };

int Cal_For[HALF_CAL] = {9,10,11,12,13,14,15,16 };
int Cal_Bac[HALF_CAL] = {8, 7, 6, 5, 4, 3, 2, 1 };

int Ref_For[REF_FCAL] = {  22, 23, 24, 25, 28, 29, 31};
int Ref_Bac[REF_FCAL] = {  21, 20, 19, 18, 27, 26, 30};

typedef vector<EVENT_PTR>  evt1; //push back 3 events 
typedef vector<evt1> evt2; // 200 zbin
typedef vector<evt2> evt3; //NCENT

class flowAsyTwist{

    public:
        flowAsyTwist();

        string filelist;
        string outfile;

        int from;
        int to;
        int nevts;
        int calib; //0 for observed, 1 for recentered, 2 for recentered+flattening
        int Zbin_Size;
        int perpartQ;
        //int addstat;

        void Init_histos();
        void Save_histos();
        void run();
        bool FillForeground(EVENT_PTR evt);
        bool FillBackground(EVENT_PTR evt0,EVENT_PTR evt1);
        int get_Q_complex();

        int get_centbin(float cent);
        int get_centPool(float cent);
        int get_zPool(float z);

    private:

        TFile *fin;
        TFile *fout;

        TChain* myChain;
        //Event information
        int RunNumber;
        int EventNumber;
        float Fcal_Et;
        float Fcal_Et_p;
        float Fcal_Et_n;
        int Centrality;
        int vx_n;
        float vx_z;
        int trk_n;
        int trk_nQ;
        float  QxTrk[3][NTYPE][NDET_TRK][NHAR];
        float  QyTrk[3][NTYPE][NDET_TRK][NHAR];
        float  QwTrk[NTYPE][NDET_TRK];
        float  QxCal[3][NDET_CAL][NHAR];
        float  QyCal[3][NDET_CAL][NHAR];
        float  QwCal[NDET_CAL];

        int centbin;
        int zbin;
     

        evt3 EventPool;

        TComplex q_Trk[3][NTYPE][NDET_TRK][NHAR];
        TComplex q_Cal[3][NDET_CAL][NHAR];

        TProfile* hfg_cmsTrk_0[NTYPE][NCENT][REF_FCAL][NHAR][4];  //numerator
        TProfile* hfg_cmsTrk_1[NTYPE][NCENT][REF_FCAL][NHAR][4];  //denominator
        TProfile* hbg_cmsTrk_0[NTYPE][NCENT][REF_FCAL][NHAR][4];  //numerator
        TProfile* hbg_cmsTrk_1[NTYPE][NCENT][REF_FCAL][NHAR][4];  //denominator

        TProfile* hfg_cmsCal_0[NCENT][REF_FCAL][NHAR][4];  //numerator
        TProfile* hfg_cmsCal_1[NCENT][REF_FCAL][NHAR][4];  //denominator
        TProfile* hbg_cmsCal_0[NCENT][REF_FCAL][NHAR][4];  //numerator
        TProfile* hbg_cmsCal_1[NCENT][REF_FCAL][NHAR][4];  //denominator

        int bad1;
        int bad2[NDET_CAL];
        TH1* hbad1;
        TH1* hbad2[NDET_CAL];
        TH1* hbad2Cal[NDET_CAL];
        //TH1* hbad3[NTYPE][NDET_TRK];

        TH1* track_check[NCENT];
        TH2* trackEt_check[NCENT];
        TH1* hcentbin;
        TH1* hcentrality;

        TH1* hq_trk[NCENT][HALF_TRK][NHAR];
        TH1* hq_cal[NCENT][REF_FCAL][NHAR];
        TH1* hpsi_trk[NCENT][HALF_TRK][NHAR];
        TH1* hpsi_cal[NCENT][REF_FCAL][NHAR];
        TH1D* hvtx_z[NCENT];
        TH1D* hvtx_zbin[NCENT];
};



