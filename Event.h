/*************************************************************************
  > File Name: Event.h
 ************************************************************************/
#include <iostream>
#include <complex> 
#include "TComplex.h"

using namespace std;

enum{ 
    NTYPE = 4,  //w_eff w_wt, w_eff wo_wt, wo_eff w_wt, wo_eff wo_wt
    NDET_CAL=32, //4*2 FCal + 9*2 Calo + 6
    NDET_TRK=14, //7*2 ID
    NHAR=5,
};

class Event{

    public:
        Event( int i, int cen, int zb);
        TComplex  getqTrk(int c, int t, int d, int h){ return qTrk[c][t][d][h]; }
        TComplex  getqCal(int c, int d, int h){ return qCal_new[c][d][h];}
        void setqTrx(TComplex qtrx[3][NTYPE][NDET_TRK][NHAR]);
        void setqCal(TComplex qcal[3][NDET_CAL][NHAR]);
        int  getCent(){ return centbin;}
        int  getZbin(){ return zbin;}
        int  getID(){ return id;}
    private:
        int id, centbin;
        int zbin;
        TComplex qTrk[3][NTYPE][NDET_TRK][NHAR];
        TComplex qCal_new[3][NDET_CAL][NHAR];
};

Event::Event(int i, int cen, int zb ):id(i),centbin(cen),zbin(zb){}

void Event::setqTrx( TComplex qtrk[3][NTYPE][NDET_TRK][NHAR]){

    for(int calib=0; calib<3; calib++){
        for(int tp=0; tp<NTYPE; tp++){
            for(int det=0; det<NDET_TRK; det++){
                for(int har=0; har<NHAR; har++){
                    qTrk[calib][tp][det][har] = qtrk[calib][tp][det][har];
                }
            }
        }
    }

}

void Event::setqCal(TComplex qcal[3][NDET_CAL][NHAR]){

    for(int calib=0; calib<3; calib++){
        for(int det=0; det<NDET_CAL; det++){
            for(int har=0; har<NHAR; har++){
                qCal_new[calib][det][har] = qcal[calib][det][har];
            }
        }
    }

}

typedef Event* EVENT_PTR;
