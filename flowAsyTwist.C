/*************************************************************************
  > File Name: flowAsyTwist.C
 ************************************************************************/
#include <fstream>
#include "flowAsyTwist.h"

using namespace std;


flowAsyTwist::flowAsyTwist():
    filelist(""),
    outfile(""),
    myChain(0){
        Zbin_Size = 4;
        calib = 2;
        perpartQ = 0;
        nevts = 40000;
        //addstat=0;
    }   

//general rule, det0 for Cal, det1 for Trk 
void flowAsyTwist::Init_histos(){

    fout = new TFile(outfile.c_str(), "recreate");
    char name[200];

    for(int tp=0; tp<NTYPE; tp++){
        for(int det0=0; det0<REF_FCAL; det0++){
            for(int ic=0; ic<NCENT; ic++){
                for(int har=0; har<NHAR; har++){
                    for(int re=0; re<4; re++){
                        sprintf(name, "hfg_cmsTrk_0_tp_%d_ic %d_cal_%d_har_%d_re_%d", tp, ic, det0, har, re);
                        hfg_cmsTrk_0[tp][ic][det0][har][re] = new TProfile(name, name, HALF_TRK, 0-0.5, HALF_TRK-0.5);
                        hfg_cmsTrk_0[tp][ic][det0][har][re] -> Sumw2();
                    }
                }
            }
        }
    }

    for(int tp=0; tp<NTYPE; tp++){
        for(int det0=0; det0<REF_FCAL; det0++){
            for(int ic=0; ic<NCENT; ic++){
                for(int har=0; har<NHAR; har++){
                    for(int re=0; re<4; re++){
                        sprintf(name, "hfg_cmsTrk_1_tp_%d_ic %d_cal_%d_har_%d_re_%d", tp, ic, det0, har, re);
                        hfg_cmsTrk_1[tp][ic][det0][har][re] = new TProfile(name, name,HALF_TRK, 0-0.5, HALF_TRK-0.5);
                        hfg_cmsTrk_1[tp][ic][det0][har][re] -> Sumw2();
                    }
                }
            }
        }
    }


    for(int tp=0; tp<NTYPE; tp++){
        for(int det0=0; det0<REF_FCAL; det0++){
            for(int ic=0; ic<NCENT; ic++){
                for(int har=0; har<NHAR; har++){
                    for(int re=0; re<4; re++){
                         sprintf(name, "hbg_cmsTrk_0_tp_%d_ic %d_cal_%d_har_%d_re_%d", tp, ic, det0, har, re);
                        hbg_cmsTrk_0[tp][ic][det0][har][re] = new TProfile(name, name,HALF_TRK, 0-0.5, HALF_TRK-0.5);
                        hbg_cmsTrk_0[tp][ic][det0][har][re] -> Sumw2();
                    }
                }
            }
        }
    }

    for(int tp=0; tp<NTYPE; tp++){
        for(int det0=0; det0<REF_FCAL; det0++){
            for(int ic=0; ic<NCENT; ic++){
                for(int har=0; har<NHAR; har++){
                    for(int re=0; re<4; re++){
                        sprintf(name, "hbg_cmsTrk_1_tp_%d_ic %d_cal_%d_har_%d_re_%d", tp, ic, det0, har, re);
                        hbg_cmsTrk_1[tp][ic][det0][har][re] = new TProfile(name, name,HALF_TRK, 0-0.5, HALF_TRK-0.5);
                        hbg_cmsTrk_1[tp][ic][det0][har][re] -> Sumw2();
                    }
                }
            }
        }
    }

    //************Cal*************************

    for(int det0=0; det0<REF_FCAL; det0++){
        for(int ic=0; ic<NCENT; ic++){
            for(int har=0; har<NHAR; har++){
                for(int re=0; re<4; re++){
                    sprintf(name, "hfg_cmsCal_0_ic_%d_cal_%d_har_%d_re_%d",  ic, det0, har, re);
                    hfg_cmsCal_0[ic][det0][har][re] = new TProfile(name, name, HALF_CAL, 0-0.5, HALF_CAL-0.5);
                    hfg_cmsCal_0[ic][det0][har][re] -> Sumw2();
                }
            }
        }
    }


    for(int det0=0; det0<REF_FCAL; det0++){
        for(int ic=0; ic<NCENT; ic++){
            for(int har=0; har<NHAR; har++){
                for(int re=0; re<4; re++){
                   sprintf(name, "hfg_cmsCal_1_ic_%d_cal_%d_har_%d_re_%d",  ic, det0, har, re);
                    hfg_cmsCal_1[ic][det0][har][re] = new TProfile(name, name, HALF_CAL, 0-0.5, HALF_CAL-0.5);
                    hfg_cmsCal_1[ic][det0][har][re] -> Sumw2();
                }
            }
        }
    }

    for(int det0=0; det0<REF_FCAL; det0++){
        for(int ic=0; ic<NCENT; ic++){
            for(int har=0; har<NHAR; har++){
                for(int re=0; re<4; re++){
                    sprintf(name, "hbg_cmsCal_0_ic_%d_cal_%d_har_%d_re_%d",  ic, det0, har, re);
                    hbg_cmsCal_0[ic][det0][har][re] = new TProfile(name, name, HALF_CAL, 0-0.5, HALF_CAL-0.5);
                    hbg_cmsCal_0[ic][det0][har][re] -> Sumw2();
                }
            }
        }
    }

    for(int det0=0; det0<REF_FCAL; det0++){
        for(int ic=0; ic<NCENT; ic++){
            for(int har=0; har<NHAR; har++){
                for(int re=0; re<4; re++){
                     sprintf(name, "hbg_cmsCal_1_ic_%d_cal_%d_har_%d_re_%d",  ic, det0, har, re);
                    hbg_cmsCal_1[ic][det0][har][re] = new TProfile(name, name, HALF_CAL, 0-0.5, HALF_CAL-0.5);
                    hbg_cmsCal_1[ic][det0][har][re] -> Sumw2();
                }
            }
        }
    }


    sprintf(name, "hbad1");
    hbad1 = new TH1D(name, name, 100, 0-0.5, 100-0.5);
    hbad1->Sumw2();

    for(int det=0; det<NDET_CAL; det++){
        sprintf(name,"hbad2_det%d", det);
        hbad2[det] = new TH1D(name, name, 100, 0-0.5, 100-0.5);
        hbad2[det]->Sumw2();
    }

    for(int det=0; det<NDET_CAL; det++){
        sprintf(name,"hbad2Cal_det%d", det);
        hbad2Cal[det] = new TH1D(name, name, 50, -100, 0);
        hbad2Cal[det]->Sumw2();
    }

    for(int ic=0; ic<NCENT; ic++){
        sprintf(name, "track_check_ic_%d",ic);
        track_check[ic] = new TH1D(name, name, 500, 0, 4000);
        track_check[ic] ->Sumw2();
    }

    for(int ic=0; ic<NCENT; ic++){
        sprintf(name, "trackEt_check_ic_%d",ic);
        trackEt_check[ic] = new TH2D(name, name, 40, 0, 20, 50, 0, 5);
        trackEt_check[ic] ->Sumw2();
    }


    hcentbin = new TH1D("hcentbin", "hcentbin", NCENT, 0-0.5, NCENT-0.5);
    hcentbin ->Sumw2();
    hcentrality = new TH1D("hcentrality", "hcentrality", 100, 0, 100);
    hcentrality ->Sumw2();

    for(int ic=0; ic<NCENT; ic++){
        for(int det1=0; det1<HALF_TRK; det1++){
            for(int ih=0; ih<NHAR; ih++){
                sprintf(name, "hq_trk_ic%d_det%d_ih%d",ic, det1, ih);
                hq_trk[ic][det1][ih] = new TH1D(name, name, 200, 0, 100);
                hq_trk[ic][det1][ih] ->Sumw2();
            }
        }
    }

    for(int ic=0; ic<NCENT; ic++){
        for(int det0=0; det0<REF_FCAL; det0++){
            for(int ih=0; ih<NHAR; ih++){
                sprintf(name, "hq_cal_ic%d_det%d_ih%d", ic, det0, ih);
                hq_cal[ic][det0][ih] = new TH1D(name, name, 200, 0, 100);
                hq_cal[ic][det0][ih] ->Sumw2();
            }
        }
    }

    for(int ic=0; ic<NCENT; ic++){
        for(int det1=0; det1<HALF_TRK; det1++){
            for(int ih=0; ih<NHAR; ih++){
                sprintf(name, "hpsi_trk_ic%d_det%d_ih%d",ic, det1, ih);
                hpsi_trk[ic][det1][ih] = new TH1D(name, name, 50, -PI, PI);
                hpsi_trk[ic][det1][ih] ->Sumw2();
            }
        }
    }

    for(int ic=0; ic<NCENT; ic++){
        for(int det0=0; det0<REF_FCAL; det0++){
            for(int ih=0; ih<NHAR; ih++){
                sprintf(name, "hpsi_cal_ic%d_det%d_ih%d", ic, det0, ih);
                hpsi_cal[ic][det0][ih] = new TH1D(name, name, 50, -PI, PI);
                hpsi_cal[ic][det0][ih] ->Sumw2();
            }
        }
    }

    for(int ic=0; ic<NCENT; ic++){
        sprintf(name, "hvtx_z_ic%d", ic);
        hvtx_z[ic] = new TH1D(name, name, 200, -100, 100);
        hvtx_z[ic] ->Sumw2();
    }

    for(int ic=0; ic<NCENT; ic++){
        sprintf(name, "hvtx_zbin_ic%d", ic);
        hvtx_zbin[ic] = new TH1D(name, name, 200,0, 100);
        hvtx_zbin[ic] ->Sumw2();
    }

    cout<<"Finish Initializing histograms"<<endl;
}

void flowAsyTwist::Save_histos(){
    fout->Write();
    cout<<"Writing to "<<outfile<<endl;
    //std::cout<<" Finish Writing TProfiles "<<endl;
    fout->Close();
}

int  flowAsyTwist::get_Q_complex(){

    for(int tp=0; tp<NTYPE; tp++){
        for(int det1=0; det1<NDET_TRK; det1++){
            if(QwTrk[tp][det1] > 0){
                for(int clib=0; clib<3; clib++){
                    for(int har=0; har<NHAR; har++){
                        float qx = 0;
                        float qy = 0;
                        if(perpartQ==0){//scalar product
                            qx = QxTrk[clib][tp][det1][har]; 
                            qy = QyTrk[clib][tp][det1][har]; 
                        }else{//per particle flow
                            qx = QxTrk[clib][tp][det1][har]/QwTrk[tp][det1]; 
                            qy = QyTrk[clib][tp][det1][har]/QwTrk[tp][det1]; 
                        }
                        q_Trk[clib][tp][det1][har] = TComplex( qx, qy);
                        if(q_Trk[clib][tp][det1][har].Im() != qy) { cout<<"Call 911"<<endl; return 1;}
                        if(q_Trk[clib][tp][det1][har].Re() != qx) { cout<<"Call 911"<<endl; return 1;}
                    }
                }
            }else{
                bad1 = 1;
            }

        }
    }

    for(int det0=0; det0<NDET_CAL; det0++){
        if( QwCal[det0] >0 ) {  
            for(int clib=0; clib<3; clib++){
                for(int har=0; har<NHAR; har++){
                    float qx=0;
                    float qy=0;
                    if(perpartQ==0){
                        qx = QxCal[clib][det0][har];
                        qy = QyCal[clib][det0][har];
                    }else{
                        qx = QxCal[clib][det0][har]/QwCal[det0];
                        qy = QyCal[clib][det0][har]/QwCal[det0];
                    }
                    q_Cal[clib][det0][har] = TComplex(qx, qy);
                    if(q_Cal[clib][det0][har].Im() != qy){ cout<<"Call 911"<<endl; return 1; }
                    if(q_Cal[clib][det0][har].Re() != qx){ cout<<"Call 911"<<endl; return 1; }
                }
            }
        }else{
            bad1 = 1;
            bad2[det0] =1;
        }
    }

    return bad1;
}

void flowAsyTwist::run(){

    myChain = new TChain("QTree");
    ifstream lis(filelist.c_str());
    int cnt=0;
    while(!lis.eof())
    {
        string filename;
        lis >> filename;

        if(cnt>=from&&cnt<to)
        {
            if(!filename.empty()) myChain->Add(filename.c_str());
            cnt++;
        }
    }
    //Fill the tree setBranchAddress
    myChain->SetBranchAddress("RunNumber",   &RunNumber);
    myChain->SetBranchAddress("EventNumber", &EventNumber);
    myChain->SetBranchAddress("Fcal_Et",     &Fcal_Et);
    myChain->SetBranchAddress("Centrality",  &Centrality );
    myChain->SetBranchAddress("Fcal_Et_p",  &Fcal_Et_p);
    myChain->SetBranchAddress("Fcal_Et_n",  &Fcal_Et_n);
    //    myChain->SetBranchAddress("vx_n",  &vx_n);
    myChain->SetBranchAddress("vx_z",  &vx_z);
    myChain->SetBranchAddress("trk_nQ",  &trk_nQ);
    myChain->SetBranchAddress("trk_n", &trk_n);
    myChain->SetBranchAddress("QxTrk", QxTrk);
    myChain->SetBranchAddress("QyTrk", QyTrk);
    myChain->SetBranchAddress("QwTrk", QwTrk);
    myChain->SetBranchAddress("QxCal", QxCal);
    myChain->SetBranchAddress("QyCal", QyCal);
    myChain->SetBranchAddress("QwCal", QwCal);

    EventPool.resize(NCENT);
    for(int ic=0; ic<NCENT; ic++){
        EventPool[ic].resize(nz);
    }

    if(calib ==0 ){cout<<"Use observed Q vector"<<endl;}
    if(calib ==1 ){cout<<"Use recentered Q vector"<<endl;}
    if(calib ==2 ){cout<<"Use recentered+flattened Q vector"<<endl;}
    if(perpartQ==0) {cout<<"Q scalar product way"<<endl;}
    if(perpartQ==1) {cout<<"Q per particle flow way"<<endl;}
    cout<<"Using Zbin_Size "<<Zbin_Size<<endl;
    //cout<<"Add more stat? "<<addstat<<endl;
    cout<<"Want to run "<<nevts<<endl;

    int Nevents = myChain->GetEntries();
    cout<<"We are reading Nevents "<<Nevents<<endl;
    for(int iev=0; iev<Nevents; iev++){
        myChain->GetEntry(iev);
        bad1 = 0;
        for(int det=0; det<NDET_CAL; det++){
            bad2[det] = 0;
        }
        centbin = get_centbin(Centrality);
        zbin = get_zPool(vx_z);
        hvtx_z[centbin] ->Fill(vx_z);
        hvtx_zbin[centbin] ->Fill(zbin);

        if(iev%20000==0) cout<<iev<<endl;

        int ntrkcut[20]={1800,1475,1200, 975, 775, 600, 450, 350,250,180,120, 70, 50, 25, 10,  2,  1,  1,  1,  1};//2010 cuts

        int cenbin0 = (int)(Centrality/5.0);
        if(trk_n<ntrkcut[cenbin0]) continue;
        if(centbin < 0) continue;
        if(zbin < 0) continue;
        track_check[centbin]->Fill(trk_n);
        trackEt_check[centbin]->Fill(trk_n*1.0, Fcal_Et);

        hcentrality->Fill(Centrality);

        int bad = get_Q_complex();

        if(bad1 ==1 ) hbad1->Fill(Centrality); 
        for(int det=0; det<NDET_CAL; det++){
            if(bad2[det] == 1){
                hbad2[det] ->Fill(Centrality);
                hbad2Cal[det] ->Fill(QwCal[det]);
            } 
        }

        if(bad==1)  continue; 
        //if(iev>=nevts){   cout<<iev<<" evts has been finished"<<endl; break;}
        Event* event0 = new Event( iev, centbin, zbin);
        hcentbin ->Fill(event0->getCent());
        event0->setqTrx(q_Trk);
        event0->setqCal(q_Cal);

        bool fsame = FillForeground( event0 );
        if(fsame != true) {cout<<"Come on, Wrong FillForeground!!!"<<endl;  delete event0; continue;}

        int KiddiePool = EventPool.at(centbin).at(zbin).size();

        bool mix;
        if(KiddiePool>0){
            for(int j=0; j<KiddiePool; j++){
                mix  =  FillBackground(event0, EventPool.at(centbin).at(zbin).at(j) );
            }
            if(mix != true) {cout<<"Come on, Wrong Background!!!"<<endl;}
        }

        if(KiddiePool<3){
         EventPool.at(centbin).at(zbin).push_back( event0 );
        }else if( KiddiePool >= 3){
            delete EventPool.at(centbin).at(zbin).at(0);
            EventPool.at(centbin).at(zbin).at(0) = 0;
            EventPool.at(centbin).at(zbin).erase( EventPool.at(centbin).at(zbin).begin() );
            EventPool.at(centbin).at(zbin).push_back( event0 );
        }
    }//end of event loop
}


bool flowAsyTwist::FillForeground( EVENT_PTR evt ){
    //int zbin = evt->getZbin();
    int cenbin = evt->getCent();
    int Calib = calib;
    TComplex Q_Trk[NTYPE][NDET_TRK][NHAR];
    TComplex Q_Cal[NDET_CAL][NHAR];
    TComplex zeroo(0, 0);
    int countzero0=0;
    int countzero1=0;

    for(int tp=0; tp<NTYPE; tp++){
        for(int det=0; det<NDET_TRK; det++){
            for(int har=0; har<NHAR; har++){
                Q_Trk[tp][det][har] = evt->getqTrk(Calib, tp, det, har);
                if( har==0 && Q_Trk[tp][det][har].Im()==0 && Q_Trk[tp][det][har].Re()==0 ) countzero0++;
            }
        }
    }
    for(int det=0; det<NDET_CAL; det++){
        for(int har=0; har<NHAR; har++){
            Q_Cal[det][har] = evt->getqCal(Calib, det, har);
            if( har==0 && Q_Cal[det][har].Im()==0 && Q_Cal[det][har].Re()==0) countzero1++;
        }
    }   

    for(int itrk=0; itrk<HALF_TRK; itrk++){
        int det1 = Trk_For[itrk];
        for(int ih=0; ih<NHAR; ih++){
            hq_trk[cenbin][itrk][ih] ->Fill(Q_Trk[0][det1][ih].Rho());
        }
    }

    for(int ical=0; ical<REF_FCAL; ical++){
        int det0 = Ref_For[ical];
        for(int ih=0; ih<NHAR; ih++){
            hq_cal[cenbin][ical][ih] ->Fill(Q_Cal[det0][ih].Rho());
        }
    }

    for(int itrk=0; itrk<HALF_TRK; itrk++){
        int det1 = Trk_For[itrk];
        for(int ih=0; ih<NHAR; ih++){
            hpsi_trk[cenbin][itrk][ih] ->Fill(Q_Trk[0][det1][ih].Theta());
        }
    }

    for(int ical=0; ical<REF_FCAL; ical++){
        int det0 = Ref_For[ical];
        for(int ih=0; ih<NHAR; ih++){
            hpsi_cal[cenbin][ical][ih] ->Fill(Q_Cal[det0][ih].Theta());
        }
    }

    if(countzero0>NTYPE*3 || countzero1>6){
        cout<<"Reallly  badddd!!!"<<evt->getID()<< " "<< endl;  
        return false;
    } 

    TComplex q[25];
    for(int j=0; j<25; j++) q[j] = TComplex(0.0, 0.0); 
    for(int har=0; har<NHAR; har++){
        //calc <q(ref)q*(ref)>  

        for(int det0=0; det0<REF_FCAL; det0++){
            for(int det1=0; det1<HALF_TRK; det1++){

                q[0] = Q_Cal[ Cal_For[det1] ][har]; //+eta
                q[1] = Q_Cal[ Cal_Bac[det1] ][har]; //-eta
                q[2] = Q_Cal[ Ref_For[det0] ][har];  //+ref
                q[3] = Q_Cal[ Ref_Bac[det0] ][har];  //-ref

                q[4] = q[1]*TComplex::Conjugate(q[2]); //-eta +ref numerator
                q[5] = q[0]*TComplex::Conjugate(q[2]); //+eta +ref
                q[6] = q[0]*TComplex::Conjugate(q[3]); //+eta -ref numerator
                q[7] = q[1]*TComplex::Conjugate(q[3]); //-eta -ref

                hfg_cmsCal_0[cenbin][det0][har][0] ->Fill(det1, q[4].Re());
                hfg_cmsCal_0[cenbin][det0][har][1] ->Fill(det1, q[4].Im());
                hfg_cmsCal_1[cenbin][det0][har][0] ->Fill(det1, q[5].Re());
                hfg_cmsCal_1[cenbin][det0][har][1] ->Fill(det1, q[5].Im());

                hfg_cmsCal_0[cenbin][det0][har][0+2] ->Fill(det1, q[6].Re());
                hfg_cmsCal_0[cenbin][det0][har][1+2] ->Fill(det1, q[6].Im());
                hfg_cmsCal_1[cenbin][det0][har][0+2] ->Fill(det1, q[7].Re());
                hfg_cmsCal_1[cenbin][det0][har][1+2] ->Fill(det1, q[7].Im());

            }
        }

    for(int tp=0; tp<NTYPE; tp++){
        for(int det0=0; det0<REF_FCAL; det0++){
                    for(int det1=0; det1<HALF_TRK; det1++){
                    q[14] = Q_Trk[tp][Trk_For[det1]][har]; //+eta
                    q[15] = Q_Trk[tp][Trk_Bac[det1]][har]; //-eta
                    q[16] = Q_Cal[ Ref_For[det0] ][har];   //+ref
                    q[17] = Q_Cal[ Ref_Bac[det0] ][har];   //-ref

                    q[18] = q[15]*TComplex::Conjugate(q[16]); //numerator
                    q[19] = q[14]*TComplex::Conjugate(q[16]);
                    q[20] = q[14]*TComplex::Conjugate(q[17]); //numerator
                    q[21] = q[15]*TComplex::Conjugate(q[17]);

                    hfg_cmsTrk_0[tp][cenbin][det0][har][0] ->Fill(det1, q[18].Re());
                    hfg_cmsTrk_0[tp][cenbin][det0][har][1] ->Fill(det1, q[18].Im());
                    hfg_cmsTrk_1[tp][cenbin][det0][har][0] ->Fill(det1, q[19].Re());
                    hfg_cmsTrk_1[tp][cenbin][det0][har][1] ->Fill(det1, q[19].Im());

                    hfg_cmsTrk_0[tp][cenbin][det0][har][0+2] ->Fill(det1, q[20].Re());
                    hfg_cmsTrk_0[tp][cenbin][det0][har][1+2] ->Fill(det1, q[20].Im());
                    hfg_cmsTrk_1[tp][cenbin][det0][har][0+2] ->Fill(det1, q[21].Re());                   
                    hfg_cmsTrk_1[tp][cenbin][det0][har][1+2] ->Fill(det1, q[21].Im());  

                }       
            }
        }

    }//end of harmonic loop
    return true;
}


bool flowAsyTwist::FillBackground( EVENT_PTR evt0, EVENT_PTR evt1){
    int cenbin = evt0->getCent();
    int Calib = calib;

    TComplex Q_Trk[4][NTYPE][NDET_TRK][NHAR];
    TComplex Q_Cal[4][NDET_CAL][NHAR];

    for(int tp=0; tp<NTYPE; tp++){
        for(int det=0; det<NDET_TRK; det++){
            for(int har=0; har<NHAR; har++){
                Q_Trk[0][tp][det][har] = evt0->getqTrk(Calib, tp, det, har);
                Q_Trk[1][tp][det][har] = evt1->getqTrk(Calib, tp, det, har);
                //     Q_Trk[2][tp][det][har] = evt2->getqTrk(Calib, tp, det, har);
                //     Q_Trk[3][tp][det][har] = evt3->getqTrk(Calib, tp, det, har);
            }
        }
    }
    for(int det=0; det<NDET_CAL; det++){
        for(int har=0; har<NHAR; har++){
            Q_Cal[0][det][har] = evt0->getqCal(Calib, det, har);
            Q_Cal[1][det][har] = evt1->getqCal(Calib, det, har);
            //    Q_Cal[2][det][har] = evt2->getqCal(Calib, det, har);
            //    Q_Cal[3][det][har] = evt3->getqCal(Calib, det, har);
        }
    }   

    TComplex q[25];
    for(int j=0; j<25; j++) q[j] = TComplex(0.0, 0.0); 

    // evt0: q(eta), evt1: q(-eta), evt2: q(ref), evt3: q(-ref)  

    for(int har=0; har<NHAR; har++){

        for(int det0=0; det0<REF_FCAL; det0++){
            for(int det1=0; det1<HALF_TRK; det1++){

                q[0] = Q_Cal[0][ Cal_For[det1] ][har]; //+eta
                q[1] = Q_Cal[0][ Cal_Bac[det1] ][har]; //-eta
                q[2] = Q_Cal[1][ Ref_For[det0] ][har];  //+ref
                q[3] = Q_Cal[1][ Ref_Bac[det0] ][har];  //-ref

                q[4] = q[1]*TComplex::Conjugate(q[2]); //-eta +ref numerator
                q[5] = q[0]*TComplex::Conjugate(q[2]); //+eta +ref
                q[6] = q[0]*TComplex::Conjugate(q[3]); //+eta -ref numerator
                q[7] = q[1]*TComplex::Conjugate(q[3]); //-eta -ref

                hfg_cmsCal_0[cenbin][det0][har][0] ->Fill(det1, q[4].Re());
                hfg_cmsCal_0[cenbin][det0][har][1] ->Fill(det1, q[4].Im());
                hfg_cmsCal_1[cenbin][det0][har][0] ->Fill(det1, q[5].Re());
                hfg_cmsCal_1[cenbin][det0][har][1] ->Fill(det1, q[5].Im());

                hfg_cmsCal_0[cenbin][det0][har][0+2] ->Fill(det1, q[6].Re());
                hfg_cmsCal_0[cenbin][det0][har][1+2] ->Fill(det1, q[6].Im());
                hfg_cmsCal_1[cenbin][det0][har][0+2] ->Fill(det1, q[7].Re());
                hfg_cmsCal_1[cenbin][det0][har][1+2] ->Fill(det1, q[7].Im());

            }
        }



        for(int det0=0; det0<REF_FCAL; det0++){
            for(int det1=0; det1<HALF_TRK; det1++){
                for(int tp=0; tp<NTYPE; tp++){
                    q[14] = Q_Trk[0][tp][Trk_For[det1]][har];  //+eta
                    q[15] = Q_Trk[0][tp][Trk_Bac[det1]][har];  //-eta
                    q[16] = Q_Cal[1][ Ref_For[det0] ][har];    //+ref
                    q[17] = Q_Cal[1][ Ref_Bac[det0] ][har];    //-ref

                    q[18] = q[15]*TComplex::Conjugate(q[16]); //numerator
                    q[19] = q[14]*TComplex::Conjugate(q[16]);

                    q[20] = q[14]*TComplex::Conjugate(q[17]);  //numerator
                    q[21] = q[15]*TComplex::Conjugate(q[17]);

                    hbg_cmsTrk_0[tp][cenbin][det0][har][0] ->Fill(det1, q[18].Re());
                    hbg_cmsTrk_0[tp][cenbin][det0][har][1] ->Fill(det1, q[18].Im());
                    hbg_cmsTrk_1[tp][cenbin][det0][har][0] ->Fill(det1, q[19].Re());
                    hbg_cmsTrk_1[tp][cenbin][det0][har][1] ->Fill(det1, q[19].Im());

                    hbg_cmsTrk_0[tp][cenbin][det0][har][0+2] ->Fill(det1, q[20].Re());
                    hbg_cmsTrk_0[tp][cenbin][det0][har][1+2] ->Fill(det1, q[20].Im());
                    hbg_cmsTrk_1[tp][cenbin][det0][har][0+2] ->Fill(det1, q[21].Re());
                    hbg_cmsTrk_1[tp][cenbin][det0][har][1+2] ->Fill(det1, q[21].Im());  

                }       
            }
        }
    }//end of harmonics
    return true;
}


int flowAsyTwist::get_centbin(float cent){
    float cent0[] = {0, 5, 10, 20, 30, 40, 50, 60};
    float cent1[] = {5, 10, 20, 30, 40, 50, 60, 70};
    int centbin = -1;
    for(int ic=0; ic<NCENT; ic++){
        if(cent>=cent0[ic] && cent<cent1[ic]) centbin = ic;
    }
    if(centbin<0 || centbin >=NCENT) return -1;
    return centbin;
}

int flowAsyTwist::get_centPool(float cent){
    float cent0[] = {0, 5, 10, 20, 30, 40, 50, 60};
    float cent1[] = {5, 10, 20, 30, 40, 50, 60, 70};
    int centbin = -1;
    for(int ic=0; ic<NCENT; ic++){
        if(cent>=cent0[ic] && cent<cent1[ic]) centbin = ic;
    }
    if(centbin<0 || centbin >=NCENT) return -1;
    return centbin;
}

int flowAsyTwist::get_zPool(float z){
    int bin = -1;
    if(fabs(z)>=100) return -1;
    bin = (z+100.0)/Zbin_Size;  // pools
    if(bin<0||bin>=nz) return -1;
    return bin;
}

