#include <TFile.h>
#include <TTreeReader.h>
#include "TreeParticles.hpp"



int run(){
    
    //TFile file_1("/user/dmeuser/master/data_johannes/v19/SinglePhoton_03Feb2017.root","read");
    //TFile file_2("/user/dmeuser/master/data/v01_danilo/SinglePhoton_03Feb2017.root","read");
    //TFile file_2("~/CMSSW_TreeWriter/CMSSW_8_0_26_patch1/src/singleEvent.root","read");
    TFile file_1("/user/dmeuser/master/data_johannes/v19/ZZ.root","read");
    TFile file_2("/user/dmeuser/master/data/ZZ.root","read");
    
    TTreeReader reader_1("TreeWriter/eventTree", &file_1);
    TTreeReader reader_2("TreeWriter/eventTree", &file_2);
    
    TTreeReaderValue<ULong64_t> evtNo_1(reader_1, "evtNo");
    TTreeReaderValue<ULong64_t> evtNo_2(reader_2, "evtNo");
    TTreeReaderValue<UInt_t> lumNo_1(reader_1, "lumNo");
    TTreeReaderValue<UInt_t> lumNo_2(reader_2, "lumNo");
    TTreeReaderValue<UInt_t> runNo_1(reader_1, "runNo");
    TTreeReaderValue<UInt_t> runNo_2(reader_2, "runNo");
    TTreeReaderValue<std::vector<tree::Photon>> photons_1(reader_1, "photons");
    TTreeReaderValue<std::vector<tree::Photon>> photons_2(reader_2, "photons");
    TTreeReaderValue<std::vector<tree::Electron>> electrons_1(reader_1, "electrons");
    TTreeReaderValue<std::vector<tree::Electron>> electrons_2(reader_2, "electrons");
    TTreeReaderValue<std::vector<tree::Jet>> jets_1(reader_1, "jets");
    TTreeReaderValue<std::vector<tree::Jet>> jets_2(reader_2, "jets");
    TTreeReaderValue<tree::MET> met_1(reader_1, "met");
    TTreeReaderValue<tree::MET> met_2(reader_2, "met");
    
    int i = 0;
    
    while (reader_2.Next()){
        while (reader_1.Next()){
            if (*evtNo_1 == *evtNo_2 && *lumNo_1 == *lumNo_2 && *runNo_1 == *runNo_2){
                if (*evtNo_1 == 972417 && *lumNo_1 == 5285 && *runNo_1 == 1){
                    //std::cout<<photons_1->at(0).pUncorrected.Pt()<<"   "<<photons_2->at(0).pUncorrected.Pt()<<std::endl;
                    std::cout<<"Photon: "<<photons_1->at(0).p.Pt()<<"   "<<photons_2->at(0).p.Pt()<<std::endl;
                    std::cout<<*runNo_1<<":"<<*lumNo_1<<":"<<*evtNo_1<<std::endl;
                    //std::cout<<"Electron: "<<electrons_1->at(0).p.Pt()<<"   "<<electrons_2->at(0).p.Pt()<<std::endl;
                    //std::cout<<"Jet: "<<jets_1->at(0).p.Pt()<<"   "<<jets_2->at(0).p.Pt()<<std::endl;
                    //std::cout<<"MET: "<<met_1->p.Pt()<<"   "<<met_2->p.Pt()<<std::endl;
                }
            }
        }
        reader_1.SetEntry(0);
        i++;
        //if (i==20) {break;}
    }
    return 0;
}

int main(){
    run();
    return 0;
}
