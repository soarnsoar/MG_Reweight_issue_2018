// -*- C++ -*-
//
// Package:    Analyzer/JHanalyzer
// Class:      JHanalyzer
// 
/**\class JHanalyzer JHanalyzer.cc Analyzer/JHanalyzer/plugins/JHanalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  JunHo Choi
//         Created:  Fri, 23 Mar 2018 03:24:18 GMT
//
//

//vector<reco::GenParticle>             "prunedGenParticles"        ""                "PAT"     

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"



using namespace edm;
using namespace reco;
using namespace std;

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"



#include <TTree.h>
#include <TFile.h>
#include <TLorentzVector.h>


//
// Class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class JHanalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit JHanalyzer(const edm::ParameterSet&);
      ~JHanalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

  //GenEventInfoProduct                   "generator"                 ""                "SIM"   

  edm::EDGetTokenT<GenParticleCollection> genParticles_Token;
  edm::EDGetTokenT<GenEventInfoProduct> genInfo_Token;


  TTree * Nphoton;
  TTree * Nphotonall;
  TTree * photoninfo;
 



  TTree * lep1;
  TTree * lep2;

  TTree * fsr_t;
  TTree * fsr_p;
  TTree * fsr_pl;
  TTree * isr_t;
  TTree * photon_t;

  TTree * fsr01_1;
  TTree * fsr01_2;

  TTree * fsr02_1;
  TTree * fsr02_2;

  TTree * fsr03_1;
  TTree * fsr03_2;

  TTree * fsr04_1;
  TTree * fsr04_2;


  TTree * truweight;

  //  int nphoton;
  
  double NPHOTON, NPHOTONALL;
  double DR1,DR2, isfsr;
 
  


  double weight;//, pu, pu_down, pu_up;

  double isr_t_px, isr_t_py, isr_t_pz,isr_t_ee; //

  double fsr_t_px, fsr_t_py, fsr_t_pz,fsr_t_ee; //default FSR (all ptls from Z/gamma*)
  double fsr_p_px, fsr_p_py, fsr_p_pz,fsr_p_ee; //only FSR photons
  double fsr_pl_px, fsr_pl_py, fsr_pl_pz,fsr_pl_ee;

  double fsr04_1_px, fsr04_1_py, fsr04_1_pz,fsr04_1_ee;
  double fsr04_2_px, fsr04_2_py, fsr04_2_pz,fsr04_2_ee;//For dressed lepton

  double fsr03_1_px, fsr03_1_py, fsr03_1_pz,fsr03_1_ee;
  double fsr03_2_px, fsr03_2_py, fsr03_2_pz,fsr03_2_ee;//For dressed lepton                                                                                   
  double fsr02_1_px, fsr02_1_py, fsr02_1_pz,fsr02_1_ee;
  double fsr02_2_px, fsr02_2_py, fsr02_2_pz,fsr02_2_ee;//For dressed lepton                                                                                   
  double fsr01_1_px, fsr01_1_py, fsr01_1_pz,fsr01_1_ee;
  double fsr01_2_px, fsr01_2_py, fsr01_2_pz,fsr01_2_ee;//For dressed lepton                                                                                   
  double photon_t_px, photon_t_py, photon_t_pz,photon_t_ee; //all status1 photons


  double lep1_px, lep1_py, lep1_pz,lep1_ee;
  double lep2_px, lep2_py, lep2_pz,lep2_ee;
  //  double promptphoton;

      // ----------member data ---------------------------
};
//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
JHanalyzer::JHanalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
 
  //vector<reco::GenParticle>             "genParticles"              ""                "SIM"     

  usesResource("TFileService");
   genParticles_Token = consumes<GenParticleCollection>(edm::InputTag("genParticles"));
   genInfo_Token = consumes<GenEventInfoProduct>(edm::InputTag("generator"));

}


JHanalyzer::~JHanalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
JHanalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  ////////////initialize/////////////
  weight=1;

  isr_t_px=0, isr_t_py=0, isr_t_pz=0,isr_t_ee=0;

  fsr_t_px=0, fsr_t_py=0, fsr_t_pz=0,fsr_t_ee=0; //default FSR (all ptls from Z/gamma*)                                                                        
  fsr_p_px=0, fsr_p_py=0, fsr_p_pz=0,fsr_p_ee=0; //only FSR photons                                                                                            
  fsr_pl_px=0, fsr_pl_py=0, fsr_pl_pz=0,fsr_pl_ee=0;

  fsr04_1_px=0, fsr04_1_py=0, fsr04_1_pz=0,fsr04_1_ee=0;
  fsr04_2_px=0, fsr04_2_py=0, fsr04_2_pz=0,fsr04_2_ee=0;//For dressed lepton                                                                                   

  fsr03_1_px=0, fsr03_1_py=0, fsr03_1_pz=0,fsr03_1_ee=0;
  fsr03_2_px=0, fsr03_2_py=0, fsr03_2_pz=0,fsr03_2_ee=0;//For dressed lepton 
                                                                                                                                                              
  fsr02_1_px=0, fsr02_1_py=0, fsr02_1_pz=0,fsr02_1_ee=0;
  fsr02_2_px=0, fsr02_2_py=0, fsr02_2_pz=0,fsr02_2_ee=0;//For dressed lepton 
                                                                                                                                                              
  fsr01_1_px=0, fsr01_1_py=0, fsr01_1_pz=0,fsr01_1_ee=0;
  fsr01_2_px=0, fsr01_2_py=0, fsr01_2_pz=0,fsr01_2_ee=0;//For dressed lepton 
                                                                                                                                                              
  photon_t_px=0, photon_t_py=0, photon_t_pz=0,photon_t_ee=0; //all status1 photons                                                                             

  lep1_px=0, lep1_py=0, lep1_pz=0,lep1_ee=0;
  lep2_px=0, lep2_py=0, lep2_pz=0,lep2_ee=0;

  //////////////////////////////////


   using namespace edm;
   Handle<reco::GenParticleCollection> genParticles;
   iEvent.getByToken(genParticles_Token, genParticles);//genParticle                                                                                         
   edm::Handle<GenEventInfoProduct> genInfo;
   iEvent.getByToken(genInfo_Token, genInfo);
   //GenEventInfoProduct                   "generator"                 ""                "SIM"   //

   weight=genInfo->weight();

   int leppid=13;
   ///////Define Hard PID/////
   vector<int> hardpid;
   hardpid.push_back(leppid);
   hardpid.push_back(-leppid);
   hardpid.push_back(23);//Z boson


   vector<int> hardindex; // to save indices of hard lepton pair && Z boson


   int lep1_i=-99;
   int lep2_i=-99;
   int nlep=0;

   int gensize= genParticles->size();
   int hardpidsize=hardpid.size();
   for(int i = 0; i < gensize; ++ i) {///scan all gen particles
     const GenParticle & p = (*genParticles)[i];
     int hardprocess = p.isHardProcess();
     int id = p.pdgId();
     int fromhard=p.statusFlags().fromHardProcess();
     int status = p.status();
     double px = p.px();
     double py = p.py();
     double pz = p.pz();
     double ee = p.energy();
     for(int j = 0; j < hardpidsize; j++){
       
       if(id==hardpid[j] && hardprocess){
	 hardindex.push_back(i);
       }       //end of if hardprocess && Z/gamma*
       
     }//end of hardpid
   

     ///////find   final prompt lepton
     if(fromhard && status==1){
       if(id==leppid){
	 lep1_px=px; lep1_py=py; lep1_pz=pz; lep1_ee=ee;
	 lep1_i=i;
	 nlep+=1;
       }//end of id==leppid
       else if(id==-leppid){
	 lep2_px=px; lep2_py=py; lep2_pz=pz; lep2_ee=ee;
	 lep2_i=i;
	 nlep+=1;
       }//end of id == -lepip
     }//end of fromhard
     
     
   }///end of gen loop
   
   
  int hardindexsize=hardindex.size();
   
   
   if(nlep!=2) return;
   ///Now we've found Z/gamma*(=llpair)
   //&& final prompt lepton
   
   TLorentzVector v1,v2;
   v1.SetPxPyPzE(lep1_px,lep1_py,lep1_pz,lep1_ee);
   v2.SetPxPyPzE(lep2_px,lep2_py,lep2_pz,lep2_ee);
   
   //Mother index//
   vector<int> gen_motherindex;
   for(int i = 0; i < gensize; ++ i) {
     const GenParticle & p = (*genParticles)[i];

     int mother = -1;
     
     for( reco::GenParticleCollection::const_iterator mit = genParticles->begin(); mit != genParticles->end(); ++mit ) {
       if( p.mother()==&(*mit) ) {
         mother = std::distance(genParticles->begin(),mit);
         break;
       }
     }
     gen_motherindex.push_back(mother);
     
   }//end of find motherindex
   //cout<<"#########################"<<endl;
   int nphoton=0;
   NPHOTONALL=0;
   // cout<<"lep1_i="<<lep1_i<<" lep2_i="<<lep2_i<<endl;
   // for(int i =0; i<hardindexsize; i++){
   //  cout<<"DYhard="<<hardindex[i]<<endl;
   // }
   //let's find particls from DY Z/gamma*
   for(int i = 0; i < gensize; ++ i) {
     const GenParticle & p = (*genParticles)[i];
     int id = p.pdgId();
     int status = p.status();
     if(id==22){

       NPHOTONALL+=1;
       if(status!=1) cout<<"non status1 photon"<<endl;
     }    
 //int prompt = p.statusFlags().isPrompt();//
     // double eta = p.eta();
     // double phi = p.phi();
     // int fromhard=p.statusFlags().fromHardProcess();//

     //     int fromhardbfFSR = p.fromHardProcessBeforeFSR();
     //int hardprocess = p.isHardProcess();//
     double px = p.px();
     double py = p.py();
     double pz = p.pz();
     double ee = p.energy();
     int mother = gen_motherindex[i];
     //   cout<<"i="<<i<<" id="<<id<<" status="<<status<<" prompt="<<prompt<<" mother="<<mother<<" hardprocess="<<hardprocess<<" fromhard="<<fromhard<<" ee="<<ee<<endl;
     TLorentzVector vfsr;
     vfsr.SetPxPyPzE(px,py,pz,ee);
     // double dR1 = vfsr.DeltaR(v1);
     //double dR2 = vfsr.DeltaR(v2);
     
     


     if(status!=1) continue;


     ////check ith particle's ancesteors////
     if(mother==-1) continue;
     if(i==lep1_i) continue;
     if(i==lep2_i) continue;
     /////////////Final particle && not prompt leptons//////////////
     if(id==22){ nphoton+=1;
       double dR1 = vfsr.DeltaR(v1);
       double dR2 = vfsr.DeltaR(v2);
       
       DR1=dR1;DR2=dR2;isfsr=0;
     


       

       photon_t_px=px; photon_t_py=py; photon_t_pz=pz; photon_t_ee=ee;
       if(dR1<0.4){
	 fsr04_1_px=px; 	 fsr04_1_py=py; 	 fsr04_1_pz=pz; 	 fsr04_1_ee=ee; 
	 // fsr04_2_px=px; 	 fsr04_2_py=py; 	 fsr04_2_pz=pz; 	 fsr04_2_ee=ee; 
       }
       else if(dR2<0.4){
	 //fsr04_1_px=px;          fsr04_1_py=py;          fsr04_1_pz=pz;          fsr04_1_ee=ee;
	 fsr04_2_px=px;          fsr04_2_py=py;          fsr04_2_pz=pz;          fsr04_2_ee=ee;
       }

       if(dR1<0.3){
	 fsr03_1_px=px;          fsr03_1_py=py;          fsr03_1_pz=pz;          fsr03_1_ee=ee;
       }
       else if(dR2<0.3){
         fsr03_2_px=px;          fsr03_2_py=py;          fsr03_2_pz=pz;          fsr03_2_ee=ee;
       }

       if(dR1<0.2){
         fsr02_1_px=px;          fsr02_1_py=py;          fsr02_1_pz=pz;          fsr02_1_ee=ee;
       }
       else if(dR2<0.2){
         fsr02_2_px=px;          fsr02_2_py=py;          fsr02_2_pz=pz;          fsr02_2_ee=ee;
       }

       if(dR1<0.1){
         fsr01_1_px=px;          fsr01_1_py=py;          fsr01_1_pz=pz;          fsr01_1_ee=ee;
       }
       else if(dR2<0.1){
         fsr01_2_px=px;          fsr01_2_py=py;          fsr01_2_pz=pz;          fsr01_2_ee=ee;
       }



     }


     bool fromhardDY=0;

     while(((*genParticles)[mother]).status() != 4){//if mother is not hardprocess parton
       
       for(int j = 0; j < hardindexsize; j++){//check for DY Z/gamma* indices
	 if(mother==hardindex[j]) fromhardDY=1;//if mother(ancsester) == DY Z/gamma*, i'th particles is from DY vtx.
       }//end of for hardindex
       mother=gen_motherindex[mother];
       if(mother==-1) break;
     }//end of while mother
     
     ///////////////////////////////////////
     if(!fromhardDY){//if it is ISR
       isr_t_px=px; isr_t_py=py; isr_t_pz=pz; isr_t_ee=ee;
       continue;  
     }
  
     ////now this particle


     /////////////Final particle && not prompt leptons && from DY vtx  => FSR//////////////   
     fsr_t_px+=px;fsr_t_py+=py;fsr_t_pz+=pz;fsr_t_ee+=ee;
     if(id==22){
       isfsr=1;
       fsr_p_px+=px;fsr_p_py+=py;fsr_p_pz+=pz;fsr_p_ee+=ee;
       fsr_pl_px+=px;fsr_pl_py+=py;fsr_pl_pz+=pz;fsr_pl_ee+=ee;
     }
     if(fabs(id)==11 || fabs(id)==13){
       //  fsr_l_px+=px;fsr_l_py+=py;fsr_l_pz+=pz;fsr_l_ee+=ee;
       fsr_pl_px+=px;fsr_pl_py+=py;fsr_pl_pz+=pz;fsr_pl_ee+=ee; 
    }
   }//End of particle loop

   lep1->Fill();
   lep2->Fill();

   fsr_t->Fill();
   fsr_p->Fill();
   fsr_pl->Fill();
   photon_t->Fill();
   isr_t->Fill();

   fsr01_1->Fill();
   fsr01_2->Fill();

   fsr02_1->Fill();
   fsr02_2->Fill();

   fsr03_1->Fill();
   fsr03_2->Fill();

   fsr04_1->Fill();
   fsr04_2->Fill();
   
   truweight->Fill();
   //photoninfo->Fill();
   //cout<<"#ofphoton="<<nphoton<<endl;
   NPHOTON=nphoton;
   Nphoton->Fill();
   Nphotonall->Fill();
   photoninfo->Fill();
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   //Handle<ExampleData> pIn;
   // iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   // ESHandle<SetupData> pSetup;
   //iSetup.get<SetupRecord>().get(pSetup);
#endif
}



/*
 fsr_t->Fill();
   fsr_p->Fill();
   fsr_pl->Fill();
   photon_t->Fill();
   isr_t->Fill();

   fsr01_1->Fill();
   fsr01_2->Fill();

   fsr02_1->Fill();
   fsr02_2->Fill();

   fsr03_1->Fill();
   fsr03_2->Fill();

   fsr04_1->Fill();
   fsr04_2->Fill();
 */

// ------------ method called once each job just before starting event loop  ------------
void 
JHanalyzer::beginJob()
{
  edm::Service<TFileService> fs;
 
  Nphoton = fs->make<TTree>("Nphoton","Nphoton");
  Nphotonall = fs->make<TTree>("Nphotonall","Nphotonall");
  photoninfo = fs->make<TTree>("photoninfo","photoninfo");


  lep1 = fs->make<TTree>("lep1","lep1");
  lep2 = fs->make<TTree>("lep2","lep2");

  fsr_t = fs->make<TTree>("fsr_t","fsr_t");
  fsr_p = fs->make<TTree>("fsr_p","fsr_p");
  fsr_pl = fs->make<TTree>("fsr_pl","fsr_pl");
  photon_t = fs->make<TTree>("photon_t","photon_t");
  isr_t = fs->make<TTree>("isr_t","isr_t");

  fsr01_1 = fs->make<TTree>("fsr01_1","fsr01_1");
  fsr01_2 = fs->make<TTree>("fsr01_2","fsr01_2");

  fsr02_1 = fs->make<TTree>("fsr02_1","fsr02_1");
  fsr02_2 = fs->make<TTree>("fsr02_2","fsr02_2");

  fsr03_1 = fs->make<TTree>("fsr03_1","fsr03_1");
  fsr03_2 = fs->make<TTree>("fsr03_2","fsr03_2");

  fsr04_1 = fs->make<TTree>("fsr04_1","fsr04_1");
  fsr04_2 = fs->make<TTree>("fsr04_2","fsr04_2");

  truweight = fs->make<TTree>("truweight","truweight");



  Nphoton->Branch("nphoton",&NPHOTON,"nphoton/D");
  Nphotonall->Branch("nphotonall",&NPHOTONALL,"nphotonall/D");
  photoninfo->Branch("isfsr",&isfsr, "isfsr/D");
  photoninfo->Branch("dR1",&DR1, "dR1/D");
  photoninfo->Branch("dR2",&DR2, "dR2/D");

  lep1->Branch("px",&lep1_px,"px/D");
  lep1->Branch("py",&lep1_py,"py/D");
  lep1->Branch("pz",&lep1_pz,"pz/D");
  lep1->Branch("ee",&lep1_ee,"ee/D");
 
  lep2->Branch("px",&lep2_px,"px/D");
  lep2->Branch("py",&lep2_py,"py/D");
  lep2->Branch("pz",&lep2_pz,"pz/D");
  lep2->Branch("ee",&lep2_ee,"ee/D");

  fsr_t->Branch("px",&fsr_t_px,"px/D");
  fsr_t->Branch("py",&fsr_t_py,"py/D");
  fsr_t->Branch("pz",&fsr_t_pz,"pz/D");
  fsr_t->Branch("ee",&fsr_t_ee,"ee/D");

  fsr_p->Branch("px",&fsr_p_px,"px/D");
  fsr_p->Branch("py",&fsr_p_py,"py/D");
  fsr_p->Branch("pz",&fsr_p_pz,"pz/D");
  fsr_p->Branch("ee",&fsr_p_ee,"ee/D");

  fsr_pl->Branch("px",&fsr_pl_px,"px/D");
  fsr_pl->Branch("py",&fsr_pl_py,"py/D");
  fsr_pl->Branch("pz",&fsr_pl_pz,"pz/D");
  fsr_pl->Branch("ee",&fsr_pl_ee,"ee/D");

  photon_t->Branch("px",&photon_t_px,"px/D");
  photon_t->Branch("py",&photon_t_py,"py/D");
  photon_t->Branch("pz",&photon_t_pz,"pz/D");
  photon_t->Branch("ee",&photon_t_ee,"ee/D");

  isr_t->Branch("px",&isr_t_px,"px/D");
  isr_t->Branch("py",&isr_t_py,"py/D");
  isr_t->Branch("pz",&isr_t_pz,"pz/D");
  isr_t->Branch("ee",&isr_t_ee,"ee/D");


  fsr01_1->Branch("px",&fsr01_1_px,"px/D");
  fsr01_1->Branch("py",&fsr01_1_py,"py/D");
  fsr01_1->Branch("pz",&fsr01_1_pz,"pz/D");
  fsr01_1->Branch("ee",&fsr01_1_ee,"ee/D");

  fsr01_2->Branch("px",&fsr01_2_px,"px/D");
  fsr01_2->Branch("py",&fsr01_2_py,"py/D");
  fsr01_2->Branch("pz",&fsr01_2_pz,"pz/D");
  fsr01_2->Branch("ee",&fsr01_2_ee,"ee/D");


  fsr02_1->Branch("px",&fsr02_1_px,"px/D");
  fsr02_1->Branch("py",&fsr02_1_py,"py/D");
  fsr02_1->Branch("pz",&fsr02_1_pz,"pz/D");
  fsr02_1->Branch("ee",&fsr02_1_ee,"ee/D");

  fsr02_2->Branch("px",&fsr02_2_px,"px/D");
  fsr02_2->Branch("py",&fsr02_2_py,"py/D");
  fsr02_2->Branch("pz",&fsr02_2_pz,"pz/D");
  fsr02_2->Branch("ee",&fsr02_2_ee,"ee/D");

  fsr03_1->Branch("px",&fsr03_1_px,"px/D");
  fsr03_1->Branch("py",&fsr03_1_py,"py/D");
  fsr03_1->Branch("pz",&fsr03_1_pz,"pz/D");
  fsr03_1->Branch("ee",&fsr03_1_ee,"ee/D");

  fsr03_2->Branch("px",&fsr03_2_px,"px/D");
  fsr03_2->Branch("py",&fsr03_2_py,"py/D");
  fsr03_2->Branch("pz",&fsr03_2_pz,"pz/D");
  fsr03_2->Branch("ee",&fsr03_2_ee,"ee/D");

  fsr04_1->Branch("px",&fsr04_1_px,"px/D");
  fsr04_1->Branch("py",&fsr04_1_py,"py/D");
  fsr04_1->Branch("pz",&fsr04_1_pz,"pz/D");
  fsr04_1->Branch("ee",&fsr04_1_ee,"ee/D");

  fsr04_2->Branch("px",&fsr04_2_px,"px/D");
  fsr04_2->Branch("py",&fsr04_2_py,"py/D");
  fsr04_2->Branch("pz",&fsr04_2_pz,"pz/D");
  fsr04_2->Branch("ee",&fsr04_2_ee,"ee/D");

  truweight->Branch("weight",&weight,"weight/D");

}

// ------------ method called once each job just after ending the event loop  ------------
void 
JHanalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JHanalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(JHanalyzer);
