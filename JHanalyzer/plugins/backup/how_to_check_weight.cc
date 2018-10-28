// -*- C++ -*-
//
// Package:    Analyzer/JHanalyzer_muon_status2223
// Class:      JHanalyzer_muon_status2223
// 
/**\class JHanalyzer_muon_status2223 JHanalyzer_muon_status2223.cc Analyzer/JHanalyzer_muon_status2223/plugins/JHanalyzer_muon_status2223.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  JunHo Choi
//         Created:  Fri, 23 Mar 2018 03:24:18 GMT
//
//
//LHEEventProduct
//externalLHEProducer
//LHERunInfoProduct
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
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/PdfInfo.h"

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

class JHanalyzer_muon_status2223 : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit JHanalyzer_muon_status2223(const edm::ParameterSet&);
      ~JHanalyzer_muon_status2223();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

  //GenEventInfoProduct                   "generator"                 ""                "SIM"   

  edm::EDGetTokenT<GenParticleCollection> genParticles_Token;
  edm::EDGetTokenT<GenEventInfoProduct> genInfo_Token;
  edm::EDGetTokenT<LHEEventProduct> LHEInfo_Token;
  //  edm::EDGetTokenT<LHERunInfoProduct> extLHEInfo_Token;  
  //   extLHEInfo_Token = consumes<LHERunInfoProduct>(edm::InputTag("externalLHEProducer"));


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

  TTree *FSRinfo;

  TTree * truweight;

  TTree *status2223;
  TTree *lheinfo;

  TTree *trupdfsclweight0;
  TTree *trupdfsclweight1;
  TTree *trupdfsclweight2;
  TTree *trupdfsclweight3;
  TTree *trupdfsclweight4;
  TTree *trupdfsclweight5;

  TTree *trupdfweight0;
  TTree *trupdfweight1;
  TTree *trupdfweight2;
  TTree *trupdfweight3;
  TTree *trupdfweight4;
  TTree *trupdfweight5;
  TTree *trupdfweight6;
  TTree *trupdfweight7;
  TTree *trupdfweight8;
  TTree *trupdfweight9;
  TTree *trupdfweight10;
  TTree *trupdfweight11;
  TTree *trupdfweight12;
  TTree *trupdfweight13;
  TTree *trupdfweight14;
  TTree *trupdfweight15;
  TTree *trupdfweight16;
  TTree *trupdfweight17;
  TTree *trupdfweight18;
  TTree *trupdfweight19;
  TTree *trupdfweight20;
  TTree *trupdfweight21;
  TTree *trupdfweight22;
  TTree *trupdfweight23;
  TTree *trupdfweight24;
  TTree *trupdfweight25;
  TTree *trupdfweight26;
  TTree *trupdfweight27;
  TTree *trupdfweight28;
  TTree *trupdfweight29;
  TTree *trupdfweight30;
  TTree *trupdfweight31;
  TTree *trupdfweight32;
  TTree *trupdfweight33;
  TTree *trupdfweight34;
  TTree *trupdfweight35;
  TTree *trupdfweight36;
  TTree *trupdfweight37;
  TTree *trupdfweight38;
  TTree *trupdfweight39;
  TTree *trupdfweight40;
  TTree *trupdfweight41;
  TTree *trupdfweight42;
  TTree *trupdfweight43;
  TTree *trupdfweight44;
  TTree *trupdfweight45;
  TTree *trupdfweight46;
  TTree *trupdfweight47;
  TTree *trupdfweight48;
  TTree *trupdfweight49;
  TTree *trupdfweight50;
  TTree *trupdfweight51;
  TTree *trupdfweight52;
  TTree *trupdfweight53;
  TTree *trupdfweight54;
  TTree *trupdfweight55;
  TTree *trupdfweight56;
  TTree *trupdfweight57;
  TTree *trupdfweight58;
  TTree *trupdfweight59;
  TTree *trupdfweight60;
  TTree *trupdfweight61;
  TTree *trupdfweight62;
  TTree *trupdfweight63;
  TTree *trupdfweight64;
  TTree *trupdfweight65;
  TTree *trupdfweight66;
  TTree *trupdfweight67;
  TTree *trupdfweight68;
  TTree *trupdfweight69;
  TTree *trupdfweight70;
  TTree *trupdfweight71;
  TTree *trupdfweight72;
  TTree *trupdfweight73;
  TTree *trupdfweight74;
  TTree *trupdfweight75;
  TTree *trupdfweight76;
  TTree *trupdfweight77;
  TTree *trupdfweight78;
  TTree *trupdfweight79;
  TTree *trupdfweight80;
  TTree *trupdfweight81;
  TTree *trupdfweight82;
  TTree *trupdfweight83;
  TTree *trupdfweight84;
  TTree *trupdfweight85;
  TTree *trupdfweight86;
  TTree *trupdfweight87;
  TTree *trupdfweight88;
  TTree *trupdfweight89;
  TTree *trupdfweight90;
  TTree *trupdfweight91;
  TTree *trupdfweight92;
  TTree *trupdfweight93;
  TTree *trupdfweight94;
  TTree *trupdfweight95;
  TTree *trupdfweight96;
  TTree *trupdfweight97;
  TTree *trupdfweight98;
  TTree *trupdfweight99;
  TTree *trupdfweight100;
  TTree *trupdfweight101;



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

  //  double boson_m,boson_pt,boson_eta,boson_phi;
  double boson_m,boson_pt;
  double NFSR,Zevent;
  double lhemass,lhept;
  double trupdfweight[102];
  double trupdfsclweight[6];
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
JHanalyzer_muon_status2223::JHanalyzer_muon_status2223(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
 
  //vector<reco::GenParticle>             "genParticles"              ""                "SIM"     

  usesResource("TFileService");//genParticles
   genParticles_Token = consumes<GenParticleCollection>(edm::InputTag("genParticles"));
   genInfo_Token = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
   LHEInfo_Token = consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"));
   //  extLHEInfo_Token = consumes<LHERunInfoProduct>(edm::InputTag("externalLHEProducer"));

}


JHanalyzer_muon_status2223::~JHanalyzer_muon_status2223()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
JHanalyzer_muon_status2223::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{


    
  edm::Handle<LHERunInfoProduct> run;
  typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
  edm::EDGetTokenT<LHEEventProduct> extLHEInfo_Token;
  //  extLHEInfo_Token = consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"));
  //iRun.getByLabel( "externalLHEProducer", run );                                                                                                                  //LHEInfo_Token                                        
  iEvent.getByToken( LHEInfo_Token, run );
  LHERunInfoProduct myLHERunInfoProduct = *(run.product());

  for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
    std::cout << iter->tag() << std::endl;
    std::vector<std::string> lines = iter->lines();
    for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
      std::cout << lines.at(iLine);
    }
  }
  
  ////////////////////////////

  /*
<weightgroup combine="envelope" name="scale_variation">
<weight id="1001"> muR=1 muF=1 </weight>->out
<weight id="1002"> muR=1 muF=2 </weight> trupdfsclweight0, w0
<weight id="1003"> muR=1 muF=0.5 </weight> trupdfsclweight1, w1 -> w3
<weight id="1004"> muR=2 muF=1 </weight> trupdfsclweight2, w2
<weight id="1005"> muR=2 muF=2 </weight> trupdfsclweight3, w3 -> w1
<weight id="1006"> muR=2 muF=0.5 </weight>->out
<weight id="1007"> muR=0.5 muF=1 </weight> trupdfsclweight4, w4
<weight id="1008"> muR=0.5 muF=2 </weight>->out
<weight id="1009"> muR=0.5 muF=0.5 </weight> trupdfsclweight5, w5


then -> 1,2 -> 2,2 -> 2,1 -> 1,0.5 -> 0.5,1 -> 0.5,0.5

2001~2102

trupdfweight0, w0
.
.
.
trupdfweight101, w101

   */


    //cout<<"!!!!!!!!!!!!!!!!!start!!!!!!!!!!!!!!!!!"<<endl;
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

  Zevent=0;
  NFSR=0;
  ////////////

//////////////////////


   using namespace edm;
   Handle<reco::GenParticleCollection> genParticles;
   iEvent.getByToken(genParticles_Token, genParticles);//genParticle                                                                                         
   edm::Handle<GenEventInfoProduct> genInfo;
   iEvent.getByToken(genInfo_Token, genInfo);
   edm::Handle<LHEEventProduct> LHEInfo;
   iEvent.getByToken(LHEInfo_Token, LHEInfo);
   //GenEventInfoProduct                   "generator"                 ""                "SIM"   //
   int lheinfoweightsize= LHEInfo->weights().size();
   
   int lheinfocommentssize = LHEInfo->comments_size();
   for (int i =0; i < lheinfocommentssize; i++){
     cout<<"comment i ="<<i<<"=" << LHEInfo->getComment(i)<<endl;
   }
   for(int i=0; i<lheinfoweightsize; i++){
     // cout<<"i="<<i<<" info="<<LHEInfo->weights()[i].id<<endl;
     
   }
   //   cout<<"before LHE weight"<<endl;
   for(int i =0; i<102;i++){
       trupdfweight[i]=LHEInfo->weights()[i+9].wgt/LHEInfo->originalXWGTUP();
   }
   /*
   <weightgroup combine="envelope" name="scale_variation">                                                                                                                                        
      <weight id="1001"> muR=1 muF=1 </weight>->out                                                                                                                                                  
      <weight id="1002"> muR=1 muF=2 </weight> trupdfsclweight0, w0                                                                                                                                  
      <weight id="1003"> muR=1 muF=0.5 </weight> trupdfsclweight1, w1                                                                                                                                
      <weight id="1004"> muR=2 muF=1 </weight> trupdfsclweight2, w2                                                                                                                                  
      <weight id="1005"> muR=2 muF=2 </weight> trupdfsclweight3, w3                                                                                                                                  
      <weight id="1006"> muR=2 muF=0.5 </weight>->out                                                                                                                                                
      <weight id="1007"> muR=0.5 muF=1 </weight> trupdfsclweight4, w4                                                                                                                                
      <weight id="1008"> muR=0.5 muF=2 </weight>->out                                                                                                                                                
      <weight id="1009"> muR=0.5 muF=0.5 </weight> trupdfsclweight5, w5  


i=0 info=1001
i=1 info=1002
i=2 info=1003
i=3 info=1004
i=4 info=1005
i=5 info=1006
i=6 info=1007
i=7 info=1008
i=8 info=1009

*/
   //   for(int i=103;i<108;i++){
   // trupdfsclweight[i-103]=LHEInfo->weights()[i].wgt/LHEInfo->originalXWGTUP();
   //}
   //   cout<<"after LHE weight"<<endl;
   
   trupdfsclweight[0]=LHEInfo->weights()[1].wgt/LHEInfo->originalXWGTUP();
   trupdfsclweight[3]=LHEInfo->weights()[2].wgt/LHEInfo->originalXWGTUP();//same order with LQNtuple
   trupdfsclweight[2]=LHEInfo->weights()[3].wgt/LHEInfo->originalXWGTUP();
   trupdfsclweight[1]=LHEInfo->weights()[4].wgt/LHEInfo->originalXWGTUP();//same order with LQNtuple
   trupdfsclweight[4]=LHEInfo->weights()[6].wgt/LHEInfo->originalXWGTUP();
   trupdfsclweight[5]=LHEInfo->weights()[8].wgt/LHEInfo->originalXWGTUP();
   //cout<<"scale5="<<trupdfsclweight[5]<<endl;
   int leppid=13;
   //LHE info///
   const lhef::HEPEUP& lheEvent = LHEInfo->hepeup();
   std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;

   Int_t nLHEParticle = lheParticles.size();
   TLorentzVector v1lhe,v2lhe,vzlhe; 
   lhemass=0;
   for( Int_t idxParticle = 0; idxParticle < nLHEParticle; ++idxParticle ){
     
     Int_t id = lheEvent.IDUP[idxParticle];
     Double_t px = lheParticles[idxParticle][0];
     Double_t py = lheParticles[idxParticle][1];
     Double_t pz = lheParticles[idxParticle][2];
     Double_t ee = lheParticles[idxParticle][3];
     Double_t mm = lheParticles[idxParticle][4];
     Int_t status = lheEvent.ISTUP[idxParticle];

     //     double px1=0,py1=0,pz1=0,ee1=0;
     //double px2=0,py2=0,pz2=0,ee2=0;
 
     //          cout<<"i="<<idxParticle<<" pid="<<id<<" status="<<status<<" px="<<px<<" py="<<py<<" pz="<<pz<<" ee="<<ee<<" mm="<<mm<<endl;
     if(status==2 && (id==23 || id == 22)){
          vzlhe.SetPxPyPzE(px,py,pz,ee);
	  lhemass=mm;
	  lhept=vzlhe.Perp();
     }
     else if(status==1 && id==leppid){
       v1lhe.SetPxPyPzE(px,py,pz,ee);
     }
     else if(status==1 && id==-leppid){
       v2lhe.SetPxPyPzE(px,py,pz,ee);
     }
   }
   if(vzlhe.E()!=0){
     //cout<<"LHEZMass=="<<vzlhe.M()<<endl;
     //cout<<"LHEZMass="<<lhemass<<endl;
   }
     else if(v1lhe.E()!=0){
       lhemass=(v1lhe+v2lhe).M();
       lhept=(v1lhe+v2lhe).Perp();
       // cout<<"LHEZMass="<<lhemass<<endl;
    //       cout<<"LHEMass="<<(v1lhe+v2lhe).M()<<endl;
     }
     if(v1lhe.E()!=0){
       //       cout<<"channel matched"<<endl;
     }


   

   ///end of LHE info//
   weight=genInfo->weight();


   ///////Define Hard PID/////
   vector<int> hardpid;
   hardpid.push_back(leppid);
   hardpid.push_back(-leppid);
   hardpid.push_back(23);//Z boson
   hardpid.push_back(22);//gamma boson


   vector<int> hardindex; // to save indices of hard lepton pair && Z boson
   vector<int> hardindex_default; // to save indices of hard lepton pair && Z boson


   int lep1_i=-99;
   int lep2_i=-99;

   int lep1_i_default=-99;
   int lep2_i_default=-99;

   int nlep=0;
   int nlep_default=0;
   int nhard=0;

   int gensize= genParticles->size();
   cout<<"gensize="<<gensize<<endl;
   int hardpidsize=hardpid.size();
   
   int lep1_i_status23=-99;
   int lep2_i_status23=-99;
   int Z_i_status22=-99;
   int Z_i_hard=-99;
   //cout<<"set hardindex(default)"<<endl;
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
     int prompt = p.statusFlags().isPrompt();//                                                                                                             
     if(hardprocess) nhard+=1;
     //cout<<"i="<<i<<" id="<<id<<" status="<<status<<" prompt="<<prompt<<" hardprocess="<<hardprocess<<" fromhard="<<fromhard<<" ee="<<ee<<endl;
     if((id==23||id==22) && status==22) Z_i_status22=i;
     if(id==leppid && status==23) lep1_i_status23=i;
     if(id==-leppid && status==23) lep2_i_status23=i;

     //     if(hardprocess) nhard+=1;
     for(int j = 0; j < hardpidsize; j++){
       
       if(id==hardpid[j] && hardprocess){
	 hardindex_default.push_back(i);
       }       //end of if hardprocess && Z/gamma*
       
     }//end of hardpid
   
     ////default method/////
     ///////find   final prompt lepton
     if(fromhard && status==1){
       if(id==leppid){
	 lep1_px=px; lep1_py=py; lep1_pz=pz; lep1_ee=ee;
	 lep1_i_default=i;
	 nlep_default+=1;
       }//end of id==leppid
       else if(id==-leppid){
	 lep2_px=px; lep2_py=py; lep2_pz=pz; lep2_ee=ee;
	 lep2_i_default=i;
	 nlep_default+=1;
       }//end of id == -lepip
     }//end of fromhard
     
     
   }///end of gen loop
   
   //
   //   if(nhard==0){
   

   //cout<<"set hardindex status2223"<<endl;
     for(int i = 0; i < gensize; ++ i) {///scan all gen particles
       const GenParticle & p = (*genParticles)[i];
       int id = p.pdgId();
       int status = p.status();
       
       for(int j = 0; j < hardpidsize; j++){

	 if(id==hardpid[j] && (status==22 || status==23)){
	   hardindex.push_back(i);
	 }       //end of if hardprocess && Z/gamma*   
	 
       }//end of hardpid  


     }//genparticle loop
     
     //   }//if nhard==0

  int hardindexsize=hardindex.size();
  int hardindexsize_default=hardindex_default.size();
  


   

   ///Now we've found Z/gamma*(=llpair)
   //&& final prompt lepton
   

   
   //Mother index//
  //cout<<"set motherindex"<<endl;
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


   //cout<<"set gen_fromhardDY"<<endl;
   //  /////////let's make hardDY//////
   vector<int> gen_fromhardDY;

   for(int i = 0; i < gensize; ++ i) {
     //    const GenParticle & p = (*genParticles)[i];

     //     int id = p.pdgId();
     // int status = p.status();
     //     int fromhard=p.statusFlags().fromHardProcess();
     int fromhardDY=0;
     int mother = gen_motherindex[i];
     if(mother<0){
       gen_fromhardDY.push_back(fromhardDY);
       continue;
     } 
     while(((*genParticles)[mother]).status() != 4){//if mother is not hardprocess parton                                                                    
       for(int j = 0; j < hardindexsize; j++){//check for DY Z/gamma* indices                                                                                
         if(mother==hardindex[j]) fromhardDY=1;//if mother(ancsester) == DY Z/gamma*, i'th particles is from DY vtx.                                         
       }//end of for hardindex                                                                                                                               
       mother=gen_motherindex[mother];
       if(mother==-1) break;
     }//end of while mother   

     gen_fromhardDY.push_back(fromhardDY);
   }

   ///let's find tau -> if its from hard DY -> veto
   //and find prompt lepton pair
   //cout<<"tau veto"<<endl;
   for(int i = 0; i < gensize; ++ i) {
     const GenParticle & p = (*genParticles)[i];
     int id = p.pdgId();
	  //     int status = p.status();
     int fromhardDY=gen_fromhardDY[i];


       if(fabs(id)==15 && fromhardDY) return; //default


   }   
   //   if(nhard==0){
   //cout<<"set DY lepton pair"<<endl;

   vector<int> vlep1_i;
   vector<int> vlep2_i;
   vector<int> vlep1_i_fake;
   vector<int> vlep2_i_fake;

   for(int i = 0; i < gensize; ++ i) {
       const GenParticle & p = (*genParticles)[i];
       int status = p.status();
       int fromhardDY=gen_fromhardDY[i];
       int id = p.pdgId();
       int mother = gen_motherindex[i];
       double px = p.px();
       double py = p.py();
       double pz = p.pz();
       double ee = p.energy();

       if(fromhardDY && (id==22) ) cout<<"################i="<<i<<" id="<<id<<" status="<<status<<" mother="<<mother<<endl;
       if(fromhardDY && (id==23||id==22)) Z_i_hard=i;
       if(status==1 && fromhardDY==1){
	 if(id==leppid){
	  

	   lep1_i=i;
   	   vlep1_i.push_back(lep1_i);
	   nlep+=1;
	   lep1_px=px; lep1_py=py; lep1_pz=pz; lep1_ee=ee;
	 }
	 else if(id==-leppid){

	   lep2_i=i;
	   vlep2_i.push_back(lep2_i);
	   nlep+=1;
	   lep2_px=px; lep2_py=py; lep2_pz=pz; lep2_ee=ee; 
	 }

       }//end of status1 fromhard

     }//end of genptlloop

     //}//end of nhard=0
     //cout<<"#########################"<<endl;
     //cout<<"lep1_i="<<lep1_i<<" lep2_i="<<lep2_i<< " nlep="<<nlep<<endl;
   int vlep1size=vlep1_i.size();
   //check if lep1_i & lep2_i is FSR leptons
   //they don't have photons OR leptons as mother in common
   //   while(nlep>=2){
     for(int l=0; l<vlep1size; l++){
       lep1_i=vlep1_i[l];
       lep2_i=vlep2_i[l];
       int mother1 = gen_motherindex[lep1_i];
       int mother2 = gen_motherindex[lep2_i];
       int mother1_pid = (*genParticles)[mother1].pdgId();
       //     int mother2_pid = (*genParticles)[mother2].pdgId();
       
       //p.pdgId
       while(mother1 > 0 && mother2 > 0){
	 
	 
	 //	 if(mother1==mother2 && (mother1_pid==22 || fabs(mother1_pid)==11 || fabs(mother1_pid)==13)){//if they are from same photon or lepton
	 if(mother1==mother2 && ( fabs(mother1_pid)==11 || fabs(mother1_pid)==13)){//if they are from same photon or lepton
	   //cout<<"FSR"<<endl;
	   vlep1_i_fake.push_back(lep1_i);
	   vlep2_i_fake.push_back(lep2_i);
	   nlep +=-2;
	   break;
	   
	   
	 }
	 mother1= gen_motherindex[mother1];
	 mother2= gen_motherindex[mother2];
	 mother1_pid = (*genParticles)[mother1].pdgId();
	 //mother2_pid = (*genParticles)[mother2].pdgId();
	 
       }//end of while mother1>0,mother2>0 //scan all mother
       
       
       //int lep1_i_fake=lep1_i;
       //int lep2_i_fake=lep2_i;
       //       vlep1_i_fake.push_back(lep1_i);
       //vlep2_i_fake.push_back(lep2_i);
     }
     //   }//end of while nlep>2
      lep1_i=-99;
     lep2_i=-99;
     nlep=0;
     int fakesize=   vlep1_i_fake.size();
     

       for(int i = 0; i < gensize; ++ i) {
	 

	 bool isfake = 0;
	 const GenParticle & p = (*genParticles)[i];
	 int status = p.status();
	 int fromhardDY=gen_fromhardDY[i];
	 int id = p.pdgId();
	 double px = p.px();
	 double py = p.py();
	 double pz = p.pz();
	 double ee = p.energy();
	 for(int l = 0; l<fakesize; l++){
	   int lep1_i_fake=vlep1_i_fake[l];
	   int lep2_i_fake=vlep2_i_fake[l];
	   if(i==lep1_i_fake){
	     isfake=1;
	     break;
	   }
	   if(i==lep2_i_fake){
	     isfake=1;
	     break;
	   }
	 }
	 if(isfake) continue;
	 if(status==1 && fromhardDY==1){
	   if(id==leppid){
	     lep1_i=i;
	     nlep+=1;
	     lep1_px=px; lep1_py=py; lep1_pz=pz; lep1_ee=ee;
	   }
	   else if(id==-leppid){
	     lep2_i=i;
	     nlep+=1;
	     lep2_px=px; lep2_py=py; lep2_pz=pz; lep2_ee=ee;
	   }
	   
	 }//end of status1 fromhard                                                                                                                             
       
     }//end of genptlloop            


   //bool isdefault=0;
   //bool isnew=0;
   //   if(nlep_default>=2) isdefault=1;
   //if(nlep>=2) isnew=1;
   if(lep1_i!=lep1_i_default){
     cout<<"lep1_i="<<lep1_i<<" lep2_i="<<lep2_i<< " nlep="<<nlep<<" fakesize="<<fakesize<<endl;

     cout<<"lep1_i_default="<<lep1_i_default<<" lep2_i_default="<<lep2_i_default<< " nlep_default="<<nlep_default<<endl;


 
  for(int i = 0; i < gensize; ++ i) {
     const GenParticle & p = (*genParticles)[i];

     int status = p.status();
     int fromhardDY=gen_fromhardDY[i];
     int id = p.pdgId();
     int fromhard=p.statusFlags().fromHardProcess();
     int hardprocess = p.isHardProcess();
     // int prompt = p.statusFlags().isPrompt();//                                                                                                            
     //     int fromhardDY = gen_fromhardDY[i];
     int mother = gen_motherindex[i];
     double ee = p.energy();

     //  cout<<"i="<<i<<" id="<<id<<" status="<<status<<" fromDY="<<fromhardDY<<" mother="<<mother<<" hardprocess="<<hardprocess<<" fromhard="<<fromhard<<" ee="<<ee<<endl;    

   }
   }//end of if defaut != statusmethod
   //cout<<"#########################"<<endl;




   if(nlep<2) return;

   TLorentzVector v1,v2;
   v1.SetPxPyPzE(lep1_px,lep1_py,lep1_pz,lep1_ee);
   v2.SetPxPyPzE(lep2_px,lep2_py,lep2_pz,lep2_ee);


   int nphoton=0;
   NPHOTONALL=0;
   // //cout<<"lep1_i="<<lep1_i<<" lep2_i="<<lep2_i<<endl;
   // for(int i =0; i<hardindexsize; i++){
   //  //cout<<"DYhard="<<hardindex[i]<<endl;
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
     int prompt = p.statusFlags().isPrompt();//
     // double eta = p.eta();
     // double phi = p.phi();
      int fromhard=p.statusFlags().fromHardProcess();//

     //     int fromhardbfFSR = p.fromHardProcessBeforeFSR();
       int hardprocess = p.isHardProcess();//
     double px = p.px();
     double py = p.py();
     double pz = p.pz();
     double ee = p.energy();
     double mm = p.mass();

     int mother = gen_motherindex[i];
     int fromhardDY = gen_fromhardDY[i];
     //     if(fabs(id)==leppid || id==23)    cout<<"i="<<i<<" id="<<id<<" status="<<status<<" prompt="<<prompt<<" mother="<<mother<<" hardprocess="<<hardprocess<<" fromhard="<<fromhard<<" ee="<<ee<<" mm="<<mm<<endl;
     
     //     cout<<"i="<<i<<" id="<<id<<" status="<<status<<" prompt="<<prompt<<" mother="<<mother<<" hardprocess="<<hardprocess<<" fromhard="<<fromhard<<" ee="<<ee<<" mm="<<mm<<endl;
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

     /*
     bool fromhardDY=0;

     while(((*genParticles)[mother]).status() != 4){//if mother is not hardprocess parton
       
       for(int j = 0; j < hardindexsize; j++){//check for DY Z/gamma* indices
	 if(mother==hardindex[j]) fromhardDY=1;//if mother(ancsester) == DY Z/gamma*, i'th particles is from DY vtx.
       }//end of for hardindex
       mother=gen_motherindex[mother];
       if(mother==-1) break;
     }//end of while mother
     */
     ///////////////////////////////////////
     if(!fromhardDY){//if it is ISR
       isr_t_px=px; isr_t_py=py; isr_t_pz=pz; isr_t_ee=ee;
      
     }
     
     ////now this particle
     //status==1, !lep1 !lep2
     else{
       /////////////Final particle && not prompt leptons && from DY vtx  => FSR//////////////   
       fsr_t_px+=px;fsr_t_py+=py;fsr_t_pz+=pz;fsr_t_ee+=ee;
       NFSR+=1;
       if(id==22){
	 isfsr=1;
	 fsr_p_px+=px;fsr_p_py+=py;fsr_p_pz+=pz;fsr_p_ee+=ee;
	 fsr_pl_px+=px;fsr_pl_py+=py;fsr_pl_pz+=pz;fsr_pl_ee+=ee;
	 
       }
       else if(fabs(id)==11 || fabs(id)==13){
	 //  fsr_l_px+=px;fsr_l_py+=py;fsr_l_pz+=pz;fsr_l_ee+=ee;
	 fsr_pl_px+=px;fsr_pl_py+=py;fsr_pl_pz+=pz;fsr_pl_ee+=ee; 
       }
     }//end of fromhardDY(FSR)
     if(id==22)   photoninfo->Fill();       
   }//End of particle loop


   if(Z_i_status22 > 0 ){
     double px,py,pz,ee;
     //     px=     (*genParticles)[Z_i_status22].px();
     //py=     (*genParticles)[Z_i_status22].py();
     //pz=     (*genParticles)[Z_i_status22].pz();
     //ee=     (*genParticles)[Z_i_status22].energy();
     px=     (*genParticles)[Z_i_hard].px();
     py=     (*genParticles)[Z_i_hard].py();
     pz=     (*genParticles)[Z_i_hard].pz();
     ee=     (*genParticles)[Z_i_hard].energy();                                                                                         

     TLorentzVector v;
     v.SetPxPyPzE(px,py,pz,ee);
     boson_m=v.M();
     // boson_phi=v.Phi();
     // boson_eta=v.Eta();
     boson_pt=v.Perp();
     Zevent=1;
   }
   else{
     double px1,py1,pz1,ee1;
     double px2,py2,pz2,ee2;

     px1= (*genParticles)[lep1_i_status23].px();
     py1= (*genParticles)[lep1_i_status23].py();
     pz1= (*genParticles)[lep1_i_status23].pz();
     ee1= (*genParticles)[lep1_i_status23].energy();

     px2= (*genParticles)[lep2_i_status23].px();
     py2= (*genParticles)[lep2_i_status23].py();
     pz2= (*genParticles)[lep2_i_status23].pz();
     ee2= (*genParticles)[lep2_i_status23].energy();

     TLorentzVector v;
     v.SetPxPyPzE(px1+px2,py1+py2,pz1+pz2,ee1+ee2);
     boson_m=v.M();
     // boson_phi=v.Phi();
     // boson_eta=v.Eta();
     boson_pt=v.Perp();
     Zevent=0;
   }
   //   cout<<"Massstatus2223="<<boson_m<<endl;
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

   status2223->Fill();
   lheinfo->Fill();
   FSRinfo->Fill();
  

   trupdfsclweight0->Fill();
   trupdfsclweight1->Fill();
   trupdfsclweight2->Fill();
   trupdfsclweight3->Fill();
   trupdfsclweight4->Fill();
   trupdfsclweight5->Fill();


   trupdfweight0->Fill();
   trupdfweight1->Fill();
   trupdfweight2->Fill();
   trupdfweight3->Fill();
   trupdfweight4->Fill();
   trupdfweight5->Fill();
   trupdfweight6->Fill();
   trupdfweight7->Fill();
   trupdfweight8->Fill();
   trupdfweight9->Fill();
   trupdfweight10->Fill();
   trupdfweight11->Fill();
   trupdfweight12->Fill();
   trupdfweight13->Fill();
   trupdfweight14->Fill();
   trupdfweight15->Fill();
   trupdfweight16->Fill();
   trupdfweight17->Fill();
   trupdfweight18->Fill();
   trupdfweight19->Fill();
   trupdfweight20->Fill();
   trupdfweight21->Fill();
   trupdfweight22->Fill();
   trupdfweight23->Fill();
   trupdfweight24->Fill();
   trupdfweight25->Fill();
   trupdfweight26->Fill();
   trupdfweight27->Fill();
   trupdfweight28->Fill();
   trupdfweight29->Fill();
   trupdfweight30->Fill();
   trupdfweight31->Fill();
   trupdfweight32->Fill();
   trupdfweight33->Fill();
   trupdfweight34->Fill();
   trupdfweight35->Fill();
   trupdfweight36->Fill();
   trupdfweight37->Fill();
   trupdfweight38->Fill();
   trupdfweight39->Fill();
   trupdfweight40->Fill();
   trupdfweight41->Fill();
   trupdfweight42->Fill();
   trupdfweight43->Fill();
   trupdfweight44->Fill();
   trupdfweight45->Fill();
   trupdfweight46->Fill();
   trupdfweight47->Fill();
   trupdfweight48->Fill();
   trupdfweight49->Fill();
   trupdfweight50->Fill();
   trupdfweight51->Fill();
   trupdfweight52->Fill();
   trupdfweight53->Fill();
   trupdfweight54->Fill();
   trupdfweight55->Fill();
   trupdfweight56->Fill();
   trupdfweight57->Fill();
   trupdfweight58->Fill();
   trupdfweight59->Fill();
   trupdfweight60->Fill();
   trupdfweight61->Fill();
   trupdfweight62->Fill();
   trupdfweight63->Fill();
   trupdfweight64->Fill();
   trupdfweight65->Fill();
   trupdfweight66->Fill();
   trupdfweight67->Fill();
   trupdfweight68->Fill();
   trupdfweight69->Fill();
   trupdfweight70->Fill();
   trupdfweight71->Fill();
   trupdfweight72->Fill();
   trupdfweight73->Fill();
   trupdfweight74->Fill();
   trupdfweight75->Fill();
   trupdfweight76->Fill();
   trupdfweight77->Fill();
   trupdfweight78->Fill();
   trupdfweight79->Fill();
   trupdfweight80->Fill();
   trupdfweight81->Fill();
   trupdfweight82->Fill();
   trupdfweight83->Fill();
   trupdfweight84->Fill();
   trupdfweight85->Fill();
   trupdfweight86->Fill();
   trupdfweight87->Fill();
   trupdfweight88->Fill();
   trupdfweight89->Fill();
   trupdfweight90->Fill();
   trupdfweight91->Fill();
   trupdfweight92->Fill();
   trupdfweight93->Fill();
   trupdfweight94->Fill();
   trupdfweight95->Fill();
   trupdfweight96->Fill();
   trupdfweight97->Fill();
   trupdfweight98->Fill();
   trupdfweight99->Fill();
   trupdfweight100->Fill();
   trupdfweight101->Fill();

   //
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
JHanalyzer_muon_status2223::beginJob()
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

  status2223 = fs->make<TTree>("status2223","status2223");
  FSRinfo = fs->make<TTree>("FSRinfo","FSRinfo");
  lheinfo = fs->make<TTree>("lheinfo","lheinfo");

  trupdfsclweight0 = fs->make<TTree>("trupdfsclweight0","trupdfsclweight0");
  trupdfsclweight1 = fs->make<TTree>("trupdfsclweight1","trupdfsclweight1");
  trupdfsclweight2 = fs->make<TTree>("trupdfsclweight2","trupdfsclweight2");
  trupdfsclweight3 = fs->make<TTree>("trupdfsclweight3","trupdfsclweight3");
  trupdfsclweight4 = fs->make<TTree>("trupdfsclweight4","trupdfsclweight4");
  trupdfsclweight5 = fs->make<TTree>("trupdfsclweight5","trupdfsclweight5");
  
  trupdfweight0 = fs->make<TTree>("trupdfweight0","trupdfweight0");
  trupdfweight1 = fs->make<TTree>("trupdfweight1","trupdfweight1");
  trupdfweight2 = fs->make<TTree>("trupdfweight2","trupdfweight2");
  trupdfweight3 = fs->make<TTree>("trupdfweight3","trupdfweight3");
  trupdfweight4 = fs->make<TTree>("trupdfweight4","trupdfweight4");
  trupdfweight5 = fs->make<TTree>("trupdfweight5","trupdfweight5");
  trupdfweight6 = fs->make<TTree>("trupdfweight6","trupdfweight6");
  trupdfweight7 = fs->make<TTree>("trupdfweight7","trupdfweight7");
  trupdfweight8 = fs->make<TTree>("trupdfweight8","trupdfweight8");
  trupdfweight9 = fs->make<TTree>("trupdfweight9","trupdfweight9");
  trupdfweight10 = fs->make<TTree>("trupdfweight10","trupdfweight10");
  trupdfweight11 = fs->make<TTree>("trupdfweight11","trupdfweight11");
  trupdfweight12 = fs->make<TTree>("trupdfweight12","trupdfweight12");
  trupdfweight13 = fs->make<TTree>("trupdfweight13","trupdfweight13");
  trupdfweight14 = fs->make<TTree>("trupdfweight14","trupdfweight14");
  trupdfweight15 = fs->make<TTree>("trupdfweight15","trupdfweight15");
  trupdfweight16 = fs->make<TTree>("trupdfweight16","trupdfweight16");
  trupdfweight17 = fs->make<TTree>("trupdfweight17","trupdfweight17");
  trupdfweight18 = fs->make<TTree>("trupdfweight18","trupdfweight18");
  trupdfweight19 = fs->make<TTree>("trupdfweight19","trupdfweight19");
  trupdfweight20 = fs->make<TTree>("trupdfweight20","trupdfweight20");
  trupdfweight21 = fs->make<TTree>("trupdfweight21","trupdfweight21");
  trupdfweight22 = fs->make<TTree>("trupdfweight22","trupdfweight22");
  trupdfweight23 = fs->make<TTree>("trupdfweight23","trupdfweight23");
  trupdfweight24 = fs->make<TTree>("trupdfweight24","trupdfweight24");
  trupdfweight25 = fs->make<TTree>("trupdfweight25","trupdfweight25");
  trupdfweight26 = fs->make<TTree>("trupdfweight26","trupdfweight26");
  trupdfweight27 = fs->make<TTree>("trupdfweight27","trupdfweight27");
  trupdfweight28 = fs->make<TTree>("trupdfweight28","trupdfweight28");
  trupdfweight29 = fs->make<TTree>("trupdfweight29","trupdfweight29");
  trupdfweight30 = fs->make<TTree>("trupdfweight30","trupdfweight30");
  trupdfweight31 = fs->make<TTree>("trupdfweight31","trupdfweight31");
  trupdfweight32 = fs->make<TTree>("trupdfweight32","trupdfweight32");
  trupdfweight33 = fs->make<TTree>("trupdfweight33","trupdfweight33");
  trupdfweight34 = fs->make<TTree>("trupdfweight34","trupdfweight34");
  trupdfweight35 = fs->make<TTree>("trupdfweight35","trupdfweight35");
  trupdfweight36 = fs->make<TTree>("trupdfweight36","trupdfweight36");
  trupdfweight37 = fs->make<TTree>("trupdfweight37","trupdfweight37");
  trupdfweight38 = fs->make<TTree>("trupdfweight38","trupdfweight38");
  trupdfweight39 = fs->make<TTree>("trupdfweight39","trupdfweight39");
  trupdfweight40 = fs->make<TTree>("trupdfweight40","trupdfweight40");
  trupdfweight41 = fs->make<TTree>("trupdfweight41","trupdfweight41");
  trupdfweight42 = fs->make<TTree>("trupdfweight42","trupdfweight42");
  trupdfweight43 = fs->make<TTree>("trupdfweight43","trupdfweight43");
  trupdfweight44 = fs->make<TTree>("trupdfweight44","trupdfweight44");
  trupdfweight45 = fs->make<TTree>("trupdfweight45","trupdfweight45");
  trupdfweight46 = fs->make<TTree>("trupdfweight46","trupdfweight46");
  trupdfweight47 = fs->make<TTree>("trupdfweight47","trupdfweight47");
  trupdfweight48 = fs->make<TTree>("trupdfweight48","trupdfweight48");
  trupdfweight49 = fs->make<TTree>("trupdfweight49","trupdfweight49");
  trupdfweight50 = fs->make<TTree>("trupdfweight50","trupdfweight50");
  trupdfweight51 = fs->make<TTree>("trupdfweight51","trupdfweight51");
  trupdfweight52 = fs->make<TTree>("trupdfweight52","trupdfweight52");
  trupdfweight53 = fs->make<TTree>("trupdfweight53","trupdfweight53");
  trupdfweight54 = fs->make<TTree>("trupdfweight54","trupdfweight54");
  trupdfweight55 = fs->make<TTree>("trupdfweight55","trupdfweight55");
  trupdfweight56 = fs->make<TTree>("trupdfweight56","trupdfweight56");
  trupdfweight57 = fs->make<TTree>("trupdfweight57","trupdfweight57");
  trupdfweight58 = fs->make<TTree>("trupdfweight58","trupdfweight58");
  trupdfweight59 = fs->make<TTree>("trupdfweight59","trupdfweight59");
  trupdfweight60 = fs->make<TTree>("trupdfweight60","trupdfweight60");
  trupdfweight61 = fs->make<TTree>("trupdfweight61","trupdfweight61");
  trupdfweight62 = fs->make<TTree>("trupdfweight62","trupdfweight62");
  trupdfweight63 = fs->make<TTree>("trupdfweight63","trupdfweight63");
  trupdfweight64 = fs->make<TTree>("trupdfweight64","trupdfweight64");
  trupdfweight65 = fs->make<TTree>("trupdfweight65","trupdfweight65");
  trupdfweight66 = fs->make<TTree>("trupdfweight66","trupdfweight66");
  trupdfweight67 = fs->make<TTree>("trupdfweight67","trupdfweight67");
  trupdfweight68 = fs->make<TTree>("trupdfweight68","trupdfweight68");
  trupdfweight69 = fs->make<TTree>("trupdfweight69","trupdfweight69");
  trupdfweight70 = fs->make<TTree>("trupdfweight70","trupdfweight70");
  trupdfweight71 = fs->make<TTree>("trupdfweight71","trupdfweight71");
  trupdfweight72 = fs->make<TTree>("trupdfweight72","trupdfweight72");
  trupdfweight73 = fs->make<TTree>("trupdfweight73","trupdfweight73");
  trupdfweight74 = fs->make<TTree>("trupdfweight74","trupdfweight74");
  trupdfweight75 = fs->make<TTree>("trupdfweight75","trupdfweight75");
  trupdfweight76 = fs->make<TTree>("trupdfweight76","trupdfweight76");
  trupdfweight77 = fs->make<TTree>("trupdfweight77","trupdfweight77");
  trupdfweight78 = fs->make<TTree>("trupdfweight78","trupdfweight78");
  trupdfweight79 = fs->make<TTree>("trupdfweight79","trupdfweight79");
  trupdfweight80 = fs->make<TTree>("trupdfweight80","trupdfweight80");
  trupdfweight81 = fs->make<TTree>("trupdfweight81","trupdfweight81");
  trupdfweight82 = fs->make<TTree>("trupdfweight82","trupdfweight82");
  trupdfweight83 = fs->make<TTree>("trupdfweight83","trupdfweight83");
  trupdfweight84 = fs->make<TTree>("trupdfweight84","trupdfweight84");
  trupdfweight85 = fs->make<TTree>("trupdfweight85","trupdfweight85");
  trupdfweight86 = fs->make<TTree>("trupdfweight86","trupdfweight86");
  trupdfweight87 = fs->make<TTree>("trupdfweight87","trupdfweight87");
  trupdfweight88 = fs->make<TTree>("trupdfweight88","trupdfweight88");
  trupdfweight89 = fs->make<TTree>("trupdfweight89","trupdfweight89");
  trupdfweight90 = fs->make<TTree>("trupdfweight90","trupdfweight90");
  trupdfweight91 = fs->make<TTree>("trupdfweight91","trupdfweight91");
  trupdfweight92 = fs->make<TTree>("trupdfweight92","trupdfweight92");
  trupdfweight93 = fs->make<TTree>("trupdfweight93","trupdfweight93");
  trupdfweight94 = fs->make<TTree>("trupdfweight94","trupdfweight94");
  trupdfweight95 = fs->make<TTree>("trupdfweight95","trupdfweight95");
  trupdfweight96 = fs->make<TTree>("trupdfweight96","trupdfweight96");
  trupdfweight97 = fs->make<TTree>("trupdfweight97","trupdfweight97");
  trupdfweight98 = fs->make<TTree>("trupdfweight98","trupdfweight98");
  trupdfweight99 = fs->make<TTree>("trupdfweight99","trupdfweight99");
  trupdfweight100 = fs->make<TTree>("trupdfweight100","trupdfweight100");
  trupdfweight101 = fs->make<TTree>("trupdfweight101","trupdfweight101");




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

  status2223->Branch("boson_m",&boson_m,"boson_m/D");
  // status2223->Branch("boson_eta",&boson_eta,"boson_eta/D");
  //status2223->Branch("boson_phi",&boson_phi,"boson_phi/D");
  status2223->Branch("boson_pt",&boson_pt,"boson_pt/D");
  lheinfo->Branch("lhemass",&lhemass,"lhemass/D");
  lheinfo->Branch("lhept",&lhept,"lhept/D");

  FSRinfo->Branch("NFSR",&NFSR,"NFSR/D");
  FSRinfo->Branch("Zevent",&Zevent,"Zevent/D");




  trupdfsclweight0->Branch("w0",&trupdfsclweight[0],"w0/D");
  trupdfsclweight1->Branch("w1",&trupdfsclweight[1],"w1/D");
  trupdfsclweight2->Branch("w2",&trupdfsclweight[2],"w2/D");
  trupdfsclweight3->Branch("w3",&trupdfsclweight[3],"w3/D");
  trupdfsclweight4->Branch("w4",&trupdfsclweight[4],"w4/D");
  trupdfsclweight5->Branch("w5",&trupdfsclweight[5],"w5/D");


  trupdfweight0->Branch("w0",&trupdfweight[0],"w0/D");
  trupdfweight1->Branch("w1",&trupdfweight[1],"w1/D");
  trupdfweight2->Branch("w2",&trupdfweight[2],"w2/D");
  trupdfweight3->Branch("w3",&trupdfweight[3],"w3/D");
  trupdfweight4->Branch("w4",&trupdfweight[4],"w4/D");
  trupdfweight5->Branch("w5",&trupdfweight[5],"w5/D");
  trupdfweight6->Branch("w6",&trupdfweight[6],"w6/D");
  trupdfweight7->Branch("w7",&trupdfweight[7],"w7/D");
  trupdfweight8->Branch("w8",&trupdfweight[8],"w8/D");
  trupdfweight9->Branch("w9",&trupdfweight[9],"w9/D");
  trupdfweight10->Branch("w10",&trupdfweight[10],"w10/D");
  trupdfweight11->Branch("w11",&trupdfweight[11],"w11/D");
  trupdfweight12->Branch("w12",&trupdfweight[12],"w12/D");
  trupdfweight13->Branch("w13",&trupdfweight[13],"w13/D");
  trupdfweight14->Branch("w14",&trupdfweight[14],"w14/D");
  trupdfweight15->Branch("w15",&trupdfweight[15],"w15/D");
  trupdfweight16->Branch("w16",&trupdfweight[16],"w16/D");
  trupdfweight17->Branch("w17",&trupdfweight[17],"w17/D");
  trupdfweight18->Branch("w18",&trupdfweight[18],"w18/D");
  trupdfweight19->Branch("w19",&trupdfweight[19],"w19/D");
  trupdfweight20->Branch("w20",&trupdfweight[20],"w20/D");
  trupdfweight21->Branch("w21",&trupdfweight[21],"w21/D");
  trupdfweight22->Branch("w22",&trupdfweight[22],"w22/D");
  trupdfweight23->Branch("w23",&trupdfweight[23],"w23/D");
  trupdfweight24->Branch("w24",&trupdfweight[24],"w24/D");
  trupdfweight25->Branch("w25",&trupdfweight[25],"w25/D");
  trupdfweight26->Branch("w26",&trupdfweight[26],"w26/D");
  trupdfweight27->Branch("w27",&trupdfweight[27],"w27/D");
  trupdfweight28->Branch("w28",&trupdfweight[28],"w28/D");
  trupdfweight29->Branch("w29",&trupdfweight[29],"w29/D");
  trupdfweight30->Branch("w30",&trupdfweight[30],"w30/D");
  trupdfweight31->Branch("w31",&trupdfweight[31],"w31/D");
  trupdfweight32->Branch("w32",&trupdfweight[32],"w32/D");
  trupdfweight33->Branch("w33",&trupdfweight[33],"w33/D");
  trupdfweight34->Branch("w34",&trupdfweight[34],"w34/D");
  trupdfweight35->Branch("w35",&trupdfweight[35],"w35/D");
  trupdfweight36->Branch("w36",&trupdfweight[36],"w36/D");
  trupdfweight37->Branch("w37",&trupdfweight[37],"w37/D");
  trupdfweight38->Branch("w38",&trupdfweight[38],"w38/D");
  trupdfweight39->Branch("w39",&trupdfweight[39],"w39/D");
  trupdfweight40->Branch("w40",&trupdfweight[40],"w40/D");
  trupdfweight41->Branch("w41",&trupdfweight[41],"w41/D");
  trupdfweight42->Branch("w42",&trupdfweight[42],"w42/D");
  trupdfweight43->Branch("w43",&trupdfweight[43],"w43/D");
  trupdfweight44->Branch("w44",&trupdfweight[44],"w44/D");
  trupdfweight45->Branch("w45",&trupdfweight[45],"w45/D");
  trupdfweight46->Branch("w46",&trupdfweight[46],"w46/D");
  trupdfweight47->Branch("w47",&trupdfweight[47],"w47/D");
  trupdfweight48->Branch("w48",&trupdfweight[48],"w48/D");
  trupdfweight49->Branch("w49",&trupdfweight[49],"w49/D");
  trupdfweight50->Branch("w50",&trupdfweight[50],"w50/D");
  trupdfweight51->Branch("w51",&trupdfweight[51],"w51/D");
  trupdfweight52->Branch("w52",&trupdfweight[52],"w52/D");
  trupdfweight53->Branch("w53",&trupdfweight[53],"w53/D");
  trupdfweight54->Branch("w54",&trupdfweight[54],"w54/D");
  trupdfweight55->Branch("w55",&trupdfweight[55],"w55/D");
  trupdfweight56->Branch("w56",&trupdfweight[56],"w56/D");
  trupdfweight57->Branch("w57",&trupdfweight[57],"w57/D");
  trupdfweight58->Branch("w58",&trupdfweight[58],"w58/D");
  trupdfweight59->Branch("w59",&trupdfweight[59],"w59/D");
  trupdfweight60->Branch("w60",&trupdfweight[60],"w60/D");
  trupdfweight61->Branch("w61",&trupdfweight[61],"w61/D");
  trupdfweight62->Branch("w62",&trupdfweight[62],"w62/D");
  trupdfweight63->Branch("w63",&trupdfweight[63],"w63/D");
  trupdfweight64->Branch("w64",&trupdfweight[64],"w64/D");
  trupdfweight65->Branch("w65",&trupdfweight[65],"w65/D");
  trupdfweight66->Branch("w66",&trupdfweight[66],"w66/D");
  trupdfweight67->Branch("w67",&trupdfweight[67],"w67/D");
  trupdfweight68->Branch("w68",&trupdfweight[68],"w68/D");
  trupdfweight69->Branch("w69",&trupdfweight[69],"w69/D");
  trupdfweight70->Branch("w70",&trupdfweight[70],"w70/D");
  trupdfweight71->Branch("w71",&trupdfweight[71],"w71/D");
  trupdfweight72->Branch("w72",&trupdfweight[72],"w72/D");
  trupdfweight73->Branch("w73",&trupdfweight[73],"w73/D");
  trupdfweight74->Branch("w74",&trupdfweight[74],"w74/D");
  trupdfweight75->Branch("w75",&trupdfweight[75],"w75/D");
  trupdfweight76->Branch("w76",&trupdfweight[76],"w76/D");
  trupdfweight77->Branch("w77",&trupdfweight[77],"w77/D");
  trupdfweight78->Branch("w78",&trupdfweight[78],"w78/D");
  trupdfweight79->Branch("w79",&trupdfweight[79],"w79/D");
  trupdfweight80->Branch("w80",&trupdfweight[80],"w80/D");
  trupdfweight81->Branch("w81",&trupdfweight[81],"w81/D");
  trupdfweight82->Branch("w82",&trupdfweight[82],"w82/D");
  trupdfweight83->Branch("w83",&trupdfweight[83],"w83/D");
  trupdfweight84->Branch("w84",&trupdfweight[84],"w84/D");
  trupdfweight85->Branch("w85",&trupdfweight[85],"w85/D");
  trupdfweight86->Branch("w86",&trupdfweight[86],"w86/D");
  trupdfweight87->Branch("w87",&trupdfweight[87],"w87/D");
  trupdfweight88->Branch("w88",&trupdfweight[88],"w88/D");
  trupdfweight89->Branch("w89",&trupdfweight[89],"w89/D");
  trupdfweight90->Branch("w90",&trupdfweight[90],"w90/D");
  trupdfweight91->Branch("w91",&trupdfweight[91],"w91/D");
  trupdfweight92->Branch("w92",&trupdfweight[92],"w92/D");
  trupdfweight93->Branch("w93",&trupdfweight[93],"w93/D");
  trupdfweight94->Branch("w94",&trupdfweight[94],"w94/D");
  trupdfweight95->Branch("w95",&trupdfweight[95],"w95/D");
  trupdfweight96->Branch("w96",&trupdfweight[96],"w96/D");
  trupdfweight97->Branch("w97",&trupdfweight[97],"w97/D");
  trupdfweight98->Branch("w98",&trupdfweight[98],"w98/D");
  trupdfweight99->Branch("w99",&trupdfweight[99],"w99/D");
  trupdfweight100->Branch("w100",&trupdfweight[100],"w100/D");
  trupdfweight101->Branch("w101",&trupdfweight[101],"w101/D");
  ////////////////////////////////////////////
  /*
  edm::Handle<LHERunInfoProduct> run; 
  typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
  edm::EDGetTokenT<LHEEventProduct> extLHEInfo_Token;
  extLHEInfo_Token = consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"));
  //iRun.getByLabel( "externalLHEProducer", run );
  // iEvent.getByToken( extLHEInfo_Token, run );
  LHERunInfoProduct myLHERunInfoProduct = *(run.product());
 
  for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
    std::cout << iter->tag() << std::endl;
    std::vector<std::string> lines = iter->lines();
    for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
      std::cout << lines.at(iLine);
    }
  }
  
  */

}

// ------------ method called once each job just after ending the event loop  ------------
void 
JHanalyzer_muon_status2223::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
JHanalyzer_muon_status2223::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(JHanalyzer_muon_status2223);
