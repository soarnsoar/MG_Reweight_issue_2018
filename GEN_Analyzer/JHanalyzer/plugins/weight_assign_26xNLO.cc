// -*- C++ -*-
//
// Package:    Analyzer/weight_assign_26xNLO
// Class:      weight_assign_26xNLO
// 
/**\class weight_assign_26xNLO weight_assign_26xNLO.cc Analyzer/weight_assign_26xNLO/plugins/weight_assign_26xNLO.cc

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
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"


using namespace edm;
using namespace reco;
using namespace std;

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Run.h"//to use edm::Run


#include <TTree.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include "GEN_Analyzer/JHanalyzer/interface/weightinfo.h"

//
// Class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.



//class weight_assign_26xNLO : public edm::one::EDAnalyzer<edm::one::SharedResources>  {

class weight_assign_26xNLO : public edm::one::EDAnalyzer<edm::one::WatchRuns,edm::one::SharedResources>  {
//class weight_assign_26xNLO : public edm::EDAnalyzer  {
   public:
      explicit weight_assign_26xNLO(const edm::ParameterSet&);
      ~weight_assign_26xNLO();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
 
  //  virtual void beginJob(edm::Run const& iEvent, edm::EventSetup const &) override;
  //  void beginRun(edm::Run const&, edm::EventSetup const&) override;//to get LHERunInfoProduct//add new method
  //void endRun(edm::Run const&, edm::EventSetup const&) override;
private:
  virtual void beginJob() override;
  //  virtual void beginJob(edm::Run const& iRun) override;
  //  virtual void beginJob(edm::Run const& iRun, edm::EventSetup const &iSetup) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  //  virtual void doBeginRun_(edm::Run const&, edm::EventSetup const&) override;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginRun(edm::Run const& iEvent, edm::EventSetup const&) override;
  //virtual void beginRun(edm::Run const& iRun, edm::EventSetup const &iSetup ) override;//to get LHERunInfoProduct//add new method
  //virtual void beginRun() override;//to get LHERunInfoProduct//add new method
  //GenEventInfoProduct                   "generator"                 ""                "SIM"   

//  void beginRun(edm::Run const& iEvent, edm::EventSetup const&) override;

  edm::EDGetTokenT<GenParticleCollection> genParticles_Token;
  edm::EDGetTokenT<GenEventInfoProduct> genInfo_Token;
  edm::EDGetTokenT<LHEEventProduct> LHEInfo_Token;
  //  edm::EDGetTokenT<LHERunInfoProduct> extLHEInfo_Token; 
  void make_weight_infos();
  

  std::vector<weightinfo> NLO26xInfo;



  TTree *DYkinematics;


  double Z_pt,Z_mass,Z_eta,Z_phi;
  double lep1_pt,lep1_eta,lep1_phi;
  double lep2_pt,lep2_eta,lep2_phi;
  double dilep_pt,dilep_mass,dilep_eta,dilep_phi;


  TTree *DYweights;
  double weight;
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
weight_assign_26xNLO::weight_assign_26xNLO(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
 
  //vector<reco::GenParticle>             "genParticles"              ""                "SIM"     

  usesResource("TFileService");
  genParticles_Token = consumes<GenParticleCollection>(edm::InputTag("genParticles"));
  genInfo_Token = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  LHEInfo_Token = consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"));
  //extLHEInfo_Token = consumes<LHERunInfoProduct>(edm::InputTag("externalLHEProducer"));  
  //  extLHEInfo_Token= consumes<LHERunInfoProduct,edm::InRun>(edm::InputTag("externalLHEProducer",""));
 
 //Let's make vector for find weights!//
  make_weight_infos();
 

}


weight_assign_26xNLO::~weight_assign_26xNLO()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
weight_assign_26xNLO::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
 
  Z_pt=0;Z_mass=0;Z_eta=0;Z_phi=0;
  lep1_pt=0;lep1_eta=0;lep1_phi=0;
  lep2_pt=0;lep2_eta=0;lep2_phi=0;
  dilep_pt=0;dilep_mass=0;dilep_eta=0;dilep_phi=0;


  edm::Handle<LHEEventProduct> LHEInfo;
  iEvent.getByToken(LHEInfo_Token, LHEInfo);

  //veto tau events//                             
  const lhef::HEPEUP& lheEvent = LHEInfo->hepeup();
  std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;
  Int_t nLHEParticle = lheParticles.size();
  for( Int_t idxParticle = 0; idxParticle < nLHEParticle; ++idxParticle ){

    Int_t id = lheEvent.IDUP[idxParticle];
    if(fabs(id)==15) return;
  }
  //////////end of veto tau/////

  ////////////initialize/////////////
  //Get weight//


  int lheinfoweightsize= LHEInfo->weights().size();
  int lheinfocommentssize = LHEInfo->comments_size();

  double w0=LHEInfo->originalXWGTUP();
  int NLO26xInfosize=NLO26xInfo.size();
  //  for (int i =0; i < lheinfoweightsize; i++){
    //cout<< LHEInfo->weights()[i].id<<endl;
    //    cout<< LHEInfo->weights()[i].wgt/w0<<endl;
  //}

  for (int i26x =0; i26x < NLO26xInfosize;i26x++){
    for (int i_lhe =0; i_lhe < lheinfoweightsize; i_lhe++){ 
      if(LHEInfo->weights()[i_lhe].id==NLO26xInfo[i26x].id()){
	NLO26xInfo[i26x].set_weight(LHEInfo->weights()[i_lhe].wgt/w0);
      }
    }
  }

  //for (int i26x =0; i26x < NLO26xInfosize;i26x++){
  //  for (int i26x =0; i26x < 9;i26x++){
  //     cout<<NLO26xInfo[i26x].id()<<"=>"<<    NLO26xInfo[i26x].weight()<<endl;
    //   cout<<    NLO26xInfo[i26x].id()<<endl;
  //}
  //  for (int i =0; i < lheinfocommentssize; i++){
    //    cout<<"comment i ="<<i<<"=" << LHEInfo->getComment(i)<<endl;
  //} 


  //weight id info

  //  edm::Handle<LHERunInfoProduct> run;
 
  //  iEvent.getByToken( extLHEInfo_Token, run );
  //iEvent.getByLabel( "externalLHEProducer", run );
  // typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
  //  LHERunInfoProduct myLHERunInfoProduct = *(run.product());



  int leppid=11;

   using namespace edm;
   Handle<reco::GenParticleCollection> genParticles;
   iEvent.getByToken(genParticles_Token, genParticles);//genParticle                                                                                         
   edm::Handle<GenEventInfoProduct> genInfo;
   iEvent.getByToken(genInfo_Token, genInfo);
   //GenEventInfoProduct                   "generator"                 ""                "SIM"   //
   weight=1;
   weight=genInfo->weight();
   //   if(weight<-99) cout<<weight<<endl;
   //   cout<<"weight="<<weight<<endl;
   vector<int> i_leptons1;
   vector<int> i_leptons2;
   vector<int> i_photons;

   int gensize= genParticles->size();
   
   for(int i = 0; i < gensize; ++ i) {///scan all gen particles
     //tau veto   
     const GenParticle & p = (*genParticles)[i];
     int id = p.pdgId();
     int status = p.status();
     if(status!=1) continue;
     //double px = p.px();
     ///double py = p.py();
     //double pz = p.pz();
     //double ee = p.energy();
     if(id==leppid) i_leptons1.push_back(i);
     else if(id== -leppid) i_leptons2.push_back(i);
     else if (id==22)  i_photons.push_back(i);
   }
   if(i_leptons1.size()<1) return;//require exact 2 leptons
   if(i_leptons2.size()<1) return;//require exact 2 leptons
   int i_lep1=i_leptons1[0];
   int i_lep2=i_leptons2[0];
   //Set 1st lep1 and lep2 as DY leptons(simply..)
   //The same as DY validation code :https://github.com/cms-sw/cmssw/blob/02d4198c0b6615287fd88e9a8ff650aea994412e/Validation/EventGenerator/plugins/DrellYanValidation.cc

   TLorentzVector v1,v2;
   v1.SetPxPyPzE((*genParticles)[i_lep1].px(), (*genParticles)[i_lep1].py(),(*genParticles)[i_lep1].pz(),(*genParticles)[i_lep1].energy());
   v2.SetPxPyPzE((*genParticles)[i_lep2].px(), (*genParticles)[i_lep2].py(),(*genParticles)[i_lep2].pz(),(*genParticles)[i_lep2].energy());

   //   if((v1+v2).M()< 60 ) return;
   
   TLorentzVector vdilep=v1+v2;
   TLorentzVector vZ=v1+v2;
   //vector<int> i_fsr;
   //photons for dressed lepton//
   int photonsize=i_photons.size();
   for(int i =0 ; i < photonsize; i++){
     const GenParticle & p = (*genParticles)[i];
     double px = p.px();
     double py = p.py();
     double pz = p.pz();                                                                                                                     
     double ee = p.energy();                                                                                                                 
     TLorentzVector vfsr;
     vfsr.SetPxPyPzE(px,py,pz,ee);
     double R1=vfsr.DeltaR(v1);
     double R2=vfsr.DeltaR(v2);
     if(R1<0.1 || R2<0.1){ //i_fsr.push_back(i);
       vZ+=vfsr;
     }
   }
   

   Z_pt=vZ.Perp(); Z_mass=vZ.M(); Z_eta=vZ.Eta(); Z_phi=vZ.Phi();
   lep1_pt=v1.Perp();lep1_eta=v1.Eta();lep1_phi=v1.Phi();
   lep2_pt=v2.Perp();lep2_eta=v2.Eta();lep2_phi=v2.Phi();
   dilep_pt=vdilep.Perp(); dilep_mass=vdilep.M(); dilep_eta=vdilep.Eta(); dilep_phi=vdilep.Phi();
   
   DYkinematics->Fill();
   DYweights->Fill();
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   //Handle<ExampleData> pIn;
   // iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   // ESHandle<SetupData> pSetup;
   //iSetup.get<SetupRecord>().get(pSetup);
#endif
}





// ------------ method called once each job just before starting event loop  ------------

void 
weight_assign_26xNLO::beginJob()
//weight_assign_26xNLO::beginJob(edm::Run const& iRun)
//weight_assign_26xNLO::beginJob(edm::Run const& iRun, edm::EventSetup const &iSetup)
{
  cout<<"begin job"<<endl;
 
  edm::Service<TFileService> fs;
  DYkinematics=fs->make<TTree>("DYkinematics","DYkinematics");
  DYkinematics->Branch("Z_pt",&Z_pt,"Z_pt/D");
  DYkinematics->Branch("Z_mass",&Z_mass,"Z_mass/D");
  DYkinematics->Branch("Z_eta",&Z_eta,"Z_eta/D");
  DYkinematics->Branch("Z_phi",&Z_phi,"Z_phi/D");

  DYkinematics->Branch("lep1_pt",&lep1_pt,"lep1_pt/D");
  DYkinematics->Branch("lep1_eta",&lep1_eta,"lep1_eta/D");
  DYkinematics->Branch("lep1_phi",&lep1_phi,"lep1_phi/D");

  DYkinematics->Branch("lep2_pt",&lep2_pt,"lep2_pt/D");
  DYkinematics->Branch("lep2_eta",&lep2_eta,"lep2_eta/D");
  DYkinematics->Branch("lep2_phi",&lep2_phi,"lep2_phi/D");

  DYkinematics->Branch("dilep_pt",&dilep_pt,"dilep_pt/D");
  DYkinematics->Branch("dilep_mass",&dilep_mass,"dilep_mass/D");
  DYkinematics->Branch("dilep_eta",&dilep_eta,"dilep_eta/D");
  DYkinematics->Branch("dilep_phi",&dilep_phi,"dilep_phi/D");

  DYweights=fs->make<TTree>("DYweights","DYweights");
  DYweights->Branch("MCweight",&weight,"MCweight/D");
  for(unsigned int i = 0;i <NLO26xInfo.size();i++){ 
  // DYweights->Branch();
    DYweights->Branch(NLO26xInfo[i].name()+"_"+NLO26xInfo[i].pdf(),&NLO26xInfo[i]._weight,NLO26xInfo[i].name()+"_"+NLO26xInfo[i].pdf()+"/D");
  }
  cout<<"end of beginjob"<<endl;
}

// ------------ method called once each job just after ending the event loop  ------------

void 
weight_assign_26xNLO::endJob() 
{
  cout<<"endjob"<<endl;
}



void 
weight_assign_26xNLO::beginRun(const Run &iEvent, EventSetup const &iSetup ){
  cout<<" beginrun"<<endl;
  cout<<"end of beginrun"<<endl;
}



void                                                                                                                                                                              
weight_assign_26xNLO::endRun(edm::Run const& iEvent, edm::EventSetup const&) 
{
  cout<<"doendrun"<<endl;

}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------


void
weight_assign_26xNLO::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//////defined by jhchoi/////
void 
weight_assign_26xNLO::make_weight_infos(){
  cout<<"make_weight_infos"<<endl;

///LO26x///
NLO26xInfo.push_back(weightinfo("1001", "306000", "scale_variation_muR_1.0_muF_1.0_306000", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1002", "306000", "scale_variation_muR_2.0_muF_1.0_306000", 2.0, 1.0));
NLO26xInfo.push_back(weightinfo("1003", "306000", "scale_variation_muR_0.5_muF_1.0_306000", 0.5, 1.0));
NLO26xInfo.push_back(weightinfo("1004", "306000", "scale_variation_muR_1.0_muF_2.0_306000", 1.0, 2.0));
NLO26xInfo.push_back(weightinfo("1005", "306000", "scale_variation_muR_2.0_muF_2.0_306000", 2.0, 2.0));
NLO26xInfo.push_back(weightinfo("1006", "306000", "scale_variation_muR_0.5_muF_2.0_306000", 0.5, 2.0));
NLO26xInfo.push_back(weightinfo("1007", "306000", "scale_variation_muR_1.0_muF_0.5_306000", 1.0, 0.5));
NLO26xInfo.push_back(weightinfo("1008", "306000", "scale_variation_muR_2.0_muF_0.5_306000", 2.0, 0.5));
NLO26xInfo.push_back(weightinfo("1009", "306000", "scale_variation_muR_0.5_muF_0.5_306000", 0.5, 0.5));
NLO26xInfo.push_back(weightinfo("1010", "306000", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306000", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1011", "306001", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306001", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1012", "306002", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306002", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1013", "306003", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306003", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1014", "306004", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306004", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1015", "306005", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306005", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1016", "306006", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306006", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1017", "306007", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306007", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1018", "306008", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306008", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1019", "306009", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306009", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1020", "306010", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306010", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1021", "306011", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306011", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1022", "306012", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306012", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1023", "306013", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306013", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1024", "306014", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306014", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1025", "306015", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306015", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1026", "306016", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306016", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1027", "306017", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306017", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1028", "306018", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306018", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1029", "306019", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306019", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1030", "306020", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306020", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1031", "306021", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306021", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1032", "306022", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306022", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1033", "306023", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306023", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1034", "306024", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306024", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1035", "306025", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306025", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1036", "306026", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306026", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1037", "306027", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306027", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1038", "306028", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306028", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1039", "306029", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306029", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1040", "306030", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306030", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1041", "306031", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306031", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1042", "306032", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306032", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1043", "306033", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306033", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1044", "306034", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306034", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1045", "306035", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306035", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1046", "306036", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306036", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1047", "306037", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306037", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1048", "306038", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306038", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1049", "306039", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306039", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1050", "306040", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306040", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1051", "306041", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306041", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1052", "306042", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306042", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1053", "306043", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306043", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1054", "306044", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306044", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1055", "306045", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306045", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1056", "306046", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306046", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1057", "306047", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306047", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1058", "306048", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306048", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1059", "306049", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306049", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1060", "306050", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306050", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1061", "306051", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306051", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1062", "306052", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306052", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1063", "306053", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306053", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1064", "306054", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306054", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1065", "306055", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306055", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1066", "306056", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306056", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1067", "306057", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306057", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1068", "306058", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306058", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1069", "306059", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306059", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1070", "306060", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306060", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1071", "306061", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306061", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1072", "306062", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306062", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1073", "306063", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306063", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1074", "306064", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306064", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1075", "306065", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306065", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1076", "306066", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306066", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1077", "306067", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306067", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1078", "306068", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306068", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1079", "306069", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306069", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1080", "306070", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306070", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1081", "306071", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306071", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1082", "306072", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306072", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1083", "306073", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306073", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1084", "306074", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306074", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1085", "306075", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306075", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1086", "306076", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306076", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1087", "306077", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306077", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1088", "306078", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306078", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1089", "306079", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306079", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1090", "306080", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306080", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1091", "306081", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306081", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1092", "306082", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306082", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1093", "306083", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306083", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1094", "306084", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306084", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1095", "306085", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306085", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1096", "306086", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306086", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1097", "306087", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306087", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1098", "306088", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306088", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1099", "306089", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306089", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1100", "306090", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306090", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1101", "306091", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306091", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1102", "306092", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306092", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1103", "306093", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306093", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1104", "306094", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306094", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1105", "306095", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306095", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1106", "306096", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306096", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1107", "306097", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306097", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1108", "306098", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306098", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1109", "306099", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306099", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1110", "306100", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306100", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1111", "306101", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306101", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1112", "306102", "PDF_variationNNPDF31_nnlo_hessian_pdfas_306102", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1113", "322500", "PDF_variation_322500", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1114", "322700", "PDF_variation_322700", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1115", "322900", "PDF_variation_322900", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1116", "323100", "PDF_variation_323100", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1117", "323300", "PDF_variation_323300", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1118", "323500", "PDF_variation_323500", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1119", "323700", "PDF_variation_323700", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1120", "323900", "PDF_variation_323900", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1121", "305800", "PDF_variationNNPDF31_nlo_hessian_pdfas_305800", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1122", "305801", "PDF_variationNNPDF31_nlo_hessian_pdfas_305801", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1123", "305802", "PDF_variationNNPDF31_nlo_hessian_pdfas_305802", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1124", "305803", "PDF_variationNNPDF31_nlo_hessian_pdfas_305803", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1125", "305804", "PDF_variationNNPDF31_nlo_hessian_pdfas_305804", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1126", "305805", "PDF_variationNNPDF31_nlo_hessian_pdfas_305805", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1127", "305806", "PDF_variationNNPDF31_nlo_hessian_pdfas_305806", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1128", "305807", "PDF_variationNNPDF31_nlo_hessian_pdfas_305807", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1129", "305808", "PDF_variationNNPDF31_nlo_hessian_pdfas_305808", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1130", "305809", "PDF_variationNNPDF31_nlo_hessian_pdfas_305809", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1131", "305810", "PDF_variationNNPDF31_nlo_hessian_pdfas_305810", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1132", "305811", "PDF_variationNNPDF31_nlo_hessian_pdfas_305811", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1133", "305812", "PDF_variationNNPDF31_nlo_hessian_pdfas_305812", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1134", "305813", "PDF_variationNNPDF31_nlo_hessian_pdfas_305813", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1135", "305814", "PDF_variationNNPDF31_nlo_hessian_pdfas_305814", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1136", "305815", "PDF_variationNNPDF31_nlo_hessian_pdfas_305815", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1137", "305816", "PDF_variationNNPDF31_nlo_hessian_pdfas_305816", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1138", "305817", "PDF_variationNNPDF31_nlo_hessian_pdfas_305817", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1139", "305818", "PDF_variationNNPDF31_nlo_hessian_pdfas_305818", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1140", "305819", "PDF_variationNNPDF31_nlo_hessian_pdfas_305819", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1141", "305820", "PDF_variationNNPDF31_nlo_hessian_pdfas_305820", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1142", "305821", "PDF_variationNNPDF31_nlo_hessian_pdfas_305821", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1143", "305822", "PDF_variationNNPDF31_nlo_hessian_pdfas_305822", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1144", "305823", "PDF_variationNNPDF31_nlo_hessian_pdfas_305823", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1145", "305824", "PDF_variationNNPDF31_nlo_hessian_pdfas_305824", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1146", "305825", "PDF_variationNNPDF31_nlo_hessian_pdfas_305825", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1147", "305826", "PDF_variationNNPDF31_nlo_hessian_pdfas_305826", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1148", "305827", "PDF_variationNNPDF31_nlo_hessian_pdfas_305827", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1149", "305828", "PDF_variationNNPDF31_nlo_hessian_pdfas_305828", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1150", "305829", "PDF_variationNNPDF31_nlo_hessian_pdfas_305829", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1151", "305830", "PDF_variationNNPDF31_nlo_hessian_pdfas_305830", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1152", "305831", "PDF_variationNNPDF31_nlo_hessian_pdfas_305831", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1153", "305832", "PDF_variationNNPDF31_nlo_hessian_pdfas_305832", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1154", "305833", "PDF_variationNNPDF31_nlo_hessian_pdfas_305833", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1155", "305834", "PDF_variationNNPDF31_nlo_hessian_pdfas_305834", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1156", "305835", "PDF_variationNNPDF31_nlo_hessian_pdfas_305835", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1157", "305836", "PDF_variationNNPDF31_nlo_hessian_pdfas_305836", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1158", "305837", "PDF_variationNNPDF31_nlo_hessian_pdfas_305837", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1159", "305838", "PDF_variationNNPDF31_nlo_hessian_pdfas_305838", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1160", "305839", "PDF_variationNNPDF31_nlo_hessian_pdfas_305839", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1161", "305840", "PDF_variationNNPDF31_nlo_hessian_pdfas_305840", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1162", "305841", "PDF_variationNNPDF31_nlo_hessian_pdfas_305841", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1163", "305842", "PDF_variationNNPDF31_nlo_hessian_pdfas_305842", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1164", "305843", "PDF_variationNNPDF31_nlo_hessian_pdfas_305843", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1165", "305844", "PDF_variationNNPDF31_nlo_hessian_pdfas_305844", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1166", "305845", "PDF_variationNNPDF31_nlo_hessian_pdfas_305845", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1167", "305846", "PDF_variationNNPDF31_nlo_hessian_pdfas_305846", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1168", "305847", "PDF_variationNNPDF31_nlo_hessian_pdfas_305847", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1169", "305848", "PDF_variationNNPDF31_nlo_hessian_pdfas_305848", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1170", "305849", "PDF_variationNNPDF31_nlo_hessian_pdfas_305849", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1171", "305850", "PDF_variationNNPDF31_nlo_hessian_pdfas_305850", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1172", "305851", "PDF_variationNNPDF31_nlo_hessian_pdfas_305851", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1173", "305852", "PDF_variationNNPDF31_nlo_hessian_pdfas_305852", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1174", "305853", "PDF_variationNNPDF31_nlo_hessian_pdfas_305853", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1175", "305854", "PDF_variationNNPDF31_nlo_hessian_pdfas_305854", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1176", "305855", "PDF_variationNNPDF31_nlo_hessian_pdfas_305855", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1177", "305856", "PDF_variationNNPDF31_nlo_hessian_pdfas_305856", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1178", "305857", "PDF_variationNNPDF31_nlo_hessian_pdfas_305857", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1179", "305858", "PDF_variationNNPDF31_nlo_hessian_pdfas_305858", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1180", "305859", "PDF_variationNNPDF31_nlo_hessian_pdfas_305859", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1181", "305860", "PDF_variationNNPDF31_nlo_hessian_pdfas_305860", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1182", "305861", "PDF_variationNNPDF31_nlo_hessian_pdfas_305861", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1183", "305862", "PDF_variationNNPDF31_nlo_hessian_pdfas_305862", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1184", "305863", "PDF_variationNNPDF31_nlo_hessian_pdfas_305863", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1185", "305864", "PDF_variationNNPDF31_nlo_hessian_pdfas_305864", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1186", "305865", "PDF_variationNNPDF31_nlo_hessian_pdfas_305865", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1187", "305866", "PDF_variationNNPDF31_nlo_hessian_pdfas_305866", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1188", "305867", "PDF_variationNNPDF31_nlo_hessian_pdfas_305867", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1189", "305868", "PDF_variationNNPDF31_nlo_hessian_pdfas_305868", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1190", "305869", "PDF_variationNNPDF31_nlo_hessian_pdfas_305869", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1191", "305870", "PDF_variationNNPDF31_nlo_hessian_pdfas_305870", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1192", "305871", "PDF_variationNNPDF31_nlo_hessian_pdfas_305871", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1193", "305872", "PDF_variationNNPDF31_nlo_hessian_pdfas_305872", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1194", "305873", "PDF_variationNNPDF31_nlo_hessian_pdfas_305873", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1195", "305874", "PDF_variationNNPDF31_nlo_hessian_pdfas_305874", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1196", "305875", "PDF_variationNNPDF31_nlo_hessian_pdfas_305875", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1197", "305876", "PDF_variationNNPDF31_nlo_hessian_pdfas_305876", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1198", "305877", "PDF_variationNNPDF31_nlo_hessian_pdfas_305877", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1199", "305878", "PDF_variationNNPDF31_nlo_hessian_pdfas_305878", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1200", "305879", "PDF_variationNNPDF31_nlo_hessian_pdfas_305879", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1201", "305880", "PDF_variationNNPDF31_nlo_hessian_pdfas_305880", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1202", "305881", "PDF_variationNNPDF31_nlo_hessian_pdfas_305881", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1203", "305882", "PDF_variationNNPDF31_nlo_hessian_pdfas_305882", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1204", "305883", "PDF_variationNNPDF31_nlo_hessian_pdfas_305883", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1205", "305884", "PDF_variationNNPDF31_nlo_hessian_pdfas_305884", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1206", "305885", "PDF_variationNNPDF31_nlo_hessian_pdfas_305885", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1207", "305886", "PDF_variationNNPDF31_nlo_hessian_pdfas_305886", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1208", "305887", "PDF_variationNNPDF31_nlo_hessian_pdfas_305887", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1209", "305888", "PDF_variationNNPDF31_nlo_hessian_pdfas_305888", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1210", "305889", "PDF_variationNNPDF31_nlo_hessian_pdfas_305889", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1211", "305890", "PDF_variationNNPDF31_nlo_hessian_pdfas_305890", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1212", "305891", "PDF_variationNNPDF31_nlo_hessian_pdfas_305891", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1213", "305892", "PDF_variationNNPDF31_nlo_hessian_pdfas_305892", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1214", "305893", "PDF_variationNNPDF31_nlo_hessian_pdfas_305893", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1215", "305894", "PDF_variationNNPDF31_nlo_hessian_pdfas_305894", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1216", "305895", "PDF_variationNNPDF31_nlo_hessian_pdfas_305895", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1217", "305896", "PDF_variationNNPDF31_nlo_hessian_pdfas_305896", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1218", "305897", "PDF_variationNNPDF31_nlo_hessian_pdfas_305897", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1219", "305898", "PDF_variationNNPDF31_nlo_hessian_pdfas_305898", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1220", "305899", "PDF_variationNNPDF31_nlo_hessian_pdfas_305899", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1221", "305900", "PDF_variationNNPDF31_nlo_hessian_pdfas_305900", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1222", "305901", "PDF_variationNNPDF31_nlo_hessian_pdfas_305901", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1223", "305902", "PDF_variationNNPDF31_nlo_hessian_pdfas_305902", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1224", "13000", "PDF_variationCT14nnlo_13000", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1225", "13001", "PDF_variationCT14nnlo_13001", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1226", "13002", "PDF_variationCT14nnlo_13002", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1227", "13003", "PDF_variationCT14nnlo_13003", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1228", "13004", "PDF_variationCT14nnlo_13004", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1229", "13005", "PDF_variationCT14nnlo_13005", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1230", "13006", "PDF_variationCT14nnlo_13006", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1231", "13007", "PDF_variationCT14nnlo_13007", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1232", "13008", "PDF_variationCT14nnlo_13008", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1233", "13009", "PDF_variationCT14nnlo_13009", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1234", "13010", "PDF_variationCT14nnlo_13010", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1235", "13011", "PDF_variationCT14nnlo_13011", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1236", "13012", "PDF_variationCT14nnlo_13012", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1237", "13013", "PDF_variationCT14nnlo_13013", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1238", "13014", "PDF_variationCT14nnlo_13014", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1239", "13015", "PDF_variationCT14nnlo_13015", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1240", "13016", "PDF_variationCT14nnlo_13016", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1241", "13017", "PDF_variationCT14nnlo_13017", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1242", "13018", "PDF_variationCT14nnlo_13018", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1243", "13019", "PDF_variationCT14nnlo_13019", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1244", "13020", "PDF_variationCT14nnlo_13020", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1245", "13021", "PDF_variationCT14nnlo_13021", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1246", "13022", "PDF_variationCT14nnlo_13022", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1247", "13023", "PDF_variationCT14nnlo_13023", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1248", "13024", "PDF_variationCT14nnlo_13024", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1249", "13025", "PDF_variationCT14nnlo_13025", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1250", "13026", "PDF_variationCT14nnlo_13026", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1251", "13027", "PDF_variationCT14nnlo_13027", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1252", "13028", "PDF_variationCT14nnlo_13028", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1253", "13029", "PDF_variationCT14nnlo_13029", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1254", "13030", "PDF_variationCT14nnlo_13030", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1255", "13031", "PDF_variationCT14nnlo_13031", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1256", "13032", "PDF_variationCT14nnlo_13032", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1257", "13033", "PDF_variationCT14nnlo_13033", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1258", "13034", "PDF_variationCT14nnlo_13034", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1259", "13035", "PDF_variationCT14nnlo_13035", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1260", "13036", "PDF_variationCT14nnlo_13036", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1261", "13037", "PDF_variationCT14nnlo_13037", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1262", "13038", "PDF_variationCT14nnlo_13038", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1263", "13039", "PDF_variationCT14nnlo_13039", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1264", "13040", "PDF_variationCT14nnlo_13040", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1265", "13041", "PDF_variationCT14nnlo_13041", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1266", "13042", "PDF_variationCT14nnlo_13042", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1267", "13043", "PDF_variationCT14nnlo_13043", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1268", "13044", "PDF_variationCT14nnlo_13044", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1269", "13045", "PDF_variationCT14nnlo_13045", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1270", "13046", "PDF_variationCT14nnlo_13046", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1271", "13047", "PDF_variationCT14nnlo_13047", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1272", "13048", "PDF_variationCT14nnlo_13048", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1273", "13049", "PDF_variationCT14nnlo_13049", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1274", "13050", "PDF_variationCT14nnlo_13050", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1275", "13051", "PDF_variationCT14nnlo_13051", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1276", "13052", "PDF_variationCT14nnlo_13052", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1277", "13053", "PDF_variationCT14nnlo_13053", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1278", "13054", "PDF_variationCT14nnlo_13054", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1279", "13055", "PDF_variationCT14nnlo_13055", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1280", "13056", "PDF_variationCT14nnlo_13056", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1281", "13065", "PDF_variation_13065", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1282", "13069", "PDF_variation_13069", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1283", "13100", "PDF_variationCT14nlo_13100", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1284", "13101", "PDF_variationCT14nlo_13101", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1285", "13102", "PDF_variationCT14nlo_13102", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1286", "13103", "PDF_variationCT14nlo_13103", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1287", "13104", "PDF_variationCT14nlo_13104", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1288", "13105", "PDF_variationCT14nlo_13105", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1289", "13106", "PDF_variationCT14nlo_13106", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1290", "13107", "PDF_variationCT14nlo_13107", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1291", "13108", "PDF_variationCT14nlo_13108", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1292", "13109", "PDF_variationCT14nlo_13109", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1293", "13110", "PDF_variationCT14nlo_13110", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1294", "13111", "PDF_variationCT14nlo_13111", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1295", "13112", "PDF_variationCT14nlo_13112", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1296", "13113", "PDF_variationCT14nlo_13113", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1297", "13114", "PDF_variationCT14nlo_13114", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1298", "13115", "PDF_variationCT14nlo_13115", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1299", "13116", "PDF_variationCT14nlo_13116", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1300", "13117", "PDF_variationCT14nlo_13117", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1301", "13118", "PDF_variationCT14nlo_13118", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1302", "13119", "PDF_variationCT14nlo_13119", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1303", "13120", "PDF_variationCT14nlo_13120", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1304", "13121", "PDF_variationCT14nlo_13121", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1305", "13122", "PDF_variationCT14nlo_13122", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1306", "13123", "PDF_variationCT14nlo_13123", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1307", "13124", "PDF_variationCT14nlo_13124", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1308", "13125", "PDF_variationCT14nlo_13125", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1309", "13126", "PDF_variationCT14nlo_13126", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1310", "13127", "PDF_variationCT14nlo_13127", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1311", "13128", "PDF_variationCT14nlo_13128", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1312", "13129", "PDF_variationCT14nlo_13129", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1313", "13130", "PDF_variationCT14nlo_13130", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1314", "13131", "PDF_variationCT14nlo_13131", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1315", "13132", "PDF_variationCT14nlo_13132", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1316", "13133", "PDF_variationCT14nlo_13133", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1317", "13134", "PDF_variationCT14nlo_13134", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1318", "13135", "PDF_variationCT14nlo_13135", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1319", "13136", "PDF_variationCT14nlo_13136", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1320", "13137", "PDF_variationCT14nlo_13137", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1321", "13138", "PDF_variationCT14nlo_13138", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1322", "13139", "PDF_variationCT14nlo_13139", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1323", "13140", "PDF_variationCT14nlo_13140", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1324", "13141", "PDF_variationCT14nlo_13141", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1325", "13142", "PDF_variationCT14nlo_13142", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1326", "13143", "PDF_variationCT14nlo_13143", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1327", "13144", "PDF_variationCT14nlo_13144", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1328", "13145", "PDF_variationCT14nlo_13145", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1329", "13146", "PDF_variationCT14nlo_13146", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1330", "13147", "PDF_variationCT14nlo_13147", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1331", "13148", "PDF_variationCT14nlo_13148", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1332", "13149", "PDF_variationCT14nlo_13149", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1333", "13150", "PDF_variationCT14nlo_13150", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1334", "13151", "PDF_variationCT14nlo_13151", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1335", "13152", "PDF_variationCT14nlo_13152", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1336", "13153", "PDF_variationCT14nlo_13153", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1337", "13154", "PDF_variationCT14nlo_13154", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1338", "13155", "PDF_variationCT14nlo_13155", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1339", "13156", "PDF_variationCT14nlo_13156", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1340", "13163", "PDF_variation_13163", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1341", "13167", "PDF_variation_13167", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1342", "13200", "PDF_variation_13200", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1343", "25200", "PDF_variationMMHT2014nlo68clas118_25200", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1344", "25201", "PDF_variationMMHT2014nlo68clas118_25201", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1345", "25202", "PDF_variationMMHT2014nlo68clas118_25202", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1346", "25203", "PDF_variationMMHT2014nlo68clas118_25203", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1347", "25204", "PDF_variationMMHT2014nlo68clas118_25204", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1348", "25205", "PDF_variationMMHT2014nlo68clas118_25205", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1349", "25206", "PDF_variationMMHT2014nlo68clas118_25206", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1350", "25207", "PDF_variationMMHT2014nlo68clas118_25207", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1351", "25208", "PDF_variationMMHT2014nlo68clas118_25208", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1352", "25209", "PDF_variationMMHT2014nlo68clas118_25209", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1353", "25210", "PDF_variationMMHT2014nlo68clas118_25210", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1354", "25211", "PDF_variationMMHT2014nlo68clas118_25211", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1355", "25212", "PDF_variationMMHT2014nlo68clas118_25212", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1356", "25213", "PDF_variationMMHT2014nlo68clas118_25213", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1357", "25214", "PDF_variationMMHT2014nlo68clas118_25214", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1358", "25215", "PDF_variationMMHT2014nlo68clas118_25215", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1359", "25216", "PDF_variationMMHT2014nlo68clas118_25216", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1360", "25217", "PDF_variationMMHT2014nlo68clas118_25217", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1361", "25218", "PDF_variationMMHT2014nlo68clas118_25218", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1362", "25219", "PDF_variationMMHT2014nlo68clas118_25219", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1363", "25220", "PDF_variationMMHT2014nlo68clas118_25220", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1364", "25221", "PDF_variationMMHT2014nlo68clas118_25221", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1365", "25222", "PDF_variationMMHT2014nlo68clas118_25222", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1366", "25223", "PDF_variationMMHT2014nlo68clas118_25223", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1367", "25224", "PDF_variationMMHT2014nlo68clas118_25224", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1368", "25225", "PDF_variationMMHT2014nlo68clas118_25225", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1369", "25226", "PDF_variationMMHT2014nlo68clas118_25226", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1370", "25227", "PDF_variationMMHT2014nlo68clas118_25227", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1371", "25228", "PDF_variationMMHT2014nlo68clas118_25228", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1372", "25229", "PDF_variationMMHT2014nlo68clas118_25229", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1373", "25230", "PDF_variationMMHT2014nlo68clas118_25230", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1374", "25231", "PDF_variationMMHT2014nlo68clas118_25231", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1375", "25232", "PDF_variationMMHT2014nlo68clas118_25232", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1376", "25233", "PDF_variationMMHT2014nlo68clas118_25233", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1377", "25234", "PDF_variationMMHT2014nlo68clas118_25234", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1378", "25235", "PDF_variationMMHT2014nlo68clas118_25235", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1379", "25236", "PDF_variationMMHT2014nlo68clas118_25236", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1380", "25237", "PDF_variationMMHT2014nlo68clas118_25237", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1381", "25238", "PDF_variationMMHT2014nlo68clas118_25238", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1382", "25239", "PDF_variationMMHT2014nlo68clas118_25239", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1383", "25240", "PDF_variationMMHT2014nlo68clas118_25240", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1384", "25241", "PDF_variationMMHT2014nlo68clas118_25241", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1385", "25242", "PDF_variationMMHT2014nlo68clas118_25242", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1386", "25243", "PDF_variationMMHT2014nlo68clas118_25243", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1387", "25244", "PDF_variationMMHT2014nlo68clas118_25244", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1388", "25245", "PDF_variationMMHT2014nlo68clas118_25245", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1389", "25246", "PDF_variationMMHT2014nlo68clas118_25246", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1390", "25247", "PDF_variationMMHT2014nlo68clas118_25247", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1391", "25248", "PDF_variationMMHT2014nlo68clas118_25248", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1392", "25249", "PDF_variationMMHT2014nlo68clas118_25249", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1393", "25250", "PDF_variationMMHT2014nlo68clas118_25250", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1394", "25300", "PDF_variationMMHT2014nnlo68cl_25300", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1395", "25301", "PDF_variationMMHT2014nnlo68cl_25301", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1396", "25302", "PDF_variationMMHT2014nnlo68cl_25302", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1397", "25303", "PDF_variationMMHT2014nnlo68cl_25303", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1398", "25304", "PDF_variationMMHT2014nnlo68cl_25304", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1399", "25305", "PDF_variationMMHT2014nnlo68cl_25305", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1400", "25306", "PDF_variationMMHT2014nnlo68cl_25306", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1401", "25307", "PDF_variationMMHT2014nnlo68cl_25307", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1402", "25308", "PDF_variationMMHT2014nnlo68cl_25308", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1403", "25309", "PDF_variationMMHT2014nnlo68cl_25309", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1404", "25310", "PDF_variationMMHT2014nnlo68cl_25310", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1405", "25311", "PDF_variationMMHT2014nnlo68cl_25311", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1406", "25312", "PDF_variationMMHT2014nnlo68cl_25312", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1407", "25313", "PDF_variationMMHT2014nnlo68cl_25313", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1408", "25314", "PDF_variationMMHT2014nnlo68cl_25314", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1409", "25315", "PDF_variationMMHT2014nnlo68cl_25315", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1410", "25316", "PDF_variationMMHT2014nnlo68cl_25316", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1411", "25317", "PDF_variationMMHT2014nnlo68cl_25317", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1412", "25318", "PDF_variationMMHT2014nnlo68cl_25318", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1413", "25319", "PDF_variationMMHT2014nnlo68cl_25319", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1414", "25320", "PDF_variationMMHT2014nnlo68cl_25320", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1415", "25321", "PDF_variationMMHT2014nnlo68cl_25321", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1416", "25322", "PDF_variationMMHT2014nnlo68cl_25322", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1417", "25323", "PDF_variationMMHT2014nnlo68cl_25323", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1418", "25324", "PDF_variationMMHT2014nnlo68cl_25324", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1419", "25325", "PDF_variationMMHT2014nnlo68cl_25325", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1420", "25326", "PDF_variationMMHT2014nnlo68cl_25326", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1421", "25327", "PDF_variationMMHT2014nnlo68cl_25327", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1422", "25328", "PDF_variationMMHT2014nnlo68cl_25328", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1423", "25329", "PDF_variationMMHT2014nnlo68cl_25329", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1424", "25330", "PDF_variationMMHT2014nnlo68cl_25330", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1425", "25331", "PDF_variationMMHT2014nnlo68cl_25331", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1426", "25332", "PDF_variationMMHT2014nnlo68cl_25332", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1427", "25333", "PDF_variationMMHT2014nnlo68cl_25333", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1428", "25334", "PDF_variationMMHT2014nnlo68cl_25334", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1429", "25335", "PDF_variationMMHT2014nnlo68cl_25335", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1430", "25336", "PDF_variationMMHT2014nnlo68cl_25336", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1431", "25337", "PDF_variationMMHT2014nnlo68cl_25337", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1432", "25338", "PDF_variationMMHT2014nnlo68cl_25338", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1433", "25339", "PDF_variationMMHT2014nnlo68cl_25339", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1434", "25340", "PDF_variationMMHT2014nnlo68cl_25340", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1435", "25341", "PDF_variationMMHT2014nnlo68cl_25341", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1436", "25342", "PDF_variationMMHT2014nnlo68cl_25342", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1437", "25343", "PDF_variationMMHT2014nnlo68cl_25343", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1438", "25344", "PDF_variationMMHT2014nnlo68cl_25344", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1439", "25345", "PDF_variationMMHT2014nnlo68cl_25345", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1440", "25346", "PDF_variationMMHT2014nnlo68cl_25346", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1441", "25347", "PDF_variationMMHT2014nnlo68cl_25347", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1442", "25348", "PDF_variationMMHT2014nnlo68cl_25348", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1443", "25349", "PDF_variationMMHT2014nnlo68cl_25349", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1444", "25350", "PDF_variationMMHT2014nnlo68cl_25350", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1445", "25000", "PDF_variation_25000", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1446", "42780", "PDF_variationABMP16als118_5_nnlo_42780", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1447", "42781", "PDF_variationABMP16als118_5_nnlo_42781", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1448", "42782", "PDF_variationABMP16als118_5_nnlo_42782", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1449", "42783", "PDF_variationABMP16als118_5_nnlo_42783", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1450", "42784", "PDF_variationABMP16als118_5_nnlo_42784", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1451", "42785", "PDF_variationABMP16als118_5_nnlo_42785", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1452", "42786", "PDF_variationABMP16als118_5_nnlo_42786", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1453", "42787", "PDF_variationABMP16als118_5_nnlo_42787", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1454", "42788", "PDF_variationABMP16als118_5_nnlo_42788", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1455", "42789", "PDF_variationABMP16als118_5_nnlo_42789", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1456", "42790", "PDF_variationABMP16als118_5_nnlo_42790", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1457", "42791", "PDF_variationABMP16als118_5_nnlo_42791", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1458", "42792", "PDF_variationABMP16als118_5_nnlo_42792", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1459", "42793", "PDF_variationABMP16als118_5_nnlo_42793", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1460", "42794", "PDF_variationABMP16als118_5_nnlo_42794", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1461", "42795", "PDF_variationABMP16als118_5_nnlo_42795", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1462", "42796", "PDF_variationABMP16als118_5_nnlo_42796", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1463", "42797", "PDF_variationABMP16als118_5_nnlo_42797", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1464", "42798", "PDF_variationABMP16als118_5_nnlo_42798", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1465", "42799", "PDF_variationABMP16als118_5_nnlo_42799", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1466", "42800", "PDF_variationABMP16als118_5_nnlo_42800", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1467", "42801", "PDF_variationABMP16als118_5_nnlo_42801", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1468", "42802", "PDF_variationABMP16als118_5_nnlo_42802", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1469", "42803", "PDF_variationABMP16als118_5_nnlo_42803", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1470", "42804", "PDF_variationABMP16als118_5_nnlo_42804", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1471", "42805", "PDF_variationABMP16als118_5_nnlo_42805", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1472", "42806", "PDF_variationABMP16als118_5_nnlo_42806", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1473", "42807", "PDF_variationABMP16als118_5_nnlo_42807", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1474", "42808", "PDF_variationABMP16als118_5_nnlo_42808", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1475", "42809", "PDF_variationABMP16als118_5_nnlo_42809", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1476", "90200", "PDF_variationPDF4LHC15_nlo_100_pdfas_90200", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1477", "90201", "PDF_variationPDF4LHC15_nlo_100_pdfas_90201", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1478", "90202", "PDF_variationPDF4LHC15_nlo_100_pdfas_90202", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1479", "90203", "PDF_variationPDF4LHC15_nlo_100_pdfas_90203", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1480", "90204", "PDF_variationPDF4LHC15_nlo_100_pdfas_90204", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1481", "90205", "PDF_variationPDF4LHC15_nlo_100_pdfas_90205", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1482", "90206", "PDF_variationPDF4LHC15_nlo_100_pdfas_90206", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1483", "90207", "PDF_variationPDF4LHC15_nlo_100_pdfas_90207", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1484", "90208", "PDF_variationPDF4LHC15_nlo_100_pdfas_90208", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1485", "90209", "PDF_variationPDF4LHC15_nlo_100_pdfas_90209", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1486", "90210", "PDF_variationPDF4LHC15_nlo_100_pdfas_90210", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1487", "90211", "PDF_variationPDF4LHC15_nlo_100_pdfas_90211", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1488", "90212", "PDF_variationPDF4LHC15_nlo_100_pdfas_90212", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1489", "90213", "PDF_variationPDF4LHC15_nlo_100_pdfas_90213", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1490", "90214", "PDF_variationPDF4LHC15_nlo_100_pdfas_90214", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1491", "90215", "PDF_variationPDF4LHC15_nlo_100_pdfas_90215", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1492", "90216", "PDF_variationPDF4LHC15_nlo_100_pdfas_90216", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1493", "90217", "PDF_variationPDF4LHC15_nlo_100_pdfas_90217", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1494", "90218", "PDF_variationPDF4LHC15_nlo_100_pdfas_90218", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1495", "90219", "PDF_variationPDF4LHC15_nlo_100_pdfas_90219", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1496", "90220", "PDF_variationPDF4LHC15_nlo_100_pdfas_90220", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1497", "90221", "PDF_variationPDF4LHC15_nlo_100_pdfas_90221", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1498", "90222", "PDF_variationPDF4LHC15_nlo_100_pdfas_90222", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1499", "90223", "PDF_variationPDF4LHC15_nlo_100_pdfas_90223", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1500", "90224", "PDF_variationPDF4LHC15_nlo_100_pdfas_90224", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1501", "90225", "PDF_variationPDF4LHC15_nlo_100_pdfas_90225", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1502", "90226", "PDF_variationPDF4LHC15_nlo_100_pdfas_90226", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1503", "90227", "PDF_variationPDF4LHC15_nlo_100_pdfas_90227", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1504", "90228", "PDF_variationPDF4LHC15_nlo_100_pdfas_90228", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1505", "90229", "PDF_variationPDF4LHC15_nlo_100_pdfas_90229", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1506", "90230", "PDF_variationPDF4LHC15_nlo_100_pdfas_90230", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1507", "90231", "PDF_variationPDF4LHC15_nlo_100_pdfas_90231", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1508", "90232", "PDF_variationPDF4LHC15_nlo_100_pdfas_90232", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1509", "90233", "PDF_variationPDF4LHC15_nlo_100_pdfas_90233", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1510", "90234", "PDF_variationPDF4LHC15_nlo_100_pdfas_90234", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1511", "90235", "PDF_variationPDF4LHC15_nlo_100_pdfas_90235", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1512", "90236", "PDF_variationPDF4LHC15_nlo_100_pdfas_90236", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1513", "90237", "PDF_variationPDF4LHC15_nlo_100_pdfas_90237", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1514", "90238", "PDF_variationPDF4LHC15_nlo_100_pdfas_90238", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1515", "90239", "PDF_variationPDF4LHC15_nlo_100_pdfas_90239", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1516", "90240", "PDF_variationPDF4LHC15_nlo_100_pdfas_90240", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1517", "90241", "PDF_variationPDF4LHC15_nlo_100_pdfas_90241", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1518", "90242", "PDF_variationPDF4LHC15_nlo_100_pdfas_90242", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1519", "90243", "PDF_variationPDF4LHC15_nlo_100_pdfas_90243", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1520", "90244", "PDF_variationPDF4LHC15_nlo_100_pdfas_90244", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1521", "90245", "PDF_variationPDF4LHC15_nlo_100_pdfas_90245", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1522", "90246", "PDF_variationPDF4LHC15_nlo_100_pdfas_90246", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1523", "90247", "PDF_variationPDF4LHC15_nlo_100_pdfas_90247", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1524", "90248", "PDF_variationPDF4LHC15_nlo_100_pdfas_90248", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1525", "90249", "PDF_variationPDF4LHC15_nlo_100_pdfas_90249", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1526", "90250", "PDF_variationPDF4LHC15_nlo_100_pdfas_90250", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1527", "90251", "PDF_variationPDF4LHC15_nlo_100_pdfas_90251", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1528", "90252", "PDF_variationPDF4LHC15_nlo_100_pdfas_90252", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1529", "90253", "PDF_variationPDF4LHC15_nlo_100_pdfas_90253", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1530", "90254", "PDF_variationPDF4LHC15_nlo_100_pdfas_90254", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1531", "90255", "PDF_variationPDF4LHC15_nlo_100_pdfas_90255", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1532", "90256", "PDF_variationPDF4LHC15_nlo_100_pdfas_90256", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1533", "90257", "PDF_variationPDF4LHC15_nlo_100_pdfas_90257", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1534", "90258", "PDF_variationPDF4LHC15_nlo_100_pdfas_90258", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1535", "90259", "PDF_variationPDF4LHC15_nlo_100_pdfas_90259", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1536", "90260", "PDF_variationPDF4LHC15_nlo_100_pdfas_90260", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1537", "90261", "PDF_variationPDF4LHC15_nlo_100_pdfas_90261", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1538", "90262", "PDF_variationPDF4LHC15_nlo_100_pdfas_90262", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1539", "90263", "PDF_variationPDF4LHC15_nlo_100_pdfas_90263", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1540", "90264", "PDF_variationPDF4LHC15_nlo_100_pdfas_90264", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1541", "90265", "PDF_variationPDF4LHC15_nlo_100_pdfas_90265", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1542", "90266", "PDF_variationPDF4LHC15_nlo_100_pdfas_90266", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1543", "90267", "PDF_variationPDF4LHC15_nlo_100_pdfas_90267", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1544", "90268", "PDF_variationPDF4LHC15_nlo_100_pdfas_90268", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1545", "90269", "PDF_variationPDF4LHC15_nlo_100_pdfas_90269", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1546", "90270", "PDF_variationPDF4LHC15_nlo_100_pdfas_90270", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1547", "90271", "PDF_variationPDF4LHC15_nlo_100_pdfas_90271", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1548", "90272", "PDF_variationPDF4LHC15_nlo_100_pdfas_90272", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1549", "90273", "PDF_variationPDF4LHC15_nlo_100_pdfas_90273", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1550", "90274", "PDF_variationPDF4LHC15_nlo_100_pdfas_90274", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1551", "90275", "PDF_variationPDF4LHC15_nlo_100_pdfas_90275", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1552", "90276", "PDF_variationPDF4LHC15_nlo_100_pdfas_90276", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1553", "90277", "PDF_variationPDF4LHC15_nlo_100_pdfas_90277", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1554", "90278", "PDF_variationPDF4LHC15_nlo_100_pdfas_90278", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1555", "90279", "PDF_variationPDF4LHC15_nlo_100_pdfas_90279", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1556", "90280", "PDF_variationPDF4LHC15_nlo_100_pdfas_90280", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1557", "90281", "PDF_variationPDF4LHC15_nlo_100_pdfas_90281", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1558", "90282", "PDF_variationPDF4LHC15_nlo_100_pdfas_90282", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1559", "90283", "PDF_variationPDF4LHC15_nlo_100_pdfas_90283", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1560", "90284", "PDF_variationPDF4LHC15_nlo_100_pdfas_90284", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1561", "90285", "PDF_variationPDF4LHC15_nlo_100_pdfas_90285", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1562", "90286", "PDF_variationPDF4LHC15_nlo_100_pdfas_90286", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1563", "90287", "PDF_variationPDF4LHC15_nlo_100_pdfas_90287", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1564", "90288", "PDF_variationPDF4LHC15_nlo_100_pdfas_90288", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1565", "90289", "PDF_variationPDF4LHC15_nlo_100_pdfas_90289", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1566", "90290", "PDF_variationPDF4LHC15_nlo_100_pdfas_90290", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1567", "90291", "PDF_variationPDF4LHC15_nlo_100_pdfas_90291", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1568", "90292", "PDF_variationPDF4LHC15_nlo_100_pdfas_90292", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1569", "90293", "PDF_variationPDF4LHC15_nlo_100_pdfas_90293", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1570", "90294", "PDF_variationPDF4LHC15_nlo_100_pdfas_90294", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1571", "90295", "PDF_variationPDF4LHC15_nlo_100_pdfas_90295", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1572", "90296", "PDF_variationPDF4LHC15_nlo_100_pdfas_90296", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1573", "90297", "PDF_variationPDF4LHC15_nlo_100_pdfas_90297", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1574", "90298", "PDF_variationPDF4LHC15_nlo_100_pdfas_90298", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1575", "90299", "PDF_variationPDF4LHC15_nlo_100_pdfas_90299", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1576", "90300", "PDF_variationPDF4LHC15_nlo_100_pdfas_90300", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1577", "90301", "PDF_variationPDF4LHC15_nlo_100_pdfas_90301", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1578", "90302", "PDF_variationPDF4LHC15_nlo_100_pdfas_90302", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1579", "91200", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91200", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1580", "91201", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91201", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1581", "91202", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91202", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1582", "91203", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91203", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1583", "91204", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91204", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1584", "91205", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91205", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1585", "91206", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91206", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1586", "91207", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91207", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1587", "91208", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91208", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1588", "91209", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91209", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1589", "91210", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91210", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1590", "91211", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91211", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1591", "91212", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91212", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1592", "91213", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91213", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1593", "91214", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91214", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1594", "91215", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91215", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1595", "91216", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91216", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1596", "91217", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91217", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1597", "91218", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91218", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1598", "91219", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91219", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1599", "91220", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91220", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1600", "91221", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91221", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1601", "91222", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91222", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1602", "91223", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91223", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1603", "91224", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91224", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1604", "91225", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91225", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1605", "91226", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91226", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1606", "91227", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91227", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1607", "91228", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91228", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1608", "91229", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91229", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1609", "91230", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91230", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1610", "91231", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91231", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1611", "91232", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91232", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1612", "91233", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91233", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1613", "91234", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91234", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1614", "91235", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91235", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1615", "91236", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91236", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1616", "91237", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91237", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1617", "91238", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91238", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1618", "91239", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91239", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1619", "91240", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91240", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1620", "91241", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91241", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1621", "91242", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91242", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1622", "91243", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91243", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1623", "91244", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91244", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1624", "91245", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91245", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1625", "91246", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91246", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1626", "91247", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91247", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1627", "91248", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91248", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1628", "91249", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91249", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1629", "91250", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91250", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1630", "91251", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91251", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1631", "91252", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91252", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1632", "91253", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91253", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1633", "91254", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91254", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1634", "91255", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91255", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1635", "91256", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91256", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1636", "91257", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91257", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1637", "91258", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91258", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1638", "91259", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91259", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1639", "91260", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91260", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1640", "91261", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91261", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1641", "91262", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91262", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1642", "91263", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91263", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1643", "91264", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91264", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1644", "91265", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91265", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1645", "91266", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91266", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1646", "91267", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91267", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1647", "91268", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91268", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1648", "91269", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91269", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1649", "91270", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91270", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1650", "91271", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91271", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1651", "91272", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91272", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1652", "91273", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91273", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1653", "91274", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91274", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1654", "91275", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91275", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1655", "91276", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91276", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1656", "91277", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91277", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1657", "91278", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91278", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1658", "91279", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91279", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1659", "91280", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91280", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1660", "91281", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91281", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1661", "91282", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91282", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1662", "91283", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91283", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1663", "91284", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91284", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1664", "91285", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91285", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1665", "91286", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91286", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1666", "91287", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91287", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1667", "91288", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91288", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1668", "91289", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91289", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1669", "91290", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91290", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1670", "91291", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91291", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1671", "91292", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91292", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1672", "91293", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91293", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1673", "91294", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91294", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1674", "91295", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91295", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1675", "91296", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91296", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1676", "91297", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91297", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1677", "91298", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91298", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1678", "91299", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91299", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1679", "91300", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91300", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1680", "91301", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91301", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1681", "91302", "PDF_variationPDF4LHC15_nnlo_100_pdfas_91302", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1682", "90400", "PDF_variationPDF4LHC15_nlo_30_pdfas_90400", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1683", "90401", "PDF_variationPDF4LHC15_nlo_30_pdfas_90401", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1684", "90402", "PDF_variationPDF4LHC15_nlo_30_pdfas_90402", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1685", "90403", "PDF_variationPDF4LHC15_nlo_30_pdfas_90403", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1686", "90404", "PDF_variationPDF4LHC15_nlo_30_pdfas_90404", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1687", "90405", "PDF_variationPDF4LHC15_nlo_30_pdfas_90405", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1688", "90406", "PDF_variationPDF4LHC15_nlo_30_pdfas_90406", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1689", "90407", "PDF_variationPDF4LHC15_nlo_30_pdfas_90407", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1690", "90408", "PDF_variationPDF4LHC15_nlo_30_pdfas_90408", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1691", "90409", "PDF_variationPDF4LHC15_nlo_30_pdfas_90409", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1692", "90410", "PDF_variationPDF4LHC15_nlo_30_pdfas_90410", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1693", "90411", "PDF_variationPDF4LHC15_nlo_30_pdfas_90411", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1694", "90412", "PDF_variationPDF4LHC15_nlo_30_pdfas_90412", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1695", "90413", "PDF_variationPDF4LHC15_nlo_30_pdfas_90413", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1696", "90414", "PDF_variationPDF4LHC15_nlo_30_pdfas_90414", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1697", "90415", "PDF_variationPDF4LHC15_nlo_30_pdfas_90415", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1698", "90416", "PDF_variationPDF4LHC15_nlo_30_pdfas_90416", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1699", "90417", "PDF_variationPDF4LHC15_nlo_30_pdfas_90417", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1700", "90418", "PDF_variationPDF4LHC15_nlo_30_pdfas_90418", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1701", "90419", "PDF_variationPDF4LHC15_nlo_30_pdfas_90419", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1702", "90420", "PDF_variationPDF4LHC15_nlo_30_pdfas_90420", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1703", "90421", "PDF_variationPDF4LHC15_nlo_30_pdfas_90421", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1704", "90422", "PDF_variationPDF4LHC15_nlo_30_pdfas_90422", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1705", "90423", "PDF_variationPDF4LHC15_nlo_30_pdfas_90423", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1706", "90424", "PDF_variationPDF4LHC15_nlo_30_pdfas_90424", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1707", "90425", "PDF_variationPDF4LHC15_nlo_30_pdfas_90425", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1708", "90426", "PDF_variationPDF4LHC15_nlo_30_pdfas_90426", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1709", "90427", "PDF_variationPDF4LHC15_nlo_30_pdfas_90427", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1710", "90428", "PDF_variationPDF4LHC15_nlo_30_pdfas_90428", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1711", "90429", "PDF_variationPDF4LHC15_nlo_30_pdfas_90429", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1712", "90430", "PDF_variationPDF4LHC15_nlo_30_pdfas_90430", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1713", "90431", "PDF_variationPDF4LHC15_nlo_30_pdfas_90431", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1714", "90432", "PDF_variationPDF4LHC15_nlo_30_pdfas_90432", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1715", "91400", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91400", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1716", "91401", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91401", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1717", "91402", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91402", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1718", "91403", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91403", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1719", "91404", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91404", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1720", "91405", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91405", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1721", "91406", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91406", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1722", "91407", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91407", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1723", "91408", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91408", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1724", "91409", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91409", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1725", "91410", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91410", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1726", "91411", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91411", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1727", "91412", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91412", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1728", "91413", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91413", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1729", "91414", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91414", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1730", "91415", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91415", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1731", "91416", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91416", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1732", "91417", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91417", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1733", "91418", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91418", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1734", "91419", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91419", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1735", "91420", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91420", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1736", "91421", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91421", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1737", "91422", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91422", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1738", "91423", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91423", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1739", "91424", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91424", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1740", "91425", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91425", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1741", "91426", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91426", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1742", "91427", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91427", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1743", "91428", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91428", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1744", "91429", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91429", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1745", "91430", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91430", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1746", "91431", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91431", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1747", "91432", "PDF_variationPDF4LHC15_nnlo_30_pdfas_91432", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1748", "61100", "PDF_variationHERAPDF20_NLO_EIG_61100", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1749", "61101", "PDF_variationHERAPDF20_NLO_EIG_61101", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1750", "61102", "PDF_variationHERAPDF20_NLO_EIG_61102", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1751", "61103", "PDF_variationHERAPDF20_NLO_EIG_61103", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1752", "61104", "PDF_variationHERAPDF20_NLO_EIG_61104", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1753", "61105", "PDF_variationHERAPDF20_NLO_EIG_61105", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1754", "61106", "PDF_variationHERAPDF20_NLO_EIG_61106", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1755", "61107", "PDF_variationHERAPDF20_NLO_EIG_61107", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1756", "61108", "PDF_variationHERAPDF20_NLO_EIG_61108", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1757", "61109", "PDF_variationHERAPDF20_NLO_EIG_61109", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1758", "61110", "PDF_variationHERAPDF20_NLO_EIG_61110", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1759", "61111", "PDF_variationHERAPDF20_NLO_EIG_61111", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1760", "61112", "PDF_variationHERAPDF20_NLO_EIG_61112", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1761", "61113", "PDF_variationHERAPDF20_NLO_EIG_61113", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1762", "61114", "PDF_variationHERAPDF20_NLO_EIG_61114", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1763", "61115", "PDF_variationHERAPDF20_NLO_EIG_61115", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1764", "61116", "PDF_variationHERAPDF20_NLO_EIG_61116", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1765", "61117", "PDF_variationHERAPDF20_NLO_EIG_61117", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1766", "61118", "PDF_variationHERAPDF20_NLO_EIG_61118", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1767", "61119", "PDF_variationHERAPDF20_NLO_EIG_61119", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1768", "61120", "PDF_variationHERAPDF20_NLO_EIG_61120", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1769", "61121", "PDF_variationHERAPDF20_NLO_EIG_61121", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1770", "61122", "PDF_variationHERAPDF20_NLO_EIG_61122", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1771", "61123", "PDF_variationHERAPDF20_NLO_EIG_61123", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1772", "61124", "PDF_variationHERAPDF20_NLO_EIG_61124", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1773", "61125", "PDF_variationHERAPDF20_NLO_EIG_61125", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1774", "61126", "PDF_variationHERAPDF20_NLO_EIG_61126", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1775", "61127", "PDF_variationHERAPDF20_NLO_EIG_61127", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1776", "61128", "PDF_variationHERAPDF20_NLO_EIG_61128", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1777", "61130", "PDF_variationHERAPDF20_NLO_VAR_61130", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1778", "61131", "PDF_variationHERAPDF20_NLO_VAR_61131", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1779", "61132", "PDF_variationHERAPDF20_NLO_VAR_61132", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1780", "61133", "PDF_variationHERAPDF20_NLO_VAR_61133", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1781", "61134", "PDF_variationHERAPDF20_NLO_VAR_61134", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1782", "61135", "PDF_variationHERAPDF20_NLO_VAR_61135", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1783", "61136", "PDF_variationHERAPDF20_NLO_VAR_61136", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1784", "61137", "PDF_variationHERAPDF20_NLO_VAR_61137", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1785", "61138", "PDF_variationHERAPDF20_NLO_VAR_61138", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1786", "61139", "PDF_variationHERAPDF20_NLO_VAR_61139", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1787", "61140", "PDF_variationHERAPDF20_NLO_VAR_61140", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1788", "61141", "PDF_variationHERAPDF20_NLO_VAR_61141", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1789", "61142", "PDF_variationHERAPDF20_NLO_VAR_61142", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1790", "61143", "PDF_variationHERAPDF20_NLO_VAR_61143", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1791", "61200", "PDF_variationHERAPDF20_NNLO_EIG_61200", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1792", "61201", "PDF_variationHERAPDF20_NNLO_EIG_61201", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1793", "61202", "PDF_variationHERAPDF20_NNLO_EIG_61202", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1794", "61203", "PDF_variationHERAPDF20_NNLO_EIG_61203", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1795", "61204", "PDF_variationHERAPDF20_NNLO_EIG_61204", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1796", "61205", "PDF_variationHERAPDF20_NNLO_EIG_61205", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1797", "61206", "PDF_variationHERAPDF20_NNLO_EIG_61206", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1798", "61207", "PDF_variationHERAPDF20_NNLO_EIG_61207", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1799", "61208", "PDF_variationHERAPDF20_NNLO_EIG_61208", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1800", "61209", "PDF_variationHERAPDF20_NNLO_EIG_61209", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1801", "61210", "PDF_variationHERAPDF20_NNLO_EIG_61210", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1802", "61211", "PDF_variationHERAPDF20_NNLO_EIG_61211", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1803", "61212", "PDF_variationHERAPDF20_NNLO_EIG_61212", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1804", "61213", "PDF_variationHERAPDF20_NNLO_EIG_61213", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1805", "61214", "PDF_variationHERAPDF20_NNLO_EIG_61214", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1806", "61215", "PDF_variationHERAPDF20_NNLO_EIG_61215", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1807", "61216", "PDF_variationHERAPDF20_NNLO_EIG_61216", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1808", "61217", "PDF_variationHERAPDF20_NNLO_EIG_61217", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1809", "61218", "PDF_variationHERAPDF20_NNLO_EIG_61218", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1810", "61219", "PDF_variationHERAPDF20_NNLO_EIG_61219", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1811", "61220", "PDF_variationHERAPDF20_NNLO_EIG_61220", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1812", "61221", "PDF_variationHERAPDF20_NNLO_EIG_61221", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1813", "61222", "PDF_variationHERAPDF20_NNLO_EIG_61222", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1814", "61223", "PDF_variationHERAPDF20_NNLO_EIG_61223", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1815", "61224", "PDF_variationHERAPDF20_NNLO_EIG_61224", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1816", "61225", "PDF_variationHERAPDF20_NNLO_EIG_61225", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1817", "61226", "PDF_variationHERAPDF20_NNLO_EIG_61226", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1818", "61227", "PDF_variationHERAPDF20_NNLO_EIG_61227", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1819", "61228", "PDF_variationHERAPDF20_NNLO_EIG_61228", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1820", "61230", "PDF_variationHERAPDF20_NNLO_VAR_61230", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1821", "61231", "PDF_variationHERAPDF20_NNLO_VAR_61231", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1822", "61232", "PDF_variationHERAPDF20_NNLO_VAR_61232", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1823", "61233", "PDF_variationHERAPDF20_NNLO_VAR_61233", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1824", "61234", "PDF_variationHERAPDF20_NNLO_VAR_61234", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1825", "61235", "PDF_variationHERAPDF20_NNLO_VAR_61235", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1826", "61236", "PDF_variationHERAPDF20_NNLO_VAR_61236", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1827", "61237", "PDF_variationHERAPDF20_NNLO_VAR_61237", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1828", "61238", "PDF_variationHERAPDF20_NNLO_VAR_61238", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1829", "61239", "PDF_variationHERAPDF20_NNLO_VAR_61239", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1830", "61240", "PDF_variationHERAPDF20_NNLO_VAR_61240", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1831", "61241", "PDF_variationHERAPDF20_NNLO_VAR_61241", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1832", "61242", "PDF_variationHERAPDF20_NNLO_VAR_61242", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1833", "61243", "PDF_variationHERAPDF20_NNLO_VAR_61243", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1834", "13400", "PDF_variationCT14qed_inc_proton_13400", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1835", "13401", "PDF_variationCT14qed_inc_proton_13401", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1836", "13402", "PDF_variationCT14qed_inc_proton_13402", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1837", "13403", "PDF_variationCT14qed_inc_proton_13403", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1838", "13404", "PDF_variationCT14qed_inc_proton_13404", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1839", "13405", "PDF_variationCT14qed_inc_proton_13405", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1840", "13406", "PDF_variationCT14qed_inc_proton_13406", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1841", "13407", "PDF_variationCT14qed_inc_proton_13407", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1842", "13408", "PDF_variationCT14qed_inc_proton_13408", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1843", "13409", "PDF_variationCT14qed_inc_proton_13409", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1844", "13410", "PDF_variationCT14qed_inc_proton_13410", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1845", "13411", "PDF_variationCT14qed_inc_proton_13411", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1846", "13412", "PDF_variationCT14qed_inc_proton_13412", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1847", "13413", "PDF_variationCT14qed_inc_proton_13413", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1848", "13414", "PDF_variationCT14qed_inc_proton_13414", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1849", "13415", "PDF_variationCT14qed_inc_proton_13415", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1850", "13416", "PDF_variationCT14qed_inc_proton_13416", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1851", "13417", "PDF_variationCT14qed_inc_proton_13417", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1852", "13418", "PDF_variationCT14qed_inc_proton_13418", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1853", "13419", "PDF_variationCT14qed_inc_proton_13419", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1854", "13420", "PDF_variationCT14qed_inc_proton_13420", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1855", "13421", "PDF_variationCT14qed_inc_proton_13421", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1856", "13422", "PDF_variationCT14qed_inc_proton_13422", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1857", "13423", "PDF_variationCT14qed_inc_proton_13423", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1858", "13424", "PDF_variationCT14qed_inc_proton_13424", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1859", "13425", "PDF_variationCT14qed_inc_proton_13425", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1860", "13426", "PDF_variationCT14qed_inc_proton_13426", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1861", "13427", "PDF_variationCT14qed_inc_proton_13427", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1862", "13428", "PDF_variationCT14qed_inc_proton_13428", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1863", "13429", "PDF_variationCT14qed_inc_proton_13429", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1864", "13430", "PDF_variationCT14qed_inc_proton_13430", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1865", "82200", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82200", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1866", "82201", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82201", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1867", "82202", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82202", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1868", "82203", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82203", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1869", "82204", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82204", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1870", "82205", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82205", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1871", "82206", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82206", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1872", "82207", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82207", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1873", "82208", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82208", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1874", "82209", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82209", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1875", "82210", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82210", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1876", "82211", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82211", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1877", "82212", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82212", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1878", "82213", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82213", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1879", "82214", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82214", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1880", "82215", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82215", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1881", "82216", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82216", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1882", "82217", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82217", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1883", "82218", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82218", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1884", "82219", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82219", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1885", "82220", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82220", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1886", "82221", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82221", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1887", "82222", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82222", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1888", "82223", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82223", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1889", "82224", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82224", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1890", "82225", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82225", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1891", "82226", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82226", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1892", "82227", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82227", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1893", "82228", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82228", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1894", "82229", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82229", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1895", "82230", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82230", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1896", "82231", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82231", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1897", "82232", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82232", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1898", "82233", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82233", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1899", "82234", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82234", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1900", "82235", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82235", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1901", "82236", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82236", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1902", "82237", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82237", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1903", "82238", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82238", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1904", "82239", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82239", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1905", "82240", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82240", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1906", "82241", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82241", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1907", "82242", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82242", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1908", "82243", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82243", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1909", "82244", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82244", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1910", "82245", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82245", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1911", "82246", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82246", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1912", "82247", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82247", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1913", "82248", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82248", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1914", "82249", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82249", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1915", "82250", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82250", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1916", "82251", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82251", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1917", "82252", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82252", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1918", "82253", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82253", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1919", "82254", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82254", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1920", "82255", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82255", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1921", "82256", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82256", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1922", "82257", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82257", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1923", "82258", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82258", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1924", "82259", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82259", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1925", "82260", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82260", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1926", "82261", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82261", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1927", "82262", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82262", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1928", "82263", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82263", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1929", "82264", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82264", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1930", "82265", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82265", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1931", "82266", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82266", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1932", "82267", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82267", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1933", "82268", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82268", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1934", "82269", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82269", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1935", "82270", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82270", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1936", "82271", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82271", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1937", "82272", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82272", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1938", "82273", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82273", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1939", "82274", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82274", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1940", "82275", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82275", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1941", "82276", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82276", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1942", "82277", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82277", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1943", "82278", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82278", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1944", "82279", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82279", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1945", "82280", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82280", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1946", "82281", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82281", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1947", "82282", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82282", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1948", "82283", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82283", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1949", "82284", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82284", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1950", "82285", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82285", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1951", "82286", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82286", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1952", "82287", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82287", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1953", "82288", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82288", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1954", "82289", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82289", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1955", "82290", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82290", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1956", "82291", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82291", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1957", "82292", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82292", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1958", "82293", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82293", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1959", "82294", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82294", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1960", "82295", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82295", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1961", "82296", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82296", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1962", "82297", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82297", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1963", "82298", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82298", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1964", "82299", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82299", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1965", "82300", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82300", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1966", "82301", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82301", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1967", "82302", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82302", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1968", "82303", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82303", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1969", "82304", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82304", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1970", "82305", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82305", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1971", "82306", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82306", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1972", "82307", "PDF_variationLUXqed17_plus_PDF4LHC15_nnlo_100_82307", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1973", "292200", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292200", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1974", "292201", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292201", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1975", "292202", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292202", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1976", "292203", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292203", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1977", "292204", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292204", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1978", "292205", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292205", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1979", "292206", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292206", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1980", "292207", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292207", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1981", "292208", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292208", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1982", "292209", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292209", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1983", "292210", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292210", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1984", "292211", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292211", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1985", "292212", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292212", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1986", "292213", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292213", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1987", "292214", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292214", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1988", "292215", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292215", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1989", "292216", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292216", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1990", "292217", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292217", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1991", "292218", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292218", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1992", "292219", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292219", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1993", "292220", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292220", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1994", "292221", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292221", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1995", "292222", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292222", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1996", "292223", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292223", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1997", "292224", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292224", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1998", "292225", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292225", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("1999", "292226", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292226", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2000", "292227", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292227", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2001", "292228", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292228", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2002", "292229", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292229", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2003", "292230", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292230", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2004", "292231", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292231", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2005", "292232", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292232", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2006", "292233", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292233", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2007", "292234", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292234", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2008", "292235", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292235", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2009", "292236", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292236", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2010", "292237", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292237", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2011", "292238", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292238", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2012", "292239", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292239", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2013", "292240", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292240", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2014", "292241", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292241", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2015", "292242", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292242", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2016", "292243", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292243", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2017", "292244", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292244", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2018", "292245", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292245", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2019", "292246", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292246", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2020", "292247", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292247", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2021", "292248", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292248", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2022", "292249", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292249", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2023", "292250", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292250", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2024", "292251", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292251", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2025", "292252", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292252", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2026", "292253", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292253", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2027", "292254", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292254", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2028", "292255", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292255", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2029", "292256", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292256", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2030", "292257", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292257", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2031", "292258", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292258", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2032", "292259", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292259", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2033", "292260", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292260", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2034", "292261", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292261", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2035", "292262", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292262", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2036", "292263", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292263", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2037", "292264", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292264", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2038", "292265", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292265", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2039", "292266", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292266", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2040", "292267", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292267", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2041", "292268", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292268", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2042", "292269", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292269", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2043", "292270", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292270", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2044", "292271", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292271", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2045", "292272", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292272", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2046", "292273", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292273", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2047", "292274", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292274", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2048", "292275", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292275", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2049", "292276", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292276", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2050", "292277", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292277", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2051", "292278", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292278", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2052", "292279", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292279", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2053", "292280", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292280", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2054", "292281", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292281", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2055", "292282", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292282", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2056", "292283", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292283", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2057", "292284", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292284", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2058", "292285", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292285", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2059", "292286", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292286", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2060", "292287", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292287", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2061", "292288", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292288", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2062", "292289", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292289", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2063", "292290", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292290", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2064", "292291", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292291", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2065", "292292", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292292", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2066", "292293", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292293", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2067", "292294", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292294", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2068", "292295", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292295", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2069", "292296", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292296", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2070", "292297", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292297", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2071", "292298", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292298", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2072", "292299", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292299", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2073", "292300", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292300", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2074", "292301", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292301", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2075", "292302", "PDF_variationNNPDF30_nlo_nf_5_pdfas_292302", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2076", "292600", "PDF_variation_292600", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2077", "315000", "PDF_variation_315000", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2078", "315200", "PDF_variation_315200", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2079", "262000", "PDF_variation_262000", 1.0, 1.0));
NLO26xInfo.push_back(weightinfo("2080", "263000", "PDF_variation_263000", 1.0, 1.0));
///End of v26x
}


//define this as a plug-in
DEFINE_FWK_MODULE(weight_assign_26xNLO);
