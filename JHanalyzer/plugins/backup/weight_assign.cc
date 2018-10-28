// -*- C++ -*-
//
// Package:    Analyzer/weight_assign
// Class:      weight_assign
// 
/**\class weight_assign weight_assign.cc Analyzer/weight_assign/plugins/weight_assign.cc

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


//
// Class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.



class weightinfo {

private:
  int _pdf;
  string _id;
  TString _name;
  double _muR,_muF,_weight;
public:
  weightinfo();
  weightinfo(string id, int pdf, TString name, double muR, double muF);
  void  set_id(string id);
  void set_pdf(int pdf);
  void set_name(TString name);
  void set_muR(double muR);
  void set_muF(double muF);
  void set_weight(double weight);
  string id();
  int pdf();
  TString name();
  double muR();
  double muF();
  double weight();
};
weightinfo::weightinfo(){
  _id="-1";
  _pdf=-1;
  _name="";
  _muR=0;
  _muF=0;
  _weight=0;
}
weightinfo::weightinfo(string id, int pdf, TString name, double muR, double muF){
  _id=id;
  _pdf=pdf;
  _name=name;
  _muR=muR;
  _muF=muF;

}
void weightinfo::set_id(string id){
  _id=id;
}
void weightinfo::set_pdf(int pdf){
  _pdf=pdf;
}
void weightinfo::set_name(TString name){
  _name=name;
}
void weightinfo::set_muR(double muR){
  _muR=muR;
}
void weightinfo::set_muF(double muF){
  _muF=muF;
}
void weightinfo::set_weight(double weight){
  _weight=weight;
}
string weightinfo::id(){
  return _id;
}
int weightinfo::pdf(){
  return _pdf;
}
TString weightinfo::name(){
  return _name;
}
double weightinfo::muR(){
  return _muR;
}
double weightinfo::muF(){
  return _muF;
}
double weightinfo::weight(){
  return _weight;
}
//class weight_assign : public edm::one::EDAnalyzer<edm::one::SharedResources>  {

class weight_assign : public edm::one::EDAnalyzer<edm::one::WatchRuns,edm::one::SharedResources>  {
//class weight_assign : public edm::EDAnalyzer  {
   public:
      explicit weight_assign(const edm::ParameterSet&);
      ~weight_assign();

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
  
  std::vector<weightinfo> LO242Info;
  std::vector<weightinfo> LO26xInfo;
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
weight_assign::weight_assign(const edm::ParameterSet& iConfig)

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
  // LO242Info.push_back(weightinfo(1080, 263000, "NNPDF30_lo_as_0130", 1, 1));

}


weight_assign::~weight_assign()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
weight_assign::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  cout<<"analyzer"<<endl;
  ////////////initialize/////////////
  //Get weight//
  edm::Handle<LHEEventProduct> LHEInfo;
  iEvent.getByToken(LHEInfo_Token, LHEInfo);
  int lheinfoweightsize= LHEInfo->weights().size();
  int lheinfocommentssize = LHEInfo->comments_size();
  cout<<"lheinfocommentssize="<<lheinfocommentssize<<endl;
  cout<<"lheinfoweightsize="<<lheinfoweightsize<<endl;
  double w0=LHEInfo->originalXWGTUP();
  int LO242Infosize=LO242Info.size();
  for (int i =0; i < lheinfoweightsize; i++){
    //cout<< LHEInfo->weights()[i].id<<endl;
    //    cout<< LHEInfo->weights()[i].wgt/w0<<endl;
  }

  for (int i242 =0; i242 < LO242Infosize;i242++){
    for (int i_lhe =0; i_lhe < lheinfoweightsize; i_lhe++){ 
    //   cout<<LO242Info[i].id()<<endl;
      //      int lheid=LHEInfo->weights()[i_lhe].id.Atoi()
      if(LHEInfo->weights()[i_lhe].id==LO242Info[i242].id()){
	//if(lheid==LO242Info[i242].id()){
	  LO242Info[i242].set_weight(LHEInfo->weights()[i_lhe].wgt/w0);
	  //cout<<"match"<<endl;
      }


    }
  }
  cout<<"LO242Infosize="<<LO242Infosize<<endl;
  for (int i242 =0; i242 < LO242Infosize;i242++){
        cout<<    LO242Info[i242].weight()<<endl;
    //   cout<<    LO242Info[i242].id()<<endl;
  }
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
   double weight=1;
   weight=genInfo->weight();
   cout<<"weight="<<weight<<endl;



   ///////Define Hard PID/////

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
weight_assign::beginJob()
//weight_assign::beginJob(edm::Run const& iRun)
//weight_assign::beginJob(edm::Run const& iRun, edm::EventSetup const &iSetup)
{
  cout<<"begin job"<<endl;
 
  edm::Service<TFileService> fs;

  /*
  edm::Handle<LHERunInfoProduct> run;
  typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;

  //  extLHEInfo_Token = consumes<LHERunInfoProduct>(edm::InputTag("externalLHEProducer"));                                                                                       

  // iRun.getByLabel( "externalLHEProducer", run );
  //iEvent.getByToken( extLHEInfo_Token, run );                                                                                                                                   
  //  iRun.getByToken( extLHEInfo_Token, run );                                                                                                                                     
  cout<<"LHERunInfo"<<endl;
  LHERunInfoProduct myLHERunInfoProduct = *(run.product());
  for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
    std::cout << iter->tag() << std::endl;
    std::vector<std::string> lines = iter->lines();
    for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
      std::cout << lines.at(iLine);
    }
  }

  */
  cout<<"end of beginjob"<<endl;
}

// ------------ method called once each job just after ending the event loop  ------------

void 
weight_assign::endJob() 
{
  cout<<"endjob"<<endl;
}



void 
weight_assign::beginRun(const Run &iEvent, EventSetup const &iSetup ){
  cout<<" beginrun"<<endl;
  cout<<"end of beginrun"<<endl;
}



void                                                                                                                                                                              
weight_assign::endRun(edm::Run const& iEvent, edm::EventSetup const&) 
{
  cout<<"doendrun"<<endl;

}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------


void
weight_assign::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//////defined by jhchoi/////
void 
weight_assign::make_weight_infos(){
  //LO242//
LO242Info.push_back(weightinfo("1", 247000, "muR_1_muF_1", 1, 1));
LO242Info.push_back(weightinfo("2", 247000, "muR_1_muF_2", 1, 2));
LO242Info.push_back(weightinfo("3", 247000, "muR_1_muF_0.5", 1, 0.5));
LO242Info.push_back(weightinfo("4", 247000, "muR_2_muF_1", 2, 1));
LO242Info.push_back(weightinfo("5", 247000, "muR_2_muF_2", 2, 2));
LO242Info.push_back(weightinfo("6", 247000, "muR_2_muF_0.5", 2, 0.5));
LO242Info.push_back(weightinfo("7", 247000, "muR_0.5_muF_1", 0.5, 1));
LO242Info.push_back(weightinfo("8", 247000, "muR_0.5_muF_2", 0.5, 2));
LO242Info.push_back(weightinfo("9", 247000, "muR_0.5_muF_0.5", 0.5, 0.5));
LO242Info.push_back(weightinfo("10", 306000, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("11", 306001, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("12", 306002, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("13", 306003, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("14", 306004, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("15", 306005, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("16", 306006, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("17", 306007, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("18", 306008, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("19", 306009, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("20", 306010, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("21", 306011, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("22", 306012, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("23", 306013, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("24", 306014, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("25", 306015, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("26", 306016, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("27", 306017, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("28", 306018, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("29", 306019, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("30", 306020, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("31", 306021, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("32", 306022, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("33", 306023, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("34", 306024, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("35", 306025, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("36", 306026, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("37", 306027, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("38", 306028, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("39", 306029, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("40", 306030, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("41", 306031, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("42", 306032, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("43", 306033, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("44", 306034, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("45", 306035, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("46", 306036, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("47", 306037, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("48", 306038, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("49", 306039, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("50", 306040, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("51", 306041, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("52", 306042, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("53", 306043, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("54", 306044, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("55", 306045, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("56", 306046, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("57", 306047, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("58", 306048, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("59", 306049, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("60", 306050, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("61", 306051, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("62", 306052, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("63", 306053, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("64", 306054, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("65", 306055, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("66", 306056, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("67", 306057, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("68", 306058, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("69", 306059, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("70", 306060, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("71", 306061, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("72", 306062, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("73", 306063, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("74", 306064, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("75", 306065, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("76", 306066, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("77", 306067, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("78", 306068, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("79", 306069, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("80", 306070, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("81", 306071, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("82", 306072, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("83", 306073, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("84", 306074, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("85", 306075, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("86", 306076, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("87", 306077, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("88", 306078, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("89", 306079, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("90", 306080, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("91", 306081, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("92", 306082, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("93", 306083, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("94", 306084, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("95", 306085, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("96", 306086, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("97", 306087, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("98", 306088, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("99", 306089, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("100", 306090, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("101", 306091, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("102", 306092, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("103", 306093, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("104", 306094, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("105", 306095, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("106", 306096, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("107", 306097, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("108", 306098, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("109", 306099, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("110", 306100, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("111", 306101, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("112", 306102, "NNPDF31_nnlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("113", 322500, "NNPDF31_nnlo_as_0108", 1, 1));
LO242Info.push_back(weightinfo("114", 322700, "NNPDF31_nnlo_as_0110", 1, 1));
LO242Info.push_back(weightinfo("115", 322900, "NNPDF31_nnlo_as_0112", 1, 1));
LO242Info.push_back(weightinfo("116", 323100, "NNPDF31_nnlo_as_0114", 1, 1));
LO242Info.push_back(weightinfo("117", 323300, "NNPDF31_nnlo_as_0117", 1, 1));
LO242Info.push_back(weightinfo("118", 323500, "NNPDF31_nnlo_as_0119", 1, 1));
LO242Info.push_back(weightinfo("119", 323700, "NNPDF31_nnlo_as_0122", 1, 1));
LO242Info.push_back(weightinfo("120", 323900, "NNPDF31_nnlo_as_0124", 1, 1));
LO242Info.push_back(weightinfo("121", 305800, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("122", 305801, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("123", 305802, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("124", 305803, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("125", 305804, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("126", 305805, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("127", 305806, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("128", 305807, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("129", 305808, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("130", 305809, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("131", 305810, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("132", 305811, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("133", 305812, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("134", 305813, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("135", 305814, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("136", 305815, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("137", 305816, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("138", 305817, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("139", 305818, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("140", 305819, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("141", 305820, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("142", 305821, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("143", 305822, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("144", 305823, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("145", 305824, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("146", 305825, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("147", 305826, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("148", 305827, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("149", 305828, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("150", 305829, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("151", 305830, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("152", 305831, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("153", 305832, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("154", 305833, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("155", 305834, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("156", 305835, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("157", 305836, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("158", 305837, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("159", 305838, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("160", 305839, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("161", 305840, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("162", 305841, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("163", 305842, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("164", 305843, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("165", 305844, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("166", 305845, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("167", 305846, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("168", 305847, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("169", 305848, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("170", 305849, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("171", 305850, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("172", 305851, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("173", 305852, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("174", 305853, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("175", 305854, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("176", 305855, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("177", 305856, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("178", 305857, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("179", 305858, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("180", 305859, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("181", 305860, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("182", 305861, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("183", 305862, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("184", 305863, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("185", 305864, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("186", 305865, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("187", 305866, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("188", 305867, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("189", 305868, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("190", 305869, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("191", 305870, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("192", 305871, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("193", 305872, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("194", 305873, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("195", 305874, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("196", 305875, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("197", 305876, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("198", 305877, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("199", 305878, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("200", 305879, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("201", 305880, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("202", 305881, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("203", 305882, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("204", 305883, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("205", 305884, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("206", 305885, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("207", 305886, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("208", 305887, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("209", 305888, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("210", 305889, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("211", 305890, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("212", 305891, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("213", 305892, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("214", 305893, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("215", 305894, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("216", 305895, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("217", 305896, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("218", 305897, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("219", 305898, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("220", 305899, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("221", 305900, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("222", 305901, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("223", 305902, "NNPDF31_nlo_hessian_pdfas", 1, 1));
LO242Info.push_back(weightinfo("224", 13000, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("225", 13001, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("226", 13002, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("227", 13003, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("228", 13004, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("229", 13005, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("230", 13006, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("231", 13007, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("232", 13008, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("233", 13009, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("234", 13010, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("235", 13011, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("236", 13012, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("237", 13013, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("238", 13014, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("239", 13015, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("240", 13016, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("241", 13017, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("242", 13018, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("243", 13019, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("244", 13020, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("245", 13021, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("246", 13022, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("247", 13023, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("248", 13024, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("249", 13025, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("250", 13026, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("251", 13027, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("252", 13028, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("253", 13029, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("254", 13030, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("255", 13031, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("256", 13032, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("257", 13033, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("258", 13034, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("259", 13035, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("260", 13036, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("261", 13037, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("262", 13038, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("263", 13039, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("264", 13040, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("265", 13041, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("266", 13042, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("267", 13043, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("268", 13044, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("269", 13045, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("270", 13046, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("271", 13047, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("272", 13048, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("273", 13049, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("274", 13050, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("275", 13051, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("276", 13052, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("277", 13053, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("278", 13054, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("279", 13055, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("280", 13056, "CT14nnlo", 1, 1));
LO242Info.push_back(weightinfo("281", 13065, "CT14nnlo_as_0116", 1, 1));
LO242Info.push_back(weightinfo("282", 13069, "CT14nnlo_as_0120", 1, 1));
LO242Info.push_back(weightinfo("283", 13100, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("284", 13101, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("285", 13102, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("286", 13103, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("287", 13104, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("288", 13105, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("289", 13106, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("290", 13107, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("291", 13108, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("292", 13109, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("293", 13110, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("294", 13111, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("295", 13112, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("296", 13113, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("297", 13114, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("298", 13115, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("299", 13116, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("300", 13117, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("301", 13118, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("302", 13119, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("303", 13120, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("304", 13121, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("305", 13122, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("306", 13123, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("307", 13124, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("308", 13125, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("309", 13126, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("310", 13127, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("311", 13128, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("312", 13129, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("313", 13130, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("314", 13131, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("315", 13132, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("316", 13133, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("317", 13134, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("318", 13135, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("319", 13136, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("320", 13137, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("321", 13138, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("322", 13139, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("323", 13140, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("324", 13141, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("325", 13142, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("326", 13143, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("327", 13144, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("328", 13145, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("329", 13146, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("330", 13147, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("331", 13148, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("332", 13149, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("333", 13150, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("334", 13151, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("335", 13152, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("336", 13153, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("337", 13154, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("338", 13155, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("339", 13156, "CT14nlo", 1, 1));
LO242Info.push_back(weightinfo("340", 13163, "CT14nlo_as_0116", 1, 1));
LO242Info.push_back(weightinfo("341", 13167, "CT14nlo_as_0120", 1, 1));
LO242Info.push_back(weightinfo("342", 13200, "CT14lo", 1, 1));
LO242Info.push_back(weightinfo("343", 25200, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("344", 25201, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("345", 25202, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("346", 25203, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("347", 25204, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("348", 25205, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("349", 25206, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("350", 25207, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("351", 25208, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("352", 25209, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("353", 25210, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("354", 25211, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("355", 25212, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("356", 25213, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("357", 25214, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("358", 25215, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("359", 25216, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("360", 25217, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("361", 25218, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("362", 25219, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("363", 25220, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("364", 25221, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("365", 25222, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("366", 25223, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("367", 25224, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("368", 25225, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("369", 25226, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("370", 25227, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("371", 25228, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("372", 25229, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("373", 25230, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("374", 25231, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("375", 25232, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("376", 25233, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("377", 25234, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("378", 25235, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("379", 25236, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("380", 25237, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("381", 25238, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("382", 25239, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("383", 25240, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("384", 25241, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("385", 25242, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("386", 25243, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("387", 25244, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("388", 25245, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("389", 25246, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("390", 25247, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("391", 25248, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("392", 25249, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("393", 25250, "MMHT2014nlo68clas118", 1, 1));
LO242Info.push_back(weightinfo("394", 25300, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("395", 25301, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("396", 25302, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("397", 25303, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("398", 25304, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("399", 25305, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("400", 25306, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("401", 25307, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("402", 25308, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("403", 25309, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("404", 25310, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("405", 25311, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("406", 25312, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("407", 25313, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("408", 25314, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("409", 25315, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("410", 25316, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("411", 25317, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("412", 25318, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("413", 25319, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("414", 25320, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("415", 25321, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("416", 25322, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("417", 25323, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("418", 25324, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("419", 25325, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("420", 25326, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("421", 25327, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("422", 25328, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("423", 25329, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("424", 25330, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("425", 25331, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("426", 25332, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("427", 25333, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("428", 25334, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("429", 25335, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("430", 25336, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("431", 25337, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("432", 25338, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("433", 25339, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("434", 25340, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("435", 25341, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("436", 25342, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("437", 25343, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("438", 25344, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("439", 25345, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("440", 25346, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("441", 25347, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("442", 25348, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("443", 25349, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("444", 25350, "MMHT2014nnlo68cl", 1, 1));
LO242Info.push_back(weightinfo("445", 25000, "MMHT2014lo68cl", 1, 1));
LO242Info.push_back(weightinfo("446", 42780, "ABMP16als118_5_nnlo", 1, 1));
LO242Info.push_back(weightinfo("447", 42781, "ABMP16als118_5_nnlo", 1, 1));
LO242Info.push_back(weightinfo("448", 42782, "ABMP16als118_5_nnlo", 1, 1));
LO242Info.push_back(weightinfo("449", 42783, "ABMP16als118_5_nnlo", 1, 1));
LO242Info.push_back(weightinfo("450", 42784, "ABMP16als118_5_nnlo", 1, 1));
LO242Info.push_back(weightinfo("451", 42785, "ABMP16als118_5_nnlo", 1, 1));
LO242Info.push_back(weightinfo("452", 42786, "ABMP16als118_5_nnlo", 1, 1));
LO242Info.push_back(weightinfo("453", 42787, "ABMP16als118_5_nnlo", 1, 1));
LO242Info.push_back(weightinfo("454", 42788, "ABMP16als118_5_nnlo", 1, 1));
LO242Info.push_back(weightinfo("455", 42789, "ABMP16als118_5_nnlo", 1, 1));
LO242Info.push_back(weightinfo("456", 42790, "ABMP16als118_5_nnlo", 1, 1));
LO242Info.push_back(weightinfo("457", 42791, "ABMP16als118_5_nnlo", 1, 1));
LO242Info.push_back(weightinfo("458", 42792, "ABMP16als118_5_nnlo", 1, 1));
LO242Info.push_back(weightinfo("459", 42793, "ABMP16als118_5_nnlo", 1, 1));
LO242Info.push_back(weightinfo("460", 42794, "ABMP16als118_5_nnlo", 1, 1));
LO242Info.push_back(weightinfo("461", 42795, "ABMP16als118_5_nnlo", 1, 1));
LO242Info.push_back(weightinfo("462", 42796, "ABMP16als118_5_nnlo", 1, 1));
LO242Info.push_back(weightinfo("463", 42797, "ABMP16als118_5_nnlo", 1, 1));
LO242Info.push_back(weightinfo("464", 42798, "ABMP16als118_5_nnlo", 1, 1));
LO242Info.push_back(weightinfo("465", 42799, "ABMP16als118_5_nnlo", 1, 1));
LO242Info.push_back(weightinfo("466", 42800, "ABMP16als118_5_nnlo", 1, 1));
LO242Info.push_back(weightinfo("467", 42801, "ABMP16als118_5_nnlo", 1, 1));
LO242Info.push_back(weightinfo("468", 42802, "ABMP16als118_5_nnlo", 1, 1));
LO242Info.push_back(weightinfo("469", 42803, "ABMP16als118_5_nnlo", 1, 1));
LO242Info.push_back(weightinfo("470", 42804, "ABMP16als118_5_nnlo", 1, 1));
LO242Info.push_back(weightinfo("471", 42805, "ABMP16als118_5_nnlo", 1, 1));
LO242Info.push_back(weightinfo("472", 42806, "ABMP16als118_5_nnlo", 1, 1));
LO242Info.push_back(weightinfo("473", 42807, "ABMP16als118_5_nnlo", 1, 1));
LO242Info.push_back(weightinfo("474", 42808, "ABMP16als118_5_nnlo", 1, 1));
LO242Info.push_back(weightinfo("475", 42809, "ABMP16als118_5_nnlo", 1, 1));
LO242Info.push_back(weightinfo("476", 90200, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("477", 90201, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("478", 90202, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("479", 90203, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("480", 90204, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("481", 90205, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("482", 90206, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("483", 90207, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("484", 90208, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("485", 90209, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("486", 90210, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("487", 90211, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("488", 90212, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("489", 90213, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("490", 90214, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("491", 90215, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("492", 90216, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("493", 90217, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("494", 90218, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("495", 90219, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("496", 90220, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("497", 90221, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("498", 90222, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("499", 90223, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("500", 90224, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("501", 90225, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("502", 90226, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("503", 90227, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("504", 90228, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("505", 90229, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("506", 90230, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("507", 90231, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("508", 90232, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("509", 90233, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("510", 90234, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("511", 90235, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("512", 90236, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("513", 90237, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("514", 90238, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("515", 90239, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("516", 90240, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("517", 90241, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("518", 90242, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("519", 90243, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("520", 90244, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("521", 90245, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("522", 90246, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("523", 90247, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("524", 90248, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("525", 90249, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("526", 90250, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("527", 90251, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("528", 90252, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("529", 90253, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("530", 90254, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("531", 90255, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("532", 90256, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("533", 90257, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("534", 90258, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("535", 90259, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("536", 90260, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("537", 90261, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("538", 90262, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("539", 90263, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("540", 90264, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("541", 90265, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("542", 90266, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("543", 90267, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("544", 90268, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("545", 90269, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("546", 90270, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("547", 90271, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("548", 90272, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("549", 90273, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("550", 90274, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("551", 90275, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("552", 90276, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("553", 90277, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("554", 90278, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("555", 90279, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("556", 90280, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("557", 90281, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("558", 90282, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("559", 90283, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("560", 90284, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("561", 90285, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("562", 90286, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("563", 90287, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("564", 90288, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("565", 90289, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("566", 90290, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("567", 90291, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("568", 90292, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("569", 90293, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("570", 90294, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("571", 90295, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("572", 90296, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("573", 90297, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("574", 90298, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("575", 90299, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("576", 90300, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("577", 90301, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("578", 90302, "PDF4LHC15_nlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("579", 91200, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("580", 91201, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("581", 91202, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("582", 91203, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("583", 91204, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("584", 91205, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("585", 91206, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("586", 91207, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("587", 91208, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("588", 91209, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("589", 91210, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("590", 91211, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("591", 91212, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("592", 91213, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("593", 91214, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("594", 91215, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("595", 91216, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("596", 91217, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("597", 91218, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("598", 91219, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("599", 91220, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("600", 91221, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("601", 91222, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("602", 91223, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("603", 91224, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("604", 91225, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("605", 91226, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("606", 91227, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("607", 91228, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("608", 91229, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("609", 91230, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("610", 91231, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("611", 91232, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("612", 91233, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("613", 91234, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("614", 91235, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("615", 91236, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("616", 91237, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("617", 91238, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("618", 91239, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("619", 91240, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("620", 91241, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("621", 91242, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("622", 91243, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("623", 91244, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("624", 91245, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("625", 91246, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("626", 91247, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("627", 91248, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("628", 91249, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("629", 91250, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("630", 91251, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("631", 91252, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("632", 91253, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("633", 91254, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("634", 91255, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("635", 91256, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("636", 91257, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("637", 91258, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("638", 91259, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("639", 91260, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("640", 91261, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("641", 91262, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("642", 91263, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("643", 91264, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("644", 91265, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("645", 91266, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("646", 91267, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("647", 91268, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("648", 91269, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("649", 91270, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("650", 91271, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("651", 91272, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("652", 91273, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("653", 91274, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("654", 91275, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("655", 91276, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("656", 91277, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("657", 91278, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("658", 91279, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("659", 91280, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("660", 91281, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("661", 91282, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("662", 91283, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("663", 91284, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("664", 91285, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("665", 91286, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("666", 91287, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("667", 91288, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("668", 91289, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("669", 91290, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("670", 91291, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("671", 91292, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("672", 91293, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("673", 91294, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("674", 91295, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("675", 91296, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("676", 91297, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("677", 91298, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("678", 91299, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("679", 91300, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("680", 91301, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("681", 91302, "PDF4LHC15_nnlo_100_pdfas", 1, 1));
LO242Info.push_back(weightinfo("682", 90400, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("683", 90401, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("684", 90402, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("685", 90403, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("686", 90404, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("687", 90405, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("688", 90406, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("689", 90407, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("690", 90408, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("691", 90409, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("692", 90410, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("693", 90411, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("694", 90412, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("695", 90413, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("696", 90414, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("697", 90415, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("698", 90416, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("699", 90417, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("700", 90418, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("701", 90419, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("702", 90420, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("703", 90421, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("704", 90422, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("705", 90423, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("706", 90424, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("707", 90425, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("708", 90426, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("709", 90427, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("710", 90428, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("711", 90429, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("712", 90430, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("713", 90431, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("714", 90432, "PDF4LHC15_nlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("715", 91400, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("716", 91401, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("717", 91402, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("718", 91403, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("719", 91404, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("720", 91405, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("721", 91406, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("722", 91407, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("723", 91408, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("724", 91409, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("725", 91410, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("726", 91411, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("727", 91412, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("728", 91413, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("729", 91414, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("730", 91415, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("731", 91416, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("732", 91417, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("733", 91418, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("734", 91419, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("735", 91420, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("736", 91421, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("737", 91422, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("738", 91423, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("739", 91424, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("740", 91425, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("741", 91426, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("742", 91427, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("743", 91428, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("744", 91429, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("745", 91430, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("746", 91431, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("747", 91432, "PDF4LHC15_nnlo_30_pdfas", 1, 1));
LO242Info.push_back(weightinfo("748", 61100, "HERAPDF20_NLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("749", 61101, "HERAPDF20_NLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("750", 61102, "HERAPDF20_NLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("751", 61103, "HERAPDF20_NLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("752", 61104, "HERAPDF20_NLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("753", 61105, "HERAPDF20_NLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("754", 61106, "HERAPDF20_NLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("755", 61107, "HERAPDF20_NLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("756", 61108, "HERAPDF20_NLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("757", 61109, "HERAPDF20_NLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("758", 61110, "HERAPDF20_NLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("759", 61111, "HERAPDF20_NLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("760", 61112, "HERAPDF20_NLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("761", 61113, "HERAPDF20_NLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("762", 61114, "HERAPDF20_NLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("763", 61115, "HERAPDF20_NLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("764", 61116, "HERAPDF20_NLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("765", 61117, "HERAPDF20_NLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("766", 61118, "HERAPDF20_NLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("767", 61119, "HERAPDF20_NLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("768", 61120, "HERAPDF20_NLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("769", 61121, "HERAPDF20_NLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("770", 61122, "HERAPDF20_NLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("771", 61123, "HERAPDF20_NLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("772", 61124, "HERAPDF20_NLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("773", 61125, "HERAPDF20_NLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("774", 61126, "HERAPDF20_NLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("775", 61127, "HERAPDF20_NLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("776", 61128, "HERAPDF20_NLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("777", 61130, "HERAPDF20_NLO_VAR", 1, 1));
LO242Info.push_back(weightinfo("778", 61131, "HERAPDF20_NLO_VAR", 1, 1));
LO242Info.push_back(weightinfo("779", 61132, "HERAPDF20_NLO_VAR", 1, 1));
LO242Info.push_back(weightinfo("780", 61133, "HERAPDF20_NLO_VAR", 1, 1));
LO242Info.push_back(weightinfo("781", 61134, "HERAPDF20_NLO_VAR", 1, 1));
LO242Info.push_back(weightinfo("782", 61135, "HERAPDF20_NLO_VAR", 1, 1));
LO242Info.push_back(weightinfo("783", 61136, "HERAPDF20_NLO_VAR", 1, 1));
LO242Info.push_back(weightinfo("784", 61137, "HERAPDF20_NLO_VAR", 1, 1));
LO242Info.push_back(weightinfo("785", 61138, "HERAPDF20_NLO_VAR", 1, 1));
LO242Info.push_back(weightinfo("786", 61139, "HERAPDF20_NLO_VAR", 1, 1));
LO242Info.push_back(weightinfo("787", 61140, "HERAPDF20_NLO_VAR", 1, 1));
LO242Info.push_back(weightinfo("788", 61141, "HERAPDF20_NLO_VAR", 1, 1));
LO242Info.push_back(weightinfo("789", 61142, "HERAPDF20_NLO_VAR", 1, 1));
LO242Info.push_back(weightinfo("790", 61143, "HERAPDF20_NLO_VAR", 1, 1));
LO242Info.push_back(weightinfo("791", 61200, "HERAPDF20_NNLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("792", 61201, "HERAPDF20_NNLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("793", 61202, "HERAPDF20_NNLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("794", 61203, "HERAPDF20_NNLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("795", 61204, "HERAPDF20_NNLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("796", 61205, "HERAPDF20_NNLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("797", 61206, "HERAPDF20_NNLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("798", 61207, "HERAPDF20_NNLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("799", 61208, "HERAPDF20_NNLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("800", 61209, "HERAPDF20_NNLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("801", 61210, "HERAPDF20_NNLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("802", 61211, "HERAPDF20_NNLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("803", 61212, "HERAPDF20_NNLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("804", 61213, "HERAPDF20_NNLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("805", 61214, "HERAPDF20_NNLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("806", 61215, "HERAPDF20_NNLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("807", 61216, "HERAPDF20_NNLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("808", 61217, "HERAPDF20_NNLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("809", 61218, "HERAPDF20_NNLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("810", 61219, "HERAPDF20_NNLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("811", 61220, "HERAPDF20_NNLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("812", 61221, "HERAPDF20_NNLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("813", 61222, "HERAPDF20_NNLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("814", 61223, "HERAPDF20_NNLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("815", 61224, "HERAPDF20_NNLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("816", 61225, "HERAPDF20_NNLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("817", 61226, "HERAPDF20_NNLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("818", 61227, "HERAPDF20_NNLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("819", 61228, "HERAPDF20_NNLO_EIG", 1, 1));
LO242Info.push_back(weightinfo("820", 61230, "HERAPDF20_NNLO_VAR", 1, 1));
LO242Info.push_back(weightinfo("821", 61231, "HERAPDF20_NNLO_VAR", 1, 1));
LO242Info.push_back(weightinfo("822", 61232, "HERAPDF20_NNLO_VAR", 1, 1));
LO242Info.push_back(weightinfo("823", 61233, "HERAPDF20_NNLO_VAR", 1, 1));
LO242Info.push_back(weightinfo("824", 61234, "HERAPDF20_NNLO_VAR", 1, 1));
LO242Info.push_back(weightinfo("825", 61235, "HERAPDF20_NNLO_VAR", 1, 1));
LO242Info.push_back(weightinfo("826", 61236, "HERAPDF20_NNLO_VAR", 1, 1));
LO242Info.push_back(weightinfo("827", 61237, "HERAPDF20_NNLO_VAR", 1, 1));
LO242Info.push_back(weightinfo("828", 61238, "HERAPDF20_NNLO_VAR", 1, 1));
LO242Info.push_back(weightinfo("829", 61239, "HERAPDF20_NNLO_VAR", 1, 1));
LO242Info.push_back(weightinfo("830", 61240, "HERAPDF20_NNLO_VAR", 1, 1));
LO242Info.push_back(weightinfo("831", 61241, "HERAPDF20_NNLO_VAR", 1, 1));
LO242Info.push_back(weightinfo("832", 61242, "HERAPDF20_NNLO_VAR", 1, 1));
LO242Info.push_back(weightinfo("833", 61243, "HERAPDF20_NNLO_VAR", 1, 1));
LO242Info.push_back(weightinfo("834", 13400, "CT14qed_inc_proton", 1, 1));
LO242Info.push_back(weightinfo("835", 13401, "CT14qed_inc_proton", 1, 1));
LO242Info.push_back(weightinfo("836", 13402, "CT14qed_inc_proton", 1, 1));
LO242Info.push_back(weightinfo("837", 13403, "CT14qed_inc_proton", 1, 1));
LO242Info.push_back(weightinfo("838", 13404, "CT14qed_inc_proton", 1, 1));
LO242Info.push_back(weightinfo("839", 13405, "CT14qed_inc_proton", 1, 1));
LO242Info.push_back(weightinfo("840", 13406, "CT14qed_inc_proton", 1, 1));
LO242Info.push_back(weightinfo("841", 13407, "CT14qed_inc_proton", 1, 1));
LO242Info.push_back(weightinfo("842", 13408, "CT14qed_inc_proton", 1, 1));
LO242Info.push_back(weightinfo("843", 13409, "CT14qed_inc_proton", 1, 1));
LO242Info.push_back(weightinfo("844", 13410, "CT14qed_inc_proton", 1, 1));
LO242Info.push_back(weightinfo("845", 13411, "CT14qed_inc_proton", 1, 1));
LO242Info.push_back(weightinfo("846", 13412, "CT14qed_inc_proton", 1, 1));
LO242Info.push_back(weightinfo("847", 13413, "CT14qed_inc_proton", 1, 1));
LO242Info.push_back(weightinfo("848", 13414, "CT14qed_inc_proton", 1, 1));
LO242Info.push_back(weightinfo("849", 13415, "CT14qed_inc_proton", 1, 1));
LO242Info.push_back(weightinfo("850", 13416, "CT14qed_inc_proton", 1, 1));
LO242Info.push_back(weightinfo("851", 13417, "CT14qed_inc_proton", 1, 1));
LO242Info.push_back(weightinfo("852", 13418, "CT14qed_inc_proton", 1, 1));
LO242Info.push_back(weightinfo("853", 13419, "CT14qed_inc_proton", 1, 1));
LO242Info.push_back(weightinfo("854", 13420, "CT14qed_inc_proton", 1, 1));
LO242Info.push_back(weightinfo("855", 13421, "CT14qed_inc_proton", 1, 1));
LO242Info.push_back(weightinfo("856", 13422, "CT14qed_inc_proton", 1, 1));
LO242Info.push_back(weightinfo("857", 13423, "CT14qed_inc_proton", 1, 1));
LO242Info.push_back(weightinfo("858", 13424, "CT14qed_inc_proton", 1, 1));
LO242Info.push_back(weightinfo("859", 13425, "CT14qed_inc_proton", 1, 1));
LO242Info.push_back(weightinfo("860", 13426, "CT14qed_inc_proton", 1, 1));
LO242Info.push_back(weightinfo("861", 13427, "CT14qed_inc_proton", 1, 1));
LO242Info.push_back(weightinfo("862", 13428, "CT14qed_inc_proton", 1, 1));
LO242Info.push_back(weightinfo("863", 13429, "CT14qed_inc_proton", 1, 1));
LO242Info.push_back(weightinfo("864", 13430, "CT14qed_inc_proton", 1, 1));
LO242Info.push_back(weightinfo("865", 82200, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("866", 82201, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("867", 82202, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("868", 82203, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("869", 82204, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("870", 82205, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("871", 82206, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("872", 82207, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("873", 82208, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("874", 82209, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("875", 82210, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("876", 82211, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("877", 82212, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("878", 82213, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("879", 82214, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("880", 82215, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("881", 82216, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("882", 82217, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("883", 82218, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("884", 82219, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("885", 82220, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("886", 82221, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("887", 82222, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("888", 82223, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("889", 82224, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("890", 82225, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("891", 82226, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("892", 82227, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("893", 82228, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("894", 82229, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("895", 82230, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("896", 82231, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("897", 82232, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("898", 82233, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("899", 82234, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("900", 82235, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("901", 82236, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("902", 82237, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("903", 82238, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("904", 82239, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("905", 82240, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("906", 82241, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("907", 82242, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("908", 82243, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("909", 82244, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("910", 82245, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("911", 82246, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("912", 82247, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("913", 82248, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("914", 82249, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("915", 82250, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("916", 82251, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("917", 82252, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("918", 82253, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("919", 82254, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("920", 82255, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("921", 82256, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("922", 82257, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("923", 82258, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("924", 82259, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("925", 82260, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("926", 82261, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("927", 82262, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("928", 82263, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("929", 82264, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("930", 82265, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("931", 82266, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("932", 82267, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("933", 82268, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("934", 82269, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("935", 82270, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("936", 82271, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("937", 82272, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("938", 82273, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("939", 82274, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("940", 82275, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("941", 82276, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("942", 82277, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("943", 82278, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("944", 82279, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("945", 82280, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("946", 82281, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("947", 82282, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("948", 82283, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("949", 82284, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("950", 82285, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("951", 82286, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("952", 82287, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("953", 82288, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("954", 82289, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("955", 82290, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("956", 82291, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("957", 82292, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("958", 82293, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("959", 82294, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("960", 82295, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("961", 82296, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("962", 82297, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("963", 82298, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("964", 82299, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("965", 82300, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("966", 82301, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("967", 82302, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("968", 82303, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("969", 82304, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("970", 82305, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("971", 82306, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("972", 82307, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1, 1));
LO242Info.push_back(weightinfo("973", 292200, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("974", 292201, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("975", 292202, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("976", 292203, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("977", 292204, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("978", 292205, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("979", 292206, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("980", 292207, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("981", 292208, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("982", 292209, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("983", 292210, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("984", 292211, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("985", 292212, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("986", 292213, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("987", 292214, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("988", 292215, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("989", 292216, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("990", 292217, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("991", 292218, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("992", 292219, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("993", 292220, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("994", 292221, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("995", 292222, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("996", 292223, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("997", 292224, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("998", 292225, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("999", 292226, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1000", 292227, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1001", 292228, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1002", 292229, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1003", 292230, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1004", 292231, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1005", 292232, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1006", 292233, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1007", 292234, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1008", 292235, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1009", 292236, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1010", 292237, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1011", 292238, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1012", 292239, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1013", 292240, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1014", 292241, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1015", 292242, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1016", 292243, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1017", 292244, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1018", 292245, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1019", 292246, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1020", 292247, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1021", 292248, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1022", 292249, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1023", 292250, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1024", 292251, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1025", 292252, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1026", 292253, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1027", 292254, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1028", 292255, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1029", 292256, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1030", 292257, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1031", 292258, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1032", 292259, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1033", 292260, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1034", 292261, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1035", 292262, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1036", 292263, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1037", 292264, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1038", 292265, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1039", 292266, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1040", 292267, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1041", 292268, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1042", 292269, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1043", 292270, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1044", 292271, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1045", 292272, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1046", 292273, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1047", 292274, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1048", 292275, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1049", 292276, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1050", 292277, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1051", 292278, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1052", 292279, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1053", 292280, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1054", 292281, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1055", 292282, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1056", 292283, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1057", 292284, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1058", 292285, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1059", 292286, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1060", 292287, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1061", 292288, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1062", 292289, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1063", 292290, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1064", 292291, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1065", 292292, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1066", 292293, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1067", 292294, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1068", 292295, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1069", 292296, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1070", 292297, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1071", 292298, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1072", 292299, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1073", 292300, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1074", 292301, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1075", 292302, "NNPDF30_nlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1076", 292600, "NNPDF30_nnlo_nf_5_pdfas", 1, 1));
LO242Info.push_back(weightinfo("1077", 315000, "NNPDF31_lo_as_0118", 1, 1));
LO242Info.push_back(weightinfo("1078", 315200, "NNPDF31_lo_as_0130", 1, 1));
LO242Info.push_back(weightinfo("1079", 262000, "NNPDF30_lo_as_0118", 1, 1));
LO242Info.push_back(weightinfo("1080", 263000, "NNPDF30_lo_as_0130", 1, 1));
/////End of LO242 setup

///LO26x///
LO26xInfo.push_back(weightinfo("1001", 247000, "muR_1.0_muF_1.0", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1006", 247000, "muR_2.0_muF_1.0", 2.0, 1.0));
LO26xInfo.push_back(weightinfo("1011", 247000, "muR_0.5_muF_1.0", 0.5, 1.0));
LO26xInfo.push_back(weightinfo("1016", 247000, "muR_1.0_muF_2.0", 1.0, 2.0));
LO26xInfo.push_back(weightinfo("1021", 247000, "muR_2.0_muF_2.0", 2.0, 2.0));
LO26xInfo.push_back(weightinfo("1026", 247000, "muR_0.5_muF_2.0", 0.5, 2.0));
LO26xInfo.push_back(weightinfo("1031", 247000, "muR_1.0_muF_0.5", 1.0, 0.5));
LO26xInfo.push_back(weightinfo("1036", 247000, "muR_2.0_muF_0.5", 2.0, 0.5));
LO26xInfo.push_back(weightinfo("1041", 247000, "muR_0.5_muF_0.5", 0.5, 0.5));
LO26xInfo.push_back(weightinfo("1046", 306000, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1047", 306001, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1048", 306002, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1049", 306003, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1050", 306004, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1051", 306005, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1052", 306006, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1053", 306007, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1054", 306008, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1055", 306009, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1056", 306010, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1057", 306011, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1058", 306012, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1059", 306013, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1060", 306014, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1061", 306015, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1062", 306016, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1063", 306017, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1064", 306018, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1065", 306019, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1066", 306020, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1067", 306021, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1068", 306022, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1069", 306023, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1070", 306024, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1071", 306025, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1072", 306026, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1073", 306027, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1074", 306028, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1075", 306029, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1076", 306030, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1077", 306031, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1078", 306032, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1079", 306033, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1080", 306034, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1081", 306035, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1082", 306036, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1083", 306037, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1084", 306038, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1085", 306039, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1086", 306040, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1087", 306041, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1088", 306042, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1089", 306043, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1090", 306044, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1091", 306045, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1092", 306046, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1093", 306047, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1094", 306048, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1095", 306049, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1096", 306050, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1097", 306051, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1098", 306052, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1099", 306053, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1100", 306054, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1101", 306055, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1102", 306056, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1103", 306057, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1104", 306058, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1105", 306059, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1106", 306060, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1107", 306061, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1108", 306062, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1109", 306063, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1110", 306064, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1111", 306065, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1112", 306066, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1113", 306067, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1114", 306068, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1115", 306069, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1116", 306070, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1117", 306071, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1118", 306072, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1119", 306073, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1120", 306074, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1121", 306075, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1122", 306076, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1123", 306077, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1124", 306078, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1125", 306079, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1126", 306080, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1127", 306081, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1128", 306082, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1129", 306083, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1130", 306084, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1131", 306085, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1132", 306086, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1133", 306087, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1134", 306088, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1135", 306089, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1136", 306090, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1137", 306091, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1138", 306092, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1139", 306093, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1140", 306094, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1141", 306095, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1142", 306096, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1143", 306097, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1144", 306098, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1145", 306099, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1146", 306100, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1147", 306101, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1148", 306102, "NNPDF31_nnlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1149", 322500, "NNPDF31_nnlo_as_0108", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1150", 322700, "NNPDF31_nnlo_as_0110", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1151", 322900, "NNPDF31_nnlo_as_0112", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1152", 323100, "NNPDF31_nnlo_as_0114", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1153", 323300, "NNPDF31_nnlo_as_0117", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1154", 323500, "NNPDF31_nnlo_as_0119", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1155", 323700, "NNPDF31_nnlo_as_0122", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1156", 323900, "NNPDF31_nnlo_as_0124", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1157", 305800, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1158", 305801, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1159", 305802, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1160", 305803, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1161", 305804, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1162", 305805, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1163", 305806, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1164", 305807, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1165", 305808, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1166", 305809, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1167", 305810, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1168", 305811, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1169", 305812, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1170", 305813, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1171", 305814, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1172", 305815, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1173", 305816, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1174", 305817, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1175", 305818, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1176", 305819, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1177", 305820, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1178", 305821, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1179", 305822, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1180", 305823, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1181", 305824, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1182", 305825, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1183", 305826, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1184", 305827, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1185", 305828, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1186", 305829, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1187", 305830, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1188", 305831, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1189", 305832, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1190", 305833, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1191", 305834, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1192", 305835, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1193", 305836, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1194", 305837, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1195", 305838, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1196", 305839, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1197", 305840, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1198", 305841, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1199", 305842, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1200", 305843, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1201", 305844, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1202", 305845, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1203", 305846, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1204", 305847, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1205", 305848, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1206", 305849, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1207", 305850, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1208", 305851, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1209", 305852, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1210", 305853, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1211", 305854, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1212", 305855, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1213", 305856, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1214", 305857, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1215", 305858, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1216", 305859, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1217", 305860, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1218", 305861, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1219", 305862, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1220", 305863, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1221", 305864, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1222", 305865, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1223", 305866, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1224", 305867, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1225", 305868, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1226", 305869, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1227", 305870, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1228", 305871, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1229", 305872, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1230", 305873, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1231", 305874, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1232", 305875, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1233", 305876, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1234", 305877, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1235", 305878, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1236", 305879, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1237", 305880, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1238", 305881, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1239", 305882, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1240", 305883, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1241", 305884, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1242", 305885, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1243", 305886, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1244", 305887, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1245", 305888, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1246", 305889, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1247", 305890, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1248", 305891, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1249", 305892, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1250", 305893, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1251", 305894, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1252", 305895, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1253", 305896, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1254", 305897, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1255", 305898, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1256", 305899, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1257", 305900, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1258", 305901, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1259", 305902, "NNPDF31_nlo_hessian_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1260", 13000, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1261", 13001, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1262", 13002, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1263", 13003, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1264", 13004, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1265", 13005, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1266", 13006, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1267", 13007, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1268", 13008, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1269", 13009, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1270", 13010, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1271", 13011, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1272", 13012, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1273", 13013, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1274", 13014, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1275", 13015, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1276", 13016, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1277", 13017, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1278", 13018, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1279", 13019, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1280", 13020, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1281", 13021, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1282", 13022, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1283", 13023, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1284", 13024, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1285", 13025, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1286", 13026, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1287", 13027, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1288", 13028, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1289", 13029, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1290", 13030, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1291", 13031, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1292", 13032, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1293", 13033, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1294", 13034, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1295", 13035, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1296", 13036, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1297", 13037, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1298", 13038, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1299", 13039, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1300", 13040, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1301", 13041, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1302", 13042, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1303", 13043, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1304", 13044, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1305", 13045, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1306", 13046, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1307", 13047, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1308", 13048, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1309", 13049, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1310", 13050, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1311", 13051, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1312", 13052, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1313", 13053, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1314", 13054, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1315", 13055, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1316", 13056, "CT14nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1317", 13065, "CT14nnlo_as_0116", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1318", 13069, "CT14nnlo_as_0120", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1319", 13100, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1320", 13101, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1321", 13102, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1322", 13103, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1323", 13104, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1324", 13105, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1325", 13106, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1326", 13107, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1327", 13108, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1328", 13109, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1329", 13110, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1330", 13111, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1331", 13112, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1332", 13113, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1333", 13114, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1334", 13115, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1335", 13116, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1336", 13117, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1337", 13118, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1338", 13119, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1339", 13120, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1340", 13121, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1341", 13122, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1342", 13123, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1343", 13124, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1344", 13125, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1345", 13126, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1346", 13127, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1347", 13128, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1348", 13129, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1349", 13130, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1350", 13131, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1351", 13132, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1352", 13133, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1353", 13134, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1354", 13135, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1355", 13136, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1356", 13137, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1357", 13138, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1358", 13139, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1359", 13140, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1360", 13141, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1361", 13142, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1362", 13143, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1363", 13144, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1364", 13145, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1365", 13146, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1366", 13147, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1367", 13148, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1368", 13149, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1369", 13150, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1370", 13151, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1371", 13152, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1372", 13153, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1373", 13154, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1374", 13155, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1375", 13156, "CT14nlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1376", 13163, "CT14nlo_as_0116", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1377", 13167, "CT14nlo_as_0120", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1378", 13200, "CT14lo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1379", 25200, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1380", 25201, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1381", 25202, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1382", 25203, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1383", 25204, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1384", 25205, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1385", 25206, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1386", 25207, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1387", 25208, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1388", 25209, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1389", 25210, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1390", 25211, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1391", 25212, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1392", 25213, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1393", 25214, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1394", 25215, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1395", 25216, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1396", 25217, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1397", 25218, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1398", 25219, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1399", 25220, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1400", 25221, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1401", 25222, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1402", 25223, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1403", 25224, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1404", 25225, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1405", 25226, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1406", 25227, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1407", 25228, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1408", 25229, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1409", 25230, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1410", 25231, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1411", 25232, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1412", 25233, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1413", 25234, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1414", 25235, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1415", 25236, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1416", 25237, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1417", 25238, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1418", 25239, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1419", 25240, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1420", 25241, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1421", 25242, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1422", 25243, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1423", 25244, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1424", 25245, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1425", 25246, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1426", 25247, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1427", 25248, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1428", 25249, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1429", 25250, "MMHT2014nlo68clas118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1430", 25300, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1431", 25301, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1432", 25302, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1433", 25303, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1434", 25304, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1435", 25305, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1436", 25306, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1437", 25307, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1438", 25308, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1439", 25309, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1440", 25310, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1441", 25311, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1442", 25312, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1443", 25313, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1444", 25314, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1445", 25315, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1446", 25316, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1447", 25317, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1448", 25318, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1449", 25319, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1450", 25320, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1451", 25321, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1452", 25322, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1453", 25323, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1454", 25324, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1455", 25325, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1456", 25326, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1457", 25327, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1458", 25328, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1459", 25329, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1460", 25330, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1461", 25331, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1462", 25332, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1463", 25333, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1464", 25334, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1465", 25335, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1466", 25336, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1467", 25337, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1468", 25338, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1469", 25339, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1470", 25340, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1471", 25341, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1472", 25342, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1473", 25343, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1474", 25344, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1475", 25345, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1476", 25346, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1477", 25347, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1478", 25348, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1479", 25349, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1480", 25350, "MMHT2014nnlo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1481", 25000, "MMHT2014lo68cl", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1482", 42780, "ABMP16als118_5_nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1483", 42781, "ABMP16als118_5_nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1484", 42782, "ABMP16als118_5_nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1485", 42783, "ABMP16als118_5_nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1486", 42784, "ABMP16als118_5_nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1487", 42785, "ABMP16als118_5_nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1488", 42786, "ABMP16als118_5_nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1489", 42787, "ABMP16als118_5_nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1490", 42788, "ABMP16als118_5_nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1491", 42789, "ABMP16als118_5_nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1492", 42790, "ABMP16als118_5_nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1493", 42791, "ABMP16als118_5_nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1494", 42792, "ABMP16als118_5_nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1495", 42793, "ABMP16als118_5_nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1496", 42794, "ABMP16als118_5_nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1497", 42795, "ABMP16als118_5_nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1498", 42796, "ABMP16als118_5_nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1499", 42797, "ABMP16als118_5_nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1500", 42798, "ABMP16als118_5_nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1501", 42799, "ABMP16als118_5_nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1502", 42800, "ABMP16als118_5_nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1503", 42801, "ABMP16als118_5_nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1504", 42802, "ABMP16als118_5_nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1505", 42803, "ABMP16als118_5_nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1506", 42804, "ABMP16als118_5_nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1507", 42805, "ABMP16als118_5_nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1508", 42806, "ABMP16als118_5_nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1509", 42807, "ABMP16als118_5_nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1510", 42808, "ABMP16als118_5_nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1511", 42809, "ABMP16als118_5_nnlo", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1512", 90200, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1513", 90201, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1514", 90202, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1515", 90203, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1516", 90204, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1517", 90205, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1518", 90206, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1519", 90207, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1520", 90208, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1521", 90209, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1522", 90210, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1523", 90211, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1524", 90212, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1525", 90213, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1526", 90214, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1527", 90215, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1528", 90216, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1529", 90217, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1530", 90218, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1531", 90219, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1532", 90220, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1533", 90221, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1534", 90222, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1535", 90223, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1536", 90224, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1537", 90225, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1538", 90226, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1539", 90227, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1540", 90228, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1541", 90229, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1542", 90230, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1543", 90231, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1544", 90232, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1545", 90233, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1546", 90234, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1547", 90235, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1548", 90236, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1549", 90237, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1550", 90238, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1551", 90239, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1552", 90240, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1553", 90241, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1554", 90242, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1555", 90243, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1556", 90244, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1557", 90245, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1558", 90246, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1559", 90247, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1560", 90248, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1561", 90249, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1562", 90250, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1563", 90251, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1564", 90252, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1565", 90253, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1566", 90254, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1567", 90255, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1568", 90256, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1569", 90257, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1570", 90258, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1571", 90259, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1572", 90260, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1573", 90261, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1574", 90262, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1575", 90263, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1576", 90264, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1577", 90265, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1578", 90266, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1579", 90267, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1580", 90268, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1581", 90269, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1582", 90270, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1583", 90271, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1584", 90272, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1585", 90273, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1586", 90274, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1587", 90275, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1588", 90276, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1589", 90277, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1590", 90278, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1591", 90279, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1592", 90280, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1593", 90281, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1594", 90282, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1595", 90283, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1596", 90284, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1597", 90285, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1598", 90286, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1599", 90287, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1600", 90288, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1601", 90289, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1602", 90290, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1603", 90291, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1604", 90292, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1605", 90293, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1606", 90294, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1607", 90295, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1608", 90296, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1609", 90297, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1610", 90298, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1611", 90299, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1612", 90300, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1613", 90301, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1614", 90302, "PDF4LHC15_nlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1615", 91200, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1616", 91201, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1617", 91202, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1618", 91203, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1619", 91204, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1620", 91205, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1621", 91206, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1622", 91207, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1623", 91208, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1624", 91209, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1625", 91210, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1626", 91211, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1627", 91212, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1628", 91213, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1629", 91214, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1630", 91215, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1631", 91216, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1632", 91217, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1633", 91218, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1634", 91219, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1635", 91220, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1636", 91221, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1637", 91222, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1638", 91223, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1639", 91224, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1640", 91225, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1641", 91226, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1642", 91227, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1643", 91228, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1644", 91229, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1645", 91230, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1646", 91231, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1647", 91232, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1648", 91233, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1649", 91234, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1650", 91235, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1651", 91236, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1652", 91237, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1653", 91238, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1654", 91239, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1655", 91240, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1656", 91241, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1657", 91242, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1658", 91243, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1659", 91244, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1660", 91245, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1661", 91246, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1662", 91247, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1663", 91248, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1664", 91249, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1665", 91250, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1666", 91251, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1667", 91252, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1668", 91253, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1669", 91254, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1670", 91255, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1671", 91256, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1672", 91257, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1673", 91258, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1674", 91259, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1675", 91260, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1676", 91261, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1677", 91262, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1678", 91263, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1679", 91264, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1680", 91265, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1681", 91266, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1682", 91267, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1683", 91268, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1684", 91269, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1685", 91270, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1686", 91271, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1687", 91272, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1688", 91273, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1689", 91274, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1690", 91275, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1691", 91276, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1692", 91277, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1693", 91278, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1694", 91279, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1695", 91280, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1696", 91281, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1697", 91282, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1698", 91283, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1699", 91284, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1700", 91285, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1701", 91286, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1702", 91287, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1703", 91288, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1704", 91289, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1705", 91290, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1706", 91291, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1707", 91292, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1708", 91293, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1709", 91294, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1710", 91295, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1711", 91296, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1712", 91297, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1713", 91298, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1714", 91299, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1715", 91300, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1716", 91301, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1717", 91302, "PDF4LHC15_nnlo_100_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1718", 90400, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1719", 90401, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1720", 90402, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1721", 90403, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1722", 90404, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1723", 90405, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1724", 90406, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1725", 90407, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1726", 90408, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1727", 90409, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1728", 90410, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1729", 90411, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1730", 90412, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1731", 90413, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1732", 90414, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1733", 90415, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1734", 90416, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1735", 90417, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1736", 90418, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1737", 90419, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1738", 90420, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1739", 90421, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1740", 90422, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1741", 90423, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1742", 90424, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1743", 90425, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1744", 90426, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1745", 90427, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1746", 90428, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1747", 90429, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1748", 90430, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1749", 90431, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1750", 90432, "PDF4LHC15_nlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1751", 91400, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1752", 91401, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1753", 91402, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1754", 91403, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1755", 91404, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1756", 91405, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1757", 91406, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1758", 91407, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1759", 91408, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1760", 91409, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1761", 91410, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1762", 91411, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1763", 91412, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1764", 91413, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1765", 91414, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1766", 91415, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1767", 91416, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1768", 91417, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1769", 91418, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1770", 91419, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1771", 91420, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1772", 91421, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1773", 91422, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1774", 91423, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1775", 91424, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1776", 91425, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1777", 91426, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1778", 91427, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1779", 91428, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1780", 91429, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1781", 91430, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1782", 91431, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1783", 91432, "PDF4LHC15_nnlo_30_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1784", 61100, "HERAPDF20_NLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1785", 61101, "HERAPDF20_NLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1786", 61102, "HERAPDF20_NLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1787", 61103, "HERAPDF20_NLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1788", 61104, "HERAPDF20_NLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1789", 61105, "HERAPDF20_NLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1790", 61106, "HERAPDF20_NLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1791", 61107, "HERAPDF20_NLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1792", 61108, "HERAPDF20_NLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1793", 61109, "HERAPDF20_NLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1794", 61110, "HERAPDF20_NLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1795", 61111, "HERAPDF20_NLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1796", 61112, "HERAPDF20_NLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1797", 61113, "HERAPDF20_NLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1798", 61114, "HERAPDF20_NLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1799", 61115, "HERAPDF20_NLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1800", 61116, "HERAPDF20_NLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1801", 61117, "HERAPDF20_NLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1802", 61118, "HERAPDF20_NLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1803", 61119, "HERAPDF20_NLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1804", 61120, "HERAPDF20_NLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1805", 61121, "HERAPDF20_NLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1806", 61122, "HERAPDF20_NLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1807", 61123, "HERAPDF20_NLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1808", 61124, "HERAPDF20_NLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1809", 61125, "HERAPDF20_NLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1810", 61126, "HERAPDF20_NLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1811", 61127, "HERAPDF20_NLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1812", 61128, "HERAPDF20_NLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1813", 61130, "HERAPDF20_NLO_VAR", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1814", 61131, "HERAPDF20_NLO_VAR", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1815", 61132, "HERAPDF20_NLO_VAR", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1816", 61133, "HERAPDF20_NLO_VAR", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1817", 61134, "HERAPDF20_NLO_VAR", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1818", 61135, "HERAPDF20_NLO_VAR", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1819", 61136, "HERAPDF20_NLO_VAR", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1820", 61137, "HERAPDF20_NLO_VAR", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1821", 61138, "HERAPDF20_NLO_VAR", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1822", 61139, "HERAPDF20_NLO_VAR", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1823", 61140, "HERAPDF20_NLO_VAR", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1824", 61141, "HERAPDF20_NLO_VAR", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1825", 61142, "HERAPDF20_NLO_VAR", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1826", 61143, "HERAPDF20_NLO_VAR", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1827", 61200, "HERAPDF20_NNLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1828", 61201, "HERAPDF20_NNLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1829", 61202, "HERAPDF20_NNLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1830", 61203, "HERAPDF20_NNLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1831", 61204, "HERAPDF20_NNLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1832", 61205, "HERAPDF20_NNLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1833", 61206, "HERAPDF20_NNLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1834", 61207, "HERAPDF20_NNLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1835", 61208, "HERAPDF20_NNLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1836", 61209, "HERAPDF20_NNLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1837", 61210, "HERAPDF20_NNLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1838", 61211, "HERAPDF20_NNLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1839", 61212, "HERAPDF20_NNLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1840", 61213, "HERAPDF20_NNLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1841", 61214, "HERAPDF20_NNLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1842", 61215, "HERAPDF20_NNLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1843", 61216, "HERAPDF20_NNLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1844", 61217, "HERAPDF20_NNLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1845", 61218, "HERAPDF20_NNLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1846", 61219, "HERAPDF20_NNLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1847", 61220, "HERAPDF20_NNLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1848", 61221, "HERAPDF20_NNLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1849", 61222, "HERAPDF20_NNLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1850", 61223, "HERAPDF20_NNLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1851", 61224, "HERAPDF20_NNLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1852", 61225, "HERAPDF20_NNLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1853", 61226, "HERAPDF20_NNLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1854", 61227, "HERAPDF20_NNLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1855", 61228, "HERAPDF20_NNLO_EIG", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1856", 61230, "HERAPDF20_NNLO_VAR", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1857", 61231, "HERAPDF20_NNLO_VAR", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1858", 61232, "HERAPDF20_NNLO_VAR", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1859", 61233, "HERAPDF20_NNLO_VAR", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1860", 61234, "HERAPDF20_NNLO_VAR", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1861", 61235, "HERAPDF20_NNLO_VAR", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1862", 61236, "HERAPDF20_NNLO_VAR", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1863", 61237, "HERAPDF20_NNLO_VAR", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1864", 61238, "HERAPDF20_NNLO_VAR", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1865", 61239, "HERAPDF20_NNLO_VAR", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1866", 61240, "HERAPDF20_NNLO_VAR", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1867", 61241, "HERAPDF20_NNLO_VAR", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1868", 61242, "HERAPDF20_NNLO_VAR", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1869", 61243, "HERAPDF20_NNLO_VAR", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1870", 13400, "CT14qed_inc_proton", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1871", 13401, "CT14qed_inc_proton", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1872", 13402, "CT14qed_inc_proton", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1873", 13403, "CT14qed_inc_proton", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1874", 13404, "CT14qed_inc_proton", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1875", 13405, "CT14qed_inc_proton", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1876", 13406, "CT14qed_inc_proton", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1877", 13407, "CT14qed_inc_proton", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1878", 13408, "CT14qed_inc_proton", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1879", 13409, "CT14qed_inc_proton", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1880", 13410, "CT14qed_inc_proton", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1881", 13411, "CT14qed_inc_proton", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1882", 13412, "CT14qed_inc_proton", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1883", 13413, "CT14qed_inc_proton", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1884", 13414, "CT14qed_inc_proton", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1885", 13415, "CT14qed_inc_proton", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1886", 13416, "CT14qed_inc_proton", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1887", 13417, "CT14qed_inc_proton", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1888", 13418, "CT14qed_inc_proton", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1889", 13419, "CT14qed_inc_proton", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1890", 13420, "CT14qed_inc_proton", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1891", 13421, "CT14qed_inc_proton", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1892", 13422, "CT14qed_inc_proton", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1893", 13423, "CT14qed_inc_proton", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1894", 13424, "CT14qed_inc_proton", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1895", 13425, "CT14qed_inc_proton", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1896", 13426, "CT14qed_inc_proton", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1897", 13427, "CT14qed_inc_proton", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1898", 13428, "CT14qed_inc_proton", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1899", 13429, "CT14qed_inc_proton", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1900", 13430, "CT14qed_inc_proton", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1901", 82200, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1902", 82201, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1903", 82202, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1904", 82203, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1905", 82204, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1906", 82205, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1907", 82206, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1908", 82207, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1909", 82208, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1910", 82209, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1911", 82210, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1912", 82211, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1913", 82212, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1914", 82213, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1915", 82214, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1916", 82215, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1917", 82216, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1918", 82217, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1919", 82218, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1920", 82219, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1921", 82220, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1922", 82221, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1923", 82222, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1924", 82223, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1925", 82224, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1926", 82225, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1927", 82226, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1928", 82227, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1929", 82228, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1930", 82229, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1931", 82230, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1932", 82231, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1933", 82232, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1934", 82233, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1935", 82234, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1936", 82235, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1937", 82236, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1938", 82237, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1939", 82238, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1940", 82239, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1941", 82240, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1942", 82241, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1943", 82242, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1944", 82243, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1945", 82244, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1946", 82245, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1947", 82246, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1948", 82247, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1949", 82248, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1950", 82249, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1951", 82250, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1952", 82251, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1953", 82252, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1954", 82253, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1955", 82254, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1956", 82255, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1957", 82256, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1958", 82257, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1959", 82258, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1960", 82259, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1961", 82260, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1962", 82261, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1963", 82262, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1964", 82263, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1965", 82264, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1966", 82265, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1967", 82266, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1968", 82267, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1969", 82268, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1970", 82269, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1971", 82270, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1972", 82271, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1973", 82272, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1974", 82273, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1975", 82274, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1976", 82275, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1977", 82276, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1978", 82277, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1979", 82278, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1980", 82279, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1981", 82280, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1982", 82281, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1983", 82282, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1984", 82283, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1985", 82284, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1986", 82285, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1987", 82286, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1988", 82287, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1989", 82288, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1990", 82289, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1991", 82290, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1992", 82291, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1993", 82292, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1994", 82293, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1995", 82294, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1996", 82295, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1997", 82296, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1998", 82297, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("1999", 82298, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2000", 82299, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2001", 82300, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2002", 82301, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2003", 82302, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2004", 82303, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2005", 82304, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2006", 82305, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2007", 82306, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2008", 82307, "LUXqed17_plus_PDF4LHC15_nnlo_100", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2009", 292200, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2010", 292201, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2011", 292202, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2012", 292203, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2013", 292204, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2014", 292205, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2015", 292206, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2016", 292207, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2017", 292208, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2018", 292209, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2019", 292210, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2020", 292211, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2021", 292212, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2022", 292213, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2023", 292214, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2024", 292215, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2025", 292216, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2026", 292217, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2027", 292218, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2028", 292219, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2029", 292220, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2030", 292221, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2031", 292222, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2032", 292223, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2033", 292224, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2034", 292225, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2035", 292226, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2036", 292227, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2037", 292228, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2038", 292229, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2039", 292230, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2040", 292231, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2041", 292232, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2042", 292233, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2043", 292234, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2044", 292235, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2045", 292236, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2046", 292237, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2047", 292238, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2048", 292239, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2049", 292240, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2050", 292241, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2051", 292242, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2052", 292243, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2053", 292244, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2054", 292245, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2055", 292246, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2056", 292247, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2057", 292248, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2058", 292249, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2059", 292250, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2060", 292251, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2061", 292252, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2062", 292253, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2063", 292254, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2064", 292255, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2065", 292256, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2066", 292257, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2067", 292258, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2068", 292259, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2069", 292260, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2070", 292261, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2071", 292262, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2072", 292263, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2073", 292264, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2074", 292265, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2075", 292266, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2076", 292267, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2077", 292268, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2078", 292269, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2079", 292270, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2080", 292271, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2081", 292272, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2082", 292273, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2083", 292274, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2084", 292275, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2085", 292276, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2086", 292277, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2087", 292278, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2088", 292279, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2089", 292280, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2090", 292281, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2091", 292282, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2092", 292283, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2093", 292284, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2094", 292285, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2095", 292286, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2096", 292287, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2097", 292288, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2098", 292289, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2099", 292290, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2100", 292291, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2101", 292292, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2102", 292293, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2103", 292294, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2104", 292295, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2105", 292296, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2106", 292297, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2107", 292298, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2108", 292299, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2109", 292300, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2110", 292301, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2111", 292302, "NNPDF30_nlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2112", 292600, "NNPDF30_nnlo_nf_5_pdfas", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2113", 315000, "NNPDF31_lo_as_0118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2114", 315200, "NNPDF31_lo_as_0130", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2115", 262000, "NNPDF30_lo_as_0118", 1.0, 1.0));
LO26xInfo.push_back(weightinfo("2116", 263000, "NNPDF30_lo_as_0130", 1.0, 1.0));
///End of v26x
}


//define this as a plug-in
DEFINE_FWK_MODULE(weight_assign);
