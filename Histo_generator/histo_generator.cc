
#include <iostream>



#include <TObjArray.h>
#include <TH1D.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <vector>
#include <TBranch.h>
#include <TSystem.h>
#include <TDirectory.h>

using namespace std;



class histo_generator{
private:
  TFile *f;

public :
  ////inner class/////
  class Name_and_X{
  public:
    TString name;
    double X;
    Name_and_X();
    Name_and_X(TString name,double X){
      this->name=name;
      this->X=X;
    }

  };
  
  ///////////////////////
  

  vector<Name_and_X> vX;
  vector<Name_and_X> vW;
  vector<vector<TH1D*>> vH;
  void run();
  void set_vX();
  void set_vW();
  void set_file_and_tree(TString filename,TString xbranchname, TString wbranchname);
  void set_histo();
  void save_result(TString tag);
  void clear_job();
  TTree *Xtree;
  TTree *Wtree;
  histo_generator(){
    vX.clear();
    vW.clear();
    vH.clear();
  }
};

void histo_generator::set_histo(){
  for(unsigned int ix=0; ix < vX.size(); ix++){
    vector<TH1D*> vHi;
    for(unsigned int iw=0; iw < vW.size(); iw++){
      TString xname=vX[ix].name;
      TString wname=vW[iw].name;
      TString histoname=xname+"_"+wname;
      int nbin=100; double xmin=0, xmax=100;
      if(xname.Contains("pt")){
	nbin=100; xmin=0; xmax=200;
      }
      else if(xname.Contains("mass")){
	nbin=60; xmin=60; xmax=120;

      }
      else if(xname.Contains("eta")){
	nbin=200; xmin=-10; xmax=10;
      }
      else if(xname.Contains("phi")){
	nbin=100; xmin=-3.5; xmax=3.5;
      }
      TH1D *h=new TH1D(histoname,histoname,nbin,xmin,xmax);
      vHi.push_back(h);
    }
    vH.push_back(vHi);//vH[ix][iw]
    vHi.clear();

  }

}

void histo_generator::set_file_and_tree(TString filename,TString xbranchname, TString wbranchname){
  f=TFile::Open(filename);
  Xtree=(TTree*)f->Get(xbranchname);
  Wtree=(TTree*)f->Get(wbranchname);
}

void histo_generator::set_vX(){
  TObjArray *ListBranch=  Xtree->TTree::GetListOfBranches();
  for(int i = 0; i<ListBranch->GetEntries(); i++){
    TBranch *br=(TBranch*)ListBranch->At(i);
    this->vX.push_back( Name_and_X(br->TBranch::GetName(),-1)  ) ;
  }
  for(unsigned int i = 0; i<vX.size(); i++){
    Xtree->SetBranchAddress(vX[i].name, &vX[i].X);
  }

}
void histo_generator::set_vW(){
  TObjArray *ListBranch=  Wtree->TTree::GetListOfBranches();
  for(int i = 0; i<ListBranch->GetEntries(); i++){
    TBranch *br=(TBranch*)ListBranch->At(i);
    this->vW.push_back( Name_and_X(br->TBranch::GetName(),-1) );
  }
  for(unsigned int i = 0; i<vW.size(); i++){
    Wtree->SetBranchAddress(vW[i].name, &vW[i].X);
  }

}


void histo_generator::run(){
  //  vector<TString,double> vtest;
  //f=TFile::Open("../OUTPUT_0.root");
  //  TTree *tree=(TTree*)f->Get("Analyzer/DYkinematics");
  //TTree *weights=(TTree*)f->Get("Analyzer/DYweights");
  
  //  double Z_pt=0;
  //double weight=0;

  //  Xtree->SetBranchAddress("Z_pt",&Z_pt);
  //Wtree->SetBranchAddress("muR_1_muF_1_247000",&weight);

  int N=Xtree->GetEntries();
  int NX=vX.size();
  int NW=vW.size();
  //  N=10;
  //muR_1_muF_1_247000
  //TH1D *h= new TH1D("Z_pt","Z_pt",200,0,200);
  for(int i = 0 ; i < N; i++){
    Xtree->GetEntry(i);
    Wtree->GetEntry(i);

    //vH[ix][iw]
    for(int ix=0; ix<NX; ix++){
      for(int iw=0;iw<NW;iw++){
	
	vH[ix][iw]->Fill(vX[ix].X, vW[iw].X);

      }
    }
  }
  

  


}
//void histo_generator::save_result(TString savedir=gDirectory->GetPath(),TString tag=""){
void histo_generator::save_result(TString tag){

  //  gSystem->Exec("mkdir -p "+savedir);

  TFile ft("OUTPUT"+tag+".root","recreate");
  



  TList *hList=new TList();
  int NX=vX.size();
  int NW=vW.size();
  
  for(int ix=0; ix<NX; ix++){
    for(int iw=0;iw<NW;iw++){
      
      //   vH[ix][iw]->Fill(vX[ix].X, vW[iw].X);
      hList->Add(vH[ix][iw]);    
    } 
  }
  hList->Write();
  ft.Write();
  ft.Close();
  delete hList;
}
void histo_generator::clear_job(){
  for(unsigned int ix=0; ix < vX.size(); ix++){
    for(unsigned int iw=0; iw < vW.size(); iw++){
      delete vH[ix][iw];
    }
  }
  vX.clear();
  vW.clear();
  vH.clear();
 
}

void run_job_i(TString tag, int i){
  //  TString dir="/cms/ldap_home/jhchoi/generator_group/reweight_issue/GEN_Analyzer/JOB_242LO_false_pdfwgt_38/OUTPUT_0.root";
  //TString tag="242LO_false_pdfwgt";
  
  cout<<tag<<"_"<<i<<endl;
  histo_generator ana;
  ana.set_file_and_tree("/cms/ldap_home/jhchoi/generator_group/reweight_issue/GEN_Analyzer/JOB_"+tag+"_"+std::to_string(i)+"/OUTPUT_0.root","Analyzer/DYkinematics","Analyzer/DYweights");
  ana.set_vX();
  ana.set_vW();
  ana.set_histo();
  ana.run();
  ana.save_result(tag+"_"+std::to_string(i));
  ana.clear_job();
  
  cout<< "end of "<<tag<<endl;
 
}
void run_job_i_batch(TString tag, int i){
  cout<<tag<<"_"<<i<<endl;
  histo_generator ana;
  ana.set_file_and_tree("/cms/ldap_home/jhchoi/generator_group/reweight_issue/GEN_Analyzer/JOB_"+tag+"_"+std::to_string(i)+"/OUTPUT_0.root","Analyzer/DYkinematics","Analyzer/DYweights");
  ana.set_vX();
  ana.set_vW();
  ana.set_histo();
  ana.run();
  
  ana.save_result(tag+"_"+std::to_string(i));
  ana.clear_job();

  cout<< "end of "<<tag<<endl;


}
void run_jobs(TString tag, int Njobs){
  for(int i =0 ; i < Njobs; i++){
    run_job_i(tag,i);
  }
}
void test(){
  TString tag = "242LO_false_pdfwgt";int i =0;
  histo_generator ana;
  ana.set_file_and_tree("/cms/ldap_home/jhchoi/generator_group/reweight_issue/GEN_Analyzer/JOB_"+tag+"_"+std::to_string(i)+"/OUTPUT_0.root","Analyzer/DYkinematics","Analyzer/DYweights");
  ana.set_vX();
  ana.set_vW();
  ana.set_histo();
  ana.run();

  ana.save_result(tag);
  ana.clear_job();

}
