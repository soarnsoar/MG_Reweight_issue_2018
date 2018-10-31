/*
scale_variation
NNPDF31_nnlo_hessian_pdfas
NNPDF31_nnlo_as_0108
NNPDF31_nnlo_as_0110
NNPDF31_nnlo_as_0112
NNPDF31_nnlo_as_0114
NNPDF31_nnlo_as_0117
NNPDF31_nnlo_as_0119
NNPDF31_nnlo_as_0122
NNPDF31_nnlo_as_0124
NNPDF31_nlo_hessian_pdfas
CT14nnlo
CT14nnlo_as_0116
CT14nnlo_as_0120
CT14nlo
CT14nlo_as_0116
CT14nlo_as_0120
CT14lo
MMHT2014nlo68clas118
MMHT2014nnlo68cl
MMHT2014lo68cl
ABMP16als118_5_nnlo
PDF4LHC15_nlo_100_pdfas
PDF4LHC15_nnlo_100_pdfas
PDF4LHC15_nlo_30_pdfas
PDF4LHC15_nnlo_30_pdfas
HERAPDF20_NLO_EIG
HERAPDF20_NLO_VAR
HERAPDF20_NNLO_EIG
HERAPDF20_NNLO_VAR
CT14qed_inc_proton
LUXqed17_plus_PDF4LHC15_nnlo_100
NNPDF30_nlo_nf_5_pdfas
NNPDF30_nnlo_nf_5_pdfas
NNPDF31_lo_as_0118
NNPDF31_lo_as_0130
NNPDF30_lo_as_0118
NNPDF30_lo_as_0130


 */
//  TH1::AddDirectory(kFALSE) ;
#include <TFile.h>
#include <TList.h>
#include <TKey.h>
#include <TDirectoryFile.h>
#include <TH1D.h>
#include <vector>
#include <TGraphAsymmErrors.h>
#include <TString.h>
#include <iostream>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TSystem.h>

using namespace std;

class plotter{

private:
  TFile *file1;
  TFile *file2;
  TList *list1;
  TList *list2;
  TKey *key;
  TKey *hkey;
  TDirectoryFile *dir;
  TList *hlist;
  TH1D* h;


  /*  TString filename1;
  TString filename2;
  TString title1;
  TString title2;
  */
public:

  //  void run();
  //vector<TString> vX;
  //vector<TString> vW;
  vector<TH1D*> vH1;
  vector<TH1D*> vH2;
  //  TH1D *envH1_high;
  // TH1D *envH2_high;//envelope histo
  // TH1D *envH1_low;
  // TH1D *envH2_low;//envelope histo                                                                                                                        
  //  TH1D  *envH1;
  // TH1D *envH2;

  TGraphAsymmErrors *envH1;
  TGraphAsymmErrors *envH2;
  TGraphAsymmErrors *envH1R;
  TGraphAsymmErrors *envH2R;
  
  TString Xname;
  TString Wname;
  TString filename1;
  TString filename2;
  TString title1;
  TString title2;
  TString alias1;
  TString alias2;

  Color_t envcolor1,color1;
  Color_t envcolor2,color2;
  
  double Gymax,Gymin;
  double legx1,legy1,legx2,legy2;
  //  void run_scale();
  plotter();
  void set_histo_vector();
  void set_envelope_histo();
  void clear();
  void set_file1(TString filename1);
  void set_file2(TString filename2);
  void set_title1(TString title1);
  void set_title2(TString title2);
  void set_alias1(TString alias1);
  void set_alias2(TString alias2);

  void set_legend();
  void run();

  TString Xtitle();
};
plotter::plotter(){
  vH1.clear();
  vH2.clear();
  envcolor1=kRed;
  color1=kOrange;
  envcolor2=kBlue;
  color2=kGreen;
  legx1=0.9;//0.9,0.85,0.45,0.7
  legy1=0.85;
  legx2=0.45;
  legy2=0.7;

}
void plotter::set_legend(){
  if(Xname.Contains("phi")){
      legx1=0.7;
      legy1=0.25;
      legx2=0.25;
      legy2=0.1;
    }
  if(Xname.Contains("eta")){
    legx1=0.9;
    legy1=0.95;
    legx2=0.45;
    legy2=0.8;

  }

  if(Xname.Contains("mass")){
    legx1=0.5;
    legy1=0.3;
    legx2=0.05;
    legy2=0.15;

  }

}
TString plotter::Xtitle(){
  TString particle="";
  TString variable="";
  TString unit="";
  if(Xname.Contains("Z")) particle="Z";
  if(Xname.Contains("dilep")) particle="ll";
  if(Xname.Contains("lep1")) particle="l^{+}";
  if(Xname.Contains("lep2")) particle="l^{-}";

  if(Xname.Contains("mass")){ variable="M"; unit="[GeV]";}
  if(Xname.Contains("pt")){ variable="P_{T}"; unit="[GeV]";}
  if(Xname.Contains("eta")) variable="#eta";
  if(Xname.Contains("phi")) variable="#phi";
  TString title=variable+"("+particle+") "+unit;
  return title;
  }
void plotter::set_file1(TString filename1){
  this->filename1=filename1;
}
void plotter::set_file2(TString filename2){
  this->filename2=filename2;
}
void plotter::set_title1(TString title1){
  this->title1=title1;
}
void plotter::set_title2(TString title2){
  this->title2=title2;
}
void plotter::set_alias1(TString alias1){
  this->alias1=alias1;
}
void plotter::set_alias2(TString alias2){
  this->alias2=alias2;
}

void plotter::clear(){

  for(size_t i = 0; i < vH1.size(); i++){
    delete vH1[i];
    vH1[i] = NULL;
  }
  for(size_t i = 0; i <vH2.size(); i++){
    delete vH2[i];
    vH2[i] = NULL;
  }
  vH1.clear();
  vH2.clear();

  /*
  delete envH1;
  envH1=NULL;
  
  delete envH2;
  envH2=NULL;
  delete envH1R;
  envH1R=NULL;

  delete envH2R;
  envH2R=NULL;
  */


  //  delete h;
  // h=NULL;
  
  //delete hkey;
  //hkey=NULL;
  //delete hlist;
  //hlist=NULL;
  //  delete dir;
  //dir=NULL;

  //  delete h;
  //  h=NULL;

  
  //  delete key;
  //key=NULL;
  //  delete file1;
  //file1=NULL;
  // delete file2;
  //file2=NULL;
  
}


void plotter::set_histo_vector(){

  //  vX.push_back("Z_pt");
  //vW.push_back("24700");



  
  //  file1=TFile::Open("OUTPUT242LO_false_pdfwgt_0.root");
  // file2=TFile::Open("OUTPUT26xLO_false_pdfwgt_0.root");
  file1=TFile::Open(filename1);
  file2=TFile::Open(filename2);


  list1= (TList*)file1->GetListOfKeys();
  list2= (TList*)file2->GetListOfKeys();


  TString object_name="";

  //TIter next(list);//
  //  //  TObject *key =;
  key=(TKey*)list1->First();
  //while( key!=list1->Last() ){
  while( 1 ){
    //    TKey *hkey=(TKey*)key->GetListOfKeys();
    TString dirname=key->GetName();
    //    cout<<"dirname="<<dirname<<endl;
    //if(hname.Contains(Xname)  &&   hname.Contains(Wname)   ){      
    if( dirname.Contains(Wname)   ){
      
      dir=(TDirectoryFile*)file1->Get(dirname);
      
      hlist=(TList*)dir->TDirectoryFile::GetListOfKeys();
      hkey=(TKey*)hlist->First();
      
      //      while( hkey!=hlist->Last() ){
      while(1 ){
	
	TString hname=hkey->GetName();
	
	if( hname.Contains(Xname)   ){
	  object_name=dirname+"/"+hname;

	  /*
	  h=(TH1D*)file1->Get(dirname+"/"+hname);
	  h->SetDirectory(0); 
	  file1->Close();
	  //	  cout<<"close file"<<endl;
	  vH1.push_back(h);
	  //cout<<"vH1.push_back(h)"<<endl;
	  */	 
	  break;
	}
	if(hkey==hlist->Last()) break;
	hkey=(TKey*)hlist->After(hkey);
      }//while hkey
      //cout<<"1"<<endl;
      //cout<<"2"<<endl;
      //delete dir;
      //dir=NULL;


      // cout<<"3"<<endl;

    }//if W match
    if(key==list1->Last()) break;
    key=(TKey*)list1->After(key);



  }//while list1
  //  h=(TH1D*)file1->Get(dirname+"/"+hname);                                                                                              
  h=(TH1D*)file1->Get(object_name);                                                                                              
  h->SetDirectory(0);                                                                                                                  
  file1->Close();                                                                                                                      
  //      cout<<"close file"<<endl;                                                                                                    
  vH1.push_back(h);       



  key=(TKey*)list2->First();
  //  while( key!=list2->Last() ){
  while( 1 ){
    //    TKey *hkey=(TKey*)key->GetListOfKeys();                                                                                              
    TString dirname=key->GetName();
    //if(hname.Contains(Xname)  &&   hname.Contains(Wname)   ){                                                                                
    if( dirname.Contains(Wname)   ){
      dir=(TDirectoryFile*)file2->Get(dirname);
      hlist=(TList*)dir->TDirectoryFile::GetListOfKeys();
      //      TH1D* h=(TH1D*)file1->Get(dirname);   
      
      hkey=(TKey*)hlist->First();
      //while( hkey!=hlist->Last() ){
      while( 1 ){
        TString hname=hkey->GetName();
        if( hname.Contains(Xname)   ){
	  /*
	  h=(TH1D*)file2->Get(dirname+"/"+hname);
          h->SetDirectory(0);
          file2->Close();
	  vH2.push_back(h);
	  cout<<"add h to f2"<<endl;
	  */
	  object_name=dirname+"/"+hname;
	  break;
        }
	if(hkey==hlist->Last()) break;
        hkey=(TKey*)hlist->After(hkey);
      }//while hkey 
      

      //      delete hkey;
      // hkey=NULL;
      //   delete hlist;
      // delete dir;
    }
    if(key==list2->Last()) break;
    key=(TKey*)list2->After(key);
  }
  //  delete key;

  //  h=(TH1D*)file2->Get(dirname+"/"+hname);                                                                                              
  h=(TH1D*)file2->Get(object_name);
  h->SetDirectory(0);                                                                                                                  
  file2->Close();                                                                                                                      
  //      cout<<"close file"<<endl;                                                                                                    
  vH2.push_back(h);       

}

  	
  ////////////////  
  /*TCanvas *c = new TCanvas();
  for(size_t i = 0; i < vH1.size(); i++){
    vH1.at(i)->Draw("E4 sames");
    vH2.at(i)->Draw("E4 sames");
  }
  c->SaveAs("temp.pdf");
  */
  

  //Now we have the set of related histograms for 1/2.
  //make envelope histo!
  //1st histo//v242
  //cout<<"clone"<<endl;
  //envH1_high=(TH1D*)vH1[0]->TH1D::Clone();
  // envH1_low=(TH1D*)vH1[0]->TH1D::Clone();
  //  envH1=(TH1D*)vH1[0]->TH1D::Clone();     
void plotter::set_envelope_histo(){ 
  double x1[vH1[0]->GetNbinsX()];
  double y1[vH1[0]->GetNbinsX()];

  double exh1[vH1[0]->GetNbinsX()];
  double exl1[vH1[0]->GetNbinsX()];

  double eyh1[vH1[0]->GetNbinsX()];
  double eyl1[vH1[0]->GetNbinsX()];


  double r1[vH1[0]->GetNbinsX()];


  double erh1[vH1[0]->GetNbinsX()];
  double erl1[vH1[0]->GetNbinsX()];

  
  Gymin=vH1[0]->GetBinContent(1)-vH1[0]->GetBinError(1);
  Gymax=vH1[0]->GetBinContent(1)+vH1[0]->GetBinError(1);
  for(size_t i = 1; i <= vH1[0]->GetNbinsX(); i++){
		     
    double ymax= vH1[0]->GetBinContent(i)+vH1[0]->GetBinError(i);
    double ymin= vH1[0]->GetBinContent(i)-vH1[0]->GetBinError(i);
    for(size_t ih=0; ih<vH1.size();ih++){
      double ymax_i= vH1[ih]->GetBinContent(i)+vH1[ih]->GetBinError(i);
      double ymin_i= vH1[ih]->GetBinContent(i)-vH1[ih]->GetBinError(i);
      if(ymax_i>ymax) ymax=ymax_i;
      if(ymin_i<ymin) ymin=ymin_i;


    }
    // envH1_high->SetBinContent(i,ymax);
    //  envH1_low->SetBinContent(i,ymin);
    //   envH1->SetBinErrorUp( ymax- vH1[0]->GetBinContent(i)  );
    //envH1->SetBinErrorDown( -ymin + vH1[0]->GetBinContent(i)  );
    x1[i-1]=vH1[0]->GetBinCenter(i);
    exh1[i-1]=0; exl1[i]=0;

    y1[i-1]=vH1[0]->GetBinContent(i);
    eyh1[i-1]=ymax- vH1[0]->GetBinContent(i);
    eyl1[i-1]=-ymin +vH1[0]->GetBinContent(i);
    //if (vH1.size()==1)cout<<eyh1[i-1]<<endl;

    r1[i-1]=1;
    if(y1[i-1]==0.){
      erh1[i-1]=0.;
      erl1[i-1]=0.;

    }
    else{
    erh1[i-1]=eyh1[i-1]/y1[i-1];
    erl1[i-1]=eyl1[i-1]/y1[i-1];
    }


    if(Gymin>ymin) Gymin=ymin;
    if(Gymax<ymax) Gymax=ymax;
  }
  envH1=new TGraphAsymmErrors(vH1[0]->GetNbinsX(),x1,y1,exl1,exh1,eyl1,eyh1);
  envH1R=new TGraphAsymmErrors(vH1[0]->GetNbinsX(),x1,r1,exl1,exh1,erl1,erh1);

  envH1->SetTitle(Xname+"_"+Wname);
  if(Wname=="24700") envH1->SetTitle(Xname+"_Scale_Variation");
  //2nd histo
  //  envH2_high=(TH1D*)vH2[0]->TH1D::Clone();
  // envH2_low=(TH1D*)vH2[0]->TH1D::Clone();
  // cout<<vH2[0]->GetNbinsX()<<endl;
  double x2[vH2[0]->GetNbinsX()];
  double y2[vH2[0]->GetNbinsX()];

  double exh2[vH2[0]->GetNbinsX()];
  double exl2[vH2[0]->GetNbinsX()];

  double eyh2[vH2[0]->GetNbinsX()];
  double eyl2[vH2[0]->GetNbinsX()];

  double r2[vH2[0]->GetNbinsX()];
  double erh2[vH2[0]->GetNbinsX()];
  double erl2[vH2[0]->GetNbinsX()];


  for(size_t i = 1; i <= vH2[0]->GetNbinsX(); i++){

    double ymax= vH2[0]->GetBinContent(i)+vH2[0]->GetBinError(i);
    double ymin= vH2[0]->GetBinContent(i)-vH2[0]->GetBinError(i);
    for(size_t ih=0; ih<vH2.size();ih++){
      double ymax_i= vH2[ih]->GetBinContent(i)+vH2[ih]->GetBinError(i);
      double ymin_i= vH2[ih]->GetBinContent(i)-vH2[ih]->GetBinError(i);
      if(ymax_i>ymax) ymax=ymax_i;
      if(ymin_i<ymin) ymin=ymin_i;
      
      
    }

    x2[i-1]=vH2[0]->GetBinCenter(i);

    exh2[i-1]=0; exl2[i]=0;

    y2[i-1]=vH2[0]->GetBinContent(i);
    eyh2[i-1]=ymax- vH2[0]->GetBinContent(i);
    eyl2[i-1]=-ymin +vH2[0]->GetBinContent(i);

    //       cout<<"x2="<<x2[i-1]<<endl;
    //cout<<"y2="<<y2[i-1]<<endl;


    if(y1[i-1]==0.){
      r2[i-1]=1;
      erh2[i-1]=0;
      erl2[i-1]=0;
 
    }
    else{
      r2[i-1]=y2[i-1]/y1[i-1];
      erh2[i-1]=(y2[i-1]+eyh2[i-1])/y1[i-1]-r2[i-1];
      erl2[i-1]=-(y2[i-1]-eyl2[i-1])/y1[i-1]+r2[i-1];
    }
    if(Gymin>ymin) Gymin=ymin;
    if(Gymax<ymax) Gymax=ymax;

    
  }
  envH2=new TGraphAsymmErrors(vH2[0]->GetNbinsX(),x2,y2,exl2,exh2,eyl2,eyh2);
  envH2R=new TGraphAsymmErrors(vH2[0]->GetNbinsX(),x2,r2,exl2,exh2,erl2,erh2);

  envH2->SetTitle(Xname+"_"+Wname);
  if(Wname=="24700") envH2->SetTitle(Xname+"_Scale_Variation");


  for(size_t i = 0; i < vH1.size(); i++){
    vH1.at(i)->SetLineColorAlpha(color1,0.35);
    vH2.at(i)->SetLineColorAlpha(color2,0.35);

    vH1.at(i)->SetMarkerStyle(color1);
    vH2.at(i)->SetMarkerStyle(color2);

    vH1.at(i)->SetMarkerColor(color1);
    vH2.at(i)->SetMarkerColor(color2);

    vH1.at(i)->SetFillColor(color1);
    vH2.at(i)->SetFillColor(color2);

  }
  
  envH1->SetLineColor(envcolor1);
  envH2->SetLineColor(envcolor2);
  
  envH1R->SetLineColor(envcolor1);
  envH2R->SetLineColor(envcolor2);


  envH1->SetMarkerColor(envcolor1);
  envH2->SetMarkerColor(envcolor2);

  envH1R->SetMarkerColor(envcolor1);
  envH2R->SetMarkerColor(envcolor2);

  envH1->SetLineWidth(2);
  envH2->SetLineWidth(2);

  envH1->SetFillColorAlpha(envcolor1,0.5);
  envH2->SetFillColorAlpha(envcolor2,0.5);

  envH1R->SetFillColorAlpha(envcolor1,0.5);
  //envH1R->SetFillColor(envcolor1);
  envH2R->SetFillColorAlpha(envcolor2,0.5);
  //envH2R->SetFillColor(envcolor2);

		   
  //  envH1->GetYaxis()->SetRangeUser(TMath::Max(Gymin,1.),Gymax*1.1);
  // envH2->GetYaxis()->SetRangeUser(TMath::Max(Gymin,1.),Gymax*1.1);
  //cout<<TMath::Max(Gymin,1.)<<endl;
  //  envH1_high->SetLineColor(envcolor1);
  // envH2_high->SetLineColor(envcolor2);

  //  envH1_low->SetLineColor(envcolor1);
  // envH2_low->SetLineColor(envcolor2);
  
  //  envH1_high->SetFillColor(0);                                                                                                                      
  //envH2_high->SetFillColor(0);
  //envH1_low->SetFillColor(0);
  //envH2_low->SetFillColor(0);


}


void plotter::run(){

  vector<TString> vX;
  vX.push_back("Z_pt");  vX.push_back("Z_mass");  vX.push_back("Z_eta");  vX.push_back("Z_phi");
  vX.push_back("dilep_pt");  vX.push_back("dilep_mass");  vX.push_back("dilep_eta");  vX.push_back("dilep_phi");
  vX.push_back("lep1_pt"); vX.push_back("lep1_eta");  vX.push_back("lep1_phi");
  vX.push_back("lep2_pt"); vX.push_back("lep2_eta");  vX.push_back("lep2_phi");
   
  vector<TString> vW;
  
  vW.push_back("24700");
  
  vW.push_back("NNPDF31_nnlo_hessian_pdfas");  
  vW.push_back("NNPDF31_nnlo_as_0108");
  vW.push_back("NNPDF31_nnlo_as_0110");
  vW.push_back("NNPDF31_nnlo_as_0112");
  vW.push_back("NNPDF31_nnlo_as_0114");
  vW.push_back("NNPDF31_nnlo_as_0117");
  vW.push_back("NNPDF31_nnlo_as_0119");
  vW.push_back("NNPDF31_nnlo_as_0122");
  vW.push_back("NNPDF31_nnlo_as_0124");
  vW.push_back("NNPDF31_nlo_hessian_pdfas");
  vW.push_back("CT14nnlo");
  vW.push_back("CT14nnlo_as_0116");
  vW.push_back("CT14nnlo_as_0120");
  vW.push_back("CT14nlo");
  vW.push_back("CT14nlo_as_0116");
  vW.push_back("CT14nlo_as_0120");
  vW.push_back("CT14lo");
  vW.push_back("MMHT2014nlo68clas118");
  vW.push_back("MMHT2014nnlo68cl");
  vW.push_back("MMHT2014lo68cl");
  vW.push_back("ABMP16als118_5_nnlo");
  vW.push_back("PDF4LHC15_nlo_100_pdfas");
  vW.push_back("PDF4LHC15_nnlo_100_pdfas");
  vW.push_back("PDF4LHC15_nlo_30_pdfas");
  vW.push_back("PDF4LHC15_nnlo_30_pdfas");
  vW.push_back("HERAPDF20_NLO_EIG");
  vW.push_back("HERAPDF20_NLO_VAR");
  vW.push_back("HERAPDF20_NNLO_EIG");
  vW.push_back("HERAPDF20_NNLO_VAR");
  vW.push_back("CT14qed_inc_proton");
  vW.push_back("LUXqed17_plus_PDF4LHC15_nnlo_100");
  vW.push_back("NNPDF30_nlo_nf_5_pdfas");
  vW.push_back("NNPDF30_nnlo_nf_5_pdfas");
  vW.push_back("NNPDF31_lo_as_0118");
  vW.push_back("NNPDF31_lo_as_0130");
  
  vW.push_back("NNPDF30_lo_as_0118");
  
  vW.push_back("NNPDF30_lo_as_0130");
  
  
  for(size_t iw=0;iw<vW.size();iw++){
    
    for(size_t ix=0; ix<vX.size();ix++){
   


      this->Xname=vX[ix];
      this->Wname=vW[iw];
      // this->Xname=Xname;
      //this->Wname=Wname;
      envcolor1=kRed;
      color1=kOrange;
      envcolor2=kBlue;
      color2=kGreen;
      //      set_title1("MGv242,LO, pdfwgt=F");
      //set_file1("OUTPUT242LO_false_pdfwgt_0.root");
      // set_title2("MGv26x,LO, pdfwgt=F");
      //set_file2("OUTPUT26xLO_false_pdfwgt_0.root");

      set_histo_vector();
      set_envelope_histo();

      if(ix==0) {
	int Nweight=vH1.size();
	std::cout<<"=====================" << Wname << ": number of weight=" << Nweight<<"=====================" << std::endl;
      }
      TCanvas *c=new TCanvas();
      c->Clear();
      c->Divide(1,2);
      c->cd(1);
      gPad->SetPad(0.05, 0.3, 0.95, 0.95);

      c->cd(2);
      gPad->SetPad(0.05,0.05,0.95,0.3);
      gPad->SetBottomMargin(0.3);
      

      c->cd(1);
      c->cd(1)->SetGrid();
      TMultiGraph *mg= new TMultiGraph();
      mg->SetTitle(Xname+"_"+Wname);
      if(Wname=="24700") mg->SetTitle(Xname+"_Scale_Variation");
      
      if(vH1.size()==1){
	mg->Add(envH1,"3L");
	mg->Add(envH2,"3L");
      }
      else{
	mg->Add(envH1,"3L");
	mg->Add(envH2,"3L");
      }
      mg->Draw("A");
      mg->GetYaxis()->SetTitle("#Events");
      mg->GetXaxis()->SetTitle("");
      
      //      envH1->Draw("SAMES E3LA");
      //      envH2->Draw("SAMES E3LA");
      //
      set_legend();
      TLegend *leg = new TLegend(legx1,legy1,legx2,legy2);
      leg->AddEntry(envH1,title1);
      leg->AddEntry(envH2,title2);
      leg->Draw();
      
      c->cd(2);
      c->cd(2)->SetGrid();

      TMultiGraph *mg2= new TMultiGraph();
      mg2->SetTitle("");
      if(vH1.size()==1){

	mg2->Add(envH1R,"3L");
	mg2->Add(envH2R,"3L");

      }
      else{
      mg2->Add(envH1R,"3L");
      mg2->Add(envH2R,"3L");
      }
      mg2->Draw("A");
      mg2->GetYaxis()->SetTitle("#frac{"+title2+"}{"+title1+"}");
      mg2->GetXaxis()->SetTitle(Xtitle());
      mg2->SetMinimum(0.9);
      mg2->SetMaximum(1.1);
      mg2 -> GetXaxis() -> SetLabelSize(0.1);
      mg2 -> GetYaxis() -> SetLabelSize(0.1);
      mg2->GetYaxis() -> SetTitleSize(0.09);
      mg2 -> GetXaxis() -> SetTitleSize(0.15);
      mg2 -> GetYaxis() -> SetTitleOffset(0.4);
      mg2 -> GetXaxis() -> SetTitleOffset(0.8);
      mg2-> GetYaxis()-> SetNdivisions(6);


      for(size_t i = 0 ; i < vH1.size(); i++){
	
	//vH1[i]->Draw("l sames");
	//vH2[i]->Draw("l sames");
      }
      gSystem->Exec("mkdir -p Envelope_plots_"+alias1+"_"+alias2+"/");
      c->SaveAs("Envelope_plots_"+alias1+"_"+alias2+"/"+Xname+"_"+Wname+".pdf");
      //      mg->GetXaxis()->SetRangeUser(TMath::Max(Gymin,1.),Gymax*1.1);
      mg->SetMinimum(TMath::Max(Gymin,1.));
      mg->SetMaximum(Gymax*1.1);
      c->cd(1)->SetLogy();
      c->SaveAs("Envelope_plots_"+alias1+"_"+alias2+"/"+Xname+"_"+Wname+"_log.pdf");
      //      c->SaveAs("Envelope_plots/"+Xname+"_"+Wname+"_log.pdf");
      clear();
      //      delete c;      
      //delete mg;
      
      //      c = NULL;
      //mg = NULL;

    }
  }
      //  TTree *tree=(TTree*)file->Get();
  //    file->ls();
  // cout<< file->GetStreamerInfoList()<<endl;
  //  #TList *list = file->GetStreamerInfoList();
  //  TList *list = file->GetList();
  //TH1D *h=(TH1D*)list->At(0);
  //  cout<<(h->GetName())<<endl;
  //  vX.clear();
  // vW.clear();
}


void Run(){
  //      set_title1("MGv242,LO, pdfwgt=F");                                                                                               
  //set_file1("OUTPUT242LO_false_pdfwgt_0.root");                                                                                          
  // set_title2("MGv26x,LO, pdfwgt=F");                                                                                                    
  //set_file2("OUTPUT26xLO_false_pdfwgt_0.root");                                                                                          

  plotter ana;
  ana.set_title1("MGv242,LO, pdfwgt=F");
  ana.set_title2("MGv26x,LO, pdfwgt=F");
 
  ana.set_file1("histosets/242LO_false_pdfwgt.root"); 
  ana.set_file2("histosets/26xLO_false_pdfwgt.root"); 
 
  ana.set_alias1("MG242LO_false_pdfwgt");
  ana.set_alias2("MG26xLO_false_pdfwgt");
  ana.run();

  ana.set_title1("MGv242,LO, pdfwgt=T");
  ana.set_title2("MGv26x,LO, pdfwgt=T");

  ana.set_file1("histosets/242LO_true_pdfwgt.root");
  ana.set_file2("histosets/26xLO_true_pdfwgt.root");

  ana.set_alias1("MG242LO_true_pdfwgt");
  ana.set_alias2("MG26xLO_true_pdfwgt");
  ana.run();
  
  ana.set_title1("MGv242,LO, pdfwgt=F");
  ana.set_title2("MGv242,LO, pdfwgt=T");

  ana.set_file1("histosets/242LO_false_pdfwgt.root");
  ana.set_file2("histosets/242LO_true_pdfwgt.root");

  ana.set_alias1("MG242LO_false_pdfwgt");
  ana.set_alias2("MG242LO_true_pdfwgt");
  ana.run();

  ana.set_title1("MGv26x,LO, pdfwgt=F");
  ana.set_title2("MGv26x,LO, pdfwgt=T");

  ana.set_file1("histosets/26xLO_false_pdfwgt.root");
  ana.set_file2("histosets/26xLO_true_pdfwgt.root");

  ana.set_alias1("MG26xLO_false_pdfwgt");
  ana.set_alias2("MG26xLO_true_pdfwgt");
  ana.run();


}
