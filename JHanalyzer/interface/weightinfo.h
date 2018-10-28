#include <TString.h>
using namespace std;
class weightinfo {

 private:
  string _pdf;
  string _id;
  TString _name;
  double _muR,_muF;//,_weight;
 public:
  weightinfo();
  weightinfo(string id, string pdf, TString name, double muR, double muF);
  void  set_id(string id);
  void set_pdf(string pdf);
  void set_name(TString name);
  void set_muR(double muR);
  void set_muF(double muF);
  void set_weight(double weight);
  string id();
  string pdf();
  TString name();
  double muR();
  double muF();
  double weight();
  double _weight;
};
