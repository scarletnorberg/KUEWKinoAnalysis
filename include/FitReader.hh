#ifndef FitReader_HH
#define FitReader_HH

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <iostream>
#include <vector>
#include <TGraphErrors.h>

#include "Category.hh"
#include "Process.hh"

using std::map;
using std::string;
using std::pair;

///////////////////////////////////////////
////////// FitReader class
///////////////////////////////////////////

class FitReader {
public:
  FitReader(const string& inputfile,
	    const string& otherfile = "",
	    const string& otherfold = "");

  virtual ~FitReader();

  void PrintCategories(bool verbose = false);
  void PrintProcesses(bool verbose = false);

  bool IsFilled(const Category&   cat,
		const Process&    proc,
		const Systematic& sys = Systematic::Default()) const;

  const TH1D* GetHistogram(const Category&   cat,
			   const Process&    proc,
			   const Systematic& sys = Systematic::Default()) const;

  // returns a 1-bin histogram added over all RISR/Mperp bins, all cats, all procs
  TH1D* GetIntegralHist(const string&       name,
			const CategoryList& cats,
			const ProcessList&  procs,
			const Systematic&   sys = Systematic::Default()) const;

  // returns a RISR/Mperp histogram added over all cats, all procs (need to have same binning or errors)
  TH1D* GetAddedHist(const string&       name,
		     const CategoryList& cats,
		     const ProcessList&  procs,
		     const Systematic&   sys = Systematic::Default()) const;

  bool IsFilled2D(const Category&   cat,
		  const Process&    proc,
		  const Systematic& sys = Systematic::Default()) const;

  const TH2D* GetHistogram2D(const Category&   cat,
			     const Process&    proc,
			     const Systematic& sys = Systematic::Default()) const;

  double Integral(const Category&   cat,
		  const Process&    proc,
		  const Systematic& sys = Systematic::Default()) const;

  bool HasSystematic(const Process& proc, const Systematic& sys) const;
  
  TCanvas* Plot1Dstack(const VS& proc,
		       const VS& lep_cat,
		       const VS& hadS_cat,
		       const VS& hadI_cat,
		       const string& canvas);

  TCanvas* Plot2D(const VS& proc,
		  const VS& lep_cat,
		  const VS& hadS_cat,
		  const VS& hadI_cat,
		  const string& canvas,
		  const string& extra = "");

  TCanvas* PlotYields(const string& can_name,
		      const VS& proc,
		      const CategoryTree& CT);

  TCanvas* Plot1Dstack(const string& can_name,
		       const VS& proc,
		       const CategoryTree& CT);

  TCanvas* Plot2D(const string& can_name,
		  const VS& proc,
		  const CategoryTree& CT);

  VS GetChannels() const;  
  const ProcessList&  GetProcesses() const;
  const CategoryList& GetCategories(const string& channel = "") const;         
  const Systematics&  GetSystematics() const;
  
private:
  mutable TFile  m_File;
  
  mutable TFile* m_FilePtr;
  string         m_FileFold;

  string m_CMSLabel;

  mutable map<Process,Systematics> m_ProcSys;
  mutable map<Process,map<Category,TH1D*> > m_ProcHist;
  mutable map<Process,map<Systematic,map<Category,pair<TH1D*,TH1D*> > > >m_ProcHistSys;
  mutable map<Process,map<Category,TH2D*> > m_ProcHist_2D;
  mutable map<Process,map<Systematic,map<Category,pair<TH2D*,TH2D*> > > >m_ProcHistSys_2D;

  ProcessList                      m_Proc;
  mutable map<string,CategoryList> m_Chan;
  mutable map<string,bool>         m_CatLabel;
  CategoryList                     m_Cat;
  Systematics                      m_Sys;
  
  ProcessBranch m_ProcBranch;
  void ReadProcesses();
  
  CategoryBranch m_CatBranch;
  void ReadCategories();

  map<string,string>          m_Title;
  map<string,int>             m_Color;
  vector<int>                 m_ColorDefault;
  vector<int>                 m_SignalColor;
  map<string,VS> m_Strings;
  void InitializeRecipes();
  string GetSignalTitle(const string& label);

  TGraphErrors* GetTotalBackground(const CategoryList& cat);

  void DrawCatTree(const CategoryTree& CT, TCanvas* can);
  void DrawMR(const FitBin& fitbin, TCanvas* can);
  
};

#endif
