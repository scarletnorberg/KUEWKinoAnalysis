#include "../include/JMETool.hh"

#include "include_jme/FactorizedJetCorrector.h"
#include "include_jme/JetCorrectorParameters.h"

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <string>

using std::string;
using std::cout;
using std::endl;

JMETool::~JMETool(){
  for(int i = 0; i < 3; i++){
    std::map<string, FormulaBinsBins*>::iterator it = m_Factors[i].begin();
    while(it != m_Factors[i].end()){
      delete it->second;
      it++;
    }
  }
}

double JMETool::GetFactor(int year, const string& name, double pT, double Eta, double A, double rho) const {
  if(year < 2016 || year > 2018)
    return 0;
  
  if(m_Factors[year-2016].count(name) == 0)
    return 0;
  
  return m_Factors[year-2016][name]->SF(Eta, pT, A, rho);
}

void JMETool::BuildMap(const std::string& JMEfolder){
  string JME[3]; // [year]
  string JME_FS[3]; // [year];

  JME[0] = "Summer16_07Aug2017_V11_MC";
  JME[1] = "Fall17_17Nov2017_V32_MC";
  JME[2] = "Autumn18_V19_MC";
  
  JME_FS[0] = "Summer16_FastSimV1_MC";
  JME_FS[1] = "Fall17_FastSimV1_MC";
  JME_FS[2] = "Autumn18_FastSimV1_MC";

  string post = "AK4PFchs.txt";

   /*
   jmrValues = {
      '2016': [1.0, 1.2, 0.8],
      '2017': [1.09, 1.14, 1.04],
      # Use 2017 values for 2018 until 2018 are released
      '2018': [1.09, 1.14, 1.04],
      'UL2017': [1.00, 1.00, 1.00],  # placeholder
        }
   jmsValues = {
      '2016': [1.00, 0.9906, 1.0094],  # nominal, down, up
      '2017': [0.982, 0.978, 0.986],
      # Use 2017 values for 2018 until 2018 are released
      '2018': [0.982, 0.978, 0.986],
      'UL2017': [1.000, 1.000, 1.000],  # placeholder
                 }
     */

  cout << "Building JME maps for all three years" << endl;
  
  for(int y = 0; y < 3; y++){
    // ParseJESUncertainty(JMEfolder+"/"+JME[y]+"/"+JME[y]+"_Uncertainty_"+post, 2016+y);
    ParseJESUncertainty(JMEfolder+"/"+JME[y]+"/"+JME[y]+"_UncertaintySources_"+post, 2016+y);
    ParseJESUncertainty(JMEfolder+"/"+JME_FS[y]+"/"+JME_FS[y]+"_Uncertainty_"+post, 2016+y, "FS");

    ParseJEC(JMEfolder+"/"+JME[y]+"/"+JME[y]+"_L1FastJet_"+post, 2016+y);
    ParseJEC(JMEfolder+"/"+JME_FS[y]+"/"+JME_FS[y]+"_L1FastJet_"+post, 2016+y, "FS");
  
    ParseJEC(JMEfolder+"/"+JME[y]+"/"+JME[y]+"_L2Relative_"+post, 2016+y);
    ParseJEC(JMEfolder+"/"+JME_FS[y]+"/"+JME_FS[y]+"_L2Relative_"+post, 2016+y, "FS");
    //cout << "...Done" << endl;
  }
}
//Reapply jec: https://github.com/cms-jet/JetMETAnalysis/blob/master/JetAnalyzers/bin/jet_apply_jec_x.cc
//https://cms-nanoaod-integration.web.cern.ch/integration/master-cmsswmaster/mc102X_doc.html varaibles
//https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/jme/jetmetUncertainties.py uncertainties 
//notneeded
//https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/jme/fatJetUncertainties.py
void JMETool::ParseJESUncertainty(const string& input, int year, const string& prefix){
  if(year < 2016 || year > 2018)
    return;
  //double jet_rawpt, jet_rawmass;
  //int Njet = jetsLVec->size();
  //for(jn>0, jn<=Njet, jn++){
  //jet_rawpt = jetsLVec.pt(jn) * (1 - jet.rawFactor);
  //jet_rawmass = jet.mass * (1 - jet.rawFactor);
  //}
  // corrector = new FactorizedJetCorrector(JetInfo::get_correction_levels(levels,L1FastJet),
  //                                        JetInfo::get_correction_tags(era,alg,levels,jecpath,L1FastJet));
  /*
  for(unsigned int ilevel=0; ilevel<levels.size(); ilevel++)
  {
  vPar.push_back(JetCorrectorParameters(string(jecpath + era + 
                                        JetInfo::get_level_tag(levels[ilevel],L1FastJet) + 
                                        jetInfo.getAlias() + 
                                        getPostfix(postfix,alg,levels[ilevel]) + ".txt")));
                                        }
                                        corrector = new FactorizedJetCorrector(vPar);
          */
  /*
  corrector->setJetPt(JRAEvt->jtpt->at(ijt));
  corrector->setJetEta(JRAEvt->jteta->at(ijt));
  corrector->setRho(JRAEvt->rho);
  JRAEvt->jtpt->at(ijt)*=jec;
  JRAEvt->jte->at(ijt) *=jec;
  */
/*
            # Get the JEC factors
            jec = jet_pt / jet_rawpt
            jecL1 = jet_pt_l1 / jet_rawpt
*/
 /*
        for jesUncertainty in self.jesUncertainties:
            jets_pt_jesUp[jesUncertainty] = []
            jets_pt_jesDown[jesUncertainty] = []
            jets_mass_jesUp[jesUncertainty] = []
            jets_mass_jesDown[jesUncertainty] = []

        if 'T1' in self.saveMETUncs:
            (met_T1_px_jerUp, met_T1_px_jerDown, met_T1_py_jerUp,
             met_T1_py_jerDown) = ({}, {}, {}, {})
        if 'T1Smear' in self.saveMETUncs:
            (met_T1Smear_px_jerUp, met_T1Smear_px_jerDown,
             met_T1Smear_py_jerUp, met_T1Smear_py_jerDown) = ({}, {}, {}, {})
        for jerID in self.splitJERIDs:
            if 'T1' in self.saveMETUncs:
                (met_T1_px_jerUp[jerID], met_T1_py_jerUp[jerID]) = (met_px,
                                                                    met_py)
                (met_T1_px_jerDown[jerID], met_T1_py_jerDown[jerID]) = (met_px,
                                                                        met_py)
            if 'T1Smear' in self.saveMETUncs:
                (met_T1Smear_px_jerUp[jerID],
                 met_T1Smear_py_jerUp[jerID]) = (met_px, met_py)
                (met_T1Smear_px_jerDown[jerID],
                 met_T1Smear_py_jerDown[jerID]) = (met_px, met_py)

        if 'T1' in self.saveMETUncs:
            (met_T1_px_jesUp, met_T1_py_jesUp) = ({}, {})
            (met_T1_px_jesDown, met_T1_py_jesDown) = ({}, {})
            for jesUncertainty in self.jesUncertainties:
                met_T1_px_jesUp[jesUncertainty] = met_px
                met_T1_py_jesUp[jesUncertainty] = met_py
                met_T1_px_jesDown[jesUncertainty] = met_px
                met_T1_py_jesDown[jesUncertainty] = met_py
        if 'T1Smear' in self.saveMETUncs:
            (met_T1Smear_px_jesUp, met_T1Smear_py_jesUp) = ({}, {})
            (met_T1Smear_px_jesDown, met_T1Smear_py_jesDown) = ({}, {})
            for jesUncertainty in self.jesUncertainties:
                met_T1Smear_px_jesUp[jesUncertainty] = met_px
                met_T1Smear_py_jesUp[jesUncertainty] = met_py
                met_T1Smear_px_jesDown[jesUncertainty] = met_px
                met_T1Smear_py_jesDown[jesUncertainty] = met_py

        # variables needed for re-applying JECs to 2017 v2 MET
        delta_x_T1Jet, delta_y_T1Jet = 0, 0
        delta_x_rawJet, delta_y_rawJet = 0, 0
for jesUncertainty in self.jesUncertainties:
                    # cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetCorUncertainties
                    # cf. https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/2000.html
                    if jesUncertainty == "HEMIssue":
                        delta = 1.
                        if iJet < nJet and jet_pt_nom > 15 and jet.jetId & 2 and jet.phi > -1.57 and jet.phi < -0.87:
                            if jet.eta > -2.5 and jet.eta < -1.3:
                                delta = 0.8
                            elif jet.eta <= -2.5 and jet.eta > -3:
                                delta = 0.65
                        jet_pt_jesUp[jesUncertainty] = jet_pt_nom
                        jet_pt_jesDown[jesUncertainty] = delta * jet_pt_nom
                        jet_mass_jesUp[jesUncertainty] = jet_mass_nom
                        jet_mass_jesDown[jesUncertainty] = delta * jet_mass_nom

                        jet_pt_jesUpT1[jesUncertainty] = jet_pt_L1L2L3
                        jet_pt_jesDownT1[jesUncertainty] = delta * \
                            jet_pt_L1L2L3

                    else:
                        self.jesUncertainty[jesUncertainty].setJetPt(
                            jet_pt_nom)
                        self.jesUncertainty[jesUncertainty].setJetEta(jet.eta)
                        delta = self.jesUncertainty[
                            jesUncertainty].getUncertainty(True)
                        jet_pt_jesUp[jesUncertainty] = jet_pt_nom * \
                            (1. + delta)
                        jet_pt_jesDown[jesUncertainty] = jet_pt_nom * \
                            (1. - delta)
                        jet_mass_jesUp[jesUncertainty] = jet_mass_nom * \
                            (1. + delta)
                        jet_mass_jesDown[jesUncertainty] = jet_mass_nom * \
                            (1. - delta)

                        # redo JES variations for T1 MET
                        self.jesUncertainty[jesUncertainty].setJetPt(
                            jet_pt_L1L2L3)
                        self.jesUncertainty[jesUncertainty].setJetEta(jet.eta)
                        delta = self.jesUncertainty[
                            jesUncertainty].getUncertainty(True)
                        jet_pt_jesUpT1[jesUncertainty] = jet_pt_L1L2L3 * \
                            (1. + delta)
                        jet_pt_jesDownT1[jesUncertainty] = jet_pt_L1L2L3 * \
                            (1. - delta)

*/
                 
   if(!ifile.is_open()){
    std::cout << "can't open JME file " << input << std::endl;
    return;
   }
   
   string line;
   size_t found;

   string label = "JESUncer_Total";
   if(prefix != "")
     label = prefix+"_"+label;
   bool start = false;
   FormulaBinsBins* Bins = nullptr;

   double ptmin, ptmax, factor;
   while(getline(ifile,line)){
     // discard comment line
     if(line.find("#") != string::npos)
       continue;

     // new name, new label
     if(line.find("[") != string::npos){
       found = line.find("[");
       line.erase(0,found+1);
       found = line.find("]");
       label = "JESUncer_"+line.substr(0,found);
       if(prefix != "")
	 label = prefix+"_"+label;
       
       continue;
     }
       
     // start new bin map
     if(line.find("{") != string::npos){
       Bins = new FormulaBinsBins();
       m_Factors[year-2016][label] = Bins;
       continue;
     }

     // new eta bin for existing bin map
     double etamin = popdouble(line);
     double etamax = popdouble(line);
     
     const FormulaBins& Bin = Bins->AddBin(etamin, etamax);

     int N = popdouble(line);

     ptmin  = popdouble(line);
     while(line.length() > 0){
       factor = popdouble(line);
       popdouble(line); // factor repeated for some reason
       if(line.length() > 0)
	 ptmax = popdouble(line);
       else
	 ptmax = ptmin + 1000000;

       Bin.AddBin(ptmin, ptmax, factor);
       
       ptmin = ptmax;
     }
   }

   ifile.close();
}

void JMETool::ParseJEC(const string& input, int year, const string& prefix){
  if(year < 2016 || year > 2018)
    return;
  
  std::ifstream ifile(input.c_str());
   if(!ifile.is_open()){
    std::cout << "can't open JME file " << input << std::endl;
    return;
   }
   
   string line;
   size_t found;

   // first line with variable and equation info
   getline(ifile,line);
   
   // remove bracket
   found = line.find("{");
   line.erase(0,found+1);

   // get how many bins specified per line
   int Nbinvar = popdouble(line);

   for(int i = 0; i < Nbinvar; i++)
     popstring(line);

   // get the number of variables in the formula
   int Nvar = popdouble(line);

   std::vector<std::string> svar;
   for(int i = 0; i < Nvar; i++)
     svar.push_back(popstring(line));
   
   string formula = popstring(line);
  
   if(Nvar == 3){
     std::vector<std::string> var;
     std::vector<std::string> sym;
     std::vector<std::string> rep;
     var.push_back("JetPt");
     var.push_back("JetA");
     var.push_back("Rho");
     sym.push_back("x");
     sym.push_back("y");
     sym.push_back("z");
     rep.push_back("q");
     rep.push_back("v");
     rep.push_back("f");

     if(formula.find("max") != string::npos)
       formula.replace(formula.find("max"), 3, "%%%");
     
     for(int i = 0; i < 3; i++){
       while(formula.find(sym[i]) != string::npos){
	 found = formula.find(sym[i]);
	 formula.replace(found, 1, rep[i]);
       }
     }

     for(int v = 0; v < 3; v++){
       for(int s = 0; s < 3; s++){
	 if(svar[s].find(var[v]) != string::npos){
	   while(formula.find(rep[s]) != string::npos){
	     found = formula.find(rep[s]);
	     formula.replace(found, 1, sym[v]);
	   }
	 }
       }
     }
   }



   if(formula.find("%%%") != string::npos)
       formula.replace(formula.find("%%%"), 3, "max");
   
   popstring(line);

   string label = popstring(line);
   found = label.find("}");
   label.erase(found,found+1);
   if(prefix != "")
     label = prefix+"_"+label;

   // get the first line to get parameter numbers
   double etamin, etamax;
   double emi, ema;
   double ptmin, ptmax;
   getline(ifile,line);
   etamin = popdouble(line);
   etamax = popdouble(line);
   if(Nbinvar == 2){
     ptmin = popdouble(line);
     ptmax = popdouble(line);
   } else {
     ptmin = 0.;
     ptmax = 1e10;
   }

   int Nparam = popdouble(line);
   
   for(int i = 0; i < Nvar*2; i++)
     popdouble(line);
   Nparam -= Nvar*2;

   std::vector<double> params;
   for(int i = 0; i < Nparam; i++)
     params.push_back(popdouble(line));
   
   FormulaBinsBins* Bins = new FormulaBinsBins(formula, Nvar, Nparam);
   m_Factors[year-2016][label] = Bins;

   const FormulaBins* Bin = &Bins->AddBin(etamin, etamax);
   Bin->AddBin(ptmin, ptmax, params);

   while(getline(ifile,line)){
     emi = popdouble(line);
     ema = popdouble(line);
     
     // new eta bin
     if(fabs(emi-etamin) + fabs(ema-etamax) > 1e-8){
       etamin = emi;
       etamax = ema;
       Bin = &Bins->AddBin(etamin, etamax);
     }
     
     if(Nbinvar == 2){
       ptmin = popdouble(line);
       ptmax = popdouble(line);
     } else {
       ptmin = 0.;
       ptmax = 1e10;
     }
     
     int N = popdouble(line);
     for(int i = 0; i < Nvar*2; i++)
       popdouble(line);
     N -= Nvar*2;

     params.clear();
     for(int i = 0; i < N; i++)
       params.push_back(popdouble(line));
     
     Bin->AddBin(ptmin, ptmax, params);
   }

   ifile.close();
}

double JMETool::popdouble(std::string& line){
  // remove leading whitespace
  while(line[0] == string(" ")[0])
    line.erase(0,1);
  
  string num;
  if(line.find(" ") != string::npos){
    size_t p = line.find(" ");
    num = line.substr(0,p);
    line.erase(0,p+1);
    while(line[0] == string(" ")[0])
      line.erase(0,1);
  } else
    num = line;

  return stod(num);
}

std::string JMETool::popstring(std::string& line){
  // remove leading whitespace
  while(line[0] == string(" ")[0])
    line.erase(0,1);
  
  string ret;
  if(line.find(" ") != string::npos){
    size_t p = line.find(" ");
    ret = line.substr(0,p);
    line.erase(0,p+1);
    while(line[0] == string(" ")[0])
      line.erase(0,1);
  } else
    ret = line;

  return ret;
}

