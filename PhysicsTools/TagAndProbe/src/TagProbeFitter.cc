#include "PhysicsTools/TagAndProbe/interface/TagProbeFitter.h"
#include <stdexcept>
#include <sys/stat.h>

#include "TROOT.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TTreeFormula.h"
#include "TH2F.h"
#include "TStyle.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooAddPdf.h"
#include "RooGlobalFunc.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooBinning.h"
#include "RooBinningCategory.h"
#include "RooMultiCategory.h"
#include "RooMappedCategory.h"
#include "RooThresholdCategory.h"
#include "Roo1DTable.h"
#include "RooMinuit.h"
#include "RooNLLVar.h"
#include "RooAbsDataStore.h"
#include "RooEfficiency.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooProdPdf.h"
#include "RooGenericPdf.h"
#include "RooExtendPdf.h"
#include "RooTrace.h"
#include "RooMsgService.h"

#include "Math/QuantFuncMathCore.h"

using namespace std;
using namespace RooFit;

TagProbeFitter::TagProbeFitter(const std::vector<std::string>& inputFileNames, string inputDirectoryName, string inputTreeName, string outputFileName, int numCPU_, bool saveWorkspace_, bool docutandcount_, bool floatShapeParameters_, const std::vector<std::string>& fixVars_):
  numCPU(numCPU_),
  saveWorkspace(saveWorkspace_),
  docutandcount(docutandcount_),
  binnedFit(false),
  massBins(0),
  weightVar(""),
  fixVars(fixVars_), 
  floatShapeParameters(floatShapeParameters_),
  quiet(false) {

  inputTree = new TChain((inputDirectoryName+"/"+inputTreeName).c_str());
  for(size_t i=0; i<inputFileNames.size(); i++){
    inputTree->Add(inputFileNames[i].c_str());
  }

  outputFile = new TFile(outputFileName.c_str(),"recreate");
  outputDirectory = outputFile->mkdir(inputDirectoryName.c_str());

  if(!floatShapeParameters && fixVars.empty()) std::cout << "TagProbeFitter: " << "You want to fix some variables but do not specify them!";
  
  gROOT->SetStyle("Plain");
  gStyle->SetTitleFillColor(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat(".2f");
}

TagProbeFitter::~TagProbeFitter(){
  //if (inputTree)
  //  delete inputTree;
  if (outputFile)
    outputFile->Close();
}

void TagProbeFitter::setQuiet(bool quiet_) { 
  quiet = quiet_; 
  if (quiet) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  } else {
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  }
}
bool TagProbeFitter::addVariable(string name, string title, double low, double hi, string units){
  variables.addClone(RooRealVar(name.c_str(), title.c_str(), low, hi, units.c_str()));
  return true;
}

bool TagProbeFitter::addCategory(string name, string title, string expression){
  RooCategory* c = (RooCategory*) parameterParser.factory(expression.c_str());
  if(!c)
    return false;

  //set the name of the category to the passed name instead of the one in the expression
  c->SetName(name.c_str());
  c->SetTitle(title.c_str());
  variables.addClone(*c);
  return true;
}

bool TagProbeFitter::addExpression(string expressionName, string title, string expression, const std::vector<string>& arguments) {
  expressionVars.push_back(make_pair(make_pair(expressionName,title), make_pair(expression,arguments)));
  return true;
}


bool TagProbeFitter::addThresholdCategory(string categoryName, string varName, double cutValue, string cutType){
  thresholdCategories.push_back(make_pair(make_pair(categoryName, cutType), make_pair(varName,cutValue)));
  return true;
}


void TagProbeFitter::addPdf(string name, vector<string>& pdfCommands){
  pdfs[name] = pdfCommands;
}

void TagProbeFitter::setBinsForMassPlots(int bins){
  massBins = bins;
}

void TagProbeFitter::setWeightVar(const std::string &var) {
  weightVar = var;
}

string TagProbeFitter::calculateEfficiency(string dirName,const std::vector<string>& effCats,const std::vector<string>& effStates, vector<string>& unbinnedVariables, map<string, vector<double> >& binnedReals, map<string, std::vector<string> >& binnedCats, vector<string>& binToPDFmap) {
  
  //go to home directory
  outputDirectory->cd();
  //make a directory corresponding to this efficiency binning
  gDirectory->mkdir(dirName.c_str())->cd();
  
  RooArgSet dataVars;
  std::vector<std::string> vars;
  
  //collect unbinned variables
  for(vector<string>::iterator v=unbinnedVariables.begin(); v!=unbinnedVariables.end(); v++){
    vars.push_back(v->c_str());
    dataVars.addClone(variables[v->c_str()], true);
    if (binnedFit && (v == unbinnedVariables.begin()))
      ((RooRealVar&)dataVars[v->c_str()]).setBins(massBins);	
  }
  
  //add weight
  if (!weightVar.empty())
    dataVars.addClone(variables[weightVar.c_str()], true);

  //add efficiency var
  for (vector<string>::const_iterator effCat = effCats.begin(); effCat != effCats.end(); ++effCat) {
    if (variables.find(effCat->c_str()) != 0)
      dataVars.addClone(variables[effCat->c_str()], true);
  }

  char aChar[100];
  //collect the binned variables and the corresponding bin categories
  RooArgSet binnedVariables;
  std::vector<std::string> binnedVariableNames;
  //RooArgSet binCategories;
  std::vector<std::vector<std::string> > binnedCategories;
  for(map<string, vector<double> >::iterator v=binnedReals.begin(); v!=binnedReals.end(); v++){
    std::vector<std::string> aString;
    //std::string name = v->first;
    binnedVariableNames.push_back(v->first);
    //if (dataVars.find(variables[v->first.c_str()]) == 0) { 
    //  cerr << "Binned variable '"<<v->first<<"' not found." << endl; 
    //  return "Error"; 
    //}
    for (unsigned int i=0; i<v->second.size()-1; i++) {
      float lowThr = v->second[i];
      float highThr = v->second[i+1];
      sprintf(aChar, "%s>%f&&%s<%f", v->first.c_str(), lowThr, v->first.c_str(), highThr);
      aString.push_back(aChar);
    }
    
    binnedCategories.push_back(aString);
    binnedVariables.addClone(variables[v->first.c_str()]);
    ((RooRealVar&)binnedVariables[v->first.c_str()]).setBinning( RooBinning(v->second.size()-1, &v->second[0]) );
    //binCategories.addClone( RooBinningCategory((v->first+"_bins").c_str(), (v->first+"_bins").c_str(), (RooRealVar&)binnedVariables[v->first.c_str()]) );
  }
  dataVars.addClone(binnedVariables, true);

  //collect the category variables and the corresponding mapped categories
  RooArgSet categories;
  RooArgSet mappedCategories;
  for(map<string, vector<string> >::iterator v=binnedCats.begin(); v!=binnedCats.end(); v++){
    TString name = v->first;
    if (variables.find(name) == 0) { cerr << "Binned category '"<<name<<"' not found." << endl; return "Error"; }
    categories.addClone(variables[name]);
    mappedCategories.addClone(RooMappedCategory(name+"_bins", name+"_bins", (RooCategory&)categories[name]));
    for(unsigned int i = 0; i<v->second.size(); i++){
      ((RooMappedCategory&)mappedCategories[name+"_bins"]).map(v->second[i].c_str(), name+"_"+TString(v->second[i].c_str()).ReplaceAll(",","_"));
    }
  }
  dataVars.addClone(categories, true);

  
  //  add all variables used in expressions
  for(vector<pair<pair<string,string>, pair<string, vector<string> > > >::const_iterator ev = expressionVars.begin(), eve = expressionVars.end(); ev != eve; ++ev){
    for (vector<string>::const_iterator it = ev->second.second.begin(), ed = ev->second.second.end(); it != ed; ++it) {
      // provided that they are real variables themselves
      if (variables.find(it->c_str())) dataVars.addClone(variables[it->c_str()], true);
    }
  }
  
  // add all real variables used in cuts
  for(vector<pair<pair<string,string>, pair<string, double> > >::const_iterator tc = thresholdCategories.begin(), tce = thresholdCategories.end(); tc != tce; ++tc){
    if (variables.find(tc->second.first.c_str())) dataVars.addClone(variables[tc->second.first.c_str()], true);
  }
   
  // Now add all expressions that are computed dynamically
  std::vector<std::string> dynamicCuts;
  std::map<std::string, std::string> dynamicExpressions;
  for(vector<pair<pair<string,string>, pair<string, vector<string> > > >::const_iterator ev = expressionVars.begin(), eve = expressionVars.end(); ev != eve; ++ev){
    for (vector<string>::const_iterator it = ev->second.second.begin(), ed = ev->second.second.end(); it != ed; ++it) {
      dynamicExpressions[ev->first.first] = ev->second.first;
    }
  }
          
  // And add all dynamic categories from thresholds 
  std::string dynamicCutString = "1";
  for(vector<pair<pair<string,string>, pair<string, double> > >::const_iterator tc = thresholdCategories.begin(), tce = thresholdCategories.end(); tc != tce; ++tc){
    std::string cutName = tc->second.first.c_str();
    if (dynamicExpressions.find(cutName) != dynamicExpressions.end())
      cutName = "("+dynamicExpressions[cutName]+")";

    if (tc->first.second == "above")
      sprintf(aChar, "%s>%f", cutName.c_str(), tc->second.second);
    else
      sprintf(aChar, "%s<%f", cutName.c_str(), tc->second.second);
    
    dynamicCuts.push_back(aChar);
    dynamicCutString += " && " + dynamicCuts.back();
  }

  // NEED TO HANDLE 1D case
  std::vector<std::string> pdfNames;
  std::vector<std::string> categorySelection;
  std::vector<std::string> catNames;

  int nCats = 0;
  for (unsigned int i=0; i<binnedCategories[0].size(); i++) {
    for (unsigned int j=0; j<binnedCategories[1].size(); j++) {

      unsigned int pdfIndex = (i*binnedCategories[1].size()+j+1)*2;

      if (pdfIndex > binToPDFmap.size())
	pdfNames.push_back(binToPDFmap[0]);
      else
	pdfNames.push_back(binToPDFmap[pdfIndex]);

      categorySelection.push_back(binnedCategories[0][i]+" && " + binnedCategories[1][j]+ " && " + dynamicCutString);
      
      sprintf(aChar, "%s_bin%d__%s_bin%d", binnedVariableNames[0].c_str(), i, binnedVariableNames[1].c_str(), j);
      catNames.push_back(aChar);
      //      std::cout << aChar << std::endl;
      nCats++;
    }
  }
  
  std::vector<TTree*> trees;
  std::vector<TTreeFormula*> formulas;

  TString rootInputFile(inputTree->GetFile()->GetName());
  TObjArray *tx = rootInputFile.Tokenize("/");
  TString basename = ((TObjString *)(tx->Last()))->String().ReplaceAll(".root", "");
  std::string treeDirectory = std::string(basename.Data())+"_"+effCats[0]+"_"+catNames.back();
  //  struct stat sb;
    //  if (stat(treeDirectory.c_str(), &sb) != 0) { // && S_ISDIR(sb.st_mode)) {   
  {
    std::map<std::string, float> treeVarsF;
    std::map<std::string, double> treeVarsD;
    std::map<std::string, int> treeVarsI;
    
    inputTree->SetBranchStatus("*", 0);
    TIterator* vit = dataVars.createIterator();
    for(RooRealVar* v = (RooRealVar*)vit->Next(); v!=0; v = (RooRealVar*)vit->Next() ){
      inputTree->SetBranchStatus(v->GetName(), 1);
      
      TLeaf* myLeaf = (TLeaf*)inputTree->GetBranch(v->GetName())->GetListOfLeaves()->At(0);
      if (strcmp(myLeaf->GetTypeName(), "Float_t") == 0)
	inputTree->SetBranchAddress(v->GetName(), &(treeVarsF[v->GetName()]));
      else if (strcmp(myLeaf->GetTypeName(), "Double_t") == 0)
	inputTree->SetBranchAddress(v->GetName(), &(treeVarsD[v->GetName()]));
      else if (strcmp(myLeaf->GetTypeName(), "Int_t") == 0)
	inputTree->SetBranchAddress(v->GetName(), &(treeVarsI[v->GetName()]));
    }
    
    mkdir(treeDirectory.c_str(), S_IRUSR | S_IWUSR | S_IXUSR);
    outputTemp = new TFile((treeDirectory+std::string("/temp.root")).c_str(), "recreate");
    for (int category=0; category<nCats; category++) {
      trees.push_back(inputTree->CloneTree(0));
      sprintf(aChar, "formula%d", category);
      formulas.push_back(new TTreeFormula(aChar, categorySelection[category].c_str(), inputTree));
    }
    
    std::cout << "Splitting original tree according to chosen categorization..." << std::endl;
    for (int z=0; z<inputTree->GetEntries(); z++) {
      if (z % 100000 == 0)
	std::cout << z << std::endl;
      inputTree->GetEntry(z);
      int chosenCat = -1;
      for (int i=0; i<nCats; i++) {
	if (formulas[i]->EvalInstance() != 0 ) {
	  chosenCat = i;
	  break;
	}
      }
      
      if (chosenCat != -1) {
	trees[chosenCat]->Fill();
      }
    }
    
    outputTemp->cd();
    for (unsigned int i=0; i<trees.size(); i++) {
      sprintf(aChar, "tree_cat%d", i);
      trees[i]->SetName(aChar);
      trees[i]->Write();
    }  
    outputTemp->Close();
  }
  outputTemp = TFile::Open((treeDirectory+std::string("/temp.root")).c_str());
  outputTemp->cd();
  trees.resize(0);
  for (int i=0; i<nCats; i++) {
    sprintf(aChar, "tree_cat%d", i);
    trees.push_back((TTree*)outputTemp->Get(aChar));
  }  

  delete inputTree;
  outputDirectory->cd();
  gDirectory->cd(dirName.c_str());

  //create the empty efficiency datasets from the binned variables
  RooRealVar efficiency("efficiency", "Efficiency", 0, 1); 
  RooDataSet fitEfficiency("fit_eff", "Efficiency from unbinned ML fit", RooArgSet(RooArgSet(binnedVariables, categories), efficiency), StoreAsymError(RooArgSet(binnedVariables, efficiency)));
  RooDataSet cntEfficiency("cnt_eff", "Efficiency from counting", RooArgSet(RooArgSet(binnedVariables, categories), efficiency), StoreAsymError(RooArgSet(binnedVariables, efficiency)));
  
  string effName;   
  for (int category=0; category<nCats; category++) {
    
    //create the dataset
    RooDataSet* data = new RooDataSet("data", "data", trees[category], dataVars, categorySelection[category].c_str(), (weightVar.empty() ? 0 : weightVar.c_str()));
    if(!floatShapeParameters) {
      //fitting whole dataset to get initial values for some parameters
      RooWorkspace* w = new RooWorkspace();
      w->import(*(data->get(0)));
      efficiency.setVal(0);//reset
      efficiency.setAsymError(0,0);
      doFitEfficiency(w, data, pdfNames[category], efficiency);
      delete w;
    }

    //get PDF name and make directory name
    TString dirName = catNames[category];
    dirName.ReplaceAll("{","").ReplaceAll("}","").ReplaceAll(";","__");
    if(pdfNames[category].size() > 0){
      dirName.Append("__").Append(pdfNames[category]);
    }
    
    cout<<"Fitting bin:  "<<dirName<<endl;
    //make a directory for each bin
    gDirectory->mkdir(dirName)->cd();
    
    //setup the efficiency category
    if (effCats.size() == 1) {
      effName = effCats.front() + "::" + effStates.front();
      RooMappedCategory efficiencyCategory("_efficiencyCategory_", "_efficiencyCategory_", (RooCategory&)dataVars[effCats.front().c_str()], "Failed");
      efficiencyCategory.map(effStates.front().c_str(), "Passed");
      data->addColumn( efficiencyCategory );
      //w->import(efficiencyCategory);
    } else {
      RooArgSet rooEffCats; 
      string multiState = "{";
      for (size_t i = 0; i < effCats.size(); ++i) {
	if (i) { multiState += ";"; effName += " && "; }
	rooEffCats.add((RooCategory &) dataVars[effCats[i].c_str()]);
	multiState += effStates[i];
	effName = effCats[i] + "::" + effStates[i];
      }
      multiState += "}";
      RooMultiCategory efficiencyMultiCategory("_efficiencyMultiCategory_", "_efficiencyMultiCategory_", rooEffCats);
      RooMappedCategory efficiencyCategory("_efficiencyCategory_", "_efficiencyCategory_", efficiencyMultiCategory, "Failed");
      efficiencyCategory.map(multiState.c_str(), "Passed");
      data->addColumn( efficiencyCategory );
    }
    RooWorkspace* w = new RooWorkspace();
    //import the data
    w->import(*(data->get(0)));
    saveDistributionsPlot(w, data);
    if(data->numEntries() > 0) {
      //set the values of binnedVariables to the mean value in this data bin
      RooArgSet meanOfVariables;
      TIterator* vit = binnedVariables.createIterator();
      for(RooRealVar* v = (RooRealVar*)vit->Next(); v!=0; v = (RooRealVar*)vit->Next() ){
	meanOfVariables.addClone(*v);
      	double mean = data->mean(*v);
      	RooBinning binning((RooBinning&)v->getBinning());
      	int ind = binning.binNumber(mean);
      	RooRealVar& newVar = (RooRealVar&)meanOfVariables[v->GetName()];
      	newVar.setVal(mean);
      	newVar.setAsymError(binning.binLow(ind)-mean, binning.binHigh(ind)-mean);
      }
      delete vit;

      //put an entry in the efficiency dataset
      //note that the category values are coming from data->get(0)
      meanOfVariables.addClone(*data->get(0), true);     
      efficiency.setVal(0); //reset
      efficiency.setAsymError(0,0);
      
      if (docutandcount) {
      	doCntEfficiency(w, data, efficiency);
      	cntEfficiency.add( RooArgSet(meanOfVariables, efficiency) );
      } else {
	doFitEfficiency(w, data, pdfNames[category], efficiency);
	fitEfficiency.add( RooArgSet(meanOfVariables, efficiency) );
      }
    }

    //save the workspace if requested
    if(saveWorkspace){
      w->Write("w");
    }
    
    //clean up
    delete data;
    delete w;

    //get back to the initial directory
    gDirectory->cd("..");
  }
  
  if (docutandcount) {
    //cntEfficiency.Write();
    gDirectory->mkdir("cnt_eff_plots")->cd();
    saveEfficiencyPlots(cntEfficiency, effName, binnedVariables, mappedCategories);
    gDirectory->cd("..");
  } else {
    //fitEfficiency.Write();
    gDirectory->mkdir("fit_eff_plots")->cd();
    saveEfficiencyPlots(fitEfficiency, effName, binnedVariables, mappedCategories);
    gDirectory->cd("..");
  }
  
  //empty string means no error
  return "";
}

void TagProbeFitter::doFitEfficiency(RooWorkspace* w, RooAbsData* data, string pdfName, RooRealVar& efficiency) {
  
  //if pdfName is empty skip the fit
  if(pdfName == ""){
    return;
  }

  //create the simultaneous pdf of name pdfName
  createPdf(w, pdfs[pdfName]);
  //set the initial values for the yields of signal and background
  setInitialValues(w, data);  
  RooFitResult* res = 0;
  
  if (binnedFit) { 
    // get variables from data, which contain also other binning or expression variables
    const RooArgSet *dataObs = data->get(0); 
    // remove everything which is not a dependency of the pdf
    RooArgSet *obs = w->pdf("simPdf")->getObservables(dataObs);
    RooDataHist *bdata = new RooDataHist("data_binned", "data_binned", *obs, *data);
    
    w->import(*bdata);
    data = bdata;
    delete obs;
  }

  double totPassing = data->sumEntries("_efficiencyCategory_==_efficiencyCategory_::Passed");
  double totFailing = data->sumEntries("_efficiencyCategory_==_efficiencyCategory_::Failed");
  
  //******* The block of code below is to make the fit converge faster.
  // ****** This part is OPTIONAL, i.e., off be default. User can activate this
  // ****** by setting the following parameters: "fixVars" and "floatShapeParameters"
  // ****** Here is the full logic:
  /////   ---> if "floatShapeParameters==true" && "fixVars is empty" : 
  ////             Just perform the fit without any of these optimizations (this is the default logic)
  /////   ---> if "floatShapeParameters==true" && "fixVars is NOT empty" : 
  ////             In each bin, fix the desired parameters and perform the fit to obtain a good starting point.
  /////            Then float these parameters and fit again.
  ////    ---> if "floatShapeParameters==false" && "fixVars is empty" : 
  ////             Do not perform any fit at all. Just print error message.
  ////    ---> if "floatShapeParameters==false" && "fixVars is NOT empty" : 
  ///              Perform a global fit to the whole sample, save the fitted values of the 
  ///              user specified parameters, and fix them for bin-by-bin fit. 
  
  if(!fixVars.empty()){
    // calculate initial values for parameters user want to fix
    if(!floatShapeParameters && fixVarValues.empty()){
      // we want to fix these parameters for each bin.
      // the following sequence fixes them, fits, releases and fits again
      // to get reasonable values. 
      // ----------------------------------------------------------------------
      // This procedure works only once with a whole dataset (without binning)
      // ----------------------------------------------------------------------
      // fix them
      varFixer(w,true);
      //do fit 
      w->pdf("simPdf")->fitTo(*data, Minimizer("migrad"), Save(true), SumW2Error(true), Extended(true), NumCPU(numCPU), PrintLevel(quiet?-1:0), PrintEvalErrors(quiet?-1:1));
      //release vars
      varFixer(w,false);
      //do fit 
      w->pdf("simPdf")->fitTo(*data, Minimizer("migrad"), Save(true), SumW2Error(true) , Extended(true), NumCPU(numCPU), PrintLevel(quiet?-1:0), PrintEvalErrors(quiet?-1:1));
      //save vars
      varSaver(w);
      // now we have a starting point. Fit will converge faster.
    }
    
    // here we can use initial values if we want (this works for each bin)
    if(!floatShapeParameters) varRestorer(w);  //restore vars
    
    // if we don't want to "floatShapeParameters" we just fix, fit, 
    //  release, and fit again. No need for global fitting above.
    //fix vars
    varFixer(w,true);
    //do fit
    res = w->pdf("simPdf")->fitTo(*data, Save(true), Hesse(false), Extended(true), NumCPU(numCPU), Minos(*w->var("efficiency")), PrintLevel(quiet?-1:0), PrintEvalErrors(quiet?-1:1), SumW2Error(true));
  } else {
    //release vars
    varFixer(w,false);
    
    //do fit
    res = w->pdf("simPdf")->fitTo(*data, Save(true), Hesse(false),Extended(true), Minos(*w->var("efficiency")), Strategy(2), PrintEvalErrors(quiet?-1:1), NumCPU(numCPU), SumW2Error(true));
  }
  
  // save everything
  res->Write("fitresults");
  w->saveSnapshot("finalState",w->components());
  saveFitPlot(w);

  //extract the efficiency parameter from the results
  RooRealVar* e = (RooRealVar*) res->floatParsFinal().find("efficiency");
  efficiency.setVal(e->getVal());
  Double_t errLo = e->getErrorLo(), errHi = e->getErrorHi();
  
  if (errLo == 0 or errHi == 0) {
    float numSignal = ((RooRealVar*) res->floatParsFinal().find("numSignalAll"))->getVal();
    float pass = numSignal*e->getVal();
    float fail = numSignal*(1-e->getVal());
    
    // Use Clopper-Pearson           
    double alpha = (1.0 - .68540158589942957)/2;
    double lo = (pass == 0) ? 0.0 : ROOT::Math::beta_quantile(   alpha, pass,   fail+1 );
    double hi = (fail == 0) ? 1.0 : ROOT::Math::beta_quantile( 1-alpha, pass+1, fail   );
    
    efficiency.setAsymError(lo-e->getVal(), hi-e->getVal());
  } else {
    efficiency.setAsymError(errLo, errHi);
  }

  if (totPassing * totFailing == 0) {
    RooRealVar* nS = (RooRealVar*) res->floatParsFinal().find("numSignalAll");
    //RooRealVar* nB = (RooRealVar*) res->floatParsFinal().find(totPassing != 0 ? "numBackgroundPass" : "numBackgroundFail");
    double cerr = ROOT::Math::beta_quantile( 1-(1.0-.68540158589942957)/2, 1, nS->getVal() ); 
    /*
    std::cout << "======================================================================================" << std::endl;
    std::cout << "======= totPassing "  << totPassing << ", totFailing " << totFailing << std::endl;
    std::cout << "======= FIT: e  "  <<  e->getVal() << ",  e Lo " << e->getErrorLo()  << ",  e Hi " <<  e->getErrorHi() << std::endl;
    std::cout << "======= FIT:nS  "  << nS->getVal() << ", nS Lo " << nS->getErrorLo() << ", nS Hi " << nS->getErrorHi() << std::endl;
    std::cout << "======= FIT:nB  "  << nB->getVal() << ", nB Lo " << nB->getErrorLo() << ", nB Hi " << nB->getErrorHi() << std::endl;
    std::cout << "======= CNT:    "  << cerr << std::endl;
    std::cout << "======================================================================================" << std::endl;
    */
    if (totPassing == 0) {
      efficiency.setVal(0);
      efficiency.setAsymError(0,cerr);
    } else {
      efficiency.setVal(1);
      efficiency.setAsymError(-cerr,0);
    }
  }
}

void TagProbeFitter::createPdf(RooWorkspace* w, vector<string>& pdfCommands){

  // create the signal and background pdfs defined by the user
  for(unsigned int i=0; i<pdfCommands.size(); i++){
    const std::string & command = pdfCommands[i];
    
    if (command.find("#import ") == 0) {
      TDirectory *here = gDirectory;
      w->import(command.substr(8).c_str());
      here->cd();
    } else if (command.find("Templated") == 0) {
      w->factory(command.c_str());
      w->Print();
      
    } else {
      TDirectory *here = gDirectory; 
      w->factory(command.c_str());
      here->cd();
    }
  } 
  
  // setup the simultaneous extended pdf
  w->factory("expr::numSignalPass('efficiency*numSignalAll', efficiency, numSignalAll[10.,10000000000])");
  w->factory("expr::numSignalFail('(1-efficiency)*numSignalAll', efficiency, numSignalAll)");
  TString sPass = "signal", sFail = "signal";
  if (w->pdf("signalPass") != 0 && w->pdf("signalFail") != 0) {
    if (w->pdf("signal") != 0) throw std::logic_error("You must either define one 'signal' PDF or two PDFs ('signalPass', 'signalFail'), not both!"); 
    sPass = "signalPass"; sFail = "signalFail";
  } else if (w->pdf("signal") != 0) {
    if (w->pdf("signalPass") != 0 || w->pdf("signalFail") != 0) {
      throw std::logic_error("You must either define one 'signal' PDF or two PDFs ('signalPass', 'signalFail'), not both!"); 
    }
  } else {
    throw std::logic_error("You must either define one 'signal' PDF or two PDFs ('signalPass', 'signalFail')");
  }
  
  w->factory("SUM::pdfPass(numSignalPass*"+sPass+", numBackgroundPass[0.,1000000000]*backgroundPass)");
  w->factory("SUM::pdfFail(numSignalFail*"+sFail+", numBackgroundFail[0.,1000000000]*backgroundFail)");
  w->factory("SIMUL::simPdf(_efficiencyCategory_, Passed=pdfPass, Failed=pdfFail)");
  // signalFractionInPassing is not used in the fit just to set the initial values
  if (w->pdf("simPdf") == 0) throw std::runtime_error("Could not create simultaneous fit pdf.");
  if(w->var("signalFractionInPassing") == 0)
    w->factory("signalFractionInPassing[0.9]");
}

void TagProbeFitter::setInitialValues(RooWorkspace* w, RooAbsData* data){
  // calculate initial values
  double signalEfficiency = w->var("efficiency")->getVal();
  double signalFractionInPassing = w->var("signalFractionInPassing")->getVal();
  //double totPassing = w->data("data")->sumEntries("_efficiencyCategory_==_efficiencyCategory_::Passed");
  //double totFailinging = w->data("data")->sumEntries("_efficiencyCategory_==_efficiencyCategory_::Failed");

  double totPassing = data->sumEntries("_efficiencyCategory_==_efficiencyCategory_::Passed");
  double totFailinging = data->sumEntries("_efficiencyCategory_==_efficiencyCategory_::Failed");

  double numSignalAll = totPassing*signalFractionInPassing/signalEfficiency;
  // check if this value is inconsistent on the failing side
  if(numSignalAll*(1-signalEfficiency) > totFailinging)
    numSignalAll = totFailinging;
  // now set the values
  w->var("numSignalAll")->setVal(numSignalAll);
  w->var("numBackgroundPass")->setVal(totPassing - numSignalAll*signalEfficiency);
  w->var("numBackgroundFail")->setVal(totFailinging -  numSignalAll*(1-signalEfficiency));
  
  if (totPassing == 0) {
    w->var("efficiency")->setVal(0.0);
    w->var("efficiency")->setAsymError(0,1);
    w->var("efficiency")->setConstant(false);
    w->var("numBackgroundPass")->setVal(0.0);
    w->var("numBackgroundPass")->setConstant(true);
    w->var("numBackgroundFail")->setConstant(false);
  } else if (totFailinging == 0) {
    w->var("efficiency")->setVal(1.0);
    w->var("efficiency")->setAsymError(-1,0);
    w->var("efficiency")->setConstant(false);
    w->var("numBackgroundPass")->setConstant(false);
    w->var("numBackgroundFail")->setVal(0.0);
    w->var("numBackgroundFail")->setConstant(true);
  } else {
    w->var("efficiency")->setConstant(false);
    w->var("numBackgroundPass")->setConstant(false);
    w->var("numBackgroundFail")->setConstant(false);
  }
  
  // if signal fraction is 1 then set the number of background events to 0.
  RooRealVar* fBkgPass = w->var("numBackgroundPass");
  RooRealVar* fBkgFail = w->var("numBackgroundFail");
  if(signalFractionInPassing==1.0) { 
    fBkgPass->setVal(0.0); 
    fBkgPass->setConstant(true); 
    fBkgFail->setVal(0.0); 
    fBkgFail->setConstant(true); 
  }
  
  // save initial state for reference
  // FIXME
  //w->saveSnapshot("initialState",w->components());
}

void TagProbeFitter::saveFitPlot(RooWorkspace* w){
  // save refferences for convenience
  RooCategory& efficiencyCategory = *w->cat("_efficiencyCategory_");
  RooAbsData* dataAll = (binnedFit ? w->data("data_binned") : w->data("data") );
  RooAbsData* dataPass = dataAll->reduce(Cut("_efficiencyCategory_==_efficiencyCategory_::Passed")); 
  RooAbsData* dataFail = dataAll->reduce(Cut("_efficiencyCategory_==_efficiencyCategory_::Failed")); 
  RooAbsPdf& pdf = *w->pdf("simPdf");
  RooArgSet *obs = pdf.getObservables(*dataAll);
  RooRealVar* mass = 0;
  TIterator* it = obs->createIterator();
  for(RooAbsArg* v = (RooAbsArg*)it->Next(); v!=0; v = (RooAbsArg*)it->Next() ){
    if(!v->InheritsFrom("RooRealVar")) continue;
    mass = (RooRealVar*)v;
    break;
  }
  if(!mass) return;
  // make a 2x2 canvas
  TCanvas canvas("fit_canvas");
  canvas.Divide(2,2);
  vector<RooPlot*> frames;
  // plot the Passing Probes
  canvas.cd(1);
  if (massBins == 0) {
      frames.push_back(mass->frame(Name("Passing"), Title("Passing Probes")));
      frames.push_back(mass->frame(Name("Failing"), Title("Failing Probes")));
      frames.push_back(mass->frame(Name("All"),     Title("All Probes")));
  } else {
      frames.push_back(mass->frame(Name("Passing"), Title("Passing Probes"), Bins(massBins)));
      frames.push_back(mass->frame(Name("Failing"), Title("Failing Probes"), Bins(massBins)));
      frames.push_back(mass->frame(Name("All"),     Title("All Probes"), Bins(massBins)));
  }
  dataPass->plotOn(frames[0]);
  pdf.plotOn(frames[0], Slice(efficiencyCategory, "Passed"), ProjWData(*dataPass), LineColor(kGreen));
  pdf.plotOn(frames[0], Slice(efficiencyCategory, "Passed"), ProjWData(*dataPass), LineColor(kGreen), Components("backgroundPass"), LineStyle(kDashed));
  frames[0]->Draw();
  // plot the Failing Probes
  canvas.cd(2);
  dataFail->plotOn(frames[1]);
  pdf.plotOn(frames[1], Slice(efficiencyCategory, "Failed"), ProjWData(*dataFail), LineColor(kRed));
  pdf.plotOn(frames[1], Slice(efficiencyCategory, "Failed"), ProjWData(*dataFail), LineColor(kRed), Components("backgroundFail"), LineStyle(kDashed));
  frames[1]->Draw();
  // plot the All Probes
  canvas.cd(3);
  dataAll->plotOn(frames.back());
  pdf.plotOn(frames.back(), ProjWData(*dataAll), LineColor(kBlue));
  pdf.plotOn(frames.back(), ProjWData(*dataAll), LineColor(kBlue), Components("backgroundPass,backgroundFail"), LineStyle(kDashed));
  frames.back()->Draw();
  // plot the Fit Results
  canvas.cd(4);
  frames.push_back(mass->frame(Name("Fit Results"), Title("Fit Results")));
  pdf.paramOn(frames.back(), dataAll, "", 0, "NELU", 0.1, 0.9, 0.9);
  // draw only the parameter box not the whole frame
  frames.back()->findObject(Form("%s_paramBox",pdf.GetName()))->Draw();
  //save and clean up
  canvas.Write();
  for (size_t i=0; i<frames.size(); i++) {
    delete frames[i];
  }
  delete dataPass;
  delete dataFail;
}

void TagProbeFitter::saveDistributionsPlot(RooWorkspace* w, RooAbsData* data) {
  
  // save pointers to datasets to manage memory
  //RooAbsData* dataAll = w->data("data");
  RooAbsData* dataPass = data->reduce(Cut("_efficiencyCategory_==_efficiencyCategory_::Passed")); 
  RooAbsData* dataFail = data->reduce(Cut("_efficiencyCategory_==_efficiencyCategory_::Failed")); 

  const RooArgSet* vars = data->get();
  vector<RooRealVar*> reals;
  TIterator* it = vars->createIterator();
  for(RooAbsArg* v = (RooAbsArg*)it->Next(); v!=0; v = (RooAbsArg*)it->Next() ){
    if(!v->InheritsFrom("RooRealVar")) continue;
    reals.push_back((RooRealVar*)v);
  }
  TCanvas canvas("distributions_canvas");
  canvas.Divide(3,reals.size());
  vector<RooPlot*> frames;
  for(unsigned int i=0; i<reals.size(); i++){
    // plot the Passing Probes
    canvas.cd(3*i+1);
    frames.push_back(reals[i]->frame(Name("Passing"), Title("Passing Probes"), Bins(100)));
    dataPass->plotOn(frames.back(), LineColor(kGreen));
    dataPass->statOn(frames.back());
    frames.back()->Draw();
    // plot the Failing Probes
    canvas.cd(3*i+2);
    frames.push_back(reals[i]->frame(Name("Failing"), Title("Failing Probes"), Bins(100)));
    dataFail->plotOn(frames.back(), LineColor(kRed));
    dataFail->statOn(frames.back());
    frames.back()->Draw();
    // plot the All Probes
    canvas.cd(3*i+3);
    frames.push_back(reals[i]->frame(Name("All"), Title("All Probes"), Bins(100)));
    data->plotOn(frames.back(), LineColor(kBlue));
    data->statOn(frames.back());
    frames.back()->Draw();
  }
  canvas.Write();
  for (size_t i=0; i<frames.size(); i++) {
    delete frames[i];
  }
  delete dataPass;
  delete dataFail;
}

void TagProbeFitter::saveEfficiencyPlots(RooDataSet& eff, const TString& effName, RooArgSet& binnedVariables, RooArgSet& mappedCategories){
  TIterator* v1it = binnedVariables.createIterator();
  for(RooRealVar* v1 = (RooRealVar*)v1it->Next(); v1!=0; v1 = (RooRealVar*)v1it->Next() ){
    RooArgSet binCategories1D;
    TIterator* v2it = binnedVariables.createIterator();
    for(RooRealVar* v2 = (RooRealVar*)v2it->Next(); v2!=0; v2 = (RooRealVar*)v2it->Next() ){
      if(v2 == v1) continue;
      binCategories1D.addClone( RooBinningCategory(TString(v2->GetName())+"_bins", TString(v2->GetName())+"_bins", *v2) );

      RooArgSet binCategories2D;
      TIterator* v3it = binnedVariables.createIterator();
      for(RooRealVar* v3 = (RooRealVar*)v3it->Next(); v3!=0; v3 = (RooRealVar*)v3it->Next() ){
        if(v3 == v1 || v3 == v2) continue;
        binCategories2D.addClone( RooBinningCategory(TString(v3->GetName())+"_bins", TString(v3->GetName())+"_bins", *v3) );
      }  
      RooMultiCategory allCats2D("allCats2D", "allCats2D", RooArgSet(binCategories2D, mappedCategories));
      if(allCats2D.numTypes()==0){
        makeEfficiencyPlot2D(eff, *v1, *v2, TString::Format("%s_%s_PLOT", v1->GetName(), v2->GetName()), "", effName);
      } else {
        RooDataSet myEff(eff);
        myEff.addColumn(allCats2D);
        TIterator* catIt = allCats2D.typeIterator();
        for(RooCatType* t = (RooCatType*)catIt->Next(); t!=0; t = (RooCatType*)catIt->Next() ){
          TString catName = t->GetName();
          if(catName.Contains("NotMapped")) 
	    continue;
          catName.ReplaceAll("{","").ReplaceAll("}","").ReplaceAll(";","_&_");
          //RooDataSet* eff_bin = (RooDataSet*) myEff.reduce( Cut(TString::Format("allCats2D==%d",t->getVal())) );       
          makeEfficiencyPlot2D(myEff, *v1, *v2, TString::Format("%s_%s_PLOT_%s",v1->GetName(), v2->GetName(), catName.Data()), catName, effName);
          //delete eff_bin;
        }        
      }
    }
    RooMultiCategory allCats1D("allCats1D", "allCats1D", RooArgSet(binCategories1D, mappedCategories));
    if(allCats1D.numTypes()==0){
      makeEfficiencyPlot1D(eff, *v1, TString::Format("%s_PLOT", v1->GetName()), "", effName);
    }else{
      RooDataSet myEff(eff);
      myEff.addColumn(allCats1D);
      TIterator* catIt = allCats1D.typeIterator();
      for(RooCatType* t = (RooCatType*)catIt->Next(); t!=0; t = (RooCatType*)catIt->Next() ){
        TString catName = t->GetName();
        if(catName.Contains("NotMapped")) 
	  continue;
        catName.ReplaceAll("{","").ReplaceAll("}","").ReplaceAll(";","_and_");
        //RooDataSet* eff_bin = (RooDataSet*) myEff.reduce( Cut(TString::Format("allCats1D==%d",t->getVal())) );
        makeEfficiencyPlot1D(myEff, *v1, TString::Format("%s_PLOT_%s", v1->GetName(), catName.Data()), catName, effName);
        //delete eff_bin;
      }
    }
  }
}

void TagProbeFitter::makeEfficiencyPlot1D(RooDataSet& eff, RooRealVar& v, const TString& plotName, const TString& plotTitle, const TString& effName){

  // FIXME x axis error bars are strange
 
  //TCanvas canvas(plotName);
  //const RooArgSet* set = eff.get();
  //RooRealVar* e = (RooRealVar*) set->find("efficiency");
  //RooPlot* p = v.frame(Name(plotName), Title(plotTitle));
  //eff.plotOnXY(p,YVar(*e));
  //p->SetYTitle(TString("Efficiency of ")+effName);
  //p->SetAxisRange(0,1,"Y");
  //p->Draw();
  //canvas.Write();
  //delete p;  

  //dataVars.Print();
  TCanvas canvas(plotName);
  const RooArgSet* set = eff.get();
  RooRealVar* e = (RooRealVar*) set->find("efficiency");
  RooRealVar* x = (RooRealVar*) set->find(v.GetName());
  
  Float_t xArray[eff.numEntries()];
  Float_t yArray[eff.numEntries()];
  //Float_t errxLowArray[eff.numEntries()];
  //Float_t errxHighArray[eff.numEntries()];
  Float_t erryLowArray[eff.numEntries()];
  Float_t erryHighArray[eff.numEntries()];
      
  for (int i=0; i<eff.numEntries(); i++) {
    eff.get(i);
    xArray[i] = x->getVal();
    yArray[i] = e->getVal();
    erryLowArray[i] = fabs(e->getErrorLo());
    erryHighArray[i] = e->getErrorHi();
    //errxLowArray[i] = fabs(x->getErrorLo());
    //erryHighArray[i] = x->getErrorHi();
  }
  canvas.cd();
  TGraphAsymmErrors* g = new TGraphAsymmErrors(eff.numEntries(), xArray, yArray, 0, 0, erryLowArray, erryHighArray);
  //TGraphAsymmErrors* g = new TGraphAsymmErrors(eff.numEntries(), xArray, yArray, errxLowArray, errxHighArray, erryLowArray, erryHighArray);

  g->SetName(plotName);
  g->SetTitle(plotTitle);
  g->GetYaxis()->SetTitle(TString("Efficiency of ")+effName);
  g->GetYaxis()->SetRangeUser(0,1);
  g->SetMarkerStyle(20);
  g->Draw("PAE");
  canvas.Write();
  delete g;  
}

void TagProbeFitter::makeEfficiencyPlot2D(RooDataSet& eff, RooRealVar& v1, RooRealVar& v2, const TString& plotName, const TString& plotTitle, const TString& effName){
  
  TCanvas canvas(plotName);
  canvas.SetRightMargin(0.15);
  TH2F* h = new TH2F(plotName, plotName, v1.getBinning().numBins(), v1.getBinning().array(), v2.getBinning().numBins(), v2.getBinning().array());
  const RooArgSet* set = eff.get();
  RooRealVar* e = (RooRealVar*) set->find("efficiency");
  RooRealVar* v1_ = (RooRealVar*) set->find(v1.GetName());
  RooRealVar* v2_ = (RooRealVar*) set->find(v2.GetName());
  h->SetTitle(TString::Format("%s;%s%s;%s%s;Efficiency of %s", plotTitle.Data(),
    v1.GetTitle(), TString(v1.getUnit()).Length()==0?"":TString::Format(" (%s)", v1.getUnit()).Data(),
    v2.GetTitle(), TString(v2.getUnit()).Length()==0?"":TString::Format(" (%s)", v2.getUnit()).Data(), effName.Data()));
  h->SetOption("colztexte");
  h->GetZaxis()->SetRangeUser(-0.001,1.001);
  h->SetStats(kFALSE);
  for(int i=0; i<eff.numEntries(); i++){
    eff.get(i);
    h->SetBinContent(h->FindBin(v1_->getVal(), v2_->getVal()), e->getVal());
    h->SetBinError(h->FindBin(v1_->getVal(), v2_->getVal()), (e->getErrorHi()-e->getErrorLo())/2.);
  }
  h->Draw();
  canvas.Write();
  delete h;
}

void TagProbeFitter::doSBSEfficiency(RooWorkspace* w, RooRealVar& efficiency){
}

void TagProbeFitter::doCntEfficiency(RooWorkspace* w, RooAbsData* data, RooRealVar& efficiency){
  int pass = data->sumEntries("_efficiencyCategory_==_efficiencyCategory_::Passed");
  int fail = data->sumEntries("_efficiencyCategory_==_efficiencyCategory_::Failed");
  double e = (pass+fail == 0) ? 0 : pass/double(pass+fail);
  std::cout << "passing: " << pass << " total: " << pass+fail << std::endl;
  std::cout << "efficiency: " << e << std::endl;
  // Use Clopper-Pearson
  double alpha = (1.0 - .68540158589942957)/2;
  double lo = (pass == 0) ? 0.0 : ROOT::Math::beta_quantile(   alpha, pass,   fail+1 );
  double hi = (fail == 0) ? 1.0 : ROOT::Math::beta_quantile( 1-alpha, pass+1, fail   );
  ////from TGraphAsymmErrors
  //double lob, hib;
  //Efficiency( pass, pass+fail, .68540158589942957, e, lob, hib );
  //std::cerr << "CNT " << pass << "/" << fail << ":  Clopper Pearson [" << lo << ", "  << hi << "], Bayes [" << lob << ", " << hib << "]" << std::endl;
  efficiency.setVal(e);
  efficiency.setAsymError(lo-e, hi-e);
}

void TagProbeFitter::varFixer(RooWorkspace* w, bool fix){
  std::vector<std::string>::const_iterator it;
  for(it=fixVars.begin(); it<fixVars.end(); it++){    
    if(w->var((*it).c_str()))
      w->var((*it).c_str())->setAttribute("Constant",fix);
    else{
      std::cout << "TagProbeFitter: " << "Can't find a variable to fix: " << *it;
    }
  }
}

void TagProbeFitter::varSaver(RooWorkspace* w){
  if(!fixVarValues.empty()){
    std::cout << "attempt to save variables more than once!" << std::endl;
    return;
  }
  std::vector<std::string>::const_iterator it;
  for(it=fixVars.begin(); it<fixVars.end(); it++){
    fixVarValues.push_back(w->var((*it).c_str())->getVal());
  }
  
}

void TagProbeFitter::varRestorer(RooWorkspace* w){
  if(fixVarValues.size()==fixVars.size())
    for(unsigned int i=0; i< fixVars.size(); i++){
      w->var(fixVars[i].c_str())->setVal(fixVarValues[i]);
    }
  else{
    std::cout << "fixVars and fixVarValues are not of the same size!" << std::endl; 
  }
}

