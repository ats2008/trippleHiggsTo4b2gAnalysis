#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

class MVATester
{
    public  :
    
    // Data members
    std::vector<string> InFileList;
    vector<string> treeNames;
    vector<string> tagNames;
    string ofileName;
    string prefix;
    TFile * outputFile;

    TFile* currentFile;
    TTree* currentTree;

    string testTrainConfig;
    
    std::map<std::string,int> Use;
    std::map<std::string,string> mvaMethordOptions;
    //std::map<string,Double_t> storageDouble;
    std::map<string,Float_t> storageDouble;
    std::map<string,Float_t> storageFloat;
    std::map<TString,TH1F*> th1fStore;

    std::vector<string> mvaTrainVars;
    std::vector<string> spectatorVars;
    
    string photonIdxMVAWeightFile;
    bool hasSetupPhotonMVA,hasWeightFiles;

    Long64_t nentries, maxEvents ;
    Int_t reportEvery;
    Bool_t initDone ;
    Double_t photonMVAValue;

    TMVA::Reader *reader ;
 // Constructor
    MVATester();
    ~MVATester();
    
    void Init( string cfgFileName);

    // Helper funtions
    void readParameters(string fname);
    
    void setUpMVA();
    void setUpTestTree(Int_t i);
    void SaveOutputs();
    void bookMVAHistos(TString tag);
    void fillMVAHistos(TString tag, Float_t );
    void closeCurrentFile();
    void analyzeAllData();
    Float_t getScore();
    void getVetorFillledFromConfigFile( fstream &cfgFile , std::vector<string> &vecList, string beginTag,string endTag, bool verbose=false);
};


MVATester::MVATester()
{
  initDone=false;
  outputFile=nullptr;
}

MVATester::~MVATester()
{
    if(outputFile) outputFile->Close();

    for (std::map<TString,TH1F *>::iterator it=th1fStore.begin() ; it!=th1fStore.end(); ++it)
    {
        //delete it->second;
    }

}


void MVATester::Init(string cfgFileName)
{
    prefix      =  "";
    ofileName   =  "output.root";
    maxEvents   =  1000;
    photonIdxMVAWeightFile="";
    
    hasWeightFiles=false;
    hasSetupPhotonMVA=false;
    reportEvery=1000;

    readParameters(cfgFileName);
    
    outputFile = new  TFile((prefix+ofileName).c_str(),"recreate");    
    
    setUpMVA();
    initDone=true;
}

void MVATester::setUpMVA()
{
    if(not hasWeightFiles)
    {
        std::cout<<"Weights are not provided exiting ! \n";
        exit(4);
    }

    std::cout<<" ******** Setting up MVA ID *********\n"
             <<"Weight file : "<<photonIdxMVAWeightFile<<"\n";
    reader =  new TMVA::Reader( "!Color:!Silent" );
    
    for(Int_t i=0;i<mvaTrainVars.size();i++)
    {
 //       std::cout<<"Adding mva vars "<<mvaTrainVars[i]<<" to "<<" reader  \n";
        storageDouble[mvaTrainVars[i]]=0.0;
        storageFloat[mvaTrainVars[i]] = storageDouble[mvaTrainVars[i]];
        reader->AddVariable(mvaTrainVars[i].c_str(), &(storageFloat[mvaTrainVars[i]]));
    }
    
  for(int i =0;i<spectatorVars.size();i++)
  {
 //       std::cout<<"Adding spectator "<<spectatorVars[i]<<" to "<<" reader \n";
        storageDouble[spectatorVars[i]]=0.0;
        storageFloat[spectatorVars[i]] = storageDouble[spectatorVars[i]];
        reader->AddSpectator(spectatorVars[i].c_str(),  &(storageFloat[spectatorVars[i]]));
  }

    reader->BookMVA("MVA", photonIdxMVAWeightFile );
    hasSetupPhotonMVA=true;
}

void MVATester::setUpTestTree( Int_t treeIdx )
{
    
    std::cout<<"Setting Current Tree to : "<<treeNames[treeIdx]<<" from "<<InFileList[treeIdx]<<"\n";
    currentFile = TFile::Open( InFileList[treeIdx].c_str() );
    currentTree =  (TTree * ) currentFile->Get(treeNames[treeIdx].c_str());

    if(not currentTree ) {
        std::cout<<"Empty tree given for MVA test Setup !! exiting .. "; std::cout<<"\n";
    }
    else {
        std::cout<<"Sucessfully set  Current Tree to : "<<treeNames[treeIdx]<<" from "<<InFileList[treeIdx]<<"\n";
    }

    for(Int_t i=0;i<spectatorVars.size();i++)
    {
        storageDouble[spectatorVars[i]]=0.0;
 //       std::cout<<" Setting spec var : "<<spectatorVars[i]<<"\n"; 
        currentTree->SetBranchAddress( spectatorVars[i].c_str(), &(storageDouble[ spectatorVars[i] ] ) ) ;
    }

    for(Int_t i=0;i<mvaTrainVars.size();i++)
    {
        storageDouble[mvaTrainVars[i]]=0.0;
 //       std::cout<<" Setting mva var : "<<mvaTrainVars[i]<<"\n"; 
        currentTree->SetBranchAddress( mvaTrainVars[i].c_str(), &(storageDouble[ mvaTrainVars[i] ] ) ) ;
    }

}

void MVATester::closeCurrentFile()
{
    currentFile->Close();
    currentTree=nullptr;
}


Float_t MVATester::getScore()
{   
    if(not hasSetupPhotonMVA)
    {
        std::cout<<" MVA not setup : "<<hasSetupPhotonMVA<<"\n";
        exit(8);
    }
    
    for(int i =0;i<mvaTrainVars.size();i++)
    {
        storageFloat[mvaTrainVars[i]] = storageDouble[mvaTrainVars[i]];
    }
    for(int i =0;i<spectatorVars.size();i++)
    {
        storageFloat[spectatorVars[i]] = storageDouble[spectatorVars[i]];
    }

    photonMVAValue = reader->EvaluateMVA("MVA");
    return photonMVAValue;
}

TString getFloatAsTstring(Float_t val,Int_t nDes=2)
{

    std::stringstream ss;
    ss << std::fixed << std::setprecision(nDes) << val;
    std::string mystring = ss.str();

    return TString(mystring.c_str());

}
void MVATester::bookMVAHistos(TString tag)
{
    TString trPrefixA=tag+"mvaVal";
    std::vector<float> edge;
    float ed=-1.0;
    for(Int_t i=0;i<9;i++)
    {
        edge.push_back(ed);
        ed+=0.25;
    }
    for(Int_t ii=0;ii<edge.size()-1;ii++)
    {       
           std::cout<<"Adding for edge : "<< edge[ii]<<" "<<edge[ii+1]<<"\n";
           auto name=trPrefixA+getFloatAsTstring(edge[ii])+"To"+getFloatAsTstring(edge[ii+1]); 
           th1fStore[name+"_mass"]= new TH1F(name+"_mass","mass", 220 , 90.0  , 200.0);
           th1fStore[name+"_mass"]->SetDirectory(0);
    }


    th1fStore[tag+"mvaScore"]= new TH1F(tag+"MVAScore","MVA Score", 220 , -1.1  , 1.1 );
    th1fStore[tag+"mvaScore"]->SetDirectory(0);
    th1fStore[tag+"SR_mvaScore"]= new TH1F(tag+"SR_MVAScore","MVA Score", 220 , -1.1  , 1.1 );
    th1fStore[tag+"SR_mvaScore"]->SetDirectory(0);
    th1fStore[tag+"CR_mvaScore"]= new TH1F(tag+"CR_MVAScore","MVA Score", 220 , -1.1  , 1.1 );
    th1fStore[tag+"CR_mvaScore"]->SetDirectory(0);

}

void MVATester::fillMVAHistos( TString tag, Float_t  mvaScore)
{
    TString trPrefixA=tag+"mvaVal";
    std::vector<float> edge;
    float ed=-1.0;
    for(Int_t i=0;i<9;i++)
    {
        edge.push_back(ed);
        ed+=0.25;
    }
     for(Int_t ii=0;ii<edge.size()-1;ii++)
    {   
        if(mvaScore > edge[ii+1])
            continue;   
        auto name=trPrefixA+getFloatAsTstring(edge[ii])+"To"+getFloatAsTstring(edge[ii+1]); 
        th1fStore[name+"_mass"]->Fill(storageDouble["CMS_hgg_mass"]);            
        break;
    }

    th1fStore[tag+"mvaScore"]->Fill(mvaScore);
    auto dr= sqrt( (storageDouble["M1jj"] -125.0)*(storageDouble["M1jj"] -125.0) + (storageDouble["M2jj"] -125.0)*(storageDouble["M2jj"] -125.0)  );
    if(dr<25.0)
    th1fStore[tag+"SR_mvaScore"]->Fill(mvaScore);  
    else if(dr < 50.0)
    th1fStore[tag+"CR_mvaScore"]->Fill(mvaScore);
}   

void MVATester::analyzeAllData()
{
    Float_t mvaVal;
    Long64_t maxEvents_;
    for(Int_t i=0;i< treeNames.size();i++)
    {
        setUpTestTree(i);
        maxEvents_=currentTree->GetEntries();
        if(maxEvents > -1)  maxEvents_ = maxEvents > maxEvents_ ? maxEvents_ : maxEvents;
        std::cout<<"Processing "<<maxEvents_<<" items from "<<tagNames[i]<<" in "<<treeNames[i]<<"\n";

        bookMVAHistos(tagNames[i]);
        auto t_start = std::chrono::high_resolution_clock::now();
        auto t_end = std::chrono::high_resolution_clock::now();

        for (Long64_t jentry=0; jentry<maxEvents_; jentry++)
        {   
            currentTree->GetEntry(jentry);
            if(jentry%reportEvery == 0 )
            {
                  t_end = std::chrono::high_resolution_clock::now();
                  std::cout<<"\t Processing Entry in event loop : "<<jentry<<" / "<<maxEvents_<<"  [ "<<100.0*jentry/maxEvents_<<"  % ]  "
                           << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
                           <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( maxEvents_ - jentry)/(1e-9 + jentry)* 0.001
                           <<std::endl;
            
            }

            mvaVal=getScore();
            fillMVAHistos(tagNames[i],mvaVal);
 //            std::cout<<"MVA VAL : "<<mvaVal<<"\n";
       }
       //closeCurrentFile();
   // break;
    }
}




void MVATester::SaveOutputs()
{    
    
    outputFile->cd();
    for (std::map<TString,TH1F *>::iterator it=th1fStore.begin() ; it!=th1fStore.end(); ++it)
    {
        std::cout<<"Writing "<<it->first<<" to file ! \n";
        auto &ahist = *(it->second); 
        std::cout<<ahist.GetName()<<"\n";
        ahist.Write();
    }
 
   // Save the output
   outputFile->Close();
}

void MVATester::readParameters(string fname)
{
    fstream cfgFile(fname,ios::in);

    Double_t aDouble;
    
    cfgFile.clear();
    cfgFile.seekg(0,ios::beg);
	bool cfgModeFlag=false;
    cfgModeFlag=false;
    std::istringstream strStream;
    std::string field;
	string line;
    Int_t tmpI;
    
	while(std::getline(cfgFile,line))
	{
	   if(line=="#PARAMS_BEG") {cfgModeFlag=true;continue;}
	   if(line=="#PARAMS_END") {cfgModeFlag=false;continue;}
	   if(not cfgModeFlag) continue;
	   if(line=="") continue;
	   if(line=="#END") break;
       strStream.clear();
       strStream.str(line);
       while (getline(strStream, field,'='))
       {
           if(field.compare("OutputFile")==0){
                 getline(strStream, field);
                 ofileName=field;
                 std::cout<<" setting ofileName = "<<ofileName<<"\n";
            }
            if(field.compare("WeightFile")==0){
                 getline(strStream, field);
                 photonIdxMVAWeightFile=field;
                 cout<<" setting photonIdxMVAWeightFile  = "<<photonIdxMVAWeightFile<<"\n";
                 hasWeightFiles=true;
            }
            if(field.compare("OutputPrefix")==0){
                 getline(strStream, field);
                 prefix=field;
                 cout<<" setting prefix = "<<prefix<<"\n";
            }
            if(field.compare("MaxEvents")==0){
                 getline(strStream, field);
                 maxEvents=std::atoi(field.c_str());
                 cout<<" setting maxEvents  = "<<maxEvents<<"\n";
            }
            if(field.compare("ReportEvery")==0){
                 getline(strStream, field);
                 reportEvery=std::atoi(field.c_str());
                 cout<<" setting reportEvery  = "<<reportEvery<<"\n";
            }
       }
    }

    getVetorFillledFromConfigFile(cfgFile, InFileList   , "#FILELIST_BEG"   , "#FILELIST_END", true);
    getVetorFillledFromConfigFile(cfgFile, treeNames    , "#TREENAMES_BEG"  , "#TREENAMES_END", true);
    getVetorFillledFromConfigFile(cfgFile, tagNames     , "#TAG_BEG"  , "#TAG_END", true);

    getVetorFillledFromConfigFile(cfgFile, mvaTrainVars   , "#MVAVARLIST_BEG", "#MVAVARLIST_END", true);
    getVetorFillledFromConfigFile(cfgFile, spectatorVars  , "#SPECTATORLIST_BEG", "#SPECTATORLIST_END", true);
    
    if(InFileList.size() != treeNames.size() )  {  std::cout<<__LINE__<<" config error !! \n"; exit(1) ;}
    if(InFileList.size() != tagNames.size()  )  {  std::cout<<__LINE__<<" config error !! \n"; exit(1) ;}
}

void MVATester::getVetorFillledFromConfigFile( fstream &cfgFile , std::vector<string> &vecList, string beginTag,string endTag, bool verbose)
{
	
    bool cfgModeFlag=false;
    cfgModeFlag=false;
    std::istringstream strStream;
    std::string field;
	string line;
    
    // getting flists
    cfgFile.clear();
	cfgFile.seekp(ios::beg);
    cfgModeFlag=false;
	int nItems(0);
    while(std::getline(cfgFile,line))
	{
	   if(line==beginTag) {cfgModeFlag=true;continue;}
	   if(line==endTag) {cfgModeFlag=false;continue;}
	   if(not cfgModeFlag) continue;
	   if(line=="") continue;
	   if(line=="#END") break;
	   vecList.push_back(line);
	   nItems++;
    }

    if(verbose)
    {
       std::cout<<" Added "<<nItems<<" between "<<beginTag<<" and "<< endTag<<"\n";
       for( auto &name : vecList)
       {
           std::cout<<"\t"<<name<<"\n";
       }
    }

}

