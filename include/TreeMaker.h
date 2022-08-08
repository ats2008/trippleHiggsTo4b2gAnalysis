#include "SCMVANtuple.h"
#include "Util.h"
#include "chrono"
#define NSTORAGE_ARRAY_MAX 5000

class genericMatrixStore
{
    public :
        genericMatrixStore(Int_t sizeX=8,Int_t sizeY=8)
        {
            matrixStore=new Float_t[sizeY*sizeX];
            nX=sizeX;
            nY=sizeY;
        }
        ~genericMatrixStore()
        {
            if(matrixStore!=nullptr)
             delete matrixStore;
        }
        void initialize(Float_t def=0.0)
        {
            for(Int_t i; i <nX;i++)
            for(Int_t j; j <nY;j++)
            {
                matrixStore[i*nY+j]=def;
            }
        }
        Float_t get(Int_t x, Int_t y)
        {
            return matrixStore[x*nY+y];
        }
        void fill(Int_t x, Int_t y, Float_t val)
        {
            matrixStore[x*nY+y]=val;
        }
        Float_t * matrixStore;
        Int_t nX,nY;
};

class TreeMaker 
{
    public  :
    
    // Data members
    std::vector<string> InFileList;
    string ofileName;
    string treeName;
    string prefix;
    TFile * outputFile ;
    TChain *treeChain ;
    SCMVANtuple ntupleRawTree;
    Double_t scAbsEtaMax,scAbsEtaMin;
    Int_t genParticlePDGID;
    Double_t drGenMatchMin;
    Double_t genParticlePtMin;
    Double_t eventGenMultiplicity;
    Bool_t genParticleIsStable;   
    // Storage Vars
        // Defenition of the storage of type Double_t
    Int_t storageIdxFilledDouble; 
    std::map<string, Float_t*> storageArrayDoubleMap;
    Double_t *storageArrayDouble;
    std::map<string, Int_t > candidateMapDouble;
        // Defenition of the storage of type Int_t
    Int_t storageIdxFilledInt;
    Int_t* storageArrayInt ;
    std::map<string, Int_t > candidateMapInt;
    std::map<string, Float_t > storageFloat;
    
    Int_t* photonSelectionCheck;
    
    // Histograms to Store
    std::map<TString,TH1F*> th1fStore;
    
    // Trees to store
    std::map<TString,TTree*> treeStore;

    // isMC
    bool isMC;
    bool doGenMatching;
    
    // book keeping vars
    bool initDone;
    std::vector<TTree*> treesToStore;
    Long64_t nentries, maxEvents ;
    Int_t reportEvery;
                    // Member functions 
    // Constructor
    TreeMaker();
    
    void Init( string cfgFileName);

    // Helper funtions
    void readParameters(string fname);
    void setupInputTree();
    void setupOutputSCTree();
    void AllocateMemory();
    
    void SaveFile();
    void setupOutPuts(bool makeTreeCopy=false);
    void SetupAnalysis(bool makeTreeCopy=false);
    
    void AddSCHistos(TString tag="");
    void AddSCTree(TString SCTreeName="SCTreeStorage");
    
    // Histogram Related Functions
    void bookHistograms();
    void fill_scHists(Int_t idx =-1);
    void fillSCVariablesToOutTree(Int_t scIDX,TString SCTreeName);
    
    void fill_genHists(Int_t idx);
    void fill_scHists(Int_t scIDX,TString tag="",Double_t dr=0.0);
    void fill_eventHists();
    void DataSCMaker();
    void genParticleSCMaker();
    void Pi0ParticleSCMaker();  
    void genParticleBMMGSCMaker();
    void genParticleBMMGSCMakerJPsiGamma();
    void AssignECALClusterVariables( Int_t scIdx);
    void AssignHCALClusterVariables( Int_t scIdx);
    void AssignPFVariables(Int_t) ;
  private :
     int   pTBins;
     float pTmin;
     float pTmax;

     int   etaBins;
     float etamin;
     float etamax;

     int   deltaRNBins;
     float deltaRMin;
     float deltaRMax;



};

TreeMaker::TreeMaker()
{    
     pTBins=50;
     pTmin=0.0;
     pTmax=50.0;

     etaBins=70;
     etamin=-3.5;
     etamax=3.5;

     deltaRNBins=300;
     deltaRMin=0.0;
     deltaRMax=3.0;


    initDone=false;
}

void TreeMaker::Init(string cfgFileName)
{
    treeName="mergedTree";
    prefix="workarea/";
    ofileName="output.root";
    InFileList.clear();
    ofileName="output.root";
    maxEvents=1000;
    isMC=false;
    doGenMatching=false;
    drGenMatchMin=0.1;
    reportEvery=10000;
    scAbsEtaMin=-1e2;
    scAbsEtaMax=1e9;
    genParticlePDGID=22;
    genParticleIsStable=true;
    genParticlePtMin=3.0;

    readParameters(cfgFileName);
    AllocateMemory();
    
    initDone=true;
}

void TreeMaker::SetupAnalysis(bool makeTreeCopy)
{

    if( not initDone) 
    {
        std::cout<<"\n Init ur vars before setupAnalysis !! \n";
        return;
    }
    
    setupInputTree();
    bookHistograms();
    setupOutPuts(makeTreeCopy);
}

void TreeMaker::setupInputTree()
{   
    if( not initDone) 
    {
        std::cout<<"\n Init ur vars before setupInputTree !! \n";
        return;
    }

    treeChain = new TChain(treeName.c_str());
    auto rslt=0;
    for(auto i=0;i<InFileList.size();i++)
    {
        rslt=treeChain->AddFile(InFileList[i].c_str(),-1);
        if(rslt!=1) exit(112);
    }
    
    nentries = treeChain->GetEntries();
    if( maxEvents < 0) maxEvents = nentries;
    maxEvents = nentries < maxEvents ? nentries: maxEvents;
    cout<<"Available total number of events "<<nentries<<" \n";
    
    //ntupleRawTree.Init(treeChain,isMC);
    ntupleRawTree.Init(treeChain);

}

void TreeMaker::setupOutPuts(bool makeTreeCopy)
{
   outputFile = new  TFile((prefix+ofileName).c_str(),"recreate");    
    
   
   if(makeTreeCopy)
   {
        setupOutputSCTree();
   }

}

void TreeMaker::AllocateMemory()
{

    // Defenition of the storage of type Double_t
    
    storageArrayDouble=new Double_t[NSTORAGE_ARRAY_MAX];
    std::cout<<"Allocated "<<sizeof(Double_t)*NSTORAGE_ARRAY_MAX/1024<<" kB of storage for "<<NSTORAGE_ARRAY_MAX<<" Doubles \n";
    storageIdxFilledDouble=0; 
    
    storageArrayInt =  new Int_t[NSTORAGE_ARRAY_MAX]; 
    std::cout<<"Allocated "<<sizeof(Int_t)*NSTORAGE_ARRAY_MAX/1024<<" kB of storage for "<<NSTORAGE_ARRAY_MAX<<" Ints\n";
    storageIdxFilledInt=0;

    candidateMapDouble["SCTreeStorage"]   = storageIdxFilledDouble ;
    storageIdxFilledDouble+=500;
}

void TreeMaker::SaveFile()
{
    outputFile->cd();

    for (std::map<TString,TH1F *>::iterator it=th1fStore.begin() ; it!=th1fStore.end(); ++it)
    {
        
        //std::cout<<"Writing "<<it->first<<" to file ! \n";
        auto &ahist = *(it->second); 
        ahist.Write();
    }
       
    for (std::map<TString,TTree *>::iterator it=treeStore.begin() ; it!=treeStore.end(); ++it)
    {
        
        //std::cout<<"Writing "<<it->first<<" to file ! \n";
        auto &atree = *(it->second); 
        atree.Write();
    }

    outputFile->Write();
    outputFile->Purge();
    outputFile->Close();

}

void TreeMaker::readParameters(string fname)
{
    fstream cfgFile(fname,ios::in);
	string line;
	bool cfgModeFlag=false;

    Double_t aDouble;
    
    cfgFile.clear();
    cfgFile.seekg(0,ios::beg);
    cfgModeFlag=false;
    std::istringstream strStream;
    std::string field;
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
            if(field.compare("OutputPrefix")==0){
                 getline(strStream, field);
                 prefix=field;
                 cout<<" setting prefix = "<<prefix<<"\n";
            }
           if(field.compare("InputTreeName")==0){
                 getline(strStream, field);
                 treeName=field;
                 cout<<" setting treeName  = "<<prefix<<"\n";
            }
            if(field.compare("ReportEvery")==0){
                 getline(strStream, field);
                 reportEvery=std::atoi(field.c_str());
                 cout<<" setting reportEvery  = "<<reportEvery<<"\n";
            }
             if(field.compare("MaxEvents")==0){
                 getline(strStream, field);
                 maxEvents=std::atoi(field.c_str());
                 cout<<" setting maxEvents  = "<<maxEvents<<"\n";
            }
            if(field.compare("DrGenMatchMin")==0){
                 getline(strStream, field);
                 drGenMatchMin=std::atof(field.c_str());
                 cout<<" setting drGenMatchMin  = "<<drGenMatchMin<<"\n";
            }
            if(field.compare("SCAbsEtaMax")==0){
                 getline(strStream, field);
                 scAbsEtaMax=std::atof(field.c_str());
                 cout<<" setting scAbsEtaMax  = "<<scAbsEtaMax<<"\n";
            }
            if(field.compare("GenParticlePtMin")==0){
                 getline(strStream, field);
                 genParticlePtMin=std::atof(field.c_str());
                 cout<<" setting genParticlePtMin  = "<<genParticlePtMin<<"\n";
            }
             if(field.compare("SCAbsEtaMin")==0){
                 getline(strStream, field);
                 scAbsEtaMin=std::atof(field.c_str());
                 cout<<" setting scAbsEtaMin  = "<<scAbsEtaMin<<"\n";
            }
            if(field.compare("GenParticlePDGID")==0){
                 getline(strStream, field);
                 genParticlePDGID=std::atoi(field.c_str());
                 cout<<" setting genParticlePDGID  = "<<genParticlePDGID<<"\n";
            }
            if(field.compare("GenParticleIsStable")==0){
                 getline(strStream, field);
                 tmpI=std::atoi(field.c_str());
                 genParticleIsStable= tmpI >0 ? 1 : 0;
                 cout<<" setting genParticleIsStable  = "<<genParticleIsStable<<"\n";
            }
            if(field.compare("IsMC")==0){
                 getline(strStream, field);
                 tmpI=std::atoi(field.c_str());
                 isMC = tmpI >0 ? 1 : 0;
                 cout<<" setting isMC  = "<<isMC<<"\n";
            }
            if(field.compare("DoGenMatching")==0){
                 getline(strStream, field);
                 tmpI=std::atoi(field.c_str());
                 doGenMatching= tmpI >0 ? 1 : 0;
                 cout<<" setting doGenMatching  = "<<doGenMatching<<"\n";
            }

       }
    }

	// getting flists
    cfgFile.clear();
	cfgFile.seekp(ios::beg);
    cfgModeFlag=false;
	while(std::getline(cfgFile,line))
	{
	   if(line=="#FILELIST_BEG") {cfgModeFlag=true;continue;}
	   if(line=="#FILELIST_END") {cfgModeFlag=false;continue;}
	   if(not cfgModeFlag) continue;
	   if(line=="") continue;
	   if(line=="#END") break;
	   InFileList.push_back(line);
	}

    std::cout<<"File List has the following files : \n";
    for( auto name : InFileList)
    {
        std::cout<<"\t"<<name<<"\n";
    }

}

#include "TreeMaker.cc"


