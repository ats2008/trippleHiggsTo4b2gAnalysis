import ROOT
import numpy as np
from array import array


class TMVAModel :
    
    def __init__(self):
            
        #  reader =  new TMVA::Reader( "!Color:!Silent" );
        self.init()
        pass

    def init(self):
        self.tmva_modelReader = ROOT.TMVA.Reader()
        self.branchList       = list()
        self.branchStore      = {}
        self.mvaName          = ""
        self.hasMVASet        = False


    def setWeightFile( self , weightFile ):
        self.weightFile=weightFile
        print(weightFile)
        self.tmva_modelReader.BookMVA( self.mvaName , self.weightFile)
    
    def addBranches( self , branchList):
        self.branchList=branchList
        for ky in self.branchList:
            self.branchStore[ky]=array('f',[0.0])
            self.tmva_modelReader.AddVariable(ky,self.branchStore[ky])

    def addSpectatorBranches( self , branchList):
        self.specBranchList=branchList
        for ky in self.specBranchList:
            self.branchStore[ky]=array('f',[0.0])
            self.tmva_modelReader.AddSpectator(ky,self.branchStore[ky])

    def setupTMVAModel(self,mvaName,weightFile,branchList,specBranchList):
        self.mvaName = mvaName
        self.addBranches(branchList)
        self.addSpectatorBranches(specBranchList)
        self.setWeightFile(weightFile)
        self.hasMVASet        = True


    def predict_(self):
        return self.tmva_modelReader.EvaluateMVA(self.mvaName)

    def predict(self,valDict):
        if not self.hasMVASet:
            print("MVA not setup ! returning nan")
            return np.nan
        
        for ky in self.branchStore:
            if ky not in valDict:
               print("Input dictionary incomplete ! returning nan")
               return np.nan
        
        for ky in self.branchStore:
            self.branchStore[ky][0]=valDict[ky]
        
        return self.predict_()


