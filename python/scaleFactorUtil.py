import ROOT ,copy
import numpy as np
import pickle
import matplotlib.pyplot as plt
import mplhep as hep

import Util as utl
from hep_ml.reweight import BinsReweighter, GBReweighter, FoldingReweighter
from hep_ml.preprocessing import IronTransformer
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import roc_curve ,roc_auc_score,auc
from sklearn.utils import shuffle
class mcScaler:
    def __init__(self):
        self.dataHist={}
        self.mcHist={}
        self.binEdges=[]
        self.scaleFactors=[]
        self.def_scaleFactor=0
        self.underflow_x={}
        self.overflow_x={}
        self.underflow_scaleFactor={}
        self.overflow_scaleFactor={}
        self.scaleFactorHist={}
        self.scaleFactorHistSet={}
        self.File=None
    
    def setSFHist(self,hist,cat='def'):
        self.scaleFactorHist[cat]=hist.Clone()
        self.scaleFactorHist[cat].SetName(cat+"_sfactorHist")
        self.underflow_x[cat]= self.scaleFactorHist[cat].GetBinLowEdge(1)
        self.overflow_x[cat]= self.scaleFactorHist[cat].GetBinLowEdge(self.scaleFactorHist[cat].GetNbinsX()+1)
        self.underflow_scaleFactor[cat] = self.scaleFactorHist[cat].GetBinContent(1)
        self.overflow_scaleFactor[cat] = self.scaleFactorHist[cat].GetBinContent(self.scaleFactorHist[cat].GetNbinsX())
        print("  Setting scale factor for category " , cat  )
        print("     Scale factogetPtDependentScaleFactorrs defined for x = [ ",self.underflow_x[cat],self.overflow_x[cat]," ]")
        print("     Under flow Scale factor     =  ",self.underflow_scaleFactor[cat])
        print("     Over flow Scale factor      = ",self.overflow_scaleFactor[cat])

        self.scaleFactorHistSet[cat]=True
    
    def setSFHistFromFile(self,fname,histNames={}):
        self.File=ROOT.TFile(fname,"READ")
        print(histNames)
        for cat in histNames:
            print("CAT : ",cat," || Setting SF hist ",histNames[cat]," from ",fname)
            hist=self.File.Get(histNames[cat])
            self.setSFHist(hist,cat)
    
    def getSFForX(self,x,cat='def'):
        if cat not in self.scaleFactorHistSet:
            print("\tWARNING ! : scaleFactor Hist not set , returning -1.0  [ cat : ",cat," available cats : ",self.scaleFactorHist.keys(),"]")
            return -1.0
        if x < self.underflow_x[cat]:
            return self.underflow_scaleFactor[cat]
        if x > self.overflow_x[cat]:
            return self.overflow_scaleFactor[cat]
        bid=self.scaleFactorHist[cat].FindBin(x)
        #print("bid = ",bid," for x ",x, " Bin edges --> ",
        #         self.scaleFactorHist.GetBinLowEdge(bid), 
        #         self.scaleFactorHist.GetBinLowEdge(bid)+self.scaleFactorHist.GetBinWidth(bid))
        return self.scaleFactorHist[cat].GetBinContent(bid)
    def getScaleFactorHist(self,cat='def'):
        return self.scaleFactorHist[cat]

    def setDataHist(self,hist,cat='def'):
        self.dataHist[cat]=hist.Clone()
        self.dataHist[cat].SetName(cat+'_dataHist')
    def setMCHist(self,hist,cat='def'):
        self.mcHist[cat]=hist.Clone()
        self.mcHist[cat].SetName(cat+'_mcHist')
    
    def getDataHist(self,cat='def'):
        return self.dataHist[cat]
    def getMCHist(self,cat='def'):
        return self.mcHist[cat]
    def Print(self):
        print("Categories  : ",self.scaleFactorHist.keys())
        
def getPtDependentScaleFactor(name="scaleFactor",
                              dataHist=None,
                              mcHists=None,
                              mcHistsToSubstract=None,
                              binEdges=None,
                              def_scaleFactor=0.0,
                              cat='def'):
    mcHistSum=mcHists[0].Clone()
    mcHistSum.Reset()
    mcHistSubsSum=mcHists[0].Clone()
    mcHistSubsSum.Reset()
    scaleFactorHist=ROOT.TH1F("","",1,0,100.0)
    if binEdges!=None:
#         binEdges.append(binEdges[-1]) # hck for fix, dont know why last bin not taken
        nBinEs=len(binEdges)
        bEdges=np.asarray(binEdges)
        scaleFactorHist=ROOT.TH1F(name,"",nBinEs-1,bEdges)
    else:
        scaleFactorHist=mcHists[0].Clone()
        scaleFactorHist.Reset()
    scaleFactorHist.SetName(name+"_scaleFactor")
    dataHTmp=scaleFactorHist.Clone()
    mcHTmp=scaleFactorHist.Clone()
    for h in mcHists:
            mcHistSum.Add(h)
    if mcHistsToSubstract:
        for h in mcHistsToSubstract:
            print("substracting hist ",h.GetName(), h.Integral())
            h=h.Clone()
            h.Scale(-1.0)
            mcHistSubsSum.Add(h)
    

    if dataHist.GetNbinsX() != mcHistSum.GetNbinsX():
        print("bin counts dont match !!")
        return None
    print("Number of bins in Data : ",dataHist.GetNbinsX())
    print("Number of bins in MC : ",mcHistSum.GetNbinsX())
    print("Width of 1st bin in Data : ",dataHist.GetBinWidth(1))
    print("Width of 1st bin in MC : ",mcHistSum.GetBinWidth(1))
    print("Number MC histograms : ",len(mcHists))
    scaler=mcScaler()        
    for i in range(1,dataHist.GetNbinsX()+1):
        x=dataHist.GetBinCenter(i)
        dValue=dataHist.GetBinContent(i)
        mcValue=mcHistSum.GetBinContent(i)
        dataHTmp.Fill(x,dValue)
        mcHTmp.Fill(x,mcValue)

    for i in range(1,dataHist.GetNbinsX()+1):
        x=dataHist.GetBinCenter(i)
        mcValue=min(mcHistSubsSum.GetBinContent(i),dataHist.GetBinContent(i))
        dataHist.Fill(x,-1*mcValue)

    for i in range(1,dataHTmp.GetNbinsX()+1):
        dValue =dataHTmp.GetBinContent(i)
        mcValue=mcHTmp.GetBinContent(i)
        if mcValue==0:
            mcValue=1e-4
            print("\tWARNING : MC value found to be 0 for [",
                  dataHTmp.GetBinLowEdge(i),",",dataHTmp.GetBinLowEdge(i+1),
                  "]", "data value = ",dValue)
        scl=dValue/mcValue
        scaleFactorHist.SetBinContent(i,scl)
        scaler.binEdges.append(dataHTmp.GetBinLowEdge(i))  
        scaler.scaleFactors.append(scl)
        print(" [ ",scaleFactorHist.GetBinLowEdge(i)," ] ", " -> " ,scaleFactorHist.GetBinContent(i) )
    scaler.binEdges.append(dataHist.GetBinLowEdge(dataHist.GetNbinsX()+1))
    scaler.setSFHist(scaleFactorHist,cat=cat)
    scaler.def_scaleFactor = 0.0
    scaler.setMCHist(mcHTmp,cat=cat)
    scaler.setDataHist(dataHTmp,cat=cat)

    return scaler


class bdtScaler:
    
    def __init__(self,lumi=1.0,balanceDataMCCounts=False):
        self.balanceDataMCCounts=balanceDataMCCounts
        self.doTransformation=False
        self.eval= True
        self.data_x= None
        self.data_w= None
        self.mc_x= None
        self.mc_w= None
        self.reweighter='GBReweighter'
        self.doFolding=False
        self.doTransformation=False
        self.reweightedWeights=None
        self.modelInitalized =False
        self.classifierInitialized =False
        self.normalizationFactor=1.0
        self.lumi=lumi

    def setData(self,X,W):
        self.data_x = X 
        self.data_w = W

    def setMC(self,X,W):
        self.mc_x = X 
        self.mc_w = W
    
    def setNFolding(n=-1):
        if n >0:
            self.doFolding=True
            self.n_folds=n
        else:
            self.doFolding=False
    def setTransformation(self,state=False):
        self.doTransformation = state
    
    def setMCX(self,X):
        self.mc_x = X 
    
    def setModel(self,n_estimators=40,min_samples_leaf=300,max_depth=3,gb_args={'subsample': 0.5}):
        
        self.n_estimators     = n_estimators
        self.min_samples_leaf = min_samples_leaf
        self.max_depth        = max_depth
        self.gb_args          = gb_args
        self.reweighter_base  = GBReweighter(n_estimators     = n_estimators,
                                             min_samples_leaf = min_samples_leaf,
                                             max_depth        = max_depth,
                                             gb_args          = gb_args)
        self.reweighter = self.reweighter_base
        if self.doFolding:
            self.reweighter = FoldingReweighter(self.reweighter_base, n_folds=self.n_folds)
        if self.doTransformation:
            self.transformer = IronTransformer()

        self.modelInitalized=True
    
    def _fitModel(self):
         
        data_x=self.data_x.T ; data_w   = self.data_w
        mc_x  =self.mc_x.T ; mc_w = self.mc_w
        
        if self.balanceDataMCCounts:
            n=min(len(mc_w),len(data_w))
            data_x,data_w = shuffle(data_x,data_w,random_state=42)
            mc_x,mc_w     = shuffle(mc_x  ,mc_w  ,random_state=42)
            data_x=data_x[:n] ; data_w = data_w[:n]
            mc_x  =mc_x[:n]     ; mc_w = mc_w[:n]

        self.reweighter.fit(original=mc_x   , original_weight=mc_w,
                            target  =data_x , target_weight  =data_w)
    
    def fit(self,data_x,data_w,mc_x,mc_w,var_list):
        if not self.modelInitalized:
            raise Exception("Model defenition is empty ")
        
        data_xp=data_x
        mc_xp  =mc_x
        if self.doTransformation:
            self.transformer.fit(data_x.T).to_numpy().T
            data_xp = transformer.transform(data_x.T).to_numpy().T
            mc_xp   = transformer.transform(mc_x.T).to_numpy().T
        self.setData(data_xp,data_w)
        self.setMC(mc_xp,mc_w)
        self.var_list=var_list
        self._fitModel()
        MC_weights  = self.predictWeight(mc_x)
        wMC_total   = np.sum(MC_weights)
        wData_Total = np.sum(self.data_w)
        self.normalizationFactor= wData_Total/wMC_total
        
        self.setTestEvent()
        self.evalTestEvent()

    def setClassifier(self,n_estimators=100, learning_rate=1.0 , max_depth=3):
        self.classifier= GradientBoostingClassifier(  n_estimators=n_estimators ,
                                                      learning_rate=learning_rate ,
                                                      max_depth=max_depth)
        self.classifierInitialized=True
    
        
    def trainClassifierModel(self,dataDict=None):
        if not self.classifierInitialized:
            raise Exception("Classifier not initialized")
        
        data_x=self.data_x ; data_w   = self.data_w
        mc_x  =self.mc_x ; mc_w = self.mc_w
        if dataDict:

            data_x=dataDict['data_x']
            mc_x=dataDict['mc_x']

            data_w=dataDict['data_w']
            mc_w=dataDict['mc_w']
        
        if self.balanceDataMCCounts:
            n=min(len(mc_w),len(data_w))
            data_x,data_w = shuffle(data_x.T,data_w,random_state=42)
            mc_x,mc_w     = shuffle(mc_x.T  ,mc_w  ,random_state=42)
            data_x=data_x[:n].T ; data_w = data_w[:n]
            mc_x  =mc_x[:n].T     ; mc_w = mc_w[:n]

        X_train=np.concatenate( [ data_x , mc_x ] , axis = -1).T
        y1=np.ones(data_x.shape[1])
        y2=np.zeros( mc_x.shape[1])  # MC is set as the positive class 
        Y_train= np.concatenate([y1,y2])
        w1= data_w# /sum(data_w)
        w2= mc_w  # /sum(mc_w)
        w_train=np.concatenate([w1,w2])
        
        self.X_train,self.Y_train,self.w_train = shuffle(X_train,Y_train,w_train, random_state=42)
        self.classifier.fit(self.X_train,self.Y_train, sample_weight = self.w_train )

    def validateClassifierModel( self , data_x, mc_x, mc_w, w_post=None):
        results={}
        if not self.classifierInitialized:
            raise Exception("Classifier is not Initialized")
        
        X=np.concatenate( [ data_x , mc_x ] , axis = -1)
        y1=np.ones(data_x.shape[1] )
        y2=np.zeros( mc_x.shape[1]  )  # MC is set as the positive class 
        Y= np.concatenate([y1,y2])
        data_w=np.ones(data_x.shape[1] )
        W=np.concatenate( [ data_w , mc_w ] )    
        
        ypred=self.classifier.predict_proba( X.T  )[:,1]
        results['data_score']= ypred[ Y > 0.5  ]
        results['mc_score']  = ypred[ Y < 0.5  ]

        roc=roc_curve( y_true=Y , y_score = ypred  )
        results['roc_raw'] ={ 'fpr': roc[0] ,'tpr': roc[1] ,'thresholds' : roc[2] }  
        #results['auc_raw'] =  roc_auc_score( y_true=Y , y_score = ypred )
        srt=np.argsort(results['roc_raw']['fpr'])
        results['auc_raw'] =  auc( results['roc_raw']['fpr'][srt], results['roc_raw']['tpr'][srt] )
        
        roc=roc_curve( y_true=Y , y_score = ypred  , sample_weight = W)
        results['roc_pre'] ={ 'fpr': roc[0] ,'tpr': roc[1] ,'thresholds' : roc[2] }  
        #results['auc_pre'] =roc_auc_score( y_true=Y , y_score = ypred , sample_weight= W)
        srt=np.argsort(results['roc_pre']['fpr'])
        results['auc_pre'] =  auc( results['roc_pre']['fpr'][srt], results['roc_pre']['tpr'][srt] )
        
        if not w_post:
            w_post = self.predictWeight(mc_x)
        Wp=np.concatenate( [ data_w , w_post ] )    
        roc  =  roc_curve( y_true=Y , y_score = ypred  , sample_weight = Wp )
        results['roc_post'] ={ 'fpr': roc[0] ,'tpr': roc[1] ,'thresholds' : roc[2] }  
        #results['auc_post'] =  roc_auc_score( y_true=Y , y_score = ypred , sample_weight= w_post)
        srt=np.argsort(results['roc_post']['fpr'])
        results['auc_post'] =  auc( results['roc_post']['fpr'][srt], results['roc_post']['tpr'][srt] )

        return results
    
    def setTestEvent(self):
        self.test_vector=self.mc_x[:,:4]
        self.test_weight=self.predictWeight(self.test_vector)

    def getTestEvents(self):    
        return self.test_vector , self.test_weight


    def evalTestEvent(self,returnEval=True):
        i=0
        wpred=self.predictWeight(self.test_vector)
        for row,w,wp in zip(self.test_vector.T , np.round(self.test_weight,3),wpred):
            i+=1
            rowp=[ i for i in np.round(row,3) ]
            print(f"{i=} | {rowp}  ,| real : { w:.3f}  eval : {wp:.3f}")
        if returnEval:
            return wpred
        return None
    def saveModel(self,f):
        if not self.modelInitalized:
            raise Exception("Model defenition is empty ")
        result={}    
        result['baseReweighter']  =self.reweighter_base
        result['reweighter']      =self.reweighter
        result['doTransformation']=self.doTransformation
        result['variables']=self.var_list
        result['haprams']={}
        for ky in ['n_estimators','min_samples_leaf','max_depth','gb_args']:
            result['haprams'][ky]=getattr(self,ky)
        
        result['classifier'] = self.classifier 

        data_x = self.data_x ;     self.data_x=None
        mc_x   = self.mc_x   ;     self.mc_x  =None
        data_w = self.data_w ;     self.data_w=None
        mc_w   = self.mc_w   ;     self.mc_w  =None

        result['model']= self
        pickle.dump(result,f)
        
        self.data_x  = data_x  
        self.mc_x    = mc_x    
        self.data_w  = data_w  
        self.mc_w    = mc_w    

    
    def getData(self):
        return self.data_x, self.data_w
      
    def getMC(self):
        return self.mc_x,self.mc_w
      
    def predictWeightWithoutLumi(self, x,lumi=None):
        if not lumi:
            lumi = self.lumi
        xp=x
        if self.doTransformation:
            xp=self.transform.transform(x.T).to_numpy().T
        return self.reweighter.predict_weights( xp.T )*self.normalizationFactor/lumi

    def predictWeight(self, x):
        xp=x
        if self.doTransformation:
            xp=self.transform.transform(x.T).to_numpy().T
        return self.reweighter.predict_weights( xp.T )*self.normalizationFactor   

    def transform(self, x):
        if doTransformation:
            return self.transform.transform(x.T).to_numpy().T
        return x

 
def plotDataMCComparison(var,x_mc,x_data,w_data,w_mc,w_mc_binned,w_mc_bdt,bins=None,saveBase=None):

    ### Plotting the Data MC of different variables after reweighting for comparison study
    varID=0
    binEdges=bins
    # binEdges=np.arange(-4.0,4.0,20)
    dolog=True
    
    f=plt.figure(figsize=(24,8))
    
    if True : # base hist
        hd =np.histogram(x_data,weights=w_data,bins=binEdges)
        hmc=np.histogram(x_mc  ,weights=w_mc  ,bins=binEdges)
    
        hData=utl.getTH1FromNumpHist(hd)
        hMC  =utl.getTH1FromNumpHist(hmc)
        hRatio=hData.Clone()
        hRatio.Divide(hMC)
    
        ax1=plt.subplot(4,3,(1,7))
        hep.histplot(hData,ax=ax1,label='Data')
        hep.histplot(hMC,ax=ax1,label='MC')
 #       ax1.set_ylim( bottom=0.9)
        #     ax1.set_xticks([])
        ax1.grid()
    
        ax2=plt.subplot(4,3,10)
        hep.histplot(hRatio,ax=ax2,label='MC')
        ax2.axhline(1.0,c='k')
    #     ax2.set_xlabel([0.0,2.0])
        ax2.grid(color='red',which='minor')
        ax1.annotate('Raw  Weights' ,(0.60,0.75),xycoords='axes fraction', fontsize=15)
        ax1.annotate(var,(0.60,0.65),xycoords='axes fraction', fontsize=15)
        ax1.legend()
        if dolog:
            ax1.semilogy()
    
    if True : # reweighted hist with binned reweighting
        hd =np.histogram(x_data,weights=w_data     ,bins=binEdges )
        hmc=np.histogram(x_mc  ,weights=w_mc_binned,bins=binEdges )
    
        hData=utl.getTH1FromNumpHist(hd)
        hMC=utl.getTH1FromNumpHist(hmc)
        hRatio=hData.Clone()
        hRatio.Divide(hMC)
    
        ax3=plt.subplot(4,3,(2,8))
#        ax3.set_ylim( bottom=0.9)
        hep.histplot(hData,ax=ax3,label='Data')
        hep.histplot(hMC,ax=ax3,label='MC')
        #     aax3.set_xticks([])
        ax3.grid()
    
        ax4=plt.subplot(4,3,11)
        hep.histplot(hRatio,ax=ax4,label='MC')
        ax4.axhline(1.0,c='k')
        ax4.set_ylim([0.0,4.0])
        ax4.grid(color='red',which='minor')
        ax3.annotate(var,(0.60,0.65),xycoords='axes fraction',fontsize=15)
        ax3.annotate('Binned Reweighting',(0.60,0.75),xycoords='axes fraction',fontsize=15)
        ax3.legend()
        if dolog:
            ax3.semilogy()
    
    if True : # reweighted hist with bdt reweighting
        hd =np.histogram(x_data,weights=w_data,bins=binEdges)
        hmc=np.histogram(x_mc  ,weights=w_mc_bdt,bins=binEdges)
    
        hData=utl.getTH1FromNumpHist(hd)
        hMC=utl.getTH1FromNumpHist(hmc)
        hRatio=hData.Clone()
        hRatio.Divide(hMC)
    
        ax3=plt.subplot(4,3,(3,9))
#        ax3.set_ylim( bottom=0.9)
        hep.histplot(hData,ax=ax3,label='Data')
        hep.histplot(hMC,ax=ax3,label='MC')
        #     aax3.set_xticks([])
        ax3.grid()
    
        ax4=plt.subplot(4,3,12)
        hep.histplot(hRatio,ax=ax4,label='MC')
        ax4.axhline(1.0,c='k')
        ax4.set_ylim([0.0,4.0])
        ax4.grid(color='red',which='minor')
        ax3.annotate(var,(0.60,0.65),xycoords='axes fraction' , fontsize=15)
        ax3.annotate('BDT Reweighting',(0.60,0.75),xycoords='axes fraction',fontsize=15)
        ax3.legend()
        if dolog:
            ax3.semilogy()
    if saveBase:
        foutname=saveBase+'/'+var+'.jpeg'
        #print("Saving file : ",foutname)
        f.savefig(foutname,bbox_inches='tight')
    plt.close(f)

       
