import matplotlib.pyplot as plt
import mplhep as hep
import uproot as urt
import ROOT
import json,sys,os,argparse
import numpy as np
import Util as utl
import pickle

from datetime import datetime

import prettytable as ptab

import scaleFactorUtil as scl
import statUtil as stat
import hep_ml as hepml
import prettytable
import classifierUtil as clsUtil

import sklearn 
from sklearn import metrics
from sklearn.metrics import log_loss
from sklearn.model_selection import train_test_split


hep.style.use("CMS")
blind='CMS_hgg_mass < 115 || CMS_hgg_mass > 135.0'

def main():
 
    parser = argparse.ArgumentParser()
    parser.add_argument('-v',"--version", help="Version of the specific derivation ",default='')
    parser.add_argument('-i',"--inputFile", help="Input File",default=None)
    parser.add_argument('-y',"--year", help="Year",default='2018')
    parser.add_argument('-f',"--fsrc", help="Input file Source",default=None)
    parser.add_argument("--test", help="do testing of model ",action='store_true')
    parser.add_argument("--exportWeights", help="export events weights ",action='store_true')
    parser.add_argument("--doSlim", help="doslim",default=False,action='store_true')
    parser.add_argument('-o',"--dest", default="workarea/tmp/classifier/",help='output data prefix')
    parser.add_argument('-c',"--clsTag", default="GB_BDT",help=' Classifier tag , selects the classifier to be trained')

    args = parser.parse_args()
    
    isTraining= True

    version = args.version
    fsrc= args.fsrc
    doIronTrnsformation=False
    inputFile  = args.inputFile
    saveOutput = True
    isTest = args.test
    exportWeights = args.exportWeights
    classifierTag = args.clsTag
    print("Clasifier Tag ",classifierTag)
    
    ## Making the variables for th reweighter training
    varForModel=[]
    with open('data/modelInfo.json') as f:
        modelDict=json.load(f)
    varForModel=modelDict['variables']    
    
    varForModel_=np.unique(varForModel)
    varForModel=[ i for i in varForModel_]
 

    outDict={}
    inputFileSource={}
    if not inputFileSource:
        if not isTest:
            setsToProcess=['test','train']
            if args.doSlim:
                inputFileSource['train'] = 'workarea/data/bdtNtuples/v8p3/slim_train/filelist.json'
                inputFileSource['test']  = 'workarea/data/bdtNtuples/v8p3/slim_train/filelist.json'
            else:
                inputFileSource['train'] = 'workarea/data/bdtNtuples/v8p3/train/filelist.json'
                inputFileSource['test']  = 'workarea/data/bdtNtuples/v8p3/test/filelist.json'
        else:
            setsToProcess=['test']
            inputFileSource['test']='workarea/data/bdtNtuples/v8p3/test/filelist.json'
    else:
        inputFileSource['train'] = fsrc.split(',')[0]
        inputFileSource['test']  = fsrc.split(',')[1]

    if isTest:
        print(" Testing the model !")
        saveOutput=False
        if not inputFile:
            print("Please provide an input sample file")
            exit(1)
        else:
            setsToProcess=['test']
    doTraining=True
    if inputFile:
        doTraining=False
        with open(inputFile, 'rb') as f:
            modelDict = pickle.load(f)
            
    fileDictGlobal={}
    for ky in inputFileSource:
        ipt=inputFileSource[ky]
        with open(ipt) as f:
            fileDictGlobal[ky]=json.load(f)
    
    prefixBase=args.dest

    ## YEARS TO PROCESS
    yearsToProcess_all=['2018','2017','2016PreVFP','2016PostVFP','run2']
    yearsToProcess=[]
    if '*' in args.year:
        yearsToProcess=yearsToProcess_all
    else:
        for yr in args.year.split(","):
            if yr not in yearsToProcess_all:
                print(yr," not in catalogue. Skipping !! ")
                continue
            yearsToProcess.append(yr)

    ## BKGS TO PROCESS
    bkgToProcess=[ 'ggBox1Bjet',
                   'ggBox2Bjet', 
                   'ggBox', 
                   'gJet20To40',
                   'gJet40ToInf'
                 ]
    
    
    outDict['varForModel']=varForModel
    
    print("Variables being used for training : ",end="\n\t")
    for i in range(len(varForModel)):
        print(varForModel[i],end = " , ")
        if i%4==0:
            print(end="\n\t")
    print()
    
    print("Outputs are stored in ",prefixBase)
    

    # We read the tree from the file and create a RDataFrame, a class that
    # allows us to interact with the data contained in the tree.
    rdataFramesGlobal={}
    for setX in setsToProcess:
        rdataFrames={}
        fileDict=fileDictGlobal[setX]
        for yr  in yearsToProcess:
            rdataFrames[yr]={'sig':{},'bkg':{},'data':{}}
            try:
                fileName=fileDict[yr]['sig']['ggHHH']
            except:
                print(yr)
            treeName = "trees/ggHHH_125_13TeV"
            rdataFrames[yr]['sig']['ggHHH']= ROOT.RDataFrame(treeName, fileName)
            print("Registering datset : ggHHH , ",yr," with tree",treeName)
            
            ky=list(fileDict[yr]['data'].keys())[0]
            fileName=fileDict[yr]['data'][ky]
            treeName = "trees/Data_13TeV_TrippleHTag_0"
            rdataFrames[yr]['data']['data']= ROOT.RDataFrame(treeName, fileName).Filter(blind)
            print("Registering datset : data , ",yr," with tree",treeName)
            
            rdataFrames[yr]['bkg']={}
            treeName = "trees/bkg_13TeV_TrippleHTag_0"
            for bkg in bkgToProcess:
                if bkg not in fileDict[yr]['bkg']:
                    print()
                    print("FILE NOT FOUND FOR : ",yr," background ",bkg)
                    print()
                    continue
                fileName = fileDict[yr]['bkg'][bkg]
                rdataFrames[yr]['bkg'][bkg]=ROOT.RDataFrame(treeName, fileName).Filter(blind)
                print("Registering datset : ",bkg," , ",yr," with tree",treeName)
        rdataFramesGlobal[setX]=rdataFrames
    
    # reading the vars from file    
    oVars=['genRecoCategory','weight','weight_bdt','lumi']
    varsToGet=[ i for i in np.unique(varForModel + oVars )]
    
    for yr in yearsToProcess:
        dataStore={}
        dataSets={}
        results={}
        print("Processing for year : ",yr)
        saveBase=prefixBase+'/'+yr+'/'
        if isTest:
            saveBase=prefixBase+'/test/'+yr+'/'

        os.system('mkdir -p '+saveBase)
        outDict['eventsStats']={}
        ridx=0
        outDictROC={}
        for setX in setsToProcess:
            print("  Processing ",setX," set")
            dataStore[yr]={}
            rdataFrames=rdataFramesGlobal[setX]
            models={}
            print("\t Converting the samples to Numpy Arrays")
            dataStore[yr]['sig']={}
            for ky in rdataFrames[yr]['sig']:
                print("\t\t sig : ",ky)
                dataStore[yr]['sig'][ky]=rdataFrames[yr]['sig'][ky].AsNumpy(varsToGet)
            dataStore[yr]['bkg']={}
            for ky in rdataFrames[yr]['bkg']:
                print("\t\t bkg : ",ky)
                dataStore[yr]['bkg'][ky]=rdataFrames[yr]['bkg'][ky].AsNumpy(varsToGet)
            dataStore[yr]['data']={}
            for ky in rdataFrames[yr]['data']:
                print("\t\t data : ",ky)
                dataStore[yr]['data'][ky]=rdataFrames[yr]['data'][ky].AsNumpy(varsToGet)
            
            data={
                    'data' :{'x' :None,'weight':None},
                    'bkg'  :{'x':None,'weight':None},
                    'sig'  :{'x':None,'weight':None}
                 }
            
            isFirst=True
            print("\t Concatinating the samples ")
            for ky in dataStore[yr]['bkg']: 
                xMC=np.stack([dataStore[yr]['bkg'][ky][k] for k in varForModel])
                print("\t\t bkg : ",ky)
                if isFirst:
                    data['bkg']['x']=xMC
                    data['bkg']['weight']=dataStore[yr]['bkg'][ky]['weight_bdt'] #*dataStore[yr]['bkg'][ky]['lumi']
                    data['bkg']['cat']=dataStore[yr]['bkg'][ky]['genRecoCategory']
                    isFirst=False
                else:
                    data['bkg']['x']=np.concatenate( [data['bkg']['x'],xMC ] , axis=-1)
                    w=dataStore[yr]['bkg'][ky]['weight_bdt']
                    data['bkg']['weight']=np.concatenate( [data['bkg']['weight'],w ] ,axis=-1)
                    data['bkg']['cat']=np.concatenate( [data['bkg']['cat'],dataStore[yr]['bkg'][ky]['genRecoCategory'] ] ,axis=-1)
            
            isFirst=True
            for ky in dataStore[yr]['sig']: 
                print("\t\t sig : ",ky)
                xMC=np.stack([dataStore[yr]['sig'][ky][k] for k in varForModel])
                if isFirst:
                    data['sig']['x']=xMC
                    data['sig']['weight']=dataStore[yr]['sig'][ky]['weight_bdt'] #*dataStore[yr]['sig'][ky]['lumi']
                    data['sig']['cat']=dataStore[yr]['sig'][ky]['genRecoCategory'] #*dataStore[yr]['sig'][ky]['lumi']
                    isFirst=False
                else:
                    data['sig']['x']=np.concatenate( [data['sig']['x'],xMC ] , axis=-1)
                    w=dataStore[yr]['sig'][ky]['weight_bdt']
                    data['sig']['weight']=np.concatenate( [data['sig']['weight'],w ] ,axis=-1)
                    data['sig']['cat']=np.concatenate( [data['sig']['cat'],dataStore[yr]['sig'][ky]['genRecoCategory'] ] ,axis=-1)
            
            ##        
            isFirst=True
            for ky in dataStore[yr]['data']:    
                print("\t\t data : ",ky)
                xData=np.stack([dataStore[yr]['data'][ky][k] for k in varForModel])
                if isFirst:
                    data['data']['x']=xData
                    data['data']['weight']=dataStore[yr]['data'][ky]['weight']
                    data['data']['cat']=dataStore[yr]['data'][ky]['genRecoCategory']
                    isFirst=False
                else:
                    data['data']['x']=np.concatenate( [data['data']['x'],xData ] , axis=-1)
                    data['data']['weight']=np.concatenate( [data['data']['weight'],dataStore[yr]['data'][ky]['weight']] ,
                                                          axis=-1)
                    data['data']['cat']=np.concatenate( [data['data']['cat'],dataStore[yr]['data'][ky]['genRecoCategory'] ] ,axis=-1)

            f,ax=plt.subplots(1,3,figsize=(15,5))
            for i,ky in enumerate(data):
                integ=np.round(np.sum(data[ky]['weight']))
                ax[i].hist(data[ky]['weight'],bins=40,label=ky+' [ '+str(integ)+' ]')
                ax[i].legend(loc=5,fontsize=10)
            f.savefig(saveBase+"/weights_"+setX+".png",bbox_inches='tight')
            plt.close(f)
            
            xVals=np.arange(0.5,6.5,1)
            c,b=np.histogram(data['sig']['cat'],bins=xVals)
            c=c/sum(c)
            b=0.5*(b[:-1]+b[1:])
            plt.bar(b,c,color='b',fill=False)
            plt.xticks(b,['Misc','All Recoed','1 Jet miss','2 Jet Miss','3 Jet Miss'])
            hep.cms.label()
            plt.ylim([0.0,0.6])
            plt.ylabel('fraction of events')
            plt.savefig(saveBase+'/sigEvtLabels_'+setX+'.png')
            plt.close()
            
            X=np.concatenate([data['sig']['x'],data['bkg']['x']],axis=-1).T
            W=np.concatenate([data['sig']['weight']/np.sum(data['sig']['weight']),
                              data['bkg']['weight']/np.sum(data['bkg']['weight'])],
                              axis=-1)

            Y=np.concatenate([np.ones( data['sig']['x'].shape[1]),
                              np.zeros(data['bkg']['x'].shape[1])
                             ],
                             axis=-1)
            C=np.concatenate([data['sig']['cat'],data['bkg']['cat']],axis=-1)
                             
            X_shuffled, Y_shuffled, W_shuffled , C_shuffled = sklearn.utils.shuffle(X, Y, W , C , random_state=42+ridx) ; ridx+=1
            catMisc = data['sig']['cat']==1
            catARec = data['sig']['cat']==2
            cat1Mis = data['sig']['cat']==3
            cat2Mis = data['sig']['cat']==4
            cat3Mis = data['sig']['cat']==5    
            print("\t\t Number of ",setX," events   : ",X_shuffled.shape[0])
            print(            "\t\t Signal events   : ",sum(Y_shuffled==1),f"[{np.sum(W_shuffled[Y_shuffled==1])}]")
            for idx in range(1,6):
                print(            "\t\t    Cat   ",idx,"    : ",sum(C_shuffled==idx),f"[{np.sum(W_shuffled[C_shuffled==idx])}]")
            print(            "\t\t Bkg.   events   : ",sum(Y_shuffled==0),f"[{np.sum(W_shuffled[Y_shuffled==0])}]")
            print() 
            dataSets[setX]={}
            dataSets[setX]['data']=data
            dataSets[setX]['x']=X_shuffled
            dataSets[setX]['w']=W_shuffled
            dataSets[setX]['y']=Y_shuffled
            dataSets[setX]['cat']=C_shuffled
            
            outDict['eventsStats'][setX]={'all' : {'w' : str(np.sum(W_shuffled) ) ,'n' : str(X_shuffled.shape[0]) } }
            outDict['eventsStats'][setX]['sig'] = {'w' : str(np.sum(W_shuffled[Y_shuffled==1])) ,'n' :  str(sum(Y_shuffled==1))} 
            outDict['eventsStats'][setX]['bkg'] = {'w' : str(np.sum(W_shuffled[Y_shuffled==0])) ,'n' :  str(sum(Y_shuffled==0))} 
            for idx in range(5):
                outDict['eventsStats'][setX]['sigCat_'+str(idx)] = {'w' : str(np.sum(W_shuffled[C_shuffled==idx]) ) ,'n' :  str(sum(C_shuffled==idx))} 

        modelParameters=modelDict['modelParameters']
        outDict['modelParameters']=modelParameters
        clf=None
        if isTraining:
            if  inputFile:
                clf=modelDict['model']
            else:
                clf,modelParameters = clsUtil.getClassifier(classifierTag,modelParameters)
            outDict['modelParameters']=modelParameters
            print("-"*80)
            print("\t\t\t Training the model !! ")
            print("-"*80)
            
            if doTraining:
                clf.fit(dataSets['train']['x'],dataSets['train']['y'])
                feature_importances=clf.clf.feature_importances_ 
                tabFeatureImportance=ptab.PrettyTable(['Variable','Importance Score' ] )
                outDict['feature_importances']={}
                for n in np.argsort(-1*feature_importances):
                    var=varForModel[n]
                    tabFeatureImportance.add_row(  [ var , feature_importances[n] ] )
                    outDict['feature_importances'][var]=str(np.round(feature_importances[n],4))
                print(tabFeatureImportance)
            if 'train' in dataSets:
                    
                y_pred=clf.predict_proba(dataSets['train']['x'])
                print("\t\t\t Done ! !! ")
                print("-"*80)
                
                results['model']=clf
                results['modelParameters']=modelParameters
                results['variables']=varForModel
 
                results['name']            = classifierTag
                results['inceptionTime']   = datetime.now().strftime("%m/%d/%Y, %H:%M:%S")

                rocAUC=sklearn.metrics.roc_auc_score(dataSets['train']['y'],y_pred[:,1],sample_weight=dataSets['train']['w'])
                results['rocCurve'] = metrics.roc_curve(dataSets['train']['y'],y_pred[:,1])
                outDict['train_auc']= str(np.round(rocAUC,5))
                outDictROC['train_fpr']=[ float(i) for i in results['rocCurve'][0] ]
                outDictROC['train_tpr']=[ float(i) for i in results['rocCurve'][1] ]
                outDictROC['train_thr']=[ float(i) for i in results['rocCurve'][2] ]
                
                plt.scatter(results['rocCurve'][1],results['rocCurve'][2],s=1,c='blue')
                plt.semilogy()
                plt.grid(which='minor',c='k')
                plt.ylabel('Background Rejection Efficiency')
                plt.xlabel('Signal Efficiency')
                plt.text(0.2,0.3,'AUC = '+str(np.round(rocAUC,4)),variant='small-caps',fontweight='bold')
                hep.cms.label('Work In Progress')
                plt.savefig(saveBase+'/fprVsTpr_training.png',bbox_inches='tight')
                plt.close()

                clsUtil.plot_mvaScores(y_pred,dataSets['train']['w'],dataSets['train']['cat'])
                clsUtil.plotTPRvsFPR(y_pred,dataSets['train']['cat'],dataSets['train']['cat'],fname=None)


        else:
            pass

        foutname=saveBase+'/'+'modelValidation.pkl'
        print("Validation outputs saved in  : ",saveBase)
        with open(foutname, 'wb') as f:
            pickle.dump(results,f)
        foutname=saveBase+'/'+'modelValidation.json'
        with open(foutname, 'w') as f:
            json.dump(outDict,f,indent=4)
        foutname=saveBase+'/'+'modelValidation_roc.json'
        with open(foutname, 'w') as f:
            json.dump(outDictROC,f,indent=4)
if __name__=='__main__':
    main( )

