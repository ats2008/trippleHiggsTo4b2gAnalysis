import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from sklearn import metrics
from sklearn.preprocessing import LabelBinarizer
from sklearn import ensemble
from xgboost import XGBClassifier

class classifier :
    def __init__(self, modelParamters):
        
        self.modelParamters = modelParamters
        self.clf=None
    
    def fit(self,X,Y,weight=None):
        raise NotImplementedError
    
    def predict_proba(self, X ):
        raise NotImplementedError
        
    def setModel(self, clf,modelParamters ):
        raise NotImplementedError
    
    def predict_probaForX(self,dataDict):
        x=np.array([ dataDict[ky] for ky in modelParamters  ] ).reshape(1,-1)
        return predict_proba(X)
        

        
class sklearn_model( classifier):
    def __init__(self, modelParamters):
        super().__init__(modelParamters)
    
    def fit(self,X,Y,sample_weight=None):
        return self.clf.fit(X,Y,sample_weight=weight)

    def predict_proba(self, X ):
        return self.clf.predict_proba( X )
        
    def setModel(self, clf,modelParamters ):
        self.clf = clf
        self.modelParamters = modelParamters

class sklearn_gbBDT( sklearn_model ):       
    def __init__(self,modelParamters):
        super().__init__(modelParamters)
        print(self.__dict__)
        print("HERE")
        self.clf = ensemble.GradientBoostingClassifier(**modelParamters)

class xgb_bdt( classifier  ):
    
    def __init__(self, modelParamters):
        super().__init__(modelParamters)
        self.clf= XGBClassifier( **self.modelParamters)
    def fit(self,X,Y,sample_weight=None):
        return self.clf.fit(X,Y,sample_weight=sample_weight)

    def predict_proba(self, X ):
        return self.clf.predict_proba( X )
 
    def setModel(self, clf,modelParamters ):
        self.clf = clf
        self.modelParamters = modelParamters

def getClassifier(tag='GB_BDT',modelParameters_={}):
    
    if tag=='GB_BDT':
        modelParameters =  {
                "n_estimators": 100,
                "max_leaf_nodes": 4,
                "max_depth": 3,
                "random_state": 2,
                "min_samples_split": 5
            }
        #modelParameters.update(modelParameters_)
        clf=sklearn_gbBDT(modelParameters)
        return clf,modelParameters
    if tag=='XGB_BDT':
        modelParameters = {
            "n_estimators": 600,
            "eta": 0.2,
            "max_depth": 4,
            "random_state": 2,
            "subsample" : 0.6,
            "verbose":True,
            'tree_method' : 'gpu_hist'
        }
        clf=xgb_bdt(modelParameters)
        return clf,modelParameters

    if tag=='XGB_BDT_v1':
        modelParameters = {
            "n_estimators": 100,
            "eta": 0.3,
            "max_depth": 2,
            "random_state": 2,
            "subsample" : 0.7,
            'tree_method' : 'gpu_hist'
        }
        clf=xgb_bdt(modelParameters)
        return clf,modelParameters
    
    if tag=='XGB_BDT_v2':
        modelParameters = {
            "n_estimators": 30,
            "eta": 0.3,
            "max_depth": 3,
            "random_state": 2,
            "subsample" : 0.7,
            'tree_method' : 'gpu_hist'
        }
        clf=xgb_bdt(modelParameters)
        return clf,modelParameters
    
    if tag=='XGB_BDT_v3':
        modelParameters = {
            "n_estimators": 1000,
            "eta": 0.2,
            "max_depth": 4,
            "random_state": 2,
            "subsample" : 0.7,
            "verbose":True,
            'tree_method' : 'gpu_hist'
        }
        clf=xgb_bdt(modelParameters)
        return clf,modelParameters
    
    if tag=='XGB_BDT_v4':
        modelParameters = {
            "n_estimators": 400,
            "eta": 0.3,
            "max_depth": 3,
            "random_state": 2,
            "subsample" : 0.7,
            'tree_method' : 'gpu_hist'
        }
        clf=xgb_bdt(modelParameters)
        return clf,modelParameters
     
    if tag=='XGB_BDT_v5':
        modelParameters = {
            "n_estimators": 400,
            "eta": 0.3,
            "max_depth": 3,
            "random_state": 2,
            "subsample" : 0.7,
            'tree_method' : 'hist'
        }
        clf=xgb_bdt(modelParameters)
        return clf,modelParameters
     
    if tag=='XGB_BDT_v5':
        modelParameters = {
            "n_estimators": 400,
            "eta": 0.3,
            "max_depth": 3,
            "random_state": 2,
            "subsample" : 0.7,
            'tree_method' : 'hist'
        }
        clf=xgb_bdt(modelParameters)
        return clf,modelParameters
  
    if tag=='XGB_BDT_v6':
        modelParameters = {
            "n_estimators": 800,
            "eta": 0.3,
            "max_depth": 3,
            "random_state": 2,
            "subsample" : 0.7,
            'tree_method' : 'hist'
        }
        clf=xgb_bdt(modelParameters)
        return clf,modelParameters
 

def plot_mvaScores(y_train_pred,W_train,C_train,y_true=None,fname=None,year=2018):    
    edges=np.linspace(0,1,50)
    edgesCentr=0.5*(edges[:-1]+edges[1:])
    
    sMask   = C_train > 0
    bMask   = C_train==0
    catBkg  = C_train==0
    catMisc = C_train==1
    catARec = C_train==2
    cat1Mis = C_train==3
    cat2Mis = C_train==4
    cat3Mis = C_train==5
    

    nMax=y_train_pred.shape[1]
    if nMax>4:
        f,ax=plt.subplots(2,3,figsize=(35,35))
        ax=np.ndarray.flatten(ax)
    elif nMax > 2:
        f,ax=plt.subplots(2,2,figsize=(25,25))
        ax=np.ndarray.flatten(ax)
    elif nMax > 0:
        f,ax=plt.subplots(1,2,figsize=(20,9))
        ax=np.ndarray.flatten(ax)
    else:
        f,ax=plt.subplots(1,1,figsize=(9,9))
        ax=np.array([ax])
    for kidx in range(nMax):
        
        count,edges=np.histogram(y_train_pred[:,kidx][catBkg],weights=W_train[catBkg],bins=edges)
        ax[kidx].scatter(edgesCentr,count/np.sum(W_train[catBkg]),label='Bkg')

        s=np.sum(W_train[sMask])
        count,edges=np.histogram(y_train_pred[:,kidx][catARec],weights=W_train[catARec],bins=edges)
        ax[kidx].scatter(edgesCentr,count/s,label='All Reco')
        count,edges=np.histogram(y_train_pred[:,kidx][cat1Mis],weights=W_train[cat1Mis],bins=edges)
        ax[kidx].scatter(edgesCentr,count/s,label='1 Jet miss')
        count,edges=np.histogram(y_train_pred[:,kidx][cat2Mis],weights=W_train[cat2Mis],bins=edges)
        ax[kidx].scatter(edgesCentr,count/s,label='2 Jet miss')
        count,edges=np.histogram(y_train_pred[:,kidx][cat3Mis],weights=W_train[cat3Mis],bins=edges)
        ax[kidx].scatter(edgesCentr,count/s,label='3 Jet miss')
        count,edges=np.histogram(y_train_pred[:,kidx][catMisc],weights=W_train[catMisc],bins=edges)

        
        ax[kidx].scatter(edgesCentr,count/s,label='misc')
        ax[kidx].set_xlabel("MVA "+str(kidx))
        ax[kidx].semilogy()
        ax[kidx].legend()
        ax[kidx].grid(which='minor',c='k')
        ax[kidx].set_ylabel('#fraction')
        ax[kidx].set_xlabel('MVA Score')
        hep.cms.label('Work In Progress',ax=ax[kidx],year=year,com=13,fontsize=20)
    
    if fname:
        f.savefig(fname,bbox_inches='tight')
    
    plt.close(f)
    
def plotTPRvsFPR(y_train_pred,y_true,W_train,year=2018,fname=None):
    
    Y_train_reduced=LabelBinarizer().fit_transform(y_true)

    f=plt.figure(figsize=(16,12))
    
    nMax=y_train_pred.shape[1]
    if Y_train_reduced.shape[1] > 1:
        for k in range(1,nMax):
            roc=metrics.roc_auc_score(Y_train_reduced[:,k],y_train_pred[:,k],sample_weight=W_train)
            fpr, tpr, thresholds=metrics.roc_curve(Y_train_reduced[:,k],y_train_pred[:,k],sample_weight=W_train)
            plt.plot(tpr,fpr,label="MVA "+str(k))
            plt.text(0.2,0.5*0.5**k,str(k)+' AUC = '+str(np.round(roc,4)),variant='small-caps',fontweight='bold')
    else:
        roc=metrics.roc_auc_score(Y_train_reduced[:,0],y_train_pred[:,1],sample_weight=W_train)
        fpr, tpr, thresholds=metrics.roc_curve(Y_train_reduced[:,0],y_train_pred[:,1],sample_weight=W_train)
        plt.plot(tpr,fpr,label="MVA ")
        plt.text(0.2,0.5,' AUC = '+str(np.round(roc,4)),variant='small-caps',fontweight='bold')

    plt.semilogy()
    plt.ylim([1e-5,1.0])
    plt.ylabel('FPR')
    plt.xlabel('TPR')
    plt.grid(which='minor',c='k')
    hep.cms.label('Work In Progress',year=year,com=13,fontsize=20)
    plt.legend()
    
    if fname:
        f.savefig(fname,bbox_inches='tight')
    
    plt.close(f)
