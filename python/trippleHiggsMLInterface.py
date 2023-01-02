import numpy as np
import pickle as pkl

def sigmoid(x):
    return 1/(1.0+np.exp(-x))

class mlScoreManger():
    
    def __init__(self, fname = None ):
        self.input_fname = fname
        self.raw_data =None
        self.nEvents=0
        if self.input_fname:
            self.load_file()
    def GetEntries(self):
        return  self.nEvents
    def load_file(self,fname=None):
        if fname:
            self.input_fname=fname
        if not self.input_fname:
            print("Input Filename not set !")
            return
        
        with open(self.input_fname,'rb') as f:
            self.raw_data=pkl.load(f)
            self.nEvents=self.raw_data['eventIndex'].shape[0]
            if self.raw_data['eventIndex'].shape[0] != np.unique(self.raw_data['eventIndex']).shape[0]:
                print("Duplicates in Event Indices !")
            
    def getEntry(self,eventIdx):
        idx=np.where(self.raw_data['eventIndex']==eventIdx)
        
        if idx[0].shape[0] != 1:
            return {'valid':False,'matches':idx[0].shape[0]}
        rslt={'valid':True,'matches':idx[0].shape[0]}
        index=idx[0][0]
        rslt['eventIndex']=self.raw_data['eventIndex'][index]
        rslt['y0']=self.raw_data['y0'][index]
        rslt['y1']=self.raw_data['y1'][index]
        rslt['label']=self.raw_data['label'][index]
        rslt['vldMask']=self.raw_data['vldMMask'][index]
        rslt['score']=self.raw_data['score'][index][:,0]
        
        return rslt    
    def printAllEventsIn(self):
        for i in self.raw_data['eventIndex']:
            print(i)

