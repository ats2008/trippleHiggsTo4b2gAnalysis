import os
listOfTags=None
tag="tthNtuplizer"
executable='python/ntupleMaker_tth.py'

tag="fitterSkimsV7_"
executable='python/ntupleMaker_analysisSkim.py'

tag="fitterSkimsV8_noPeakingMVA"
executable='python/ntupleMaker_analysisSkimNoPeaking.py'

#tag="jetIDStudy"
#executable='python/ntupleMaker_h24bseletion.py'

tag="cutBased_variableSearch"
executable='python/flashGGAnalyzer.py'
listOfTags=['data2018A','data2018B','data2018C','data2018D','signal']

tag="cutBased_ntupleSkim"
executable='python/ntupleMaker_fggCutBased.py'

listOfTags=['data2018A','data2018B','data2018C','data2018D','signal']
qcdSamples=['ggM80Inc','ggM80Jbox1bjet','ggM80Jbox2bjet','gJet20To40','gJet40ToInf']
doubleHTags=['gluGluToHH','vbfHH','WToQQHHTo2B2G','ZToBBHHTo2B2G','ggZTo2BHHTo2B2G','ttHHTo2B2G']
singleHTags=['gluGluToH','vbfH','THQ','ttHJetToGG']

for head in [qcdSamples,singleHTags,doubleHTags]:
    for ky in head:
        listOfTags.append(ky)

doFID=True
allDataFiles={
    'data2018A'   :   "workarea/data/flashGGNtuples/RunIISummer19/2018/output_EGamma_alesauva-UL2018_0-10_6_4-v0-Run2018A-12Nov2019_UL2018-v2-981b04a73c9458401b9ffd78fdd24189_USER.root",
    'data2018B'   :   "workarea/data/flashGGNtuples/RunIISummer19/2018/output_EGamma_alesauva-UL2018_0-10_6_4-v0-Run2018B-12Nov2019_UL2018-v2-981b04a73c9458401b9ffd78fdd24189_USER.root",
    'data2018C'   :   "workarea/data/flashGGNtuples/RunIISummer19/2018/output_EGamma_alesauva-UL2018_0-10_6_4-v0-Run2018C-12Nov2019_UL2018-v2-981b04a73c9458401b9ffd78fdd24189_USER.root",
    'data2018D'   :   "workarea/data/flashGGNtuples/RunIISummer19/2018/output_EGamma_alesauva-UL2018_0-10_6_4-v0-Run2018D-12Nov2019_UL2018-v4-981b04a73c9458401b9ffd78fdd24189_USER.root",
    'signal'        :   "workarea/data/flashGGNtuples/RunIISummer19/2018/c3_1_c4_1_HHHto4b2gamma_UL17.root",
    'ggM80Inc'      :   "workarea/data/flashGGNtuples/RunIISummer19/2018/DiPhotonJetsBox_MGG-80toInf_13TeV-sherpa_v1.root",
    'ggM80Jbox1bjet':   "workarea/data/flashGGNtuples/RunIISummer19/2018/output_DiPhotonJetsBox1BJet_MGG-80toInf_13TeV-sherpa.root",
    'ggM80Jbox2bjet':   "workarea/data/flashGGNtuples/RunIISummer19/2018/output_DiPhotonJetsBox2BJets_MGG-80toInf_13TeV-sherpa.root",
    'ttgJets'       :   "workarea/data/flashGGNtuples/RunIISummer19/2018/output_TGJets_TuneCP5_13TeV-amcatnlo-madspin-pythia8.root",
    'ttgg'          :   "workarea/data/flashGGNtuples/RunIISummer19/2018/output_TTGG_TuneCP5_13TeV-amcatnlo-pythia8.root",
    'tGJets'        :   "workarea/data/flashGGNtuples/RunIISummer19/2018/output_TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8.root",
    'gJet20To40'    :   "workarea/data/flashGGNtuples/RunIISummer19/2018/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8.root",
    'gJet40ToInf'   :   "workarea/data/flashGGNtuples/RunIISummer19/2018/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8.root",
    "ttHJetToGG"    :   "workarea/data/flashGGNtuples/RunIISummer19/2018/ttHJetToGG_M125_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8.root",
    "ggZTo2BHHTo2B2G"    :   "workarea/data/flashGGNtuples/RunIISummer19/2018/ggZTo2BHHTo2B2G_UL17.root",
    "ttHHTo2B2G"         :   "workarea/data/flashGGNtuples/RunIISummer19/2018/ttHHTo2B2G_UL17.root",
    "ttWToQQHTo2G"       :   "workarea/data/flashGGNtuples/RunIISummer19/2018/ttWToQQHTo2G_UL17.root",
    "WToQQHHTo2B2G"      :   "workarea/data/flashGGNtuples/RunIISummer19/2018/WToQQHHTo2B2G_UL17.root",
    "ZToBBHHTo2B2G"      :   "workarea/data/flashGGNtuples/RunIISummer19/2018/ZToBBHHTo2B2G_UL17.root",
    "gluGluToHH"         :   "workarea/data/flashGGNtuples/RunIISummer19/2018/GluGluToHHTo2B2G_node_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8.root",
    "vbfHH"              :   "workarea/data/flashGGNtuples/RunIISummer19/2018/VBFHHTo2B2G_CV_1_C2V_1_C3_1_dipoleRecoilOn-TuneCP5_PSweights_13TeV-madgraph-pythia8.root",
    "vhTogg"             :   "workarea/data/flashGGNtuples/RunIISummer19/2018/VHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8.root",
    "gluGluToH"          :   "workarea/data/flashGGNtuples/RunIISummer19/2018/GluGluHToGG_M-125_TuneCP5_13TeV-powheg-pythia8.root",
    "vbfH"               :   "workarea/data/flashGGNtuples/RunIISummer19/2018/VBFHToGG_M125_TuneCP5_13TeV-amcatnlo-pythia8_storeWeights.root",
    "THQ"                :   "workarea/data/flashGGNtuples/RunIISummer19/2018/THQ_ctcvcp_HToGG_M125_TuneCP5_13TeV-madgraph-pythia8.root"  
}
if listOfTags==None:
    listOfTags=list(allDataFiles.keys())

params={
    'data2018A'      : { 'processID':'DATA' ,'TREE':'tagsDumper/trees/Data_13TeV_TrippleHTag_0', 'doPtReWeighting':'0', 'weightScale': 1.0 ,'fID':{'proc' : 'Data', 'mass' : None , 'sqrts' :'13TeV' }  } ,
    'data2018B'      : { 'processID':'DATA' ,'TREE':'tagsDumper/trees/Data_13TeV_TrippleHTag_0', 'doPtReWeighting':'0', 'weightScale': 1.0 ,'fID':{'proc' : 'Data', 'mass' : None , 'sqrts' :'13TeV' }  } ,
    'data2018C'      : { 'processID':'DATA' ,'TREE':'tagsDumper/trees/Data_13TeV_TrippleHTag_0', 'doPtReWeighting':'0', 'weightScale': 1.0 ,'fID':{'proc' : 'Data', 'mass' : None , 'sqrts' :'13TeV' }  } ,
    'data2018D'      : { 'processID':'DATA' ,'TREE':'tagsDumper/trees/Data_13TeV_TrippleHTag_0', 'doPtReWeighting':'0', 'weightScale': 1.0 ,'fID':{'proc' : 'Data', 'mass' : None , 'sqrts' :'13TeV' }  } ,
    'ggM80Inc'       : { 'processID':'MC'   ,'TREE':'tagsDumper/trees/bkg_13TeV_TrippleHTag_0' , 'doPtReWeighting':'1', 'weightScale': 1.0 } ,
    'ggM80Jbox1bjet' : { 'processID':'MC'   ,'TREE':'tagsDumper/trees/bkg_13TeV_TrippleHTag_0' , 'doPtReWeighting':'1', 'weightScale': 1.0 } ,
    'ggM80Jbox2bjet' : { 'processID':'MC'   ,'TREE':'tagsDumper/trees/bkg_13TeV_TrippleHTag_0' , 'doPtReWeighting':'1', 'weightScale': 1.0 } ,
    'gJet20To40'     : { 'processID':'MC'   ,'TREE':'tagsDumper/trees/bkg_13TeV_TrippleHTag_0' , 'doPtReWeighting':'0', 'weightScale': 1.0 } ,
    'gJet40ToInf'    : { 'processID':'MC'   ,'TREE':'tagsDumper/trees/bkg_13TeV_TrippleHTag_0' , 'doPtReWeighting':'0', 'weightScale': 1.0} ,
    
    'ttgJets'        : { 'processID':'MC'   ,'TREE':'tagsDumper/trees/bkg_13TeV_TrippleHTag_0' , 'doPtReWeighting':'0', 'weightScale': 1.0 } ,
    'ttgg'           : { 'processID':'MC'   ,'TREE':'tagsDumper/trees/bkg_13TeV_TrippleHTag_0' , 'doPtReWeighting':'0', 'weightScale': 1.0 } ,
    'tGJets'         : { 'processID':'MC'   ,'TREE':'tagsDumper/trees/bkg_13TeV_TrippleHTag_0' , 'doPtReWeighting':'0', 'weightScale': 1.0 } ,
    
    
    "gluGluToH"      : { 'processID':'MC'   ,'TREE':'tagsDumper/trees/bkg_13TeV_TrippleHTag_0' , 'doPtReWeighting':'0', 'weightScale': 1.0,'fID':{'proc' : 'ggH'  , 'mass' : '125', 'sqrts' :'13TeV' } } ,
    "vhTogg"         : { 'processID':'MC'   ,'TREE':'tagsDumper/trees/bkg_13TeV_TrippleHTag_0' , 'doPtReWeighting':'0', 'weightScale': 1.0,'fID':{'proc' : 'vH'   , 'mass' : '125', 'sqrts' :'13TeV' } } ,
    "vbfH"           : { 'processID':'MC'   ,'TREE':'tagsDumper/trees/bkg_13TeV_TrippleHTag_0' , 'doPtReWeighting':'0', 'weightScale': 1.0,'fID':{'proc' : 'vbfH' , 'mass' : '125', 'sqrts' :'13TeV' } } ,
    "THQ"            : { 'processID':'MC'   ,'TREE':'tagsDumper/trees/bkg_13TeV_TrippleHTag_0' , 'doPtReWeighting':'0', 'weightScale': 1.0,'fID':{'proc' : 'tHq'  , 'mass' : '125', 'sqrts' :'13TeV' } } ,
    "ttHJetToGG"     : { 'processID':'MC'   ,'TREE':'tagsDumper/trees/bkg_13TeV_TrippleHTag_0' , 'doPtReWeighting':'0', 'weightScale': 1.0,'fID':{'proc' : 'ttH'  , 'mass' : '125', 'sqrts' :'13TeV' } } ,
    
    "gluGluToHH"     : { 'processID':'MC'   ,'TREE':'tagsDumper/trees/bkg_13TeV_TrippleHTag_0' , 'doPtReWeighting':'0', 'resetWeight' : 1.0 , 'weightScale': 3.28e-7,'fID':{'proc' : 'vHH'  , 'mass' : '125', 'sqrts' :'13TeV' } } ,
    "vbfHH"          : { 'processID':'MC'   ,'TREE':'tagsDumper/trees/bkg_13TeV_TrippleHTag_0' , 'doPtReWeighting':'0', 'resetWeight' : 1.0 , 'weightScale': 1.83e-8,'fID':{'proc' : 'vbfHH', 'mass' : '125', 'sqrts' :'13TeV' } } ,
    "WToQQHHTo2B2G"  : { 'processID':'MC'   ,'TREE':'tagsDumper/trees/MC_13TeV_TrippleHTag_0'  , 'doPtReWeighting':'0', 'resetWeight' : 1.0 , 'weightScale': 3.58e-9,'fID':{'proc' : 'vHH'  , 'mass' : '125', 'sqrts' :'13TeV' } } ,
    "ZToBBHHTo2B2G"  : { 'processID':'MC'   ,'TREE':'tagsDumper/trees/MC_13TeV_TrippleHTag_0'  , 'doPtReWeighting':'0', 'resetWeight' : 1.0 , 'weightScale': 4.34e-10,'fID':{'proc' : 'vHH'  , 'mass' : '125', 'sqrts' :'13TeV' } } ,
    "ggZTo2BHHTo2B2G": { 'processID':'MC'   ,'TREE':'tagsDumper/trees/MC_13TeV_TrippleHTag_0'  , 'doPtReWeighting':'0', 'resetWeight' : 1.0 , 'weightScale': 6.77e-11,'fID':{'proc' : 'zHH'  , 'mass' : '125', 'sqrts' :'13TeV' } } ,
    "ttHHTo2B2G"     : { 'processID':'MC'   ,'TREE':'tagsDumper/trees/MC_13TeV_TrippleHTag_0'  , 'doPtReWeighting':'0', 'resetWeight' : 1.0 , 'weightScale': 8.20e-9 ,'fID':{'proc' : 'ttHH' , 'mass' : '125', 'sqrts' :'13TeV' } } ,
    "ttWToQQHTo2G"   : { 'processID':'MC'   ,'TREE':'tagsDumper/trees/MC_13TeV_TrippleHTag_0'  , 'doPtReWeighting':'0', 'resetWeight' : 1.0 , 'weightScale': 8.11e-9 ,'fID':{'proc' : 'ttWH' , 'mass' : '125', 'sqrts' :'13TeV' } } ,
    'signal'         : { 'processID':'MC'   ,'TREE':'tagsDumper/trees/MC_13TeV_TrippleHTag_0'  , 'doPtReWeighting':'0', 'resetWeight' : 1.0 , 'weightScale': 3.09e-10,'fID':{'proc' : 'ggHHH', 'mass' : '125', 'sqrts' :'13TeV' } } ,
}

for ky in ['gluGluToH','vhTogg','vbfH','THQ','ttHJetToGG'] :
    params[ky]['fID']['proc']='singleH'

for ky in ['gluGluToHH','WToQQHHTo2B2G','ZToBBHHTo2B2G','ggZTo2BHHTo2B2G','vbfHH','ttWToQQHTo2G','ttHHTo2B2G'] :
    params[ky]['fID']['proc']='doubleH'

base='analysis_cfgs/'+tag
obase='workarea/batch/'+tag



if not os.path.exists(base):
    os.system('mkdir -p '+base)
if not os.path.exists(obase):
    os.system('mkdir -p '+obase)

template="""
#HEADER_BEG
analysis_cfgs/mvaHeader.cfg
#HEADER_END

#PARAMS_BEG
OutpuFileName=@@OBASE/@@OFNAME
MaxEvents=-1
WeightScale=@@WeightScale
treeName=@@TREENAME
outTreeName=@@OTREE
processID=@@PROCESSID
resetWeight=@@resetWeight
doBjetConting=False
maxNBjetsFromMC=0
minNBjetsFromMC=0
doPtReWeighting=@@DoPtReWeighting
pTReweitingFile=metaData/pTScaleFactorFile.root
pTReweitingHistName=pTDependentpT_BB_scaleFactor,pTDependentpT_BE_scaleFactor,pTDependentpT_EB_scaleFactor,pTDependentpT_EE_scaleFactor
pTReweitingHistCatagories=BB,BE,EB,EE
#PARAMS_END

#FNAMES_BEG
@@IFNAME
#FNAMES_END

"""

os.system('mkdir -p '+base)
subName=base+'/subAll.sh'
fSub=open(subName,'w')
#for ky in allDataFiles:
for ky in listOfTags:
    tmp=template.replace('@@IFNAME',allDataFiles[ky])
    tmp=tmp.replace('@@TAG',tag)
    tmp=tmp.replace('@@PROCESSID',params[ky]['processID'])
    tmp=tmp.replace('@@TREENAME',params[ky]['TREE'])
    tmp=tmp.replace('@@WeightScale',str(params[ky]['weightScale']))
    ofname=allDataFiles[ky].split('/')[-1]
    
    tree=params[ky]['TREE'].split('/')[-1]
    if doFID:
        if 'fID' in params[ky]:
            tree=params[ky]['fID']['proc']
            if 'mass' in params[ky]['fID'] :
                if params[ky]['fID']['mass']:
                    tree+='_'+params[ky]['fID']['mass']
            if 'sqrts' in params[ky]['fID']:
                if params[ky]['fID']['sqrts']:
                    tree+='_'+params[ky]['fID']['sqrts']
    resetWeight=-1e5
    if 'resetWeight' in params[ky]:
        resetWeight=params[ky]['resetWeight']
    tmp=tmp.replace('@@resetWeight',str(resetWeight))
    tmp=tmp.replace('@@OTREE',tree)
    tmp=tmp.replace('@@OTREE',tree)
    tmp=tmp.replace('@@OBASE',obase)
    tmp=tmp.replace('@@OFNAME',ofname)
    tmp=tmp.replace('@@DoPtReWeighting',params[ky]['doPtReWeighting'])
    cfgName=base+'/'+ky+'.cfg'
    f= open(cfgName,'w')
    f.write(tmp)
    fSub.write( 'python -W ignore '+executable+' '+cfgName+' > logs/' +ky+ '.log &\n'  )

fSub.close()

print("sub at ",subName)


""" 
    withe correcyted weights
    "gluGluToHH"     : { 'processID':'MC'   ,'TREE':'tagsDumper/trees/bkg_13TeV_TrippleHTag_0' , 'resetWeight' : 1.0 , 'weightScale': 31.05,'fID':{'proc' : 'vHH'  , 'mass' : '125', 'sqrts' :'13TeV' } } ,
    "vbfHH"          : { 'processID':'MC'   ,'TREE':'tagsDumper/trees/bkg_13TeV_TrippleHTag_0' , 'resetWeight' : 1.0 , 'weightScale': 1.726,'fID':{'proc' : 'vbfHH', 'mass' : '125', 'sqrts' :'13TeV' } } ,
    "WToQQHHTo2B2G"  : { 'processID':'MC'   ,'TREE':'tagsDumper/trees/MC_13TeV_TrippleHTag_0'  , 'resetWeight' : 1.0 , 'weightScale': 7.78e-4,'fID':{'proc' : 'vHH'  , 'mass' : '125', 'sqrts' :'13TeV' } } ,
    "ZToBBHHTo2B2G"  : { 'processID':'MC'   ,'TREE':'tagsDumper/trees/MC_13TeV_TrippleHTag_0'  , 'resetWeight' : 1.0 , 'weightScale': 9.63e-5,'fID':{'proc' : 'vHH'  , 'mass' : '125', 'sqrts' :'13TeV' } } ,
    "ggZTo2BHHTo2B2G": { 'processID':'MC'   ,'TREE':'tagsDumper/trees/MC_13TeV_TrippleHTag_0'  , 'resetWeight' : 1.0 , 'weightScale': 1.50e-5,'fID':{'proc' : 'zHH'  , 'mass' : '125', 'sqrts' :'13TeV' } } ,
    "ttHHTo2B2G"     : { 'processID':'MC'   ,'TREE':'tagsDumper/trees/MC_13TeV_TrippleHTag_0'  , 'resetWeight' : 1.0 , 'weightScale': 1.78e-3,'fID':{'proc' : 'ttHH' , 'mass' : '125', 'sqrts' :'13TeV' } } ,
    "ttWToQQHTo2G"   : { 'processID':'MC'   ,'TREE':'tagsDumper/trees/MC_13TeV_TrippleHTag_0'  , 'resetWeight' : 1.0 , 'weightScale': 1.76e-3,'fID':{'proc' : 'ttWH' , 'mass' : '125', 'sqrts' :'13TeV' } } ,
    'signal'         : { 'processID':'MC'   ,'TREE':'tagsDumper/trees/MC_13TeV_TrippleHTag_0'  , 'resetWeight' : 1.0 , 'weightScale': 1.14e-4,'fID':{'proc' : 'ggHHH', 'mass' : '125', 'sqrts' :'13TeV' } } ,

"""

