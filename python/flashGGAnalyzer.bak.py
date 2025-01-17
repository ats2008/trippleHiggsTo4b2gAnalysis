from __future__ import print_function
import ROOT 
import numpy as np
from trippleHiggsSelector import *
from Util  import *

import os,sys

cfgFileName=''
if len(sys.argv) <2:
    print("Usage\n\t ~$python recoAnalyzer.py <configFile>\n")
    exit(1)
else:
    cfgFileName=sys.argv[1]
maxEvtsSuperSeeder=-1
if len(sys.argv) >2:
    maxEvtsSuperSeeder=int(sys.argv[2])
print("Loading cfg file ",cfgFileName)
f=open(cfgFileName,'r')
cfgTxt=f.readlines()
f.close()

allFnames=getListOfStringsFromConfigs(cfgTxt,"#FNAMES_BEG","#FNAMES_END")
foutName=getValueFromConfigs(cfgTxt,"OutpuFileName")

maxEvents=-1
tmp_=getValueFromConfigs(cfgTxt,"MaxEvents")
if tmp_!='':
    maxEvents=int(tmp_)

for i in allFnames:
    print(" file : ",i)
if(maxEvtsSuperSeeder > 0):
    maxEvents=maxEvtsSuperSeeder
print("maxevents : ",maxEvents)

histStore = getRecoHistos()

nNoHtoGammaGamma=0
nNoHHto4B=0
totalEvents=0

isMC = True
isMC = False

histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_tightID"]= histStore["hhhTo4b2gamma"]["h1MassVsh2Mass"].Clone()
histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_tightID"].SetName("h1MassVsh2Mass_tightID")

histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_mediumID"]= histStore["hhhTo4b2gamma"]["h1MassVsh2Mass"].Clone()
histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_mediumID"].SetName("h1MassVsh2Mass_mediumID")

histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_looseID"]= histStore["hhhTo4b2gamma"]["h1MassVsh2Mass"].Clone()
histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_looseID"].SetName("h1MassVsh2Mass_looseID")

histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_looseID"]= histStore["hhhTo4b2gamma"]["h1MassVsh2Mass"].Clone()
histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_looseID"].SetName("h1MassVsh2Mass_looseID")

histStore["h1Tobb"]["mass_tightID" ] =histStore["h1Tobb"]["mass" ].Clone()
histStore["h2Tobb"]["mass_tightID" ] =histStore["h2Tobb"]["mass" ].Clone()
histStore["h1Tobb"]["mass_tightID" ].SetName("mass_tightID" )
histStore["h2Tobb"]["mass_tightID" ].SetName("mass_tightID" )

histStore["h1Tobb"]["mass_mediumID" ] =histStore["h1Tobb"]["mass" ].Clone()
histStore["h2Tobb"]["mass_mediumID" ] =histStore["h2Tobb"]["mass" ].Clone()
histStore["h1Tobb"]["mass_mediumID" ].SetName("mass_mediumID" )
histStore["h2Tobb"]["mass_mediumID" ].SetName("mass_mediumID")

histStore["h1Tobb"]["mass_looseID" ] =histStore["h1Tobb"]["mass" ].Clone()
histStore["h2Tobb"]["mass_looseID" ] =histStore["h2Tobb"]["mass" ].Clone()
histStore["h1Tobb"]["mass_looseID" ].SetName("mass_looseID" )
histStore["h2Tobb"]["mass_looseID" ].SetName("mass_looseID" )

for fname in allFnames:
    
    simFile = ROOT.TFile(fname,'READ')
    eTree=simFile.Get('tagsDumper/trees/MC_13TeV_TrippleHTag_0')
    if not eTree:
        eTree=simFile.Get('tagsDumper/trees/Data_13TeV_TrippleHTag_0')
    maxEvents_ = eTree.GetEntries()
    if(maxEvents >0  and (totalEvents+maxEvents_) > maxEvents):
        maxEvents_= (maxEvents - totalEvents)
    totalEvents+=maxEvents_
    
    for i in range(maxEvents_):
        eTree.GetEntry(i)
        if(i%1000==0):
            print("Doing i = ",i," / ",maxEvents_,
                  " n No HiggsCand : ", nNoHtoGammaGamma,
                  " n No 4BCands  : ", nNoHHto4B)
            print(" gg mass , h1MassVsh2Mass  : ",eTree.CMS_hgg_mass ,"( ",eTree.M1jj , eTree.M2jj,")" )
            
        histStore['events']['nEvents'].Fill("nEvents",1);
        histStore["diPhotons"]["mass"].Fill( eTree.CMS_hgg_mass  )
        
        histStore["hhhTo4b2gamma"]["h1MassVsh2Mass"].Fill(eTree.M1jj , eTree.M2jj )
        histStore["h1Tobb"]["mass"].Fill( eTree.M1jj )
        histStore["h2Tobb"]["mass"].Fill( eTree.M2jj )
        continue

        if(eTree.h1LeadingJet_DeepFlavour > 0.7100 and eTree.h1SubleadingJet_DeepFlavour > 0.7100 
            and  eTree.h2LeadingJet_DeepFlavour > 0.7100 and eTree.h2SubleadingJet_DeepFlavour > 0.7100 ):
            histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_tightID"].Fill(eTree.M1jj , eTree.M2jj )
            histStore["h1Tobb"]["mass_tightID" ].Fill(eTree.M1jj)
            histStore["h2Tobb"]["mass_tightID" ].Fill(eTree.M2jj)

        if(eTree.h1LeadingJet_DeepFlavour > 0.2783 and eTree.h1SubleadingJet_DeepFlavour > 0.2783 
            and  eTree.h2LeadingJet_DeepFlavour > 0.2783 and eTree.h2SubleadingJet_DeepFlavour > 0.2783 ):
            histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_mediumID"].Fill(eTree.M1jj , eTree.M2jj )
            histStore["h1Tobb"]["mass_mediumID" ].Fill(eTree.M1jj)
            histStore["h2Tobb"]["mass_mediumID" ].Fill(eTree.M2jj)

        
        if(eTree.h1LeadingJet_DeepFlavour > 0.0490 and eTree.h1SubleadingJet_DeepFlavour > 0.0490 
            and  eTree.h2LeadingJet_DeepFlavour > 0.0490 and eTree.h2SubleadingJet_DeepFlavour > 0.0490 ):
            histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_looseID"].Fill(eTree.M1jj , eTree.M2jj )
            histStore["h1Tobb"]["mass_looseID" ].Fill(eTree.M1jj)
            histStore["h2Tobb"]["mass_looseID" ].Fill(eTree.M2jj)

        
        
        ###               SKIPPING          ###
        continue
        
        if( not triggerSelection(eTree)):        continue
        histStore['events']['nEvents'].Fill("nTriggered",1);

        diPhotons_=getDiPhotons(eTree)
       
        if(not diPhotons_['valid']):
            nNoHtoGammaGamma+=1
            #print("!! no h->yy cand ")
            continue
        histStore['events']['nEvents'].Fill("nDiPhotonSelected",1);
       
        idx=diPhotons_['diPhotonCand']
        if idx < 0:
            continue
        histStore['events']['nEvents'].Fill("nDiPhotonCands",1);
        histStore['events']['nEvents'].Fill("nDiPhotonCandsInWindow120to130",1);
       
        bJets=getBJetParis(eTree,diPhotons_['diPhotons'][idx]['index'])
        if(not bJets['valid']): 
            #print("!! no hh->4b cand ")
            nNoHHto4B+=1
            continue
        
        quad=bJets['bJetQuad']
        histStore['events']['nEvents'].Fill("nBJetQuads",1);
        if quad['r_HH'] < 25.0 :
            histStore['events']['nEvents'].Fill("nDiPhotonCandEventsIn25GevMbbSignalWindow",1);
        
        diPhotons=diPhotons_['diPhotons']

        histStore["diPhotons"]["pT"].Fill( diPhotons[idx]["pT"])
        histStore["diPhotons"]["eta"].Fill( diPhotons[idx]["eta"])
        histStore["diPhotons"]["y"].Fill( diPhotons[idx]["y"])
        histStore["diPhotons"]["phi"].Fill( diPhotons[idx]["phi"])
        histStore["diPhotons"]["g1g2_Dr"].Fill( diPhotons[idx]["g1g2_Dr"])
        histStore["diPhotons"]["detIdx"].Fill( diPhotons[idx]["detIdx"])
        g1,g2=diPhotons[idx]['index']
        histStore["diPhotons"]["gamma1_pT"].Fill( eTree.photons_pt[g1])
        histStore["diPhotons"]["gamma1_eta"].Fill( eTree.photons_eta[g1])
        histStore["diPhotons"]["gamma1_phi"].Fill( eTree.photons_phi[g1])
        histStore["diPhotons"]["gamma2_pT"].Fill( eTree.photons_pt[g2])
        histStore["diPhotons"]["gamma2_eta"].Fill( eTree.photons_eta[g2])
        histStore["diPhotons"]["gamma2_phi"].Fill( eTree.photons_phi[g2])
    
        i1,i2,i3,i4=quad['idxs']
        idxs=[]
        if isMC:
            #print("!! Got Quad r_HH : ",quad['r_HH'])
            #print("!! Got Quad DH1_min : ",quad['D_HH'])
            isValid,hasHH,hasLeadHiggs,hasSubLeadHiggs=checkGenMatches(eTree,quad['idxs'])
            idxs=getBestGetMatchesGlobal(eTree)
            if(quad['r_HH'] < 25.0):
                hasRecoedH1=True
                hasRecoedH2=True
                hasRecoedHH=True
                if idxs[0] < 0 or idxs[1] <0 :
                    hasRecoedH1=False
                if idxs[2] < 0 or idxs[3] <0 :
                    hasRecoedH2=False
                hasRecoedHH = hasRecoedH1 and hasRecoedH2
                #print("!! Gen Match Result : ",isValid," : ",hasHH,hasLeadHiggs,hasSubLeadHiggs)
                if hasRecoedH1:
                    histStore["events"]["nHiggsMatch"].Fill("hasAllRecoedH1",1)      
                if hasRecoedH2:
                    histStore["events"]["nHiggsMatch"].Fill("hasAllRecoedH2",1)      
                if hasRecoedHH:
                    histStore["events"]["nHiggsMatch"].Fill("hasAllRecoedHH",1)      
                if isValid :
                    histStore["events"]["nHiggsMatch"].Fill("AllInSR",1)      
                if isValid and hasHH:
                    histStore["events"]["nHiggsMatch"].Fill("HH",1)      
                if isValid and hasLeadHiggs:
                    histStore["events"]["nHiggsMatch"].Fill("H1",1)      
                if isValid and hasSubLeadHiggs:
                    histStore["events"]["nHiggsMatch"].Fill("H2",1)      
        
        histStore["h1Tobb"]["mass"].Fill( quad["m1"] )
        histStore["h1Tobb"]["mass_preReg"].Fill( quad["m1_preReg"] )
        histStore["h1Tobb"]["pT"].Fill( quad["p4_h1"].Pt() )
        histStore["h1Tobb"]["eta"].Fill( quad["p4_h1"].Eta() )
        histStore["h1Tobb"]["y"].Fill( quad["p4_h1"].Rapidity() )
        histStore["h1Tobb"]["phi"].Fill( quad["p4_h1"].Phi() )
        
        histStore["h2Tobb"]["mass"].Fill( quad["m2"] )
        histStore["h2Tobb"]["mass_preReg"].Fill( quad["m2_preReg"] )
        histStore["h2Tobb"]["pT"].Fill( quad["p4_h2"].Pt() )
        histStore["h2Tobb"]["eta"].Fill( quad["p4_h2"].Eta() )
        histStore["h2Tobb"]["y"].Fill( quad["p4_h2"].Rapidity() )
        histStore["h2Tobb"]["phi"].Fill( quad["p4_h2"].Phi() )

        for htag,hJets in zip(['1','2'] , [[i1,i2],[i3,i4]])
            for jtag, idx in zip(['1','2'] , hJets): 
                sidx=str(idx)
                histStore["h"+htag+"Tobb"]["bJet"+jtag+"_pT_preReg"].Fill( getattr(eTree,'jets_'+sidx+'_pt'[i1]) )
                histStore["h"+htag+"Tobb"]["bJet"+jtag+"_pT"].Fill(        getattr(eTree,'jets_'+sidx+'_pt'*eTree.jets_bJetRegCorr[i1]) )
                histStore["h"+htag+"Tobb"]["bJet"+jtag+"_eta"].Fill(       getattr(eTree,'jets_'+sidx+'_eta'))
                histStore["h"+htag+"Tobb"]["bJet"+jtag+"_phi"].Fill(       getattr(eTree,'jets_'+sidx+'_phi'))
                histStore["h"+htag+"Tobb"]["bJet"+jtag+"_deepCSVScore"].Fill( getattr(eTree,'jets_'+sidx+'deepCSVScore') );
                histStore["h"+htag+"Tobb"]["bJet"+jtag+"_deepJetScore"].Fill( getattr(eTree,'jets_'+sidx+'deepJetScore') );
                histStore["h"+htag+"Tobb"]["bJet"+jtag+"_hFlavour"].Fill( getattr(eTree,'jets_'+sidx+'_flavour' ) );
                histStore["h"+htag+"Tobb"]["bJet"+jtag+"_pFlavour"].Fill( getattr(eTree,'jets_'+sidx+'_pFlavour') );
        
        histStore["h1Tobb"]["b1b2_Dr"].Fill(deltaR(eTree.jets_eta[i1],eTree.jets_phi[i1],eTree.jets_eta[i2],eTree.jets_phi[i2]) )
        histStore["h2Tobb"]["b1b2_Dr"].Fill(deltaR(eTree.jets_eta[i3],eTree.jets_phi[i3],eTree.jets_eta[i4],eTree.jets_phi[i4]) )
        
        #hhTo4b
        histStore["hhTo4b"]["mass"].Fill( quad["mass"] )
        histStore["hhTo4b"]["mass_preReg"].Fill( quad["mass_preReg"] )
        histStore["hhTo4b"]["pT"].Fill( quad["pT"] )
        histStore["hhTo4b"]["eta"].Fill( quad["eta"] )
        histStore["hhTo4b"]["y"].Fill( quad["y"] )
        histStore["hhTo4b"]["phi"].Fill( quad["phi"] )
        
        #hhhTo4b2gamma
        allDauP4 = (quad['p4_h1']+quad['p4_h2']+diPhotons[idx]['p4'])
        allDauP4_preReg = (quad['p4_h1_preReg']+quad['p4_h2_preReg']+diPhotons[idx]['p4'])
        

        histStore["hhhTo4b2gamma"]["mass"].Fill( allDauP4.M() )
        histStore["hhhTo4b2gamma"]["mass_preReg"].Fill( allDauP4_preReg.M()  )
        histStore["hhhTo4b2gamma"]["pT"].Fill( allDauP4.Pt() )
        histStore["hhhTo4b2gamma"]["eta"].Fill( allDauP4.Eta() )
        histStore["hhhTo4b2gamma"]["y"].Fill( allDauP4.Rapidity() )
        histStore["hhhTo4b2gamma"]["phi"].Fill( allDauP4.Phi() )
        histStore["hhhTo4b2gamma"]["h1MassVsh2Mass"].Fill(quad['m1'],quad['m2'] )
        histStore["hhhTo4b2gamma"]["h1MassVsh2Mass_preReg"].Fill(quad['m1_preReg'],quad['m2_preReg'] )
        
        mx0 =  allDauP4.M() - ( quad['p4_h1'].M()+quad['p4_h2'].M()+diPhotons[idx]['p4'].M() - 375.0 )
        mx1 = (quad['p4_h1']+quad['p4_h2']).M() - ( quad['p4_h1'].M()+quad['p4_h2'].M() - 250.0 ) 
        mx2 = (quad['p4_h1']+diPhotons[idx]['p4']).M()-( quad['p4_h1'].M()+diPhotons[idx]['p4'].M() - 250.0 ) 
        mx3 = (quad['p4_h2']+diPhotons[idx]['p4']).M()-( quad['p4_h2'].M()+diPhotons[idx]['p4'].M() - 250.0 ) 

        histStore["hhhTo4b2gamma"]["massX0"].Fill( mx0 )
        histStore["hhhTo4b2gamma"]["massX1"].Fill( mx1 )
        histStore["hhhTo4b2gamma"]["massX2"].Fill( mx2 )
        histStore["hhhTo4b2gamma"]["massX3"].Fill( mx3 )

        histStore["hhhTo4b2gamma"]["massX0X1"].Fill( mx0,mx1  )
        histStore["hhhTo4b2gamma"]["massX2X3"].Fill( mx2,mx3  )
        histStore["hhhTo4b2gamma"]["massX0X2"].Fill( mx0,mx2  )
        histStore["hhhTo4b2gamma"]["massX1X3"].Fill( mx1,mx3  )
        histStore["hhhTo4b2gamma"]["massX0X3"].Fill( mx0,mx3  )
        histStore["hhhTo4b2gamma"]["massX1X2"].Fill( mx1,mx2  )
        
        diPhoton=diPhotons[idx]
        bjetQuad=quad
        fillTrippleHVariables(eTree,histStore,bjetQuad,diPhoton);

        jetPts=np.array(eTree.jets_pt)
        mask= jetPts>25
        jetEta=np.array(eTree.jets_eta)
        mask=np.logical_and(mask , abs(jetEta) < 2.5 )
        nGoodJets=sum(mask)
        histStore['hhhTo4b2gamma']['eventJetMultiplicity'].Fill(nGoodJets)
    simFile.Close()           
saveTheDictionary(histStore,foutName)
print(" File written out  : ",foutName)

