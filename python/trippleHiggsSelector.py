from __future__ import print_function
import ROOT 
import numpy as np
import itertools as itrTools
from Util import *
from trippleHiggsUtils import *
"""
from inspect import currentframe

def get_lineNumber():
    cf = currentframe()
    return cf.f_back.f_lineno
"""

def getNBsFromCand(eTree):
    bCount=0
    if  eTree.h1LeadingJet_hflav==5:
        bCount+=1
    if  eTree.h1SubleadingJet_hflav==5:
        bCount+=1
    if  eTree.h2LeadingJet_hflav==5:
        bCount+=1
    if  eTree.h2SubleadingJet_hflav==5:
        bCount+=1
    return bCount

def triggerSelection(eTree):
    if(eTree.nPhotons < 2): return False
    if(eTree.photons_pt[0] < 30) :return False;
    if(eTree.photons_pt[1] < 22): return False;
    pho1LV=ROOT.TLorentzVector()
    pho2LV=ROOT.TLorentzVector()
    pho1LV.SetPtEtaPhiE(eTree.photons_pt[0],eTree.photons_eta[0],
                        eTree.photons_phi[0],eTree.photons_e[0])
    pho2LV.SetPtEtaPhiE(eTree.photons_pt[1],eTree.photons_eta[1],
                        eTree.photons_phi[1],eTree.photons_e[1])
    if((pho1LV+pho2LV).M() < 90 ): return False
    return True
 
def getHiggsDauP4s(eTree,pdgId):
    genPt =[1e4 for i in range(4)]
    genEta=[1e4 for i in range(4)]
    genPhi=[1e4 for i in range(4)]
    genE  =[1e4 for i in range(4)]
    genHPt=[-1e4 for i in range(4)]
    k=0
    #print(eTree.gen_H1_dau1_pdgId,eTree.gen_H2_dau1_pdgId,eTree.gen_H3_dau1_pdgId)
    if abs(eTree.gen_H1_dau1_pdgId)==pdgId:
        genPt[k]=eTree.gen_H1_dau1_pt
        genEta[k]=eTree.gen_H1_dau1_eta
        genPhi[k]=eTree.gen_H1_dau1_phi
        genE[k]=eTree.gen_H1_dau1_e
        k+=1
        genPt[k]=eTree.gen_H1_dau2_pt
        genEta[k]=eTree.gen_H1_dau2_eta
        genPhi[k]=eTree.gen_H1_dau2_phi
        genE[k]=eTree.gen_H1_dau2_e
        genHPt[k/2]=eTree.gen_H1_pt
        k+=1
    if abs(eTree.gen_H2_dau1_pdgId)==pdgId:
        genPt[k]=eTree.gen_H2_dau1_pt
        genEta[k]=eTree.gen_H2_dau1_eta
        genPhi[k]=eTree.gen_H2_dau1_phi
        genE[k]=eTree.gen_H2_dau1_e
        k+=1
        genPt[k]=eTree.gen_H2_dau2_pt
        genEta[k]=eTree.gen_H2_dau2_eta
        genPhi[k]=eTree.gen_H2_dau2_phi
        genE[k]=eTree.gen_H2_dau2_e
        genHPt[k/2]=eTree.gen_H2_pt
        k+=1

    if abs(eTree.gen_H3_dau1_pdgId)==pdgId:
        genPt[k]=eTree.gen_H3_dau1_pt
        genEta[k]=eTree.gen_H3_dau1_eta
        genPhi[k]=eTree.gen_H3_dau1_phi
        genE[k]=eTree.gen_H3_dau1_e
        k+=1
        genPt[k]=eTree.gen_H3_dau2_pt
        genEta[k]=eTree.gen_H3_dau2_eta
        genPhi[k]=eTree.gen_H3_dau2_phi
        genE[k]=eTree.gen_H3_dau2_e
        genHPt[k/2]=eTree.gen_H3_pt
        k+=1
    allP4=[]
    for i in range(k):
      lv=ROOT.TLorentzVector()
      lv.SetPtEtaPhiE(genPt[i],genEta[i],genPhi[i],genE[i])
      allP4.append(lv)
    
    order=[0,1]
    if allP4[order[0]].Pt() < allP4[order[1]].Pt():
        order=[1,0]
    if len(allP4) > 2:
        #print("\n in getDau : ",end= "")
        #for ii in range(4):
        #    print(ii," ) ",np.round(allP4[ii].Pt(),2 ), " | ",end="")
        if genHPt[1] > genHPt[0]:
            order=[2,3,order[0],order[1]]
        else:
            order=[order[0],order[1],2,3]
        if allP4[order[0]].Pt() < allP4[order[1]].Pt():
            t=order[1] ; order[1]=order[0] ; order[0]=t;
        #print("\tallP4[",order[2]",].Pt() < allP4[",order[3],"].Pt() ",allP4[order[2]].Pt(),"  < "allP4[order[3]].Pt())
        if allP4[order[2]].Pt() < allP4[order[3]].Pt():
            t=order[2] ; order[2]=order[3] ; order[3]=t;
           
        #print("  order : ",order ," ) \n\t",end="")
        #for ii in range(4):
        #    print(ii," ) ",np.round(allP4[order[ii]].Pt(),2 ), " | ",end="")

    return [allP4[i] for i in order ]
    
def checkGenMatches(eTree,jetIndices):
    idxs=[-1 for i in range(4)]
    genPt=[1e4 for i in range(4)]
    genEta=[1e4 for i in range(4)]
    genPhi=[1e4 for i in range(4)]
    idxs=[-1 for i in range(4)]
    #print("Jet Indices : ",jetIndices)
    jetRegCorr=np.array(eTree.jets_bJetRegCorr)[jetIndices]
    jetPt=np.array(eTree.jets_pt)[jetIndices]
    jetEta=np.array(eTree.jets_eta)[jetIndices]
    jetPhi=np.array(eTree.jets_phi)[jetIndices]
    jetPt =np.array(eTree.jets_pt)[jetIndices]

    k=0
    #print("H1  Dau 1" , eTree.gen_H1_dau1_pt ,eTree.gen_H1_dau1_eta  , eTree.gen_H1_dau1_phi ,eTree.gen_H1_dau1_pdgId )
    #print("H1  Dau 1" , eTree.gen_H1_dau2_pt ,eTree.gen_H1_dau2_eta  , eTree.gen_H1_dau1_phi ,eTree.gen_H1_dau1_pdgId )
    #print("H2  Dau 2" , eTree.gen_H2_dau1_pt ,eTree.gen_H2_dau1_eta  , eTree.gen_H2_dau1_phi ,eTree.gen_H2_dau1_pdgId )
    #print("H3  Dau 2" , eTree.gen_H2_dau2_pt ,eTree.gen_H2_dau2_eta  , eTree.gen_H2_dau2_phi ,eTree.gen_H2_dau2_pdgId )
    #print("H3  Dau 2" , eTree.gen_H3_dau1_pt ,eTree.gen_H3_dau1_eta  , eTree.gen_H3_dau1_phi ,eTree.gen_H3_dau1_pdgId )
    #print("H3  Dau 2" , eTree.gen_H3_dau2_pt ,eTree.gen_H3_dau2_eta  , eTree.gen_H3_dau2_phi ,eTree.gen_H3_dau2_pdgId )

    if(abs(eTree.gen_H1_dau1_eta) > 2.5 ): return False,False,False,False
    if(abs(eTree.gen_H1_dau2_eta) > 2.5 ): return False,False,False,False
    if(abs(eTree.gen_H2_dau1_eta) > 2.5 ): return False,False,False,False
    if(abs(eTree.gen_H2_dau2_eta) > 2.5 ): return False,False,False,False
    if(abs(eTree.gen_H3_dau1_eta) > 2.5 ): return False,False,False,False
    if(abs(eTree.gen_H3_dau2_eta) > 2.5 ): return False,False,False,False

    #for i in range(eTree.nJets):
    #    print("J : ",i," --> ",eTree.jets_pt[i]*eTree.jets_bJetRegCorr[i] ,eTree.jets_eta[i]  , eTree.jets_phi[i] ," || ", eTree.jets_deepJetScore[i])

    if abs(eTree.gen_H1_dau1_pdgId)==5:
        genPt[k]=eTree.gen_H1_dau1_pt
        genEta[k]=eTree.gen_H1_dau1_eta
        genPhi[k]=eTree.gen_H1_dau1_phi
        k+=1
        genPt[k]=eTree.gen_H1_dau2_pt
        genEta[k]=eTree.gen_H1_dau2_eta
        genPhi[k]=eTree.gen_H1_dau2_phi
        k+=1
    if abs(eTree.gen_H2_dau1_pdgId)==5:
        genPt[k]=eTree.gen_H2_dau1_pt
        genEta[k]=eTree.gen_H2_dau1_eta
        genPhi[k]=eTree.gen_H2_dau1_phi
        k+=1
        genPt[k]=eTree.gen_H2_dau2_pt
        genEta[k]=eTree.gen_H2_dau2_eta
        genPhi[k]=eTree.gen_H2_dau2_phi
        k+=1

    if abs(eTree.gen_H3_dau1_pdgId)==5:
        genPt[k]=eTree.gen_H3_dau1_pt
        genEta[k]=eTree.gen_H3_dau1_eta
        genPhi[k]=eTree.gen_H3_dau1_phi
        k+=1
        genPt[k]=eTree.gen_H3_dau2_pt
        genEta[k]=eTree.gen_H3_dau2_eta
        genPhi[k]=eTree.gen_H3_dau2_phi
        k+=1
    if k > 4:
        print(" GEN GONE WRONG , more than 4 b as higgs daughters !! \n")
        exit(1)
    
    possiblities  =[
                      [0,1,2,3],
                      [1,0,2,3],
                      [0,1,3,2],
                      [1,0,3,2],
                      [2,3,0,1],
                      [2,3,1,0],
                      [3,2,0,1],
                      [3,2,1,0]
                   ]

    isGenMatch=False
    hA_matched=False
    hB_matched=False
 #   for i in range(4):
 #       print("!! \tGen : ",genEta[i],genPhi[i]," | Reco : ",jetEta[i],jetPhi[i],jetPt[i])
    for pos in possiblities:
        i,j,k,l=pos
        mat1=False
        mat2=False
        drA=deltaR(genEta[i],genPhi[i],jetEta[0],jetPhi[0])
        drB=deltaR(genEta[j],genPhi[j],jetEta[1],jetPhi[1])
        
        isPtMatch=False
        if abs(jetPt[0]*jetRegCorr[0]/genPt[i] - 1.0) < 0.2 and  abs(jetPt[1]*jetRegCorr[1]/genPt[j] - 1.0) < 0.2 :
            isPtMatch=True
        isPtMatch=True
        if(drA < 0.4) and (drB<0.4)  and isPtMatch :
            hA_matched = True
            mat1=True
   #     print("!!\t\t Drs ",i,j," -- > ",drA,drB,mat1)
        
        drA=deltaR(genEta[k],genPhi[k],jetEta[2],jetPhi[2])
        drB=deltaR(genEta[l],genPhi[l],jetEta[3],jetPhi[3])

        isPtMatch=False
        if abs(jetPt[2]*jetRegCorr[2]/genPt[k] - 1.0) < 0.2 and  abs(jetPt[3]*jetRegCorr[3]/genPt[l] - 1.0) < 0.2 :
            isPtMatch=True
        isPtMatch=True
        if(drA < 0.4 and drB<0.4  and isPtMatch):
            hB_matched = True
            mat2=True

    #    print("!!\t\t Drs ",k,l," -- > ",drA,drB,mat2)
        if mat1 and mat2:
            isGenMatch=True
            break
    return True,isGenMatch,hA_matched,hB_matched


def getBestGetMatchesGlobalFGG(eTree,jetMask=[],drMax=0.4):
#def getBestGetMatchesGlobalFGG(eTree,drMax=0.4,etaCut=2.5,pTMin=25.0):
    #print("In getBestGetMatchesGlobalFGG")

    idxs  =[-1 for i in range(4)]
    
    jetPt=[]
    jetEta=[]
    jetPhi=[]
    jetMass=[]
    fggIdx=[]
    jetBtagScore=[]
    flav=[]
    N_JET_MAX=8

    if jetMask==[]:
        jetMask=[True for i in range(N_JET_MAX) ]
    for i in range(N_JET_MAX):
        if not jetMask[i]:
            continue
            
        jetPt.append(getattr(eTree,'jet_'+str(i)+'_pt'))
        jetEta.append(getattr(eTree,'jet_'+str(i)+'_eta'))
        jetPhi.append(getattr(eTree,'jet_'+str(i)+'_phi'))
        jetMass.append(getattr(eTree,'jet_'+str(i)+'_mass'))
        jetBtagScore.append(getattr(eTree , 'jet_'+str(i)+'_deepJetScore') ) 
        flav.append(getattr(eTree,'jet_'+str(i)+'_flavour'))
        fggIdx.append(i)

    scoreorder=np.argsort(np.array(jetBtagScore)*-1)
    bJetScoreRankMap={ k : n for n,k in enumerate(scoreorder)}
    jetPt =np.array(jetPt)
    jetEta=np.array(jetEta)
    jetPhi=np.array(jetPhi)
    allGammas=getHiggsDauP4s(eTree,22)
    allBquarks=getHiggsDauP4s(eTree,5)
        
    genPt=[]
    genEta=[]
    genPhi=[]
    for genb in allBquarks:
        genPt.append( genb.Pt() )
        genEta.append(genb.Eta())
        genPhi.append(genb.Phi())

    genPt =np.array(genPt)
    genEta=np.array(genEta)
    genPhi=np.array(genPhi)
    
    drMins=[]
    drMinsPostFlavourMatch=[]
    jetMassMatch=[]
    deepJetScores=[]
    allP4s=[]
    hFlavour=[]
    hFlavourAll=[]

    jetIdLoose=[]
    jetIdTight=[]
    jetIdTight2017=[]
    jetIdTight2018=[]
    jetPUScore=[]

    drMaxs=[drMax for i in range(4)]
    #print("jet Eta ",np.round(jetEta,3))
    #print("jet Phi ",np.round(jetPhi,3))
    if len(jetPt) > 0 :
        for i in range(4):
            drs=deltaR(genEta[i],genPhi[i],jetEta,jetPhi)
            idxDrsWise=np.argsort(drs)
            fggIdxDrWise= [ fggIdx[k] for k in idxDrsWise]
            genB=allBquarks[i]
    #    print()
    #    print("New Jet @gen ",genEta[i],genPhi[i]," [ ",allBquarks[i].Eta(),allBquarks[i].Phi(),"] ")
    #    print(np.round(drs,3))
    #    print(idxDrsWise, " | ",fggIdxDrWise)
    #    print("gen b ",i,"  : ",genB.Pt(),genB.Eta(),genB.Phi())
    #    for n,idx in enumerate(idxDrsWise):
    #        print("\t",idx," : pt ",jetPt[idx]," , eta ",jetEta[idx]," , phi ",jetPhi[idx]," flav : ",flav[idx], " : dr ",drs[idx] ," btag ",jetBtagScore[idx],bJetScoreRankMap[idx])
            drMins.append(drs[idxDrsWise[0]])
            drMinsPostFlavourMatch.append(drs[idxDrsWise[0]])

            jetMassMatch.append(jetMass[ idxDrsWise[0] ])
            p4=ROOT.TLorentzVector()
            p4.SetPtEtaPhiM( jetPt[idxDrsWise[0]] , jetEta[idxDrsWise[0]], jetPhi[idxDrsWise[0]],jetMass[idxDrsWise[0]] )
            allP4s.append(p4)
            fggIdx0=fggIdx[idxDrsWise[0]]
            hFlavourAll.append(getattr(eTree,'jet_'+str( fggIdx0 )+'_flavour'))
            deepJetScores.append(getattr(eTree,'jet_'+str( fggIdx0 )+'_deepCSVScore'))
            jetIdLoose.append(getattr(eTree,'jet_'+str( fggIdx0 )+'_isLoose'))
            jetIdTight.append(getattr(eTree,'jet_'+str( fggIdx0 )+'_isTight'))
            jetIdTight2017.append(getattr(eTree,'jet_'+str( fggIdx0 )+'_isTight2017'))
            jetIdTight2018.append(getattr(eTree,'jet_'+str( fggIdx0 )+'_isTight2018'))
            jetPUScore.append(getattr(eTree,'jet_'+str( fggIdx0 )+'_puJetIdMVA'))
            hFlavour.append(getattr(eTree,'jet_'+str( fggIdx0 )+'_flavour'))
            drMin_=1e9
            drMax = drMaxs[i]
            for fgg_idx,idx in zip(fggIdxDrWise,idxDrsWise) :
                #print("   ",fgg_idx," | ",idx, " mathing flav : ",getattr(eTree,'jet_'+str( fgg_idx )+'_flavour')  )
                if abs(getattr(eTree,'jet_'+str( fgg_idx )+'_flavour')) <4.80:
                    continue
                if drs[idx]<drMinsPostFlavourMatch[-1]:
                    #print("   flavour match found with dr ",drs[idx])
                    drMinsPostFlavourMatch[-1]=drs[idx]
                    hFlavourAll[-1]=getattr(eTree,'jet_'+str( fgg_idx )+'_flavour')
                if drs[idx] > drMax : #and abs(jetPt[idx]*jetRegCorr[idx]/genPt[i] - 1.0) < 0.2:
                    break
                if idxs[i] > -1:
                   if drs[idx] > drMin_:
                        continue
                idxs[i]=idx
                drMin_=drs[idx]
                hFlavour[-1]=getattr(eTree,'jet_'+str( fgg_idx )+'_flavour')
                deepJetScores[-1]=getattr(eTree,'jet_'+str( fgg_idx )+'_deepCSVScore')
                jetIdLoose[-1] = getattr(eTree,'jet_'+str( fgg_idx )+'_isLoose')
                jetIdTight[-1] = getattr(eTree,'jet_'+str( fgg_idx )+'_isTight')
                jetIdTight2017[-1] = getattr(eTree,'jet_'+str( fgg_idx )+'_isTight2017')
                jetIdTight2018[-1] = getattr(eTree,'jet_'+str( fgg_idx )+'_isTight2018')
                jetPUScore[-1] = getattr(eTree,'jet_'+str( fgg_idx )+'_puJetIdMVA')
                allP4s[-1].SetPtEtaPhiM( jetPt[idx] , jetEta[idx], jetPhi[idx],jetMass[idx] ) ;allP4s[-1]=p4;
                jetMassMatch[-1]=jetMass[idx]
                #print("    Setting matched idx as fgg : ",fgg_idx,"  , idx ",idx )
    #print("Allocated idxs ",idxs,"  jetPts : ",jetPt)
    idx=copy.deepcopy(idxs)
    #if idx[1] >-1 :
    #    if (idx[0] >-1) and (jetPt[idx[0]] > jetPt[idx[1]]):
    #        pass
    #    else:
    #        t=idx[0]
    #        print("\t FALSE jetPt[",idx[0],"] > jetPt[",idx[1],"] ", jetPt[idx[0]] ,">", jetPt[idx[1]] ," ==> ", jetPt[idx[0]] > jetPt[idx[1]])
    #        idx[0]=idx[1]
    #        idx[1]=t

    #if idx[3] >-1 :
    #    if (idx[2] >-1) and (jetPt[idx[2]] > jetPt[idx[3]]):
    #        pass
    #    else:
    #        print("\t  jetPt[",idx[0],"] > jetPt[",idx[1],"] ", jetPt[idx[0]] ,">", jetPt[idx[1]] ," ==> ", jetPt[idx[0]] > jetPt[idx[1]])
    #        t=idx[2]
    #        idx[2]=idx[3]
    #        idx[3]=t
    fgg_idxs_matched=[]
    bjetScoreRanks=[]
    for k in idx:
        if k >-1:
            fgg_idxs_matched.append(fggIdx[k])
            bjetScoreRanks.append(bJetScoreRankMap[k])
        else:
            fgg_idxs_matched.append(-1)
            bjetScoreRanks.append(-1)
     
    #print("-->  idxs ",idx,"-->",fgg_idxs_matched," btag --> ",bjetScoreRanks)
    return {
               'idxs':idx , 'fgg_idxs' : fgg_idxs_matched,'bJetScoreRanks': bjetScoreRanks,
               'drMins' : drMins , 'drMinsPostFlavourMatch' : drMinsPostFlavourMatch,
               'jetMass':jetMassMatch , 'p4s': allP4s ,
               'deepJetScores' : deepJetScores , 'hFlavour':hFlavour,'hFlavourAll':hFlavourAll,
               'isLoose': jetIdLoose , 'isTight' : jetIdTight ,
               'isTight2017': jetIdTight2017 , 'isTight2018' : jetIdTight2018 ,
            }

def getBestGetMatchesGlobal(eTree):

    genEta=[1e4 for i in range(4)]
    genPhi=[1e4 for i in range(4)]
    genPt=[1e4 for i in range(4)]
    idxs=[-1 for i in range(4)]

    jetRegCorr=np.array(eTree.jets_bJetRegCorr)
    jetPt=np.array(eTree.jets_pt)
    jetEta=np.array(eTree.jets_eta)
    jetPhi=np.array(eTree.jets_phi)
    
    mask= jetPt>25
    jetEta=np.array(eTree.jets_eta)
    mask=np.logical_and(mask , abs(jetEta) < 2.5 )
    
    jetRegCorr=jetRegCorr[mask]
    jetPt     =jetPt[mask]     
    jetEta    =jetEta[mask]    
    jetPhi    =jetPhi[mask]    
    

    
    allGammas=getHiggsDauP4s(eTree,22)
    allBquarks=getHiggsDauP4s(eTree,5)
    
    
    genPt=[]
    genEta=[]
    genPhi=[]
    for genb in allBquarks:
        genPt.append( genb.Pt())
        genEta.append(genb.Eta())
        genPhi.append(genb.Phi())
    genPt =np.array(genPt)
    genEta=np.array(genEta)
    genPhi=np.array(genPhi)
    drMins=[]
    for i in range(4):
        drs=deltaR(genEta[i],genPhi[i],jetEta,jetPhi)
        idxDrsWise=np.argsort(drs)
        #genB=allBquarks[i]
        #print("gen b ",i,"  : ",genB.Pt(),genB.Eta(),genB.Phi())
        #for idx in idxDrsWise:
        #    print("\t",idx," : ",jetPt[idx]," , ",jetEta[idx]," , ",jetPhi[idx])
        drMins.append(drs[idxDrsWise[0]])

        for idx in idxDrsWise :
            if drs[idx] < 0.4 : #and abs(jetPt[idx]*jetRegCorr[idx]/genPt[i] - 1.0) < 0.2:
                idxs[i]=idx
            else:
                break
    
    idx=copy.deepcopy(idxs)
    if idx[1] >0 :
        if idx[0] >0 and jetPt[idx[0]] > jetPt[idx[1]]:
            pass
        else:
                t=idx[0]
                idx[0]=idx[1]
                idx[1]=t

    if idx[3] >0 :
        if idx[2] >0 and jetPt[idx[2]] > jetPt[idx[3]]:
            pass
        else:
                t=idx[2]
                idx[2]=idx[3]
                idx[3]=t
     
    return {'idxs':idxs,'drMins' : drMins }


def getCorrectedJetP4(p4_raw, corrrectionFactor):
    bLV=ROOT.TLorentzVector()
    bLV.SetPtEtaPhiE(p4_raw.Pt()*corrrectionFactor,p4_raw.Eta(),p4_raw.Phi(),p4_raw.E()*corrrectionFactor)
    return bLV
    
def getDiPhotons(eTree):
    result={'valid':False , 'nDiPhotons':0,'diPhotonCand':-1,'diPhotons':[]}
    if (eTree.nPhotons < 2):
        return result
    diPhoPtMax=-1e3
    for i in range(eTree.nPhotons):
        pho1LV=ROOT.TLorentzVector()
        pho1LV.SetPtEtaPhiE(eTree.photons_pt[i],eTree.photons_eta[i],
                        eTree.photons_phi[i],eTree.photons_e[i])
        if(abs(eTree.photons_eta[i]) > 2.5 ): continue
        if((abs(eTree.photons_eta[i]) > 1.44  ) and (abs(eTree.photons_eta[i]) < 1.567) ): continue

        for j in range(i+1,eTree.nPhotons):
            if(abs(eTree.photons_eta[j]) > 2.5 ): continue
            if((abs(eTree.photons_eta[j]) > 1.44  ) and  (abs(eTree.photons_eta[j]) < 1.567) ): continue
            pho2LV=ROOT.TLorentzVector()
            pho2LV.SetPtEtaPhiE(eTree.photons_pt[j],eTree.photons_eta[j],
                        eTree.photons_phi[j],eTree.photons_e[j])
            diPho_p4=pho1LV+pho2LV
            m=diPho_p4.M()
            if( m < 100.0 ) : continue
            if( m > 180.0 ) : continue
            if(pho1LV.Pt()/m <1/3.0) : continue
            if(pho2LV.Pt()/m <1/4.0) : continue
            
            diPho={'index':(i,j)}
            if pho2LV.Pt() > pho1LV.Pt():
                diPho['index']=(j,i)
            diPho['mass']=m
            diPho['p4']=diPho_p4
            diPho['pT']=diPho_p4.Pt()
            diPho['eta']=diPho_p4.Eta()
            diPho['y']=diPho_p4.Rapidity()
            diPho['phi']=diPho_p4.Phi()
            diPho['g1g2_Dr']=pho2LV.DeltaR(pho1LV)
            detId=0
            if(abs(pho1LV.Eta()) < 1.4 and abs(pho2LV.Eta()) < 1.4 ): detId=0
            elif(abs(pho1LV.Eta()) < 1.4 and abs(pho2LV.Eta()) > 1.4 ) : detId=1
            elif(abs(pho1LV.Eta()) > 1.4 and abs(pho2LV.Eta()) < 1.4 ):  detId=2
            elif(abs(pho1LV.Eta()) > 1.4 and abs(pho2LV.Eta()) > 1.4 ) : detId=3
            else:    detId=4
            diPho['detIdx'] = detId
            result['valid']=True
            
            if diPhoPtMax < diPho_p4.Pt():
                diPhoPtMax=diPho_p4.Pt()
                result['diPhotonCand']=result['nDiPhotons']
                
            result['diPhotons'].append(diPho)
            result['nDiPhotons']+=1
    return result

def getBJetParisFGG(eTree,mask=[],doMisclassificationCorrection=False):
    result={'isValid':False , 'allBJetQuads':[],'bJetQuad':[],'nJetQuads':0,'fail':'Pass'}
    X0OverY0=1.05
    jets_pt=[]
    jets_eta=[]
    jets_phi=[]
    jets_mass=[]
    jets_deepJetScore=[]
    jets_fgg_index=[]
    N_JET_MAX=8 

    if len(mask)==0:
        mask=np.zeros(N_JET_MAX , dtype=bool)
    for i in range(N_JET_MAX):
        if (getattr(eTree,'jet_'+str(i)+'_isValid') < 0.25):
            continue
        if not mask[i]:
            continue

        jets_pt.append(getattr(eTree,'jet_'+str(i)+'_pt'))
        jets_eta.append(getattr(eTree,'jet_'+str(i)+'_eta'))
        jets_phi.append(getattr(eTree,'jet_'+str(i)+'_phi'))
        jets_mass.append(getattr(eTree,'jet_'+str(i)+'_mass'))
        jets_deepJetScore.append( getattr(eTree,'jet_'+str( i )+'_deepJetScore') )
        jets_fgg_index.append( i  )
    #jets_pt =np.array(jets_pt)
    #jets_eta=np.array(jets_eta)
    #jets_phi=np.array(jets_phi)
    
    photons_eta=[ eTree.leadingPhoton_eta ,eTree.subleadingPhoton_eta  ]
    photons_phi=[ eTree.leadingPhoton_eta ,eTree.subleadingPhoton_phi  ]
    nJets=len(jets_pt)
    a=X0OverY0
    b=np.sqrt(1.0+X0OverY0*X0OverY0)
    
    if(nJets < 4 ) : 
        #print(" REQUIRES ATLEAST 4 JETS FOR MAKING CADIDATES ! pairing failed for event id ",eTree.event, "[ :( ]")
        result['fail']='nValidJets'
        return result 
    
    
    jetPts=np.array(jets_pt)
    jetEta=np.array(jets_eta)
    jetPhi=np.array(jets_phi)
    jet_deepJetScore=np.array(jets_deepJetScore)
    jetScoreOrder=np.argsort(jet_deepJetScore*-1.0)
    
    mask = jetPts>0

    #print("\t Score of jets under consideratio n : ",jet_deepJetScore)
    #print("\t Sorted ordering of jets under consideratio n : ",jetScoreOrder)
    #print("\t Sorted ordering of jets under consideratio n : ",jet_deepJetScore[jetScoreOrder])

    nMax =  4
    mask[jetScoreOrder[nMax:]]=False
    goodJets=np.where(mask)[0]
    if(len(goodJets) < 4 ): 
        print("SANITY CHECK FAILS !! att FGG BJet Pair maker ")
        result['fail']='goodJets'
        return result
    allJetCombinations=itrTools.combinations(goodJets,4)
    jetsIdx=np.arange(0,4,1)
    jetLVs=[]
    for i in range(4):
        jetLVs.append(ROOT.TLorentzVector())
    possiblilities_=[[0,1,2,3],[0,2,1,3],[0,3,1,2]]
    
    jetCobinationScores=[]
    metric_min=1e9
    etric_min_idx=-1
 #   print("Good jets  : ",goodJets)
    for jetCombination in allJetCombinations:
        # setting 4 jet values
 #       print("Doing Jet Combination : ",jetCombination)
        for j in range(4):
            ii=jetCombination[j]
            jetLVs[j].SetPtEtaPhiM(jets_pt[ii],jets_eta[ii],
                           jets_phi[ii],jets_mass[ii])
        # scanning combination in 4 jets
        for combi in possiblilities_:
#            print("  Permutation  : ",combi)
            quad_={}
            quad_['idxs']=[jetCombination[i] for i in combi]

            if jets_pt[quad_['idxs'][0]] < jets_pt[quad_['idxs'][1]]:
                t=quad_['idxs'][1] ; quad_['idxs'][1]=quad_['idxs'][0] ;quad_['idxs'][0]=t
                t=combi[1] ; combi[1]=combi[0] ;combi[0]=t

            if jets_pt[quad_['idxs'][2]] < jets_pt[quad_['idxs'][3]]:
                t=quad_['idxs'][3] ; quad_['idxs'][3]=quad_['idxs'][2] ;quad_['idxs'][2]=t
                t=combi[3] ; combi[3]=combi[2] ;combi[2]=t
            quad_['fgg_idxs']  = [ jets_fgg_index[ idx ] for idx in quad_['idxs'] ]         
            correctedLV0=getCorrectedJetP4(jetLVs[combi[0]],1.0)
            correctedLV1=getCorrectedJetP4(jetLVs[combi[1]],1.0)
            correctedLV2=getCorrectedJetP4(jetLVs[combi[2]],1.0)
            correctedLV3=getCorrectedJetP4(jetLVs[combi[3]],1.0)
            
            p4_h1_preReg=jetLVs[combi[0]]+jetLVs[combi[1]]
            p4_h2_preReg=jetLVs[combi[2]]+jetLVs[combi[3]]
            p4_h1=correctedLV0 + correctedLV1
            p4_h2=correctedLV2 + correctedLV3

            if p4_h2.Pt() > p4_h1.Pt():
                t=p4_h1_preReg
                p4_h1_preReg=p4_h2_preReg
                p4_h2_preReg=t

                t=p4_h1
                p4_h1=p4_h2
                p4_h2=t

            quad_['p4_h1_preReg']=p4_h1_preReg
            quad_['p4_h2_preReg']=p4_h2_preReg
            quad_['p4_h1']=p4_h1
            quad_['p4_h2']=p4_h2
            quad_['m1_preReg']=p4_h1_preReg.M()
            quad_['m2_preReg']=p4_h2_preReg.M()
            quad_['m1']=p4_h1.M()
            quad_['m2']=p4_h2.M()
            quad_['mass']=(p4_h1+p4_h2).M()
            quad_['mass_preReg']=(p4_h1_preReg+p4_h2_preReg).M()
            quad_['pT']=(p4_h1+p4_h2).Pt()
            quad_['y']=(p4_h1+p4_h2).Rapidity()
            quad_['eta']=(p4_h1+p4_h2).Eta()
            quad_['phi']=(p4_h1+p4_h2).Phi()
            quad_['r_HH']=np.sqrt( (p4_h1.M() - 125.0)*(p4_h1.M() - 125.0) + 
                                   (p4_h2.M() - 125.0)*(p4_h2.M() - 125.0))
            quad_['D_HH']=np.abs( p4_h1.M() - a*p4_h2.M() )/b     
     #       print("!! D_HH for combi : ",combi," : ",quad_['D_HH']," r_HH : ", quad_['r_HH'])
            #if quad_['r_HH'] < metric_min:
            if quad_['D_HH'] < metric_min:
                metric_min_idx=len(result['allBJetQuads'])
                metric_min=quad_['D_HH']
            result['allBJetQuads'].append(quad_)
            result['isValid']=True
            result['nJetQuads']+=1
    if doMisclassificationCorrection:
        idx2ndMin=-1
        metri2ndMin=1e10
        for idx in range(len(result['allBJetQuads'])):
            if idx==metric_min_idx:
                continue
            if result['allBJetQuads'][idx]['D_HH'] < metri2ndMin:
                metri2ndMin=result['allBJetQuads'][idx]['D_HH']
                idx2ndMin=idx
        if( idx2ndMin < 0 ):
            print("Problem !! no second quad !! setting idx2ndMin = metric_min_idx")
            idx2ndMin=metric_min_idx
        print("idx2ndMin, metric_min_idx",idx2ndMin,metric_min_idx)
        if abs(result['allBJetQuads'][idx2ndMin]['D_HH']-result['allBJetQuads'][metric_min_idx]['D_HH']) <30.0:
            if result['allBJetQuads'][idx2ndMin]['pT'] > result['allBJetQuads'][metric_min_idx]['pT']:
                metric_min_idx=idx2ndMin
    
    result['bJetQuad']=result['allBJetQuads'][metric_min_idx]
    return result           

def getBJetParis_wrapper(eTree,methord='default',**kwargs):

    if methord.lower()=='default':
        return getBJetParisFGG(eTree,**kwargs)

    if methord.lower()=='mha':
        return getBJetParisFGG_MHA(eTree,**kwargs)



def getBJetParisFGG_MHA( eTree, mask=None, mlScoreTag='', threshold=-1e3, doMisclassificationCorrection=True,returnOnlyJets=False):
    result={'isValid':False , 'jetsValid' : False,
            'allBJetQuads':[],'bJetQuad':[],
            'allJetsSelected':[],
            'nJetQuads':0}
    X0OverY0=1.05
    jets_pt=[]
    jets_eta=[]
    jets_phi=[]
    jets_mass=[]
    jets_deepJetScore=[]
    jets_mlScore=[]
    jets_fgg_index=[]
    N_JET_MAX=8
    nJets=0
    if len(mask)==0:
        mask=np.zeros(N_JET_MAX , dtype=bool)
    for i in range(N_JET_MAX):

        if (getattr(eTree,'jet_'+str(i)+'_isValid') < 0.25):
            continue
        if not mask[i]:
            continue

        jets_pt.append(getattr(eTree,'jet_'+str(i)+'_pt'))
        jets_eta.append(getattr(eTree,'jet_'+str(i)+'_eta'))
        jets_phi.append(getattr(eTree,'jet_'+str(i)+'_phi'))
        jets_mass.append(getattr(eTree,'jet_'+str(i)+'_mass'))
        jets_deepJetScore.append( getattr(eTree,'jet_'+str( i )+'_deepCSVScore') )
        jets_mlScore.append(getattr(eTree,'jet_'+str(i)+mlScoreTag+'_score'))
        jets_fgg_index.append( i  )
        nJets+=1
    #jets_pt =np.array(jets_pt)
    #jets_eta=np.array(jets_eta)
    #jets_phi=np.array(jets_phi)
    
    if(nJets < 4 ) : 
        result['fail']='nValidJets'
        return result 

    photons_eta=[ eTree.leadingPhoton_eta ,eTree.subleadingPhoton_eta  ]
    photons_phi=[ eTree.leadingPhoton_eta ,eTree.subleadingPhoton_phi  ]

    nJets=len(jets_pt)
    a=X0OverY0
    b=np.sqrt(1.0+X0OverY0*X0OverY0 )
    
    # Add the masks for jet selection
    jetPts=np.array(jets_pt)
    jetEta=np.array(jets_eta)
    jet_deepJetScore=np.array(jets_deepJetScore)
    jets_mlScore=np.array(jets_mlScore)
    jetScoreOrder=np.argsort(jets_mlScore*-1.0)
    maskToProcess=np.ones(len(jetPts),dtype=bool)

#    print("\t Score of jets under consideratio n : ",jets_mlScore)
#    print("\t Sorted ordering of jets under consideratio n : ",jetScoreOrder)
#    print("\t Sorted ordering of jets under consideratio n : ",jets_mlScore[jetScoreOrder])

    
    nMax= 4
    for jetIDx in jetScoreOrder[:nMax]:
        result['allJetsSelected'].append( jets_fgg_index[jetIDx] )

    if returnOnlyJets:
        return result

    #print(jets_mlScore)
    #print(maskToProcess)
    #print(jetScoreOrder)

    maskToProcess[jetScoreOrder[nMax:]]=False
    goodJets=np.where(maskToProcess)[0]
    result['jetsValid']=True

    #print(mask)
    #print(jetScoreOrder)
    #print(goodJets)
    allJetCombinations=itrTools.combinations(goodJets,4)

    jetsIdx=np.arange(0,4,1)
    jetLVs=[]
    for i in range(4):
        jetLVs.append(ROOT.TLorentzVector())
    possiblilities_=[[0,1,2,3],[0,2,1,3],[0,3,1,2]]
    
    jetCobinationScores=[]
    metric_min=1e9
    metric_min_idx=-1
#   print("Good jets  : ",goodJets)
    for jetCombination in allJetCombinations:
        # setting 4 jet values
        mlScore=[-2.0 for ll in range(4)]
        for j in range(4):
            ii=jetCombination[j]
            jetLVs[j].SetPtEtaPhiM(jets_pt[ii],jets_eta[ii],
                           jets_phi[ii],jets_mass[ii])
            mlScore[j] = jets_mlScore[ii]
        
        # print("Doing Jet Combination : ",jetCombination)
        # scanning combination in 4 jets

        for combi in possiblilities_:
            quad_={}
            quad_['idxs']=[jetCombination[i] for i in combi]

            if jets_pt[quad_['idxs'][0]] < jets_pt[quad_['idxs'][1]]:
                t=quad_['idxs'][1] ; quad_['idxs'][1]=quad_['idxs'][0] ;quad_['idxs'][0]=t
                t=combi[1] ; combi[1]=combi[0] ;combi[0]=t

            if jets_pt[quad_['idxs'][2]] < jets_pt[quad_['idxs'][3]]:
                t=quad_['idxs'][3] ; quad_['idxs'][3]=quad_['idxs'][2] ;quad_['idxs'][2]=t
                t=combi[3] ; combi[3]=combi[2] ;combi[2]=t
            
            quad_['fgg_idxs']  = [ jets_fgg_index[ idx ] for idx in quad_['idxs'] ]         
            quad_['mlScore'] =np.array([ mlScore[combi[0]] , mlScore[combi[1]] ,mlScore[combi[2]] , mlScore[combi[3]] ])
            quad_['n_goodJets'] = sum(quad_['mlScore'] > threshold)
            correctedLV0=getCorrectedJetP4(jetLVs[combi[0]],1.0)
            correctedLV1=getCorrectedJetP4(jetLVs[combi[1]],1.0)
            correctedLV2=getCorrectedJetP4(jetLVs[combi[2]],1.0)
            correctedLV3=getCorrectedJetP4(jetLVs[combi[3]],1.0)
            
            p4_h1_preReg=jetLVs[combi[0]]+jetLVs[combi[1]]
            p4_h2_preReg=jetLVs[combi[2]]+jetLVs[combi[3]]
            p4_h1=correctedLV0 + correctedLV1
            p4_h2=correctedLV2 + correctedLV3

            if p4_h2.Pt() > p4_h1.Pt():
                t=p4_h1_preReg
                p4_h1_preReg=p4_h2_preReg
                p4_h2_preReg=t

                t=p4_h1
                p4_h1=p4_h2
                p4_h2=t

            quad_['p4_h1_preReg']=p4_h1_preReg
            quad_['p4_h2_preReg']=p4_h2_preReg
            quad_['p4_h1']=p4_h1
            quad_['p4_h2']=p4_h2
            quad_['m1_preReg']=p4_h1_preReg.M()
            quad_['m2_preReg']=p4_h2_preReg.M()
            quad_['m1']=p4_h1.M()
            quad_['m2']=p4_h2.M()
            quad_['mass']=(p4_h1+p4_h2).M()
            quad_['mass_preReg']=(p4_h1_preReg+p4_h2_preReg).M()
            quad_['pT']=(p4_h1+p4_h2).Pt()
            quad_['y']=(p4_h1+p4_h2).Rapidity()
            quad_['eta']=(p4_h1+p4_h2).Eta()
            quad_['phi']=(p4_h1+p4_h2).Phi()
            quad_['r_HH']=np.sqrt( (p4_h1.M() - 125.0)*(p4_h1.M() - 125.0) + 
                                   (p4_h2.M() - 125.0)*(p4_h2.M() - 125.0))
            quad_['D_HH']=np.abs( p4_h1.M() - a*p4_h2.M() )/b     
     #      print("!! D_HH for combi : ",combi," : ",quad_['D_HH']," r_HH : ", quad_['r_HH'])
            #if quad_['r_HH'] < metric_min:
            if quad_['D_HH'] < metric_min:
                metric_min_idx=len(result['allBJetQuads'])
                metric_min=quad_['D_HH']
            result['allBJetQuads'].append(quad_)
            result['isValid']=True
            result['nJetQuads']+=1
    if doMisclassificationCorrection:
        idx2ndMin=-1
        metri2ndMin=1e10
        for idx in range(len(result['allBJetQuads'])):
            if idx==metric_min_idx:
                continue
            if result['allBJetQuads'][idx]['D_HH'] < metri2ndMin:
                metri2ndMin=result['allBJetQuads'][idx]['D_HH']
                idx2ndMin=idx
        if( idx2ndMin < 0 ):
            print("Problem !! no second quad !! setting idx2ndMin = metric_min_idx")
            idx2ndMin=metric_min_idx
        if abs(result['allBJetQuads'][idx2ndMin]['D_HH']-result['allBJetQuads'][metric_min_idx]['D_HH']) <30.0:
            if result['allBJetQuads'][idx2ndMin]['pT'] > result['allBJetQuads'][metric_min_idx]['pT']:
                metric_min_idx=idx2ndMin
    
    result['bJetQuad']=result['allBJetQuads'][metric_min_idx]
    return result           



def getBJetParisFGG_MHA_v0( eTree, etaCut=2.5, pTMin=25.0, mlScoreTag='', threshold=-1e3,nGoodMLJets=0):
    result={'isValid':False , 'allBJetQuads':[],'bJetQuad':[],'nJetQuads':0}
    X0OverY0=1.05
    jets_pt=[]
    jets_eta=[]
    jets_phi=[]
    jets_mass=[]
    jets_deepJetScore=[]
    jets_mlScore=[]
    jets_fgg_index=[]
    N_JET_MAX=8
    
    for i in range(N_JET_MAX):
        if (getattr(eTree,'jet_'+str(i)+'_isValid') < 0.25):
            continue
            
        jets_pt.append(getattr(eTree,'jet_'+str(i)+'_pt'))
        jets_eta.append(getattr(eTree,'jet_'+str(i)+'_eta'))
        jets_phi.append(getattr(eTree,'jet_'+str(i)+'_phi'))
        jets_mass.append(getattr(eTree,'jet_'+str(i)+'_mass'))
        jets_deepJetScore.append( getattr(eTree,'jet_'+str( i )+'_deepCSVScore') )
        jets_mlScore.append(getattr(eTree,'jet_'+str(i)+mlScoreTag+'_score'))
        jets_fgg_index.append( i  )
    
    #jets_pt =np.array(jets_pt)
    #jets_eta=np.array(jets_eta)
    #jets_phi=np.array(jets_phi)
    photons_eta=[ eTree.leadingPhoton_eta ,eTree.subleadingPhoton_eta  ]
    photons_phi=[ eTree.leadingPhoton_eta ,eTree.subleadingPhoton_phi  ]
    nJets=len(jets_pt)
    a=X0OverY0
    b=np.sqrt(1.0+X0OverY0*X0OverY0)
    if(nJets < 4 ) : 
        result['fail']='nValidJets'
        return result 
    if(jets_pt[3] < pTMin ) :
        result['fail']='pTMin'
        return result 
    
    
    # Add the masks for jet selection
    jetPts=np.array(jets_pt)
    mask= jetPts>pTMin
    n=sum(mask)
    if n<4:
        result['fail']='nJetsPt'
        return result
    
    jetEta=np.array(jets_eta)
    mask=np.logical_and(mask , abs(jetEta) < etaCut )
    n=sum(mask)
    if n<4:
        result['fail']='nJetsEta'
        return result
    
    #jetPhi=np.array(jets_phi)
    #dR=deltaR(jetEta,jetPhi,photons_eta[0],photons_phi[0])
    #mask=np.logical_and(  mask , dR > -0.4 )
    #dr=deltaR(jetEta,jetPhi,photons_eta[1],photons_phi[1])
    #mask=np.logical_and(  mask , dR > -0.4 )
    
    jet_deepJetScore=np.array(jets_deepJetScore)
    jet_deepJetScore[np.logical_not(mask)] = -1e3

    jets_mlScore=np.array(jets_mlScore)
    #print()
    #print("pre mask :",mask)
    #print("jet scores ",jets_mlScore)
    #print("post mask :",mask)
    m=jets_mlScore > threshold
    n=sum(m)
    if n<nGoodMLJets:
        result['fail']='nJetsMLScore'
        return result
    
    jetScoreOrder=np.argsort(jets_mlScore*-1.0)

    #print("\t Score of jets under consideratio n : ",jet_deepJetScore)
    #print("\t Sorted ordering of jets under consideratio n : ",jetScoreOrder)
    #print("\t Sorted ordering of jets under consideratio n : ",jet_deepJetScore[jetScoreOrder])

    nMax= 4
    mask[jetScoreOrder[nMax:]]=False
    
    goodJets=np.where(mask)[0]
    if(len(goodJets) < 3 ): 
        result['fail']='goodJets'
        return result
    allJetCombinations=itrTools.combinations(goodJets,4)
    jetsIdx=np.arange(0,4,1)
    jetLVs=[]
    for i in range(4):
        jetLVs.append(ROOT.TLorentzVector())
    possiblilities_=[[0,1,2,3],[0,2,1,3],[0,3,1,2]]
    
    jetCobinationScores=[]
    metric_min=1e9
    metric_min_idx=-1
#     print("Good jets  : ",goodJets)
    for jetCombination in allJetCombinations:
        # setting 4 jet values
        mlScore=[-2.0 for ll in range(4)]
        for j in range(4):
            ii=jetCombination[j]
            jetLVs[j].SetPtEtaPhiM(jets_pt[ii],jets_eta[ii],
                           jets_phi[ii],jets_mass[ii])
            mlScore[j] = jets_mlScore[ii]

#         print("Doing Jet Combination : ",jetCombination)
        # scanning combination in 4 jets
        for combi in possiblilities_:
            quad_={}
            quad_['idxs']=[jetCombination[i] for i in combi]

            if jets_pt[quad_['idxs'][0]] < jets_pt[quad_['idxs'][1]]:
                t=quad_['idxs'][1] ; quad_['idxs'][1]=quad_['idxs'][0] ;quad_['idxs'][0]=t
                t=combi[1] ; combi[1]=combi[0] ;combi[0]=t

            if jets_pt[quad_['idxs'][2]] < jets_pt[quad_['idxs'][3]]:
                t=quad_['idxs'][3] ; quad_['idxs'][3]=quad_['idxs'][2] ;quad_['idxs'][2]=t
                t=combi[3] ; combi[3]=combi[2] ;combi[2]=t
            quad_['fgg_idxs']  = [ jets_fgg_index[ idx ] for idx in quad_['idxs'] ]         
            quad_['mlScore'] =np.array([ mlScore[combi[0]] , mlScore[combi[1]] ,mlScore[combi[2]] , mlScore[combi[3]] ])
            quad_['n_goodJets'] = sum(quad_['mlScore'] > threshold)
            correctedLV0=getCorrectedJetP4(jetLVs[combi[0]],getattr(eTree,'jets_'+quad_['idxs'][0]+'_bJetRegCorr'))
            correctedLV1=getCorrectedJetP4(jetLVs[combi[1]],getattr(eTree,'jets_'+quad_['idxs'][1]+'_bJetRegCorr'))
            correctedLV2=getCorrectedJetP4(jetLVs[combi[2]],getattr(eTree,'jets_'+quad_['idxs'][2]+'_bJetRegCorr'))
            correctedLV3=getCorrectedJetP4(jetLVs[combi[3]],getattr(eTree,'jets_'+quad_['idxs'][3]+'_bJetRegCorr'))
            
            p4_h1_preReg=jetLVs[combi[0]]+jetLVs[combi[1]]
            p4_h2_preReg=jetLVs[combi[2]]+jetLVs[combi[3]]
            p4_h1=correctedLV0 + correctedLV1
            p4_h2=correctedLV2 + correctedLV3

            if p4_h2.Pt() > p4_h1.Pt():
                t=p4_h1_preReg
                p4_h1_preReg=p4_h2_preReg
                p4_h2_preReg=t

                t=p4_h1
                p4_h1=p4_h2
                p4_h2=t

            quad_['p4_h1_preReg']=p4_h1_preReg
            quad_['p4_h2_preReg']=p4_h2_preReg
            quad_['p4_h1']=p4_h1
            quad_['p4_h2']=p4_h2
            quad_['m1_preReg']=p4_h1_preReg.M()
            quad_['m2_preReg']=p4_h2_preReg.M()
            quad_['m1']=p4_h1.M()
            quad_['m2']=p4_h2.M()
            quad_['mass']=(p4_h1+p4_h2).M()
            quad_['mass_preReg']=(p4_h1_preReg+p4_h2_preReg).M()
            quad_['pT']=(p4_h1+p4_h2).Pt()
            quad_['y']=(p4_h1+p4_h2).Rapidity()
            quad_['eta']=(p4_h1+p4_h2).Eta()
            quad_['phi']=(p4_h1+p4_h2).Phi()
            quad_['r_HH']=np.sqrt( (p4_h1.M() - 125.0)*(p4_h1.M() - 125.0) + 
                                   (p4_h2.M() - 125.0)*(p4_h2.M() - 125.0))
            quad_['D_HH']=np.abs( p4_h1.M() - a*p4_h2.M() )/b     
     #      print("!! D_HH for combi : ",combi," : ",quad_['D_HH']," r_HH : ", quad_['r_HH'])
            #if quad_['r_HH'] < metric_min:
            if quad_['D_HH'] < metric_min:
                metric_min_idx=len(result['allBJetQuads'])
                metric_min=quad_['D_HH']
            result['allBJetQuads'].append(quad_)
            result['isValid']=True
            result['nJetQuads']+=1
    if True:
        idx2ndMin=-1
        metri2ndMin=1e10
        for idx in range(len(result['allBJetQuads'])):
            if idx==metric_min_idx:
                continue
            if result['allBJetQuads'][idx]['D_HH'] < metri2ndMin:
                metri2ndMin=result['allBJetQuads'][idx]['D_HH']
                idx2ndMin=idx
        if( idx2ndMin < 0 ):
            print("Problem !! no second quad !! setting idx2ndMin = metric_min_idx")
            idx2ndMin=metric_min_idx
        if abs(result['allBJetQuads'][idx2ndMin]['D_HH']-result['allBJetQuads'][metric_min_idx]['D_HH']) <30.0:
            if result['allBJetQuads'][idx2ndMin]['pT'] > result['allBJetQuads'][metric_min_idx]['pT']:
                metric_min_idx=idx2ndMin
    
    result['bJetQuad']=result['allBJetQuads'][metric_min_idx]
    return result           

def getBJetParis(eTree,HggCandidateIndexes):
    result={'valid':False , 'allBJetQuads':[],'bJetQuad':[],'nJetQuads':0}
    X0OverY0=1.05
    
    a=X0OverY0
    b=np.sqrt(1.0+X0OverY0*X0OverY0)
#     print("\t\t nj :",eTree.nJets)
    if(eTree.nJets < 4 ) : return result 
#     print("\t\t pt[3] :",eTree.jets_pt[3])
    
    if(eTree.jets_pt[3] < 25.0 ) : return result 
    
    # Add the masks for jet selection
    jetPts=np.array(eTree.jets_pt)
    mask= jetPts>25
    jetEta=np.array(eTree.jets_eta)
    mask=np.logical_and(mask , abs(jetEta) < 2.5 )
    
    jetPhi=np.array(eTree.jets_phi)
    dR=deltaR(jetEta,jetPhi,eTree.photons_eta[HggCandidateIndexes[0]],eTree.photons_phi[HggCandidateIndexes[0]])
    mask=np.logical_and(mask , dR > 0.4 )
    dr=deltaR(jetEta,jetPhi,eTree.photons_eta[HggCandidateIndexes[1]],eTree.photons_phi[HggCandidateIndexes[1]])
    mask=np.logical_and(mask , dR > 0.4 )
    
    jet_deepJetScore=np.array(eTree.jets_deepJetScore)
    jet_deepJetScore[np.logical_not(mask)] = -1e3
    jetScoreOrder=np.argsort(jet_deepJetScore*-1.0)

    #print("\t Score of jets under consideratio n : ",jet_deepJetScore)
    #print("\t Sorted ordering of jets under consideratio n : ",jetScoreOrder)
    #print("\t Sorted ordering of jets under consideratio n : ",jet_deepJetScore[jetScoreOrder])
    nMax= 4
    mask[jetScoreOrder[nMax:]]=False
    
    goodJets=np.where(mask)[0]
    if(len(goodJets) < 4 ): return result
    allJetCombinations=itrTools.combinations(goodJets,4)
    jetsIdx=np.arange(0,4,1)
    jetLVs=[]
    for i in range(4):
        jetLVs.append(ROOT.TLorentzVector())
    possiblilities_=[[0,1,2,3],[0,2,1,3],[0,3,1,2]]
    
    jetCobinationScores=[]
    metric_min=1e9
    metric_min_idx=-1
#     print("Good jets  : ",goodJets)
    for jetCombination in allJetCombinations:
        # setting 4 jet values
        for j in range(4):
            ii=jetCombination[j]
            jetLVs[j].SetPtEtaPhiM(eTree.jets_pt[ii],eTree.jets_eta[ii],
                            eTree.jets_phi[ii],eTree.jets_mass[ii])
#         print("Doing Jet Combination : ",jetCombination)
        # scanning combination in 4 jets
        for combi in possiblilities_:
            quad_={}
            quad_['idxs']=[jetCombination[i] for i in combi]

            if eTree.jets_pt[quad_['idxs'][0]] < eTree.jets_pt[quad_['idxs'][1]]:
                t=quad_['idxs'][1] ; quad_['idxs'][1]=quad_['idxs'][0] ;quad_['idxs'][0]=t
                t=combi[1] ; combi[1]=combi[0] ;combi[0]=t

            if eTree.jets_pt[quad_['idxs'][2]] < eTree.jets_pt[quad_['idxs'][3]]:
                t=quad_['idxs'][3] ; quad_['idxs'][3]=quad_['idxs'][2] ;quad_['idxs'][2]=t
                t=combi[3] ; combi[3]=combi[2] ;combi[2]=t
            
            correctedLV0=getCorrectedJetP4(jetLVs[combi[0]],eTree.jets_bJetRegCorr[quad_['idxs'][0]])
            correctedLV1=getCorrectedJetP4(jetLVs[combi[1]],eTree.jets_bJetRegCorr[quad_['idxs'][1]])
            correctedLV2=getCorrectedJetP4(jetLVs[combi[2]],eTree.jets_bJetRegCorr[quad_['idxs'][2]])
            correctedLV3=getCorrectedJetP4(jetLVs[combi[3]],eTree.jets_bJetRegCorr[quad_['idxs'][3]])
            
            p4_h1_preReg=jetLVs[combi[0]]+jetLVs[combi[1]]
            p4_h2_preReg=jetLVs[combi[2]]+jetLVs[combi[3]]
            p4_h1=correctedLV0 + correctedLV1
            p4_h2=correctedLV2 + correctedLV3

            if p4_h2.Pt() > p4_h1.Pt():
                t=p4_h1_preReg
                p4_h1_preReg=p4_h2_preReg
                p4_h2_preReg=t

                t=p4_h1
                p4_h1=p4_h2
                p4_h2=t

            quad_['p4_h1_preReg']=p4_h1_preReg
            quad_['p4_h2_preReg']=p4_h2_preReg
            quad_['p4_h1']=p4_h1
            quad_['p4_h2']=p4_h2
            quad_['m1_preReg']=p4_h1_preReg.M()
            quad_['m2_preReg']=p4_h2_preReg.M()
            quad_['m1']=p4_h1.M()
            quad_['m2']=p4_h2.M()
            quad_['mass']=(p4_h1+p4_h2).M()
            quad_['mass_preReg']=(p4_h1_preReg+p4_h2_preReg).M()
            quad_['pT']=(p4_h1+p4_h2).Pt()
            quad_['y']=(p4_h1+p4_h2).Rapidity()
            quad_['eta']=(p4_h1+p4_h2).Eta()
            quad_['phi']=(p4_h1+p4_h2).Phi()
            quad_['r_HH']=np.sqrt( (p4_h1.M() - 125.0)*(p4_h1.M() - 125.0) + 
                                   (p4_h2.M() - 125.0)*(p4_h2.M() - 125.0))
            quad_['D_HH']=np.abs( p4_h1.M() - a*p4_h2.M() )/b     
     #       print("!! D_HH for combi : ",combi," : ",quad_['D_HH']," r_HH : ", quad_['r_HH'])
            #if quad_['r_HH'] < metric_min:
            if quad_['D_HH'] < metric_min:
                metric_min_idx=len(result['allBJetQuads'])
                metric_min=quad_['D_HH']
            result['allBJetQuads'].append(quad_)
            result['valid']=True
            result['nJetQuads']+=1
    if True:
        idx2ndMin=-1
        metri2ndMin=1e10
        for idx in range(len(result['allBJetQuads'])):
            if idx==metric_min_idx:
                continue
            if result['allBJetQuads'][idx]['D_HH'] < metri2ndMin:
                metri2ndMin=result['allBJetQuads'][idx]['D_HH']
                idx2ndMin=idx
        if( idx2ndMin < 0 ):
            print("Problem !! no second quad !! setting idx2ndMin = metric_min_idx")
            idx2ndMin=metric_min_idx
        if abs(result['allBJetQuads'][idx2ndMin]['D_HH']-result['allBJetQuads'][metric_min_idx]['D_HH']) <30.0:
            if result['allBJetQuads'][idx2ndMin]['pT'] > result['allBJetQuads'][metric_min_idx]['pT']:
                metric_min_idx=idx2ndMin
    
    result['bJetQuad']=result['allBJetQuads'][metric_min_idx]
    return result           


def fillTrippleHVariables(eTree,histStore,quad,diPhoton):
    
    g1,g2=diPhoton['index']
    j1,j2,k1,k2=quad['idxs']
    LVStore={}
    LVStore['j1LV']=ROOT.TLorentzVector();
    LVStore['j2LV']=ROOT.TLorentzVector();
    LVStore['k1LV']=ROOT.TLorentzVector();
    LVStore['k2LV']=ROOT.TLorentzVector();
    
    LVStore['g1LV']=ROOT.TLorentzVector()
    LVStore['g1LV'].SetPtEtaPhiE(eTree.photons_pt[g1],eTree.photons_eta[g1],eTree.photons_phi[g1],eTree.photons_e[g1])

    LVStore['g2LV']=ROOT.TLorentzVector()
    LVStore['g2LV'].SetPtEtaPhiE(eTree.photons_pt[g2],eTree.photons_eta[g2],eTree.photons_phi[g2],eTree.photons_e[g2])

    LVStore['j1LV'].SetPtEtaPhiM(eTree.jets_pt[j1],eTree.jets_eta[j1],eTree.jets_phi[j1],eTree.jets_mass[j1]);
    LVStore['j1LV']=getCorrectedJetP4(LVStore['j1LV'],eTree.jets_bJetRegCorr[j1])
    LVStore['j2LV'].SetPtEtaPhiM(eTree.jets_pt[j2],eTree.jets_eta[j2],eTree.jets_phi[j2],eTree.jets_mass[j2]);
    LVStore['j2LV']=getCorrectedJetP4(LVStore['j2LV'],eTree.jets_bJetRegCorr[j2])
    LVStore['k1LV'].SetPtEtaPhiM(eTree.jets_pt[k1],eTree.jets_eta[k1],eTree.jets_phi[k1],eTree.jets_mass[k1]);
    LVStore['k1LV']=getCorrectedJetP4(LVStore['k1LV'],eTree.jets_bJetRegCorr[k1])
    LVStore['k2LV'].SetPtEtaPhiM(eTree.jets_pt[k2],eTree.jets_eta[k2],eTree.jets_phi[k2],eTree.jets_mass[k2]);
    LVStore['k2LV']=getCorrectedJetP4(LVStore['k2LV'],eTree.jets_bJetRegCorr[k2])
    
    LVStore['HggLV']=diPhoton['p4']
    LVStore['H1bbLV']=quad['p4_h1']
    LVStore['H2bbLV']=quad['p4_h2']
    
    
    LVStore['H1LV']=LVStore['HggLV']
    LVStore['H2LV']=LVStore['HggLV']
    LVStore['H3LV']=LVStore['HggLV']
    
    if LVStore['HggLV'].Pt() < LVStore['H1bbLV'].Pt():
        LVStore['H1LV']=LVStore['H1bbLV']
        if LVStore['HggLV'].Pt() < LVStore['H2bbLV'].Pt():
            LVStore['H2LV']=LVStore['H2bbLV']
        else:
            LVStore['H3LV']=LVStore['H2bbLV']
    else:
        LVStore['H2LV']=LVStore['H1bbLV']
        LVStore['H3LV']=LVStore['H2bbLV']
    
    LVStore['HHHLV']=LVStore['H1LV']+LVStore['H2LV']+LVStore['H3LV']
    
    hgg_=ROOT.TLorentzVector(LVStore['HggLV']) ;  hgg_.Boost(-1.0*LVStore['HHHLV'].BoostVector())
    histStore['vars_v1']['HHHCosThetaHgg'].Fill( hgg_.CosTheta())
    lv=ROOT.TLorentzVector(LVStore['g1LV']) ;    lv.Boost(-1*LVStore['HggLV'].BoostVector())
    histStore['vars_v1']['HggCosThetaLeadGamma'].Fill(lv.CosTheta())
    lv=ROOT.TLorentzVector(LVStore['j1LV']) ;    lv.Boost(-1*LVStore['H1bbLV'].BoostVector())
    histStore['vars_v1']['H1bbCosThetaLeadJet'].Fill(lv.CosTheta())
    lv=ROOT.TLorentzVector(LVStore['k1LV']) ;    lv.Boost(-1*LVStore['H2bbLV'].BoostVector())
    
    histStore['vars_v1']['H2bbCosThetaLeadJet'].Fill(lv.CosTheta())
    fillKinematicVarsFromLV(LVStore,histStore["kinematicVars"])

def getSelectedJetCollectionMaskPt(eTree ,jetMask=[],pTMin=25.0):
    if jetMask==[]:
        for i in range(8):
            jetMask.append(True)
    for i in range(8):
        if not jetMask[i]:
            continue
        if abs(getattr(eTree,'jet_'+str(i)+'_isValid') )  < 0.5:
            jetMask[i]=False
        if abs(getattr(eTree,'jet_'+str(i)+'_pt') )  < pTMin:
            jetMask[i]=False
    return jetMask

def getSelectedJetCollectionMaskEta(eTree ,jetMask=[],etaMax=2.5):
    if jetMask==[]:
        for i in range(8):
            jetMask.append(True)
    for i in range(8):
        if not jetMask[i]:
            continue
        if abs(getattr(eTree,'jet_'+str(i)+'_isValid') )  < 0.5:
            jetMask[i]=False
        if abs(getattr(eTree,'jet_'+str(i)+'_eta') ) > etaMax:
            jetMask[i]=False
    return jetMask

def getSelectedJetCollectionMaskBadBJet(eTree ,jetMask=[]):
    if jetMask==[]:
        for i in range(8):
            jetMask.append(True)
    for i in range(8):
        if not jetMask[i]:
            continue
        if abs(getattr(eTree,'jet_'+str(i)+'_isValid') )  < 0.5:
            jetMask[i]=False
        if abs(getattr(eTree,'jet_'+str(i)+'_deepJetScore') ) > 1.05:
            jetMask[i]=False
    return jetMask

def getSelectedJetCollectionMaskOverLap(eTree ,jetMask=[] , overlapRemovalDRMax=0.4) :
    if jetMask==[]:
        for i in range(8):
            jetMask.append(True)
    for i in range(8):
        if not jetMask[i]:
            continue
        if abs(getattr(eTree,'jet_'+str(i)+'_isValid') )  < 0.5:
            jetMask[i]=False
        dr1 =  deltaR(getattr(eTree,'jet_'+str(i)+'_eta') , getattr(eTree,'jet_'+str(i)+'_phi') ,eTree.leadingPhoton_eta,eTree.leadingPhoton_phi )
        dr2 =  deltaR(getattr(eTree,'jet_'+str(i)+'_eta') , getattr(eTree,'jet_'+str(i)+'_phi') ,eTree.subleadingPhoton_eta,eTree.subleadingPhoton_phi )
        if dr1  < overlapRemovalDRMax :
            jetMask[i]=False
        if dr2 < overlapRemovalDRMax :
            jetMask[i]=False
    return np.array(jetMask)      

def vetoOverCountings(nBs,dataTag):
    if dataTag=='ggJets':
        if nBs>0:
            return True
    if dataTag=='ggJets1b':
        if nBs>1:
            return True

    return False
