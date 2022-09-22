import FWCore.ParameterSet.Config as cms
#process = cms.Process("Validation")
process = cms.Process("ntuplizer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v15')
#process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2018_realistic_v15')
process.load("PhysicsTools.PatAlgos.slimming.slimmedAddPileupInfo_cfi")
process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
"file:/afs/cern.ch/work/m/mstamenk/public/HHH/HHH4b2y_RunIISummer20UL17MINIAODSIM.root",

                    ),
                duplicateCheckMode=cms.untracked.string("noDuplicateCheck"),
                            )

########## b Jet Energy Corrections ################
# Modules from NanoAOD
from PhysicsTools.NanoAOD.jets_cff import bJetVars,bjetNN
#### update modules parameters
# create variables for regression
process.bJetVars = bJetVars.clone(
    src = cms.InputTag("slimmedJets"),
)
# add variables to the jet collection
# (using only bJetVars, which are variables relevant to regression)
process.updatedJetsWithUserData = cms.EDProducer("PATJetUserDataEmbedder",
    src = cms.InputTag("slimmedJets"),
    userFloats = cms.PSet(
         leadTrackPt = cms.InputTag("bJetVars:leadTrackPt"),
         leptonPtRel = cms.InputTag("bJetVars:leptonPtRel"),
         leptonPtRatio = cms.InputTag("bJetVars:leptonPtRatio"),
         leptonPtRelInv = cms.InputTag("bJetVars:leptonPtRelInv"),
         leptonPtRelv0 = cms.InputTag("bJetVars:leptonPtRelv0"),
         leptonPtRatiov0 = cms.InputTag("bJetVars:leptonPtRatiov0"),
         leptonPtRelInvv0 = cms.InputTag("bJetVars:leptonPtRelInvv0"),
         leptonDeltaR = cms.InputTag("bJetVars:leptonDeltaR"),
         leptonPt = cms.InputTag("bJetVars:leptonPt"),
         vtxPt = cms.InputTag("bJetVars:vtxPt"),
         vtxMass = cms.InputTag("bJetVars:vtxMass"),
         vtx3dL = cms.InputTag("bJetVars:vtx3dL"),
         vtx3deL = cms.InputTag("bJetVars:vtx3deL"),
         ptD = cms.InputTag("bJetVars:ptD"),
         genPtwNu = cms.InputTag("bJetVars:genPtwNu"),
#        pileupJetIdDiscriminant = cms.InputTag("pileupJetIdUpdated:fullDiscriminant"),

         ),
    userInts = cms.PSet(
        vtxNtrk = cms.InputTag("bJetVars:vtxNtrk"),
        leptonPdgId = cms.InputTag("bJetVars:leptonPdgId"),
#       pileupJetId = cms.InputTag("pileupJetIdUpdated:fullId"),
    ),
)
# obtain the b jet regression corrections
process.bjetNN = bjetNN.clone(
    src = cms.InputTag("updatedJetsWithUserData"),
)

# add the b jet regression corrections to the jet collection
process.slimmedJetsbReg = cms.EDProducer("PATJetUserDataEmbedder",
     src = cms.InputTag("updatedJetsWithUserData"),
     userFloats = cms.PSet(
         bJetRegCorr = cms.InputTag("bjetNN:corr"),
         bJetRegRes = cms.InputTag("bjetNN:res"),
#         pileupJetIdDiscriminant = cms.InputTag("pileupJetIdUpdated:fullDiscriminant"),
         ),
#     userInts = cms.PSet(
#         pileupJetId = cms.InputTag("pileupJetIdUpdated:fullId"),
#         ),
)

process.ntuplizer = cms.EDAnalyzer(
       "genNtuple",
       label_pileUp = cms.untracked.InputTag("slimmedAddPileupInfo"),
       pruned     = cms.InputTag("prunedGenParticles"),
       photons = cms.InputTag("slimmedPhotons"),
       recJet = cms.InputTag("slimmedJetsbReg"),
       genJet = cms.InputTag("slimmedGenJets"),
    )


process.TFileService = cms.Service("TFileService", fileName = cms.string('pu_trees_VH.root'))

process.p = cms.Path(
     process.bJetVars*
     process.updatedJetsWithUserData*
     process.bjetNN*
     process.slimmedJetsbReg*
     process.ntuplizer
    )


process.schedule = cms.Schedule(process.p)
#process.output_step = cms.EndPath(process.out)
#process.schedule.extend([process.output_step])

process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000))


f=open("flist.cfg")
txt=f.readlines()
f.close()
maxEvents=int(txt.pop(0)[:-1])
fOutName=txt.pop(0)[:-1]
process.source.fileNames=cms.untracked.vstring()
for l in txt:
    process.source.fileNames.append(l[:-1])
process.maxEvents.input = maxEvents
process.TFileService.fileName =  fOutName
