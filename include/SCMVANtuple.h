//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Feb 14 01:31:57 2022 by ROOT version 6.14/09
// from TTree EventTree/Event data
// found on file: bsToJPsiGamma_run2UL2018_0.root
//////////////////////////////////////////////////////////

#ifndef SCMVANtuple_h
#define SCMVANtuple_h
#define SCMVANtuple_cxx

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class SCMVANtuple {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxscMinDrWithGsfElectornSC = 1;
   static constexpr Int_t kMaxscFoundGsfMatch = 1;

   // Declaration of leaf types
   UInt_t          run;
   ULong64_t       event;
   UInt_t          lumis;
   Bool_t          isData;
   Int_t           gen_nBs;
   vector<double>  *gen_Bs_pt;
   vector<double>  *gen_Bs_energy;
   vector<double>  *gen_Bs_eta;
   vector<double>  *gen_Bs_phi;
   vector<double>  *gen_Bs_pz;
   vector<double>  *gen_Bs_pdgId;
   Int_t           gen_nBsMuonM;
   vector<double>  *gen_BsMuonM_pt;
   vector<double>  *gen_BsMuonM_eta;
   vector<double>  *gen_BsMuonM_phi;
   Int_t           gen_nBsMuonP;
   vector<double>  *gen_BsMuonP_pt;
   vector<double>  *gen_BsMuonP_eta;
   vector<double>  *gen_BsMuonP_phi;
   Int_t           gen_nBsPhoton;
   vector<double>  *gen_BsPhoton_pt;
   vector<double>  *gen_BsPhoton_energy;
   vector<double>  *gen_BsPhoton_eta;
   vector<double>  *gen_BsPhoton_phi;
   Int_t           nMC;
   vector<int>     *mcPID;
   vector<int>     *mcStatus;
   vector<float>   *mcVtx_x;
   vector<float>   *mcVtx_y;
   vector<float>   *mcVtx_z;
   vector<float>   *mcPt;
   vector<float>   *mcEta;
   vector<float>   *mcPhi;
   vector<float>   *mcE;
   vector<float>   *mcEt;
   vector<float>   *mcMass;
   Int_t           nSC;
   vector<float>   *scE;
   vector<float>   *scEt;
   vector<float>   *scRawE;
   vector<float>   *scEta;
   vector<float>   *scPhi;
   vector<float>   *scX;
   vector<float>   *scY;
   vector<float>   *scZ;
   vector<float>   *scEtaWidth;
   vector<float>   *scPhiWidth;
   vector<float>   *scRawEt;
   vector<float>   *scMinDrWithGsfElectornSC_;
   vector<bool>    *scFoundGsfMatch_;
   vector<float>   *scE5x5;
   vector<float>   *scE2x2Ratio;
   vector<float>   *scE3x3Ratio;
   vector<float>   *scEMaxRatio;
   vector<float>   *scE2ndRatio;
   vector<float>   *scETopRatio;
   vector<float>   *scERightRatio;
   vector<float>   *scEBottomRatio;
   vector<float>   *scELeftRatio;
   vector<float>   *scE2x5MaxRatio;
   vector<float>   *scE2x5TopRatio;
   vector<float>   *scE2x5RightRatio;
   vector<float>   *scE2x5BottomRatio;
   vector<float>   *scE2x5LeftRatio;
   vector<float>   *scSwissCross;
   vector<float>   *scR9;
   vector<float>   *scSigmaIetaIeta;
   vector<float>   *scSigmaIetaIphi;
   vector<float>   *scSigmaIphiIphi;
   vector<float>   *scFull5x5_e5x5;
   vector<float>   *scFull5x5_e2x2Ratio;
   vector<float>   *scFull5x5_e3x3Ratio;
   vector<float>   *scFull5x5_eMaxRatio;
   vector<float>   *scFull5x5_e2ndRatio;
   vector<float>   *scFull5x5_eTopRatio;
   vector<float>   *scFull5x5_eRightRatio;
   vector<float>   *scFull5x5_eBottomRatio;
   vector<float>   *scFull5x5_eLeftRatio;
   vector<float>   *scFull5x5_e2x5MaxRatio;
   vector<float>   *scFull5x5_e2x5TopRatio;
   vector<float>   *scFull5x5_e2x5RightRatio;
   vector<float>   *scFull5x5_e2x5BottomRatio;
   vector<float>   *scFull5x5_e2x5LeftRatio;
   vector<float>   *scFull5x5_swissCross;
   vector<float>   *scFull5x5_r9;
   vector<float>   *scFull5x5_sigmaIetaIeta;
   vector<float>   *scFull5x5_sigmaIetaIphi;
   vector<float>   *scFull5x5_sigmaIphiIphi;
   vector<float>   *scNHcalRecHitInDIEta5IPhi5;
   vector<float>   *scEFromHcalRecHitInDIEta5IPhi5;
   vector<float>   *scNHcalRecHitInDIEta2IPhi2;
   vector<float>   *scEFromHcalRecHitInDIEta2IPhi2;
   vector<float>   *scPFChIso1;
   vector<float>   *scPFChIso2;
   vector<float>   *scPFChIso3;
   vector<float>   *scPFChIso4;
   vector<float>   *scPFChIso5;
   vector<float>   *scPFPhoIso1;
   vector<float>   *scPFPhoIso2;
   vector<float>   *scPFPhoIso3;
   vector<float>   *scPFPhoIso4;
   vector<float>   *scPFPhoIso5;
   vector<float>   *scPFNeuIso1;
   vector<float>   *scPFNeuIso2;
   vector<float>   *scPFNeuIso3;
   vector<float>   *scPFNeuIso4;
   vector<float>   *scPFNeuIso5;
   Int_t           nPFCandidates;
   Float_t         pf_ecalEnergy[3445];   //[nPFCandidates]
   Float_t         pf_ecalRawEnergy[3445];   //[nPFCandidates]
   Float_t         pf_hcalEnergy[3445];   //[nPFCandidates]
   Float_t         pf_hcalRawEnergy[3445];   //[nPFCandidates]
   Float_t         pf_HoE[3445];   //[nPFCandidates]
   Float_t         pf_mvaIso[3445];   //[nPFCandidates]
   Float_t         pf_vertexX[3445];   //[nPFCandidates]
   Float_t         pf_vertexY[3445];   //[nPFCandidates]
   Float_t         pf_vertexZ[3445];   //[nPFCandidates]
   Float_t         pf_ecalEntryEta[3445];   //[nPFCandidates]
   Float_t         pf_ecalEntryPhi[3445];   //[nPFCandidates]
   Float_t         pf_ecalEntryX[3445];   //[nPFCandidates]
   Float_t         pf_ecalEntryY[3445];   //[nPFCandidates]
   Float_t         pf_ecalEntryZ[3445];   //[nPFCandidates]
   Float_t         pf_eta[3445];   //[nPFCandidates]
   Float_t         pf_phi[3445];   //[nPFCandidates]
   Float_t         pf_pt[3445];   //[nPFCandidates]
   Float_t         pf_id[3445];   //[nPFCandidates]
   Float_t         pf_mass[3445];   //[nPFCandidates]
   Int_t           nHCALClusters;
   Int_t           nHCALClusterEnergyMatrix;
   Float_t         clusterHCAL_energy[278];   //[nHCALClusters]
   Float_t         clusterHCAL_correctedEnergy[278];   //[nHCALClusters]
   Float_t         clusterHCAL_time[278];   //[nHCALClusters]
   Float_t         clusterHCAL_x[278];   //[nHCALClusters]
   Float_t         clusterHCAL_y[278];   //[nHCALClusters]
   Float_t         clusterHCAL_z[278];   //[nHCALClusters]
   Float_t         clusterHCAL_pt[278];   //[nHCALClusters]
   Float_t         clusterHCAL_eta[278];   //[nHCALClusters]
   Float_t         clusterHCAL_phi[278];   //[nHCALClusters]
   Float_t         clusterHCAL_size[278];   //[nHCALClusters]
   Int_t           nECALClusters;
   Int_t           nECALClusterEnergyMatrix;
   Float_t         clusterECAL_energy[627];   //[nECALClusters]
   Float_t         clusterECAL_correctedEnergy[627];   //[nECALClusters]
   Float_t         clusterECAL_time[627];   //[nECALClusters]
   Float_t         clusterECAL_eta[627];   //[nECALClusters]
   Float_t         clusterECAL_phi[627];   //[nECALClusters]
   Float_t         clusterECAL_x[627];   //[nECALClusters]
   Float_t         clusterECAL_y[627];   //[nECALClusters]
   Float_t         clusterECAL_z[627];   //[nECALClusters]
   Float_t         clusterECAL_pt[627];   //[nECALClusters]
   Float_t         clusterECAL_size[627];   //[nECALClusters]
   Float_t         clusterECAL_energyMatrix[15675];   //[nECALClusterEnergyMatrix]
   Int_t           nGeneralTracks;
   Float_t         generalTracks_outer_x[2926];   //[nGeneralTracks]
   Float_t         generalTracks_outer_y[2926];   //[nGeneralTracks]
   Float_t         generalTracks_outer_z[2926];   //[nGeneralTracks]
   Float_t         generalTracks_outer_px[2926];   //[nGeneralTracks]
   Float_t         generalTracks_outer_py[2926];   //[nGeneralTracks]
   Float_t         generalTracks_outer_pz[2926];   //[nGeneralTracks]
   Float_t         generalTracks_charge[2926];   //[nGeneralTracks]
   Int_t           nGenBs;
   Float_t         genBs_nDaus[1];   //[nGenBs]
   Float_t         genBs_pt[1];   //[nGenBs]
   Float_t         genBs_eta[1];   //[nGenBs]
   Float_t         genBs_phi[1];   //[nGenBs]
   Float_t         genBs_mass[1];   //[nGenBs]
   Int_t           nGenBsJpsi;
   Float_t         genBsJpsi_pt[1746];   //[nGenBsJpsi]
   Float_t         genBsJpsi_eta[1746];   //[nGenBsJpsi]
   Float_t         genBsJpsi_phi[1746];   //[nGenBsJpsi]
   Float_t         genBsJpsi_mass[1746];   //[nGenBsJpsi]
   Int_t           nGenBsMuP;
   Float_t         genBsMuP_pt[1];   //[nGenBsMuP]
   Float_t         genBsMuP_eta[1];   //[nGenBsMuP]
   Float_t         genBsMuP_phi[1];   //[nGenBsMuP]
   Float_t         genBsMuP_mass[1];   //[nGenBsMuP]
   Int_t           nGenBsMuM;
   Float_t         genBsMuM_pt[1];   //[nGenBsMuM]
   Float_t         genBsMuM_eta[1];   //[nGenBsMuM]
   Float_t         genBsMuM_phi[1];   //[nGenBsMuM]
   Float_t         genBsMuM_mass[1];   //[nGenBsMuM]
   Int_t           nGenBsPho;
   Float_t         genBsPho_pt[1];   //[nGenBsPho]
   Float_t         genBsPho_eta[1];   //[nGenBsPho]
   Float_t         genBsPho_phi[1];   //[nGenBsPho]
   Float_t         genBsPho_mass[1];   //[nGenBsPho]
   Int_t           nGenBsDaus;
   Float_t         genBsDaus_pt[2];   //[nGenBsDaus]
   Float_t         genBsDaus_eta[2];   //[nGenBsDaus]
   Float_t         genBsDaus_phi[2];   //[nGenBsDaus]
   Float_t         genBsDaus_mass[2];   //[nGenBsDaus]
   Float_t         genBsDaus_pdgid[2];   //[nGenBsDaus]
   Int_t           nPrimaryVertex;
   Float_t         primaryVertex_isFake[61];   //[nPrimaryVertex]
   Float_t         primaryVertex_x[61];   //[nPrimaryVertex]
   Float_t         primaryVertex_y[61];   //[nPrimaryVertex]
   Float_t         primaryVertex_z[61];   //[nPrimaryVertex]
   Float_t         primaryVertex_t[61];   //[nPrimaryVertex]
   Float_t         primaryVertex_covXX[61];   //[nPrimaryVertex]
   Float_t         primaryVertex_covXY[61];   //[nPrimaryVertex]
   Float_t         primaryVertex_covXZ[61];   //[nPrimaryVertex]
   Float_t         primaryVertex_covYY[61];   //[nPrimaryVertex]
   Float_t         primaryVertex_covYZ[61];   //[nPrimaryVertex]
   Float_t         primaryVertex_covZZ[61];   //[nPrimaryVertex]
   Float_t         primaryVertex_x_error[61];   //[nPrimaryVertex]
   Float_t         primaryVertex_y_error[61];   //[nPrimaryVertex]
   Float_t         primaryVertex_z_error[61];   //[nPrimaryVertex]
   Float_t         primaryVertex_t_error[61];   //[nPrimaryVertex]
   Float_t         primaryVertex_ntracks[61];   //[nPrimaryVertex]
   Float_t         primaryVertex_ndof[61];   //[nPrimaryVertex]
   Float_t         primaryVertex_chi2[61];   //[nPrimaryVertex]
   Float_t         primaryVertex_normalizedChi2[61];   //[nPrimaryVertex]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_gen_nBs;   //!
   TBranch        *b_gen_Bs_pt;   //!
   TBranch        *b_gen_Bs_energy;   //!
   TBranch        *b_gen_Bs_eta;   //!
   TBranch        *b_gen_Bs_phi;   //!
   TBranch        *b_gen_Bs_pz;   //!
   TBranch        *b_gen_Bs_pdgId;   //!
   TBranch        *b_gen_nBsMuonM;   //!
   TBranch        *b_gen_BsMuonM_pt;   //!
   TBranch        *b_gen_BsMuonM_eta;   //!
   TBranch        *b_gen_BsMuonM_phi;   //!
   TBranch        *b_gen_nBsMuonP;   //!
   TBranch        *b_gen_BsMuonP_pt;   //!
   TBranch        *b_gen_BsMuonP_eta;   //!
   TBranch        *b_gen_BsMuonP_phi;   //!
   TBranch        *b_gen_nBsPhoton;   //!
   TBranch        *b_gen_BsPhoton_pt;   //!
   TBranch        *b_gen_BsPhoton_energy;   //!
   TBranch        *b_gen_BsPhoton_eta;   //!
   TBranch        *b_gen_BsPhoton_phi;   //!
   TBranch        *b_nMC;   //!
   TBranch        *b_mcPID;   //!
   TBranch        *b_mcStatus;   //!
   TBranch        *b_mcVtx_x;   //!
   TBranch        *b_mcVtx_y;   //!
   TBranch        *b_mcVtx_z;   //!
   TBranch        *b_mcPt;   //!
   TBranch        *b_mcEta;   //!
   TBranch        *b_mcPhi;   //!
   TBranch        *b_mcE;   //!
   TBranch        *b_mcEt;   //!
   TBranch        *b_mcMass;   //!
   TBranch        *b_nSC;   //!
   TBranch        *b_scE;   //!
   TBranch        *b_scEt;   //!
   TBranch        *b_scRawE;   //!
   TBranch        *b_scEta;   //!
   TBranch        *b_scPhi;   //!
   TBranch        *b_scX;   //!
   TBranch        *b_scY;   //!
   TBranch        *b_scZ;   //!
   TBranch        *b_scEtaWidth;   //!
   TBranch        *b_scPhiWidth;   //!
   TBranch        *b_scRawEt;   //!
   TBranch        *b_scMinDrWithGsfElectornSC_;   //!
   TBranch        *b_scFoundGsfMatch_;   //!
   TBranch        *b_scE5x5;   //!
   TBranch        *b_scE2x2Ratio;   //!
   TBranch        *b_scE3x3Ratio;   //!
   TBranch        *b_scEMaxRatio;   //!
   TBranch        *b_scE2ndRatio;   //!
   TBranch        *b_scETopRatio;   //!
   TBranch        *b_scERightRatio;   //!
   TBranch        *b_scEBottomRatio;   //!
   TBranch        *b_scELeftRatio;   //!
   TBranch        *b_scE2x5MaxRatio;   //!
   TBranch        *b_scE2x5TopRatio;   //!
   TBranch        *b_scE2x5RightRatio;   //!
   TBranch        *b_scE2x5BottomRatio;   //!
   TBranch        *b_scE2x5LeftRatio;   //!
   TBranch        *b_scSwissCross;   //!
   TBranch        *b_scR9;   //!
   TBranch        *b_scSigmaIetaIeta;   //!
   TBranch        *b_scSigmaIetaIphi;   //!
   TBranch        *b_scSigmaIphiIphi;   //!
   TBranch        *b_scFull5x5_e5x5;   //!
   TBranch        *b_scFull5x5_e2x2Ratio;   //!
   TBranch        *b_scFull5x5_e3x3Ratio;   //!
   TBranch        *b_scFull5x5_eMaxRatio;   //!
   TBranch        *b_scFull5x5_e2ndRatio;   //!
   TBranch        *b_scFull5x5_eTopRatio;   //!
   TBranch        *b_scFull5x5_eRightRatio;   //!
   TBranch        *b_scFull5x5_eBottomRatio;   //!
   TBranch        *b_scFull5x5_eLeftRatio;   //!
   TBranch        *b_scFull5x5_e2x5MaxRatio;   //!
   TBranch        *b_scFull5x5_e2x5TopRatio;   //!
   TBranch        *b_scFull5x5_e2x5RightRatio;   //!
   TBranch        *b_scFull5x5_e2x5BottomRatio;   //!
   TBranch        *b_scFull5x5_e2x5LeftRatio;   //!
   TBranch        *b_scFull5x5_swissCross;   //!
   TBranch        *b_scFull5x5_r9;   //!
   TBranch        *b_scFull5x5_sigmaIetaIeta;   //!
   TBranch        *b_scFull5x5_sigmaIetaIphi;   //!
   TBranch        *b_scFull5x5_sigmaIphiIphi;   //!
   TBranch        *b_scNHcalRecHitInDIEta5IPhi5;   //!
   TBranch        *b_scEFromHcalRecHitInDIEta5IPhi5;   //!
   TBranch        *b_scNHcalRecHitInDIEta2IPhi2;   //!
   TBranch        *b_scEFromHcalRecHitInDIEta2IPhi2;   //!
   TBranch        *b_scPFChIso1;   //!
   TBranch        *b_scPFChIso2;   //!
   TBranch        *b_scPFChIso3;   //!
   TBranch        *b_scPFChIso4;   //!
   TBranch        *b_scPFChIso5;   //!
   TBranch        *b_scPFPhoIso1;   //!
   TBranch        *b_scPFPhoIso2;   //!
   TBranch        *b_scPFPhoIso3;   //!
   TBranch        *b_scPFPhoIso4;   //!
   TBranch        *b_scPFPhoIso5;   //!
   TBranch        *b_scPFNeuIso1;   //!
   TBranch        *b_scPFNeuIso2;   //!
   TBranch        *b_scPFNeuIso3;   //!
   TBranch        *b_scPFNeuIso4;   //!
   TBranch        *b_scPFNeuIso5;   //!
   TBranch        *b_nPFCandidates;   //!
   TBranch        *b_pf_ecalEnergy;   //!
   TBranch        *b_pf_ecalRawEnergy;   //!
   TBranch        *b_pf_hcalEnergy;   //!
   TBranch        *b_pf_hcalRawEnergy;   //!
   TBranch        *b_pf_HoE;   //!
   TBranch        *b_pf_mvaIso;   //!
   TBranch        *b_pf_vertexX;   //!
   TBranch        *b_pf_vertexY;   //!
   TBranch        *b_pf_vertexZ;   //!
   TBranch        *b_pf_ecalEntryEta;   //!
   TBranch        *b_pf_ecalEntryPhi;   //!
   TBranch        *b_pf_ecalEntryX;   //!
   TBranch        *b_pf_ecalEntryY;   //!
   TBranch        *b_pf_ecalEntryZ;   //!
   TBranch        *b_pf_eta;   //!
   TBranch        *b_pf_phi;   //!
   TBranch        *b_pf_pt;   //!
   TBranch        *b_pf_id;   //!
   TBranch        *b_pf_mass;   //!
   TBranch        *b_nHCALClusters;   //!
   TBranch        *b_nHCALClusterEnergyMatrix;   //!
   TBranch        *b_clusterHCAL_energy;   //!
   TBranch        *b_clusterHCAL_correctedEnergy;   //!
   TBranch        *b_clusterHCAL_time;   //!
   TBranch        *b_clusterHCAL_x;   //!
   TBranch        *b_clusterHCAL_y;   //!
   TBranch        *b_clusterHCAL_z;   //!
   TBranch        *b_clusterHCAL_pt;   //!
   TBranch        *b_clusterHCAL_eta;   //!
   TBranch        *b_clusterHCAL_phi;   //!
   TBranch        *b_clusterHCAL_size;   //!
   TBranch        *b_nECALClusters;   //!
   TBranch        *b_nECALClusterEnergyMatrix;   //!
   TBranch        *b_clusterECAL_energy;   //!
   TBranch        *b_clusterECAL_correctedEnergy;   //!
   TBranch        *b_clusterECAL_time;   //!
   TBranch        *b_clusterECAL_eta;   //!
   TBranch        *b_clusterECAL_phi;   //!
   TBranch        *b_clusterECAL_x;   //!
   TBranch        *b_clusterECAL_y;   //!
   TBranch        *b_clusterECAL_z;   //!
   TBranch        *b_clusterECAL_pt;   //!
   TBranch        *b_clusterECAL_size;   //!
   TBranch        *b_clusterECAL_energyMatrix;   //!
   TBranch        *b_nGeneralTracks;   //!
   TBranch        *b_generalTracks_outer_x;   //!
   TBranch        *b_generalTracks_outer_y;   //!
   TBranch        *b_generalTracks_outer_z;   //!
   TBranch        *b_generalTracks_outer_px;   //!
   TBranch        *b_generalTracks_outer_py;   //!
   TBranch        *b_generalTracks_outer_pz;   //!
   TBranch        *b_generalTracks_charge;   //!
   TBranch        *b_nGenBs;   //!
   TBranch        *b_genBs_nDaus;   //!
   TBranch        *b_genBs_pt;   //!
   TBranch        *b_genBs_eta;   //!
   TBranch        *b_genBs_phi;   //!
   TBranch        *b_genBs_mass;   //!
   TBranch        *b_nGenBsJpsi;   //!
   TBranch        *b_genBsJpsi_pt;   //!
   TBranch        *b_genBsJpsi_eta;   //!
   TBranch        *b_genBsJpsi_phi;   //!
   TBranch        *b_genBsJpsi_mass;   //!
   TBranch        *b_nGenBsMuP;   //!
   TBranch        *b_genBsMuP_pt;   //!
   TBranch        *b_genBsMuP_eta;   //!
   TBranch        *b_genBsMuP_phi;   //!
   TBranch        *b_genBsMuP_mass;   //!
   TBranch        *b_nGenBsMuM;   //!
   TBranch        *b_genBsMuM_pt;   //!
   TBranch        *b_genBsMuM_eta;   //!
   TBranch        *b_genBsMuM_phi;   //!
   TBranch        *b_genBsMuM_mass;   //!
   TBranch        *b_nGenBsPho;   //!
   TBranch        *b_genBsPho_pt;   //!
   TBranch        *b_genBsPho_eta;   //!
   TBranch        *b_genBsPho_phi;   //!
   TBranch        *b_genBsPho_mass;   //!
   TBranch        *b_nGenBsDaus;   //!
   TBranch        *b_genBsDaus_pt;   //!
   TBranch        *b_genBsDaus_eta;   //!
   TBranch        *b_genBsDaus_phi;   //!
   TBranch        *b_genBsDaus_mass;   //!
   TBranch        *b_genBsDaus_pdgid;   //!
   TBranch        *b_nPrimaryVertex;   //!
   TBranch        *b_primaryVertex_isFake;   //!
   TBranch        *b_primaryVertex_x;   //!
   TBranch        *b_primaryVertex_y;   //!
   TBranch        *b_primaryVertex_z;   //!
   TBranch        *b_primaryVertex_t;   //!
   TBranch        *b_primaryVertex_covXX;   //!
   TBranch        *b_primaryVertex_covXY;   //!
   TBranch        *b_primaryVertex_covXZ;   //!
   TBranch        *b_primaryVertex_covYY;   //!
   TBranch        *b_primaryVertex_covYZ;   //!
   TBranch        *b_primaryVertex_covZZ;   //!
   TBranch        *b_primaryVertex_x_error;   //!
   TBranch        *b_primaryVertex_y_error;   //!
   TBranch        *b_primaryVertex_z_error;   //!
   TBranch        *b_primaryVertex_t_error;   //!
   TBranch        *b_primaryVertex_ntracks;   //!
   TBranch        *b_primaryVertex_ndof;   //!
   TBranch        *b_primaryVertex_chi2;   //!
   TBranch        *b_primaryVertex_normalizedChi2;   //!

   SCMVANtuple(TTree *tree=0);
   virtual ~SCMVANtuple();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
 //  virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef SCMVANtuple_cxx
SCMVANtuple::SCMVANtuple(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
        std::cout<<"Initializing empty tree !! \n";
   }
   else {
   Init(tree);
   }

}

SCMVANtuple::~SCMVANtuple()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t SCMVANtuple::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t SCMVANtuple::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void SCMVANtuple::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   gen_Bs_pt = 0;
   gen_Bs_energy = 0;
   gen_Bs_eta = 0;
   gen_Bs_phi = 0;
   gen_Bs_pz = 0;
   gen_Bs_pdgId = 0;
   gen_BsMuonM_pt = 0;
   gen_BsMuonM_eta = 0;
   gen_BsMuonM_phi = 0;
   gen_BsMuonP_pt = 0;
   gen_BsMuonP_eta = 0;
   gen_BsMuonP_phi = 0;
   gen_BsPhoton_pt = 0;
   gen_BsPhoton_energy = 0;
   gen_BsPhoton_eta = 0;
   gen_BsPhoton_phi = 0;
   mcPID = 0;
   mcStatus = 0;
   mcVtx_x = 0;
   mcVtx_y = 0;
   mcVtx_z = 0;
   mcPt = 0;
   mcEta = 0;
   mcPhi = 0;
   mcE = 0;
   mcEt = 0;
   mcMass = 0;
   scE = 0;
   scEt = 0;
   scRawE = 0;
   scEta = 0;
   scPhi = 0;
   scX = 0;
   scY = 0;
   scZ = 0;
   scEtaWidth = 0;
   scPhiWidth = 0;
   scRawEt = 0;
   scMinDrWithGsfElectornSC_ = 0;
   scFoundGsfMatch_ = 0;
   scE5x5 = 0;
   scE2x2Ratio = 0;
   scE3x3Ratio = 0;
   scEMaxRatio = 0;
   scE2ndRatio = 0;
   scETopRatio = 0;
   scERightRatio = 0;
   scEBottomRatio = 0;
   scELeftRatio = 0;
   scE2x5MaxRatio = 0;
   scE2x5TopRatio = 0;
   scE2x5RightRatio = 0;
   scE2x5BottomRatio = 0;
   scE2x5LeftRatio = 0;
   scSwissCross = 0;
   scR9 = 0;
   scSigmaIetaIeta = 0;
   scSigmaIetaIphi = 0;
   scSigmaIphiIphi = 0;
   scFull5x5_e5x5 = 0;
   scFull5x5_e2x2Ratio = 0;
   scFull5x5_e3x3Ratio = 0;
   scFull5x5_eMaxRatio = 0;
   scFull5x5_e2ndRatio = 0;
   scFull5x5_eTopRatio = 0;
   scFull5x5_eRightRatio = 0;
   scFull5x5_eBottomRatio = 0;
   scFull5x5_eLeftRatio = 0;
   scFull5x5_e2x5MaxRatio = 0;
   scFull5x5_e2x5TopRatio = 0;
   scFull5x5_e2x5RightRatio = 0;
   scFull5x5_e2x5BottomRatio = 0;
   scFull5x5_e2x5LeftRatio = 0;
   scFull5x5_swissCross = 0;
   scFull5x5_r9 = 0;
   scFull5x5_sigmaIetaIeta = 0;
   scFull5x5_sigmaIetaIphi = 0;
   scFull5x5_sigmaIphiIphi = 0;
   scNHcalRecHitInDIEta5IPhi5 = 0;
   scEFromHcalRecHitInDIEta5IPhi5 = 0;
   scNHcalRecHitInDIEta2IPhi2 = 0;
   scEFromHcalRecHitInDIEta2IPhi2 = 0;
   scPFChIso1 = 0;
   scPFChIso2 = 0;
   scPFChIso3 = 0;
   scPFChIso4 = 0;
   scPFChIso5 = 0;
   scPFPhoIso1 = 0;
   scPFPhoIso2 = 0;
   scPFPhoIso3 = 0;
   scPFPhoIso4 = 0;
   scPFPhoIso5 = 0;
   scPFNeuIso1 = 0;
   scPFNeuIso2 = 0;
   scPFNeuIso3 = 0;
   scPFNeuIso4 = 0;
   scPFNeuIso5 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("gen_nBs", &gen_nBs, &b_gen_nBs);
   fChain->SetBranchAddress("gen_Bs_pt", &gen_Bs_pt, &b_gen_Bs_pt);
   fChain->SetBranchAddress("gen_Bs_energy", &gen_Bs_energy, &b_gen_Bs_energy);
   fChain->SetBranchAddress("gen_Bs_eta", &gen_Bs_eta, &b_gen_Bs_eta);
   fChain->SetBranchAddress("gen_Bs_phi", &gen_Bs_phi, &b_gen_Bs_phi);
   fChain->SetBranchAddress("gen_Bs_pz", &gen_Bs_pz, &b_gen_Bs_pz);
   fChain->SetBranchAddress("gen_Bs_pdgId", &gen_Bs_pdgId, &b_gen_Bs_pdgId);
   fChain->SetBranchAddress("gen_nBsMuonM", &gen_nBsMuonM, &b_gen_nBsMuonM);
   fChain->SetBranchAddress("gen_BsMuonM_pt", &gen_BsMuonM_pt, &b_gen_BsMuonM_pt);
   fChain->SetBranchAddress("gen_BsMuonM_eta", &gen_BsMuonM_eta, &b_gen_BsMuonM_eta);
   fChain->SetBranchAddress("gen_BsMuonM_phi", &gen_BsMuonM_phi, &b_gen_BsMuonM_phi);
   fChain->SetBranchAddress("gen_nBsMuonP", &gen_nBsMuonP, &b_gen_nBsMuonP);
   fChain->SetBranchAddress("gen_BsMuonP_pt", &gen_BsMuonP_pt, &b_gen_BsMuonP_pt);
   fChain->SetBranchAddress("gen_BsMuonP_eta", &gen_BsMuonP_eta, &b_gen_BsMuonP_eta);
   fChain->SetBranchAddress("gen_BsMuonP_phi", &gen_BsMuonP_phi, &b_gen_BsMuonP_phi);
   fChain->SetBranchAddress("gen_nBsPhoton", &gen_nBsPhoton, &b_gen_nBsPhoton);
   fChain->SetBranchAddress("gen_BsPhoton_pt", &gen_BsPhoton_pt, &b_gen_BsPhoton_pt);
   fChain->SetBranchAddress("gen_BsPhoton_energy", &gen_BsPhoton_energy, &b_gen_BsPhoton_energy);
   fChain->SetBranchAddress("gen_BsPhoton_eta", &gen_BsPhoton_eta, &b_gen_BsPhoton_eta);
   fChain->SetBranchAddress("gen_BsPhoton_phi", &gen_BsPhoton_phi, &b_gen_BsPhoton_phi);
   fChain->SetBranchAddress("nMC", &nMC, &b_nMC);
   fChain->SetBranchAddress("mcPID", &mcPID, &b_mcPID);
   fChain->SetBranchAddress("mcStatus", &mcStatus, &b_mcStatus);
   fChain->SetBranchAddress("mcVtx_x", &mcVtx_x, &b_mcVtx_x);
   fChain->SetBranchAddress("mcVtx_y", &mcVtx_y, &b_mcVtx_y);
   fChain->SetBranchAddress("mcVtx_z", &mcVtx_z, &b_mcVtx_z);
   fChain->SetBranchAddress("mcPt", &mcPt, &b_mcPt);
   fChain->SetBranchAddress("mcEta", &mcEta, &b_mcEta);
   fChain->SetBranchAddress("mcPhi", &mcPhi, &b_mcPhi);
   fChain->SetBranchAddress("mcE", &mcE, &b_mcE);
   fChain->SetBranchAddress("mcEt", &mcEt, &b_mcEt);
   fChain->SetBranchAddress("mcMass", &mcMass, &b_mcMass);
   fChain->SetBranchAddress("nSC", &nSC, &b_nSC);
   fChain->SetBranchAddress("scE", &scE, &b_scE);
   fChain->SetBranchAddress("scEt", &scEt, &b_scEt);
   fChain->SetBranchAddress("scRawE", &scRawE, &b_scRawE);
   fChain->SetBranchAddress("scEta", &scEta, &b_scEta);
   fChain->SetBranchAddress("scPhi", &scPhi, &b_scPhi);
   fChain->SetBranchAddress("scX", &scX, &b_scX);
   fChain->SetBranchAddress("scY", &scY, &b_scY);
   fChain->SetBranchAddress("scZ", &scZ, &b_scZ);
   fChain->SetBranchAddress("scEtaWidth", &scEtaWidth, &b_scEtaWidth);
   fChain->SetBranchAddress("scPhiWidth", &scPhiWidth, &b_scPhiWidth);
   fChain->SetBranchAddress("scRawEt", &scRawEt, &b_scRawEt);
   fChain->SetBranchAddress("scMinDrWithGsfElectornSC_", &scMinDrWithGsfElectornSC_, &b_scMinDrWithGsfElectornSC_);
   fChain->SetBranchAddress("scFoundGsfMatch_", &scFoundGsfMatch_, &b_scFoundGsfMatch_);
   fChain->SetBranchAddress("scE5x5", &scE5x5, &b_scE5x5);
   fChain->SetBranchAddress("scE2x2Ratio", &scE2x2Ratio, &b_scE2x2Ratio);
   fChain->SetBranchAddress("scE3x3Ratio", &scE3x3Ratio, &b_scE3x3Ratio);
   fChain->SetBranchAddress("scEMaxRatio", &scEMaxRatio, &b_scEMaxRatio);
   fChain->SetBranchAddress("scE2ndRatio", &scE2ndRatio, &b_scE2ndRatio);
   fChain->SetBranchAddress("scETopRatio", &scETopRatio, &b_scETopRatio);
   fChain->SetBranchAddress("scERightRatio", &scERightRatio, &b_scERightRatio);
   fChain->SetBranchAddress("scEBottomRatio", &scEBottomRatio, &b_scEBottomRatio);
   fChain->SetBranchAddress("scELeftRatio", &scELeftRatio, &b_scELeftRatio);
   fChain->SetBranchAddress("scE2x5MaxRatio", &scE2x5MaxRatio, &b_scE2x5MaxRatio);
   fChain->SetBranchAddress("scE2x5TopRatio", &scE2x5TopRatio, &b_scE2x5TopRatio);
   fChain->SetBranchAddress("scE2x5RightRatio", &scE2x5RightRatio, &b_scE2x5RightRatio);
   fChain->SetBranchAddress("scE2x5BottomRatio", &scE2x5BottomRatio, &b_scE2x5BottomRatio);
   fChain->SetBranchAddress("scE2x5LeftRatio", &scE2x5LeftRatio, &b_scE2x5LeftRatio);
   fChain->SetBranchAddress("scSwissCross", &scSwissCross, &b_scSwissCross);
   fChain->SetBranchAddress("scR9", &scR9, &b_scR9);
   fChain->SetBranchAddress("scSigmaIetaIeta", &scSigmaIetaIeta, &b_scSigmaIetaIeta);
   fChain->SetBranchAddress("scSigmaIetaIphi", &scSigmaIetaIphi, &b_scSigmaIetaIphi);
   fChain->SetBranchAddress("scSigmaIphiIphi", &scSigmaIphiIphi, &b_scSigmaIphiIphi);
   fChain->SetBranchAddress("scFull5x5_e5x5", &scFull5x5_e5x5, &b_scFull5x5_e5x5);
   fChain->SetBranchAddress("scFull5x5_e2x2Ratio", &scFull5x5_e2x2Ratio, &b_scFull5x5_e2x2Ratio);
   fChain->SetBranchAddress("scFull5x5_e3x3Ratio", &scFull5x5_e3x3Ratio, &b_scFull5x5_e3x3Ratio);
   fChain->SetBranchAddress("scFull5x5_eMaxRatio", &scFull5x5_eMaxRatio, &b_scFull5x5_eMaxRatio);
   fChain->SetBranchAddress("scFull5x5_e2ndRatio", &scFull5x5_e2ndRatio, &b_scFull5x5_e2ndRatio);
   fChain->SetBranchAddress("scFull5x5_eTopRatio", &scFull5x5_eTopRatio, &b_scFull5x5_eTopRatio);
   fChain->SetBranchAddress("scFull5x5_eRightRatio", &scFull5x5_eRightRatio, &b_scFull5x5_eRightRatio);
   fChain->SetBranchAddress("scFull5x5_eBottomRatio", &scFull5x5_eBottomRatio, &b_scFull5x5_eBottomRatio);
   fChain->SetBranchAddress("scFull5x5_eLeftRatio", &scFull5x5_eLeftRatio, &b_scFull5x5_eLeftRatio);
   fChain->SetBranchAddress("scFull5x5_e2x5MaxRatio", &scFull5x5_e2x5MaxRatio, &b_scFull5x5_e2x5MaxRatio);
   fChain->SetBranchAddress("scFull5x5_e2x5TopRatio", &scFull5x5_e2x5TopRatio, &b_scFull5x5_e2x5TopRatio);
   fChain->SetBranchAddress("scFull5x5_e2x5RightRatio", &scFull5x5_e2x5RightRatio, &b_scFull5x5_e2x5RightRatio);
   fChain->SetBranchAddress("scFull5x5_e2x5BottomRatio", &scFull5x5_e2x5BottomRatio, &b_scFull5x5_e2x5BottomRatio);
   fChain->SetBranchAddress("scFull5x5_e2x5LeftRatio", &scFull5x5_e2x5LeftRatio, &b_scFull5x5_e2x5LeftRatio);
   fChain->SetBranchAddress("scFull5x5_swissCross", &scFull5x5_swissCross, &b_scFull5x5_swissCross);
   fChain->SetBranchAddress("scFull5x5_r9", &scFull5x5_r9, &b_scFull5x5_r9);
   fChain->SetBranchAddress("scFull5x5_sigmaIetaIeta", &scFull5x5_sigmaIetaIeta, &b_scFull5x5_sigmaIetaIeta);
   fChain->SetBranchAddress("scFull5x5_sigmaIetaIphi", &scFull5x5_sigmaIetaIphi, &b_scFull5x5_sigmaIetaIphi);
   fChain->SetBranchAddress("scFull5x5_sigmaIphiIphi", &scFull5x5_sigmaIphiIphi, &b_scFull5x5_sigmaIphiIphi);
   fChain->SetBranchAddress("scNHcalRecHitInDIEta5IPhi5", &scNHcalRecHitInDIEta5IPhi5, &b_scNHcalRecHitInDIEta5IPhi5);
   fChain->SetBranchAddress("scEFromHcalRecHitInDIEta5IPhi5", &scEFromHcalRecHitInDIEta5IPhi5, &b_scEFromHcalRecHitInDIEta5IPhi5);
   fChain->SetBranchAddress("scNHcalRecHitInDIEta2IPhi2", &scNHcalRecHitInDIEta2IPhi2, &b_scNHcalRecHitInDIEta2IPhi2);
   fChain->SetBranchAddress("scEFromHcalRecHitInDIEta2IPhi2", &scEFromHcalRecHitInDIEta2IPhi2, &b_scEFromHcalRecHitInDIEta2IPhi2);
   fChain->SetBranchAddress("scPFChIso1", &scPFChIso1, &b_scPFChIso1);
   fChain->SetBranchAddress("scPFChIso2", &scPFChIso2, &b_scPFChIso2);
   fChain->SetBranchAddress("scPFChIso3", &scPFChIso3, &b_scPFChIso3);
   fChain->SetBranchAddress("scPFChIso4", &scPFChIso4, &b_scPFChIso4);
   fChain->SetBranchAddress("scPFChIso5", &scPFChIso5, &b_scPFChIso5);
   fChain->SetBranchAddress("scPFPhoIso1", &scPFPhoIso1, &b_scPFPhoIso1);
   fChain->SetBranchAddress("scPFPhoIso2", &scPFPhoIso2, &b_scPFPhoIso2);
   fChain->SetBranchAddress("scPFPhoIso3", &scPFPhoIso3, &b_scPFPhoIso3);
   fChain->SetBranchAddress("scPFPhoIso4", &scPFPhoIso4, &b_scPFPhoIso4);
   fChain->SetBranchAddress("scPFPhoIso5", &scPFPhoIso5, &b_scPFPhoIso5);
   fChain->SetBranchAddress("scPFNeuIso1", &scPFNeuIso1, &b_scPFNeuIso1);
   fChain->SetBranchAddress("scPFNeuIso2", &scPFNeuIso2, &b_scPFNeuIso2);
   fChain->SetBranchAddress("scPFNeuIso3", &scPFNeuIso3, &b_scPFNeuIso3);
   fChain->SetBranchAddress("scPFNeuIso4", &scPFNeuIso4, &b_scPFNeuIso4);
   fChain->SetBranchAddress("scPFNeuIso5", &scPFNeuIso5, &b_scPFNeuIso5);
   fChain->SetBranchAddress("nPFCandidates", &nPFCandidates, &b_nPFCandidates);
   fChain->SetBranchAddress("pf_ecalEnergy", pf_ecalEnergy, &b_pf_ecalEnergy);
   fChain->SetBranchAddress("pf_ecalRawEnergy", pf_ecalRawEnergy, &b_pf_ecalRawEnergy);
   fChain->SetBranchAddress("pf_hcalEnergy", pf_hcalEnergy, &b_pf_hcalEnergy);
   fChain->SetBranchAddress("pf_hcalRawEnergy", pf_hcalRawEnergy, &b_pf_hcalRawEnergy);
   fChain->SetBranchAddress("pf_HoE", pf_HoE, &b_pf_HoE);
   fChain->SetBranchAddress("pf_mvaIso", pf_mvaIso, &b_pf_mvaIso);
   fChain->SetBranchAddress("pf_vertexX", pf_vertexX, &b_pf_vertexX);
   fChain->SetBranchAddress("pf_vertexY", pf_vertexY, &b_pf_vertexY);
   fChain->SetBranchAddress("pf_vertexZ", pf_vertexZ, &b_pf_vertexZ);
   fChain->SetBranchAddress("pf_ecalEntryEta", pf_ecalEntryEta, &b_pf_ecalEntryEta);
   fChain->SetBranchAddress("pf_ecalEntryPhi", pf_ecalEntryPhi, &b_pf_ecalEntryPhi);
   fChain->SetBranchAddress("pf_ecalEntryX", pf_ecalEntryX, &b_pf_ecalEntryX);
   fChain->SetBranchAddress("pf_ecalEntryY", pf_ecalEntryY, &b_pf_ecalEntryY);
   fChain->SetBranchAddress("pf_ecalEntryZ", pf_ecalEntryZ, &b_pf_ecalEntryZ);
   fChain->SetBranchAddress("pf_eta", pf_eta, &b_pf_eta);
   fChain->SetBranchAddress("pf_phi", pf_phi, &b_pf_phi);
   fChain->SetBranchAddress("pf_pt", pf_pt, &b_pf_pt);
   fChain->SetBranchAddress("pf_id", pf_id, &b_pf_id);
   fChain->SetBranchAddress("pf_mass", pf_mass, &b_pf_mass);
   fChain->SetBranchAddress("nHCALClusters", &nHCALClusters, &b_nHCALClusters);
   fChain->SetBranchAddress("nHCALClusterEnergyMatrix", &nHCALClusterEnergyMatrix, &b_nHCALClusterEnergyMatrix);
   fChain->SetBranchAddress("clusterHCAL_energy", clusterHCAL_energy, &b_clusterHCAL_energy);
   fChain->SetBranchAddress("clusterHCAL_correctedEnergy", clusterHCAL_correctedEnergy, &b_clusterHCAL_correctedEnergy);
   fChain->SetBranchAddress("clusterHCAL_time", clusterHCAL_time, &b_clusterHCAL_time);
   fChain->SetBranchAddress("clusterHCAL_x", clusterHCAL_x, &b_clusterHCAL_x);
   fChain->SetBranchAddress("clusterHCAL_y", clusterHCAL_y, &b_clusterHCAL_y);
   fChain->SetBranchAddress("clusterHCAL_z", clusterHCAL_z, &b_clusterHCAL_z);
   fChain->SetBranchAddress("clusterHCAL_pt", clusterHCAL_pt, &b_clusterHCAL_pt);
   fChain->SetBranchAddress("clusterHCAL_eta", clusterHCAL_eta, &b_clusterHCAL_eta);
   fChain->SetBranchAddress("clusterHCAL_phi", clusterHCAL_phi, &b_clusterHCAL_phi);
   fChain->SetBranchAddress("clusterHCAL_size", clusterHCAL_size, &b_clusterHCAL_size);
   fChain->SetBranchAddress("nECALClusters", &nECALClusters, &b_nECALClusters);
   fChain->SetBranchAddress("nECALClusterEnergyMatrix", &nECALClusterEnergyMatrix, &b_nECALClusterEnergyMatrix);
   fChain->SetBranchAddress("clusterECAL_energy", clusterECAL_energy, &b_clusterECAL_energy);
   fChain->SetBranchAddress("clusterECAL_correctedEnergy", clusterECAL_correctedEnergy, &b_clusterECAL_correctedEnergy);
   fChain->SetBranchAddress("clusterECAL_time", clusterECAL_time, &b_clusterECAL_time);
   fChain->SetBranchAddress("clusterECAL_eta", clusterECAL_eta, &b_clusterECAL_eta);
   fChain->SetBranchAddress("clusterECAL_phi", clusterECAL_phi, &b_clusterECAL_phi);
   fChain->SetBranchAddress("clusterECAL_x", clusterECAL_x, &b_clusterECAL_x);
   fChain->SetBranchAddress("clusterECAL_y", clusterECAL_y, &b_clusterECAL_y);
   fChain->SetBranchAddress("clusterECAL_z", clusterECAL_z, &b_clusterECAL_z);
   fChain->SetBranchAddress("clusterECAL_pt", clusterECAL_pt, &b_clusterECAL_pt);
   fChain->SetBranchAddress("clusterECAL_size", clusterECAL_size, &b_clusterECAL_size);
   fChain->SetBranchAddress("clusterECAL_energyMatrix", clusterECAL_energyMatrix, &b_clusterECAL_energyMatrix);
   fChain->SetBranchAddress("nGeneralTracks", &nGeneralTracks, &b_nGeneralTracks);
   fChain->SetBranchAddress("generalTracks_outer_x", generalTracks_outer_x, &b_generalTracks_outer_x);
   fChain->SetBranchAddress("generalTracks_outer_y", generalTracks_outer_y, &b_generalTracks_outer_y);
   fChain->SetBranchAddress("generalTracks_outer_z", generalTracks_outer_z, &b_generalTracks_outer_z);
   fChain->SetBranchAddress("generalTracks_outer_px", generalTracks_outer_px, &b_generalTracks_outer_px);
   fChain->SetBranchAddress("generalTracks_outer_py", generalTracks_outer_py, &b_generalTracks_outer_py);
   fChain->SetBranchAddress("generalTracks_outer_pz", generalTracks_outer_pz, &b_generalTracks_outer_pz);
   fChain->SetBranchAddress("generalTracks_charge", generalTracks_charge, &b_generalTracks_charge);
   fChain->SetBranchAddress("nGenBs", &nGenBs, &b_nGenBs);
   fChain->SetBranchAddress("genBs_nDaus", genBs_nDaus, &b_genBs_nDaus);
   fChain->SetBranchAddress("genBs_pt", genBs_pt, &b_genBs_pt);
   fChain->SetBranchAddress("genBs_eta", genBs_eta, &b_genBs_eta);
   fChain->SetBranchAddress("genBs_phi", genBs_phi, &b_genBs_phi);
   fChain->SetBranchAddress("genBs_mass", genBs_mass, &b_genBs_mass);
   fChain->SetBranchAddress("nGenBsJpsi", &nGenBsJpsi, &b_nGenBsJpsi);
   fChain->SetBranchAddress("genBsJpsi_pt", genBsJpsi_pt, &b_genBsJpsi_pt);
   fChain->SetBranchAddress("genBsJpsi_eta", genBsJpsi_eta, &b_genBsJpsi_eta);
   fChain->SetBranchAddress("genBsJpsi_phi", genBsJpsi_phi, &b_genBsJpsi_phi);
   fChain->SetBranchAddress("genBsJpsi_mass", genBsJpsi_mass, &b_genBsJpsi_mass);
   fChain->SetBranchAddress("nGenBsMuP", &nGenBsMuP, &b_nGenBsMuP);
   fChain->SetBranchAddress("genBsMuP_pt", &genBsMuP_pt, &b_genBsMuP_pt);
   fChain->SetBranchAddress("genBsMuP_eta", &genBsMuP_eta, &b_genBsMuP_eta);
   fChain->SetBranchAddress("genBsMuP_phi", &genBsMuP_phi, &b_genBsMuP_phi);
   fChain->SetBranchAddress("genBsMuP_mass", &genBsMuP_mass, &b_genBsMuP_mass);
   fChain->SetBranchAddress("nGenBsMuM", &nGenBsMuM, &b_nGenBsMuM);
   fChain->SetBranchAddress("genBsMuM_pt", &genBsMuM_pt, &b_genBsMuM_pt);
   fChain->SetBranchAddress("genBsMuM_eta", &genBsMuM_eta, &b_genBsMuM_eta);
   fChain->SetBranchAddress("genBsMuM_phi", &genBsMuM_phi, &b_genBsMuM_phi);
   fChain->SetBranchAddress("genBsMuM_mass", &genBsMuM_mass, &b_genBsMuM_mass);
   fChain->SetBranchAddress("nGenBsPho", &nGenBsPho, &b_nGenBsPho);
   fChain->SetBranchAddress("genBsPho_pt", genBsPho_pt, &b_genBsPho_pt);
   fChain->SetBranchAddress("genBsPho_eta", genBsPho_eta, &b_genBsPho_eta);
   fChain->SetBranchAddress("genBsPho_phi", genBsPho_phi, &b_genBsPho_phi);
   fChain->SetBranchAddress("genBsPho_mass", genBsPho_mass, &b_genBsPho_mass);
   fChain->SetBranchAddress("nGenBsDaus", &nGenBsDaus, &b_nGenBsDaus);
   fChain->SetBranchAddress("genBsDaus_pt", genBsDaus_pt, &b_genBsDaus_pt);
   fChain->SetBranchAddress("genBsDaus_eta", genBsDaus_eta, &b_genBsDaus_eta);
   fChain->SetBranchAddress("genBsDaus_phi", genBsDaus_phi, &b_genBsDaus_phi);
   fChain->SetBranchAddress("genBsDaus_mass", genBsDaus_mass, &b_genBsDaus_mass);
   fChain->SetBranchAddress("genBsDaus_pdgid", genBsDaus_pdgid, &b_genBsDaus_pdgid);
   fChain->SetBranchAddress("nPrimaryVertex", &nPrimaryVertex, &b_nPrimaryVertex);
   fChain->SetBranchAddress("primaryVertex_isFake", primaryVertex_isFake, &b_primaryVertex_isFake);
   fChain->SetBranchAddress("primaryVertex_x", primaryVertex_x, &b_primaryVertex_x);
   fChain->SetBranchAddress("primaryVertex_y", primaryVertex_y, &b_primaryVertex_y);
   fChain->SetBranchAddress("primaryVertex_z", primaryVertex_z, &b_primaryVertex_z);
   fChain->SetBranchAddress("primaryVertex_t", primaryVertex_t, &b_primaryVertex_t);
   fChain->SetBranchAddress("primaryVertex_covXX", primaryVertex_covXX, &b_primaryVertex_covXX);
   fChain->SetBranchAddress("primaryVertex_covXY", primaryVertex_covXY, &b_primaryVertex_covXY);
   fChain->SetBranchAddress("primaryVertex_covXZ", primaryVertex_covXZ, &b_primaryVertex_covXZ);
   fChain->SetBranchAddress("primaryVertex_covYY", primaryVertex_covYY, &b_primaryVertex_covYY);
   fChain->SetBranchAddress("primaryVertex_covYZ", primaryVertex_covYZ, &b_primaryVertex_covYZ);
   fChain->SetBranchAddress("primaryVertex_covZZ", primaryVertex_covZZ, &b_primaryVertex_covZZ);
   fChain->SetBranchAddress("primaryVertex_x_error", primaryVertex_x_error, &b_primaryVertex_x_error);
   fChain->SetBranchAddress("primaryVertex_y_error", primaryVertex_y_error, &b_primaryVertex_y_error);
   fChain->SetBranchAddress("primaryVertex_z_error", primaryVertex_z_error, &b_primaryVertex_z_error);
   fChain->SetBranchAddress("primaryVertex_t_error", primaryVertex_t_error, &b_primaryVertex_t_error);
   fChain->SetBranchAddress("primaryVertex_ntracks", primaryVertex_ntracks, &b_primaryVertex_ntracks);
   fChain->SetBranchAddress("primaryVertex_ndof", primaryVertex_ndof, &b_primaryVertex_ndof);
   fChain->SetBranchAddress("primaryVertex_chi2", primaryVertex_chi2, &b_primaryVertex_chi2);
   fChain->SetBranchAddress("primaryVertex_normalizedChi2", primaryVertex_normalizedChi2, &b_primaryVertex_normalizedChi2);
   Notify();
}

Bool_t SCMVANtuple::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void SCMVANtuple::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t SCMVANtuple::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef SCMVANtuple_cxx
