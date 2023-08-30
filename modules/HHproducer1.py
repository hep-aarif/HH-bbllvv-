import ROOT
from ROOT import TLorentzVector
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection 
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

import math
import os
import numpy as np
from numpy import sign
from numpy import argsort

class HHproducer1(Module):
  def __init__(self , year):
    self.year = year
  def beginJob(self):
    pass
  def endJob(self):
    pass
  def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    self.out = wrappedOutputTree
    self.out.branch("HLT_passEle32WPTight", "I")
    self.out.branch("lhe_nlepton", "I")
    self.out.branch("nHad_tau", "I")
    self.out.branch("met_user","F")
    self.out.branch("met_phi_user","F")

    self.out.branch("TightFatJet_id","I",lenVar="nTightFatJet")
    self.out.branch("TightFatJet_initid","I",lenVar="nTightFatJet")
    self.out.branch("TightFatJet_bbvsQCD","F",lenVar="nTightFatJet")
    self.out.branch("TightFatJet_ccvsQCD","F",lenVar="nTightFatJet")
    self.out.branch("TightFatJet_qqvsQCD","F",lenVar="nTightFatJet")
    self.out.branch("TightFatJet_drl1","F",lenVar="nTightFatJet")
    self.out.branch("TightFatJet_drl2","F",lenVar="nTightFatJet")
    self.out.branch("TightFatJet_pt","F",lenVar="nTightFatJet")
    self.out.branch("TightFatJet_eta","F",lenVar="nTightFatJet")
    self.out.branch("TightFatJet_phi","F",lenVar="nTightFatJet")
    self.out.branch("TightFatJet_mass","F",lenVar="nTightFatJet")
    self.out.branch("TightFatJet_dRinit","F",lenVar="nTightFatJet")
    self.out.branch("H_1AK8", "B")
    self.out.branch("H_1AK8_case", "I")
    self.out.branch("H_2AK4", "B")
    self.out.branch("H_2AK4_case", "I")
    self.out.branch("H_AK8Jet_id", "I")
    self.out.branch("H_AK8Jet_initid", "I")
    self.out.branch("H_AK8Jet_pt", "F")
    self.out.branch("H_AK8Jet_eta", "F")
    self.out.branch("H_AK8Jet_phi", "F")
    self.out.branch("H_AK8Jet_mass", "F")
    self.out.branch("H_AK8Jet_PNmass", "F")
    self.out.branch("H_AK8Jet_SDmass", "F")
    self.out.branch("H_AK8Jet_drl1", "F")
    self.out.branch("H_AK8Jet_drl2", "F")

    self.out.branch("n_bjet_DeepB_M", "I")
    self.out.branch("n_bjet_DeepB_L", "I")
    self.out.branch("n_tight_nob", "I")
    self.out.branch("HT", "F")
    self.out.branch("TightAK4Jet_pt","F",lenVar="nTightAK4Jet")
    self.out.branch("TightAK4Jet_eta","F",lenVar="nTightAK4Jet")
    self.out.branch("TightAK4Jet_phi","F",lenVar="nTightAK4Jet")
    self.out.branch("TightAK4Jet_mass","F",lenVar="nTightAK4Jet")
    self.out.branch("TightAK4Jet_id","I",lenVar="nTightAK4Jet")
    self.out.branch("TightAK4Jet_nob_id","I",lenVar="nTightAK4Jet")
    self.out.branch("TightAK4Jet_b_DeepCSVmedium_id","I",lenVar="nTightAK4Jet")
    self.out.branch("TightAK4Jet_b_DeepCSVloose_id","I",lenVar="nTightAK4Jet")
    self.out.branch("TightAK4Jet_drl1","F",lenVar="nTightAK4Jet")
    self.out.branch("TightAK4Jet_drl2","F",lenVar="nTightAK4Jet")
    self.out.branch("TightAK4Jet_dRinit","F",lenVar="nTightAK4Jet")
    self.out.branch("Had_tau_id","I",lenVar="nTau")

    self.out.branch("h_j1_pt", "F")
    self.out.branch("h_j1_eta", "F")
    self.out.branch("h_j1_phi", "F")
    self.out.branch("h_j1_mass", "F")
    self.out.branch("h_j1_id", "I")
    self.out.branch("h_j1_drl1", "F")
    self.out.branch("h_j1_drl2", "F")
    self.out.branch("h_j2_pt", "F")
    self.out.branch("h_j2_eta", "F")
    self.out.branch("h_j2_phi", "F")
    self.out.branch("h_j2_mass", "F")
    self.out.branch("h_j2_id", "I")
    self.out.branch("h_j2_drl1", "F")
    self.out.branch("h_j2_drl2", "F")
    self.out.branch("h_mjj", "F")
    self.out.branch("h_detajj", "F")
    self.out.branch("h_dRjj", "F")
    self.out.branch("h_dphijj", "F")

    self.is_mc = bool(inputTree.GetBranch("GenJet_pt"))
    self.is_lhe = bool(inputTree.GetBranch("nLHEPart"))

  def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    pass

  def analyze(self, event):
    
    # PV selection
    if (event.PV_npvsGood<1): return False
    # if abs(event.GEN_zj1_pdgid)==5: return False

    # trigger selection
    # special action for 2017 single ele HLT, https://twiki.cern.ch/twiki/bin/viewauth/CMS/Egamma2017DataRecommendations#Single_Electron_Triggers
    # TrigObj_filterBits, 1024 = 1e (32_L1DoubleEG_AND_L1SingleEGOr)
    HLT_passEle32WPTight=0
    if self.year=="2017":
      trgobjs=Collection(event, 'TrigObj')
      if event.HLT_Ele32_WPTight_Gsf_L1DoubleEG==1:
	for iobj in range(0,event.nTrigObj):
	  if trgobjs[iobj].id==11 and (trgobjs[iobj].filterBits & (1<<10))== (1<<10):
	    HLT_passEle32WPTight=1

    self.out.fillBranch("HLT_passEle32WPTight",HLT_passEle32WPTight)

    lhe_nlepton=0
    if self.is_lhe:
      lheparticle = Collection(event, 'LHEPart')
      for ilhe in range(0, event.nLHEPart):
        if lheparticle[ilhe].status==1 and (abs(lheparticle[ilhe].pdgId)==11 or abs(lheparticle[ilhe].pdgId)==13 or abs(lheparticle[ilhe].pdgId)==15):
          lhe_nlepton=lhe_nlepton+1

    self.out.fillBranch("lhe_nlepton", lhe_nlepton)

    l1v4_tmp=TLorentzVector()
    l2v4_tmp=TLorentzVector()
    l1v4_tmp.SetPtEtaPhiM(event.l1_pt,event.l1_eta,event.l1_phi,event.l1_mass)
    l2v4_tmp.SetPtEtaPhiM(event.l2_pt,event.l2_eta,event.l2_phi,event.l2_mass)

    tau_v4_temp=TLorentzVector()
    taus = Collection(event, 'Tau')
    nHad_tau=0
    Had_tau_id=[]
    for itau in range(0, event.nTau):
      tau_v4_temp.SetPtEtaPhiM(taus[itau].pt, taus[itau].eta, taus[itau].phi, taus[itau].mass)
      pass_tau_lep_Dr=1
      if taus[itau].pt>20 and abs(taus[itau].eta)<2.3 and taus[itau].idDecayModeOldDMs and taus[itau].idDeepTau2017v2p1VSe>=4 and taus[itau].idDeepTau2017v2p1VSjet>=4 and taus[itau].idDeepTau2017v2p1VSmu>=1:
        if tau_v4_temp.DeltaR(l1v4_tmp)<0.4:pass_tau_lep_Dr=0
        if tau_v4_temp.DeltaR(l2v4_tmp)<0.4:pass_tau_lep_Dr=0
	if pass_tau_lep_Dr:
	  nHad_tau=nHad_tau+1
	  Had_tau_id.append(itau)
    self.out.fillBranch("nHad_tau", nHad_tau)

    met_user=-99
    met_phi_user=-99

    if self.is_mc:
      met_user=event.MET_T1Smear_pt
      met_phi_user=event.MET_T1Smear_phi
    else:
      met_user=event.MET_T1_pt
      met_phi_user=event.MET_T1_phi

    self.out.fillBranch("met_user",met_user)
    self.out.fillBranch("met_phi_user",met_phi_user)


    fatjet_v4_temp=TLorentzVector()
    # Fat jet selection
    # due to the potential overlap between jet and lepton, Tight ID is used. (*10>=2)
    # fatjets = Collection(event, 'FatJetNoVlep')
    # fatjets_init = Collection(event, 'FatJet')
    fatjets = Collection(event, 'FatJet')
    l1v4_tmp.SetPtEtaPhiM(event.l1_pt,event.l1_eta,event.l1_phi,event.l1_mass)
    l2v4_tmp.SetPtEtaPhiM(event.l2_pt,event.l2_eta,event.l2_phi,event.l2_mass)

    TightFatJet_id = []
    TightFatJet_initid = []
    TightFatJet_bbvsQCD = []
    TightFatJet_ccvsQCD = []
    TightFatJet_qqvsQCD = []
    TightFatJet_drl1 = []
    TightFatJet_drl2 = []
    TightFatJet_pt = []
    TightFatJet_eta = []
    TightFatJet_phi = []
    TightFatJet_mass = []
    TightFatJet_dRinit = []
    TightFatJet_v4_all = []
    # local jet id, with dr(lepton)>0.8 or <0.8
    TightFatJet_drLa08_id = []
    TightFatJet_drSm08_id = []

    nFatjet_MP=-1

    for ijet in range(0, event.nFatJet):
      if abs(fatjets[ijet].eta)>2.4 :continue
      if not (fatjets[ijet].jetId>1 and fatjets[ijet].pt>300):continue
      # default jet mass is FatJet_mass
      fatjet_v4_temp.SetPtEtaPhiM(fatjets[ijet].pt_nom,fatjets[ijet].eta,fatjets[ijet].phi,fatjets[ijet].mass_nom)
      # order of jet in LepSubtraction FatJet collection
      TightFatJet_id.append(ijet)
      #TightFatJet_dRinit.append(fatjets[ijet].dRinit)
      # order of jet in initial NanoAOD FatJet collection
      #TightFatJet_initid.append(fatjets[ijet].id)
      TightFatJet_pt.append(fatjets[ijet].pt_nom)
      TightFatJet_eta.append(fatjets[ijet].eta)
      TightFatJet_phi.append(fatjets[ijet].phi)
      TightFatJet_mass.append(fatjets[ijet].mass_nom)
      TightFatJet_bbvsQCD.append(fatjets[ijet].particleNetMD_Xbb/(fatjets[ijet].particleNetMD_Xbb + fatjets[ijet].particleNetMD_QCD))
      TightFatJet_ccvsQCD.append(fatjets[ijet].particleNetMD_Xcc/(fatjets[ijet].particleNetMD_Xcc + fatjets[ijet].particleNetMD_QCD))
      TightFatJet_qqvsQCD.append((fatjets[ijet].particleNetMD_Xcc + fatjets[ijet].particleNetMD_Xqq)/(fatjets[ijet].particleNetMD_Xcc + fatjets[ijet].particleNetMD_Xqq + fatjets[ijet].particleNetMD_QCD))
      # https://github.com/cms-sw/cmssw/blob/6d2f66057131baacc2fcbdd203588c41c885b42c/RecoBTag/ONNXRuntime/python/pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags_cfi.py#L74
#      fatJets_qqvsQCD.append((fatjets[ijet].deepTagMD_ZvsQCD)
      TightFatJet_drl1.append(l1v4_tmp.DeltaR(fatjet_v4_temp))
      TightFatJet_drl2.append(l2v4_tmp.DeltaR(fatjet_v4_temp))
      TightFatJet_v4_all.append(fatjet_v4_temp.Clone())


    nFatjet_MP=len(TightFatJet_id)
    if nFatjet_MP>0:
      for ij_tmp in range(nFatjet_MP):
        # e.g., if TightFatJet_id=[0,2,3], it mean the 1th,3rd,4th jet in FatJetNoVLep pass the MP
        if TightFatJet_drl1[ij_tmp]>0.8 and TightFatJet_drl1[ij_tmp]>0.8:
          TightFatJet_drLa08_id.append(ij_tmp)
        else:
          TightFatJet_drSm08_id.append(ij_tmp)

    # check the fat jet property
    # for ParticleNetMD, WP for bbVsQCD, in BTV-22-001, HP,MP,LP correspond to signal eff 40%, 60% and 80%
    # 2016apv pre-VFP HP:0.9883, MP:0.9737, LP:0.9088
    # 2016 post-VFP HP:0.9883, MP:0.9735, LP:0.9137
    # 2017 HP:0.9870, MP:0.9714, LP:0.9105
    # 2018 HP:0.9880, MP:0.9734, LP:0.9172
    Hbb_th=-1
    if self.year=="2016apv":Hbb_threshold=0.9088
    if self.year=="2016":Hbb_threshold=0.9137
    if self.year=="2017":Hbb_threshold=0.9105
    if self.year=="2018":Hbb_threshold=0.9172

    # for ParticleNetMD, WP for ccVsQCD, in BTV-22-001, HP,MP,LP correspond to signal eff 40%, 60% and 80%
    # 2016apv pre-VFP HP:0.9909, MP:0.9751, LP:0.9252
    # 2016 post-VFP HP:0.9905, MP:0.9743, LP:0.9252
    # 2017 HP:0.9909, MP:0.9765, LP:0.9347
    # 2018 HP:0.9917, MP:0.9777 LP:0.9368
    Hcc_th=-1
    if self.year=="2016apv":Hcc_threshold=0.9252
    if self.year=="2016":Hcc_threshold=0.9252
    if self.year=="2017":Hcc_threshold=0.9252
    if self.year=="2018":Hcc_threshold=0.9368
    
    
    # two categories (for H and Z decay), 2AK4: Two AK4 jets, 1AK8: One AK8 Jet
    H_1AK8=False
    H_1AK8_case=-1  
    H_2AK4=False
    H_2AK4_case=-1  

    H_AK8Jet_id=-1
    H_AK8Jet_initid=-1
    H_AK8Jet_pt=-99
    H_AK8Jet_eta=-99
    H_AK8Jet_phi=-99
    #FatJet_mass
    H_AK8Jet_mass=-99
    #FatJet_particleNet_mass
    H_AK8Jet_PNmass=-99
    #FatJet_msoftdrop
    H_AK8Jet_SDmass=-99
    H_AK8Jet_drl1=-99
    H_AK8Jet_drl2=-99


    # no fat jet passing tight ID
    if nFatjet_MP==0: 
      H_2AK4=True
    
    # only one fat jet passing tight ID
    # if it pass the MP then it's a Higgs->bb fat jet
    if nFatjet_MP==1:
      H_1AK8=True
      H_1AK8_case=1
      if TightFatJet_drl1[0]>0.8 and TightFatJet_drl2[0]>0.8:
        if TightFatJet_bbvsQCD[0]>Hbb_threshold and TightFatJet_ccvsQCD[0]<Hcc_threshold:
          H_1AK8_case=0
          H_AK8Jet_id=TightFatJet_id[0]
          #H_AK8Jet_initid=TightFatJet_initid[0]
          H_AK8Jet_pt=TightFatJet_v4_all[0].Pt()
          H_AK8Jet_eta=TightFatJet_v4_all[0].Eta()
          H_AK8Jet_phi=TightFatJet_v4_all[0].Phi()
          H_AK8Jet_mass=TightFatJet_v4_all[0].M()
          #H_AK8Jet_PNmass=event.FatJet_particleNet_mass[H_AK8Jet_initid]
          #H_AK8Jet_SDmass=event.FatJet_msoftdrop[H_AK8Jet_initid]
          H_AK8Jet_drl1=TightFatJet_drl1[0]
          H_AK8Jet_drl2=TightFatJet_drl2[0]        
      else:
        H_1AK8=False
    else:
        H_1AK8=False



    self.out.fillBranch("TightFatJet_id", TightFatJet_id)
    self.out.fillBranch("TightFatJet_bbvsQCD", TightFatJet_bbvsQCD)
    self.out.fillBranch("TightFatJet_ccvsQCD", TightFatJet_ccvsQCD)
    self.out.fillBranch("TightFatJet_qqvsQCD", TightFatJet_qqvsQCD)
    self.out.fillBranch("TightFatJet_drl1", TightFatJet_drl1)
    self.out.fillBranch("TightFatJet_drl2", TightFatJet_drl2)
    self.out.fillBranch("TightFatJet_pt", TightFatJet_pt)
    self.out.fillBranch("TightFatJet_eta", TightFatJet_eta)
    self.out.fillBranch("TightFatJet_phi", TightFatJet_phi)
    self.out.fillBranch("TightFatJet_mass", TightFatJet_mass)
    self.out.fillBranch("TightFatJet_dRinit", TightFatJet_dRinit)
    self.out.fillBranch("H_1AK8",H_1AK8)
    self.out.fillBranch("H_1AK8_case",H_1AK8_case)
    self.out.fillBranch("H_AK8Jet_id",H_AK8Jet_id)
    self.out.fillBranch("H_AK8Jet_initid",H_AK8Jet_initid)
    self.out.fillBranch("H_AK8Jet_pt",H_AK8Jet_pt)
    self.out.fillBranch("H_AK8Jet_eta",H_AK8Jet_eta)
    self.out.fillBranch("H_AK8Jet_phi",H_AK8Jet_phi)
    self.out.fillBranch("H_AK8Jet_mass",H_AK8Jet_mass)
    self.out.fillBranch("H_AK8Jet_PNmass",H_AK8Jet_PNmass)
    self.out.fillBranch("H_AK8Jet_SDmass",H_AK8Jet_SDmass)
    self.out.fillBranch("H_AK8Jet_drl1",H_AK8Jet_drl1)
    self.out.fillBranch("H_AK8Jet_drl2",H_AK8Jet_drl2)

    # https://btv-wiki.docs.cern.ch/ScaleFactors/#sf-campaigns
    # tight PF jets (ak4), *10>=2
    # medium B-tag WP
    # DeepFlavor=(nanoaod btagDeepFlavB) 
    # 2016apv loose: 0.0508, medium: 0.2598, tight: 0.6502
    # 2016 loose: 0.0480, medium: 0.2489, tight: 0.6377
    # 2017 loose: 0.0532, medium: 0.3040, tight: 0.7476
    # 2018 loose: 0.0490, medium: 0.2783, tight: 0.7100

    # jets = Collection(event, 'JetNoVlep')
    # jets_init = Collection(event, 'Jet')
    jets = Collection(event, 'Jet')


    TightAK4Jet_id = []
    TightAK4Jet_pt = []
    TightAK4Jet_eta = []
    TightAK4Jet_phi = []
    TightAK4Jet_mass = []
    TightAK4Jet_drl1 = []
    TightAK4Jet_drl2 = []
    #TightAK4Jet_dRinit = []

    TightAK4Jet_nob_id = []
    TightAK4Jet_b_DeepCSVmedium_id = []
    TightAK4Jet_b_DeepCSVloose_id = []
    nobjet_v4_passdrlep_id = []
    nobjet_v4_faildrlep_id = []

    jet_v4_all = []
    bjet_v4_all = []
    nobjet_v4_all = []
    nobjet_v4_passdrlep = []
    nobjet_v4_faildrlep = []

    # deepflavB, default is 2016apv WP
    medium_Bcut = 0.2598
    loose_Bcut = 0.0508
    if self.year=="2016":
      medium_Bcut = 0.2489
      loose_Bcut= 0.048
    if self.year=="2017":
      medium_Bcut = 0.3040
      loose_Bcut= 0.0532
    if self.year=="2018":
      medium_Bcut = 0.2783
      loose_Bcut= 0.0490

    jet_v4_temp=TLorentzVector()

    for ijet in range(0, event.nJet):

      if abs(jets[ijet].eta)>2.4:continue

      jet_is_tau=0
      if nHad_tau>0:
        for ita in Had_tau_id:
          if ijet==event.Tau_jetIdx[ita]:jet_is_tau=1
      if jet_is_tau:continue

      # pass tight ID, i.e., not require dR(jet,lep) here
      if not (jets[ijet].jetId>1 and jets[ijet].pt_nom>30):continue 
      jet_v4_temp.SetPtEtaPhiM(jets[ijet].pt_nom,jets[ijet].eta,jets[ijet].phi,jets[ijet].mass_nom)
      l1v4_tmp.SetPtEtaPhiM(event.l1_pt,event.l1_eta,event.l1_phi,event.l1_mass)
      l2v4_tmp.SetPtEtaPhiM(event.l2_pt,event.l2_eta,event.l2_phi,event.l2_mass)

      # remove overlap with FatJet
      if H_AK8Jet_id>-1:
        fatjetv4_tmp_=TLorentzVector()
        fatjetv4_tmp_.SetPtEtaPhiM(H_AK8Jet_pt,H_AK8Jet_eta,H_AK8Jet_phi,H_AK8Jet_mass)
        if jet_v4_temp.DeltaR(fatjetv4_tmp_)<0.8:continue

      TightAK4Jet_id.append(ijet)
      TightAK4Jet_pt.append(jets[ijet].pt_nom)
      TightAK4Jet_eta.append(jets[ijet].eta)
      TightAK4Jet_phi.append(jets[ijet].phi)
      TightAK4Jet_mass.append(jets[ijet].mass_nom)
      TightAK4Jet_drl1.append(l1v4_tmp.DeltaR(jet_v4_temp))
      TightAK4Jet_drl2.append(l2v4_tmp.DeltaR(jet_v4_temp))
      #TightAK4Jet_dRinit.append(jets[ijet].dRinit)
      jet_v4_all.append(jet_v4_temp.Clone())
 
      # only use flavor tag infor for those jet with dR(lep)>0.4, i.e., the lepton subtraction will change the flavor tag infor and can't be used anymore after lepton subtraction
      if l1v4_tmp.DeltaR(jet_v4_temp)>0.4 and l2v4_tmp.DeltaR(jet_v4_temp)>0.4:
        if not jets[ijet].btagDeepFlavB > loose_Bcut:
          nobjet_v4_passdrlep_id.append(ijet)
          nobjet_v4_passdrlep.append(jet_v4_temp.Clone())
        else:
          TightAK4Jet_b_DeepCSVloose_id.append(ijet)
          bjet_v4_all.append(jet_v4_temp.Clone())
          if jets[ijet].btagDeepFlavB > medium_Bcut:
            TightAK4Jet_b_DeepCSVmedium_id.append(ijet)
      else:
        nobjet_v4_faildrlep_id.append(ijet)
        nobjet_v4_faildrlep.append(jet_v4_temp.Clone())
 
    HT=0
    for ijet in TightAK4Jet_id:
      HT=HT + jets[ijet].pt_nom
    self.out.fillBranch("HT",HT)

    TightAK4Jet_nob_id = [x for x in TightAK4Jet_id if (x not in TightAK4Jet_b_DeepCSVloose_id)]
    nobjet_v4_all = [x for x in jet_v4_all if x not in bjet_v4_all]

    n_tight_jet=len(TightAK4Jet_id)

    n_bjet_DeepB_M = len(TightAK4Jet_b_DeepCSVmedium_id)
    n_bjet_DeepB_L = len(TightAK4Jet_b_DeepCSVloose_id)
    n_tight_nob = len(TightAK4Jet_nob_id)
    self.out.fillBranch("n_bjet_DeepB_M",n_bjet_DeepB_M)
    self.out.fillBranch("n_bjet_DeepB_L",n_bjet_DeepB_L)
    self.out.fillBranch("n_tight_nob",n_tight_nob)

    Had_tau_id.extend(np.zeros(event.nTau-len(Had_tau_id),int)-1)
    self.out.fillBranch("Had_tau_id", Had_tau_id)
    
#      pass_jet_lep_Dr=1
#      for ilep in range(0,len(tightLeptons)):
#	if jet_v4_temp.DeltaR(tightLeptons[ilep])<0.4:pass_jet_lep_Dr=0
#
#      if not (pass_jet_lep_Dr>0):continue

    h_j1_pt=-99
    h_j1_eta=-99
    h_j1_phi=-99
    h_j1_mass=-99
    h_j1_id=-99
    h_j1_drl1=-99
    h_j1_drl2=-99
    h_j2_pt=-99
    h_j2_eta=-99
    h_j2_phi=-99
    h_j2_mass=-99
    h_j2_id=-99
    h_j2_drl1=-99
    h_j2_drl2=-99
    h_mjj=-99
    h_detajj=-99
    h_dRjj=-99
    h_dphijj=-99

    jet_v4_temp.SetPtEtaPhiM(jets[ijet].pt_nom,jets[ijet].eta,jets[ijet].phi,jets[ijet].mass_nom)
    l1v4_tmp.SetPtEtaPhiM(event.l1_pt,event.l1_eta,event.l1_phi,event.l1_mass)
    l2v4_tmp.SetPtEtaPhiM(event.l2_pt,event.l2_eta,event.l2_phi,event.l2_mass)





    # 2 resolved jets category
    if H_2AK4:
      # at least 2 good jets
      if n_tight_jet<2:return False
      # if only one b-tag jet, must be medium ID 
      if n_bjet_DeepB_L<1:return False
      if n_bjet_DeepB_L==1 and n_bjet_DeepB_M<1:return False

      hbb_mass_threshold=99
      hbb_v4_temp=-1
  
      if n_bjet_DeepB_L==1:
        if len(nobjet_v4_passdrlep)>0:
          for ij in range(0,len(nobjet_v4_passdrlep)):
            if abs((bjet_v4_all[0]+nobjet_v4_passdrlep[ij]).M()-125.)<hbb_mass_threshold:
              hbb_mass_threshold=abs((bjet_v4_all[0]+nobjet_v4_passdrlep[ij]).M()-125.)
              hbb_v4_temp=ij
    
          if bjet_v4_all[0].Pt()>nobjet_v4_passdrlep[hbb_v4_temp].Pt():
            H_2AK4_case=0
            h_j1_pt=bjet_v4_all[0].Pt()
            h_j1_eta=bjet_v4_all[0].Eta()
            h_j1_phi=bjet_v4_all[0].Phi()
            h_j1_mass=bjet_v4_all[0].M()
            h_j1_id=TightAK4Jet_b_DeepCSVloose_id[0]
            h_j1_drl1=jets[h_j1_id].DeltaR(l1v4_tmp)
            h_j1_drl2=jets[h_j1_id].DeltaR(l2v4_tmp)
            h_j2_pt=nobjet_v4_passdrlep[hbb_v4_temp].Pt()
            h_j2_eta=nobjet_v4_passdrlep[hbb_v4_temp].Eta()
            h_j2_phi=nobjet_v4_all[hbb_v4_temp].Phi()
            h_j2_mass=nobjet_v4_all[hbb_v4_temp].M()
            h_j2_id=nobjet_v4_passdrlep_id[hbb_v4_temp]
            h_j2_drl1=jets[h_j2_id].DeltaR(l1v4_tmp)
            h_j2_drl2=jets[h_j2_id].DeltaR(l2v4_tmp)
          else:
            H_2AK4_case=1
            h_j1_pt=nobjet_v4_all[hbb_v4_temp].Pt()
            h_j1_eta=nobjet_v4_all[hbb_v4_temp].Eta()
            h_j1_phi=nobjet_v4_all[hbb_v4_temp].Phi()
            h_j1_mass=nobjet_v4_all[hbb_v4_temp].M()
            h_j1_id=nobjet_v4_passdrlep_id[hbb_v4_temp]
            h_j1_drl1=jets[h_j1_id].DeltaR(l1v4_tmp)
            h_j1_drl2=jets[h_j1_id].DeltaR(l2v4_tmp)
            h_j2_pt=bjet_v4_all[0].Pt()
            h_j2_eta=bjet_v4_all[0].Eta()
            h_j2_phi=bjet_v4_all[0].Phi()
            h_j2_mass=bjet_v4_all[0].M()
            h_j2_id=TightAK4Jet_b_DeepCSVloose_id[0]
            h_j2_drl1=jets[h_j2_id].DeltaR(l2v4_tmp)
            h_j2_drl2=jets[h_j2_id].DeltaR(l2v4_tmp)

          h_mjj=(bjet_v4_all[0]+nobjet_v4_passdrlep[hbb_v4_temp]).M()
          h_detajj=abs(h_j1_eta-h_j2_eta)
          h_dRjj=bjet_v4_all[0].DeltaR(nobjet_v4_passdrlep[hbb_v4_temp])
          h_dphijj=bjet_v4_all[0].DeltaPhi(nobjet_v4_passdrlep[hbb_v4_temp])
          
      if n_bjet_DeepB_L==2:
        H_2AK4_case=2
        h_j1_pt=bjet_v4_all[0].Pt()
        h_j1_eta=bjet_v4_all[0].Eta()
        h_j1_phi=bjet_v4_all[0].Phi()
        h_j1_mass=bjet_v4_all[0].M()
        h_j1_id=TightAK4Jet_b_DeepCSVloose_id[0]
        h_j1_drl1=jets[h_j1_id].DeltaR(l1v4_tmp)
        h_j1_drl2=jets[h_j1_id].DeltaR(l2v4_tmp)
        h_j2_pt=bjet_v4_all[1].Pt()
        h_j2_eta=bjet_v4_all[1].Eta()
        h_j2_phi=bjet_v4_all[1].Phi()
        h_j2_mass=bjet_v4_all[1].M()
        h_j2_id=TightAK4Jet_b_DeepCSVloose_id[1]
        h_j2_drl1=jets[h_j2_id].DeltaR(l1v4_tmp)
        h_j2_drl2=jets[h_j2_id].DeltaR(l2v4_tmp)
        h_mjj=(bjet_v4_all[0]+bjet_v4_all[1]).M()
        h_detajj=abs(h_j1_eta-h_j2_eta)
        h_dRjj=bjet_v4_all[0].DeltaR(bjet_v4_all[1])
        h_dphijj=bjet_v4_all[0].DeltaPhi(bjet_v4_all[1])
          
      hbb_id_item_temp=[]
      hbb_mass_temp=[]
      hbb_min_item_temp=-99
      if n_bjet_DeepB_L>2:
        for b1_id_temp,b1_v4_temp in enumerate(bjet_v4_all):
          for b2_id_temp,b2_v4_temp in enumerate(bjet_v4_all):
            if b1_id_temp<b2_id_temp:
              hbb_id_item_temp.append((b1_id_temp,b2_id_temp))
              hbb_mass_temp.append(abs((b1_v4_temp+b2_v4_temp).M()-125.))
        hbb_min_item_temp=hbb_mass_temp.index(min(hbb_mass_temp))
        h_j1_v4_temp=bjet_v4_all[hbb_id_item_temp[hbb_min_item_temp][0]]
        h_j2_v4_temp=bjet_v4_all[hbb_id_item_temp[hbb_min_item_temp][1]]
    
        H_2AK4_case=3
        h_j1_pt=h_j1_v4_temp.Pt()
        h_j1_eta=h_j1_v4_temp.Eta()
        h_j1_phi=h_j1_v4_temp.Phi()
        h_j1_mass=h_j1_v4_temp.M()
        h_j1_id=TightAK4Jet_b_DeepCSVloose_id[hbb_id_item_temp[hbb_min_item_temp][0]]
        h_j1_drl1=jets[h_j1_id].DeltaR(l1v4_tmp)
        h_j1_drl2=jets[h_j1_id].DeltaR(l2v4_tmp)
        h_j2_pt=h_j2_v4_temp.Pt()
        h_j2_eta=h_j2_v4_temp.Eta()
        h_j2_phi=h_j2_v4_temp.Phi()
        h_j2_mass=h_j2_v4_temp.M()
        h_j2_id=TightAK4Jet_b_DeepCSVloose_id[hbb_id_item_temp[hbb_min_item_temp][1]]
        h_j2_drl1=jets[h_j2_id].DeltaR(l1v4_tmp)
        h_j2_drl2=jets[h_j2_id].DeltaR(l2v4_tmp)
        h_mjj=(h_j1_v4_temp+h_j2_v4_temp).M()
        h_detajj=abs(h_j1_eta-h_j2_eta)
        h_dRjj=h_j1_v4_temp.DeltaR(h_j2_v4_temp)
        h_dphijj=h_j1_v4_temp.DeltaPhi(h_j2_v4_temp)


    self.out.fillBranch("h_j1_pt",h_j1_pt)
    self.out.fillBranch("h_j1_eta",h_j1_eta)
    self.out.fillBranch("h_j1_phi",h_j1_phi)
    self.out.fillBranch("h_j1_mass",h_j1_mass)
    self.out.fillBranch("h_j1_id",h_j1_id)
    self.out.fillBranch("h_j1_drl1",h_j1_drl1)
    self.out.fillBranch("h_j1_drl2",h_j1_drl2)
    self.out.fillBranch("h_j2_pt",h_j2_pt)
    self.out.fillBranch("h_j2_eta",h_j2_eta)
    self.out.fillBranch("h_j2_phi",h_j2_phi)
    self.out.fillBranch("h_j2_mass",h_j2_mass)
    self.out.fillBranch("h_j2_id",h_j2_id)
    self.out.fillBranch("h_j2_drl1",h_j2_drl1)
    self.out.fillBranch("h_j2_drl2",h_j2_drl2)
    self.out.fillBranch("h_mjj",h_mjj)
    self.out.fillBranch("h_detajj",h_detajj)
    self.out.fillBranch("h_dRjj",h_dRjj)
    self.out.fillBranch("h_dphijj",h_dphijj)

    TightAK4Jet_nob_id.extend(np.zeros(n_tight_jet-n_tight_nob,int)-1)
    TightAK4Jet_b_DeepCSVmedium_id.extend(np.zeros(n_tight_jet-n_bjet_DeepB_M,int)-1)
    TightAK4Jet_b_DeepCSVloose_id.extend(np.zeros(n_tight_jet-n_bjet_DeepB_L,int)-1)

    self.out.fillBranch("TightAK4Jet_id",TightAK4Jet_id)
    self.out.fillBranch("TightAK4Jet_nob_id",TightAK4Jet_nob_id)
    self.out.fillBranch("TightAK4Jet_b_DeepCSVmedium_id",TightAK4Jet_b_DeepCSVmedium_id)
    self.out.fillBranch("TightAK4Jet_b_DeepCSVloose_id",TightAK4Jet_b_DeepCSVloose_id)
    self.out.fillBranch("TightAK4Jet_drl1",TightAK4Jet_drl1)
    self.out.fillBranch("TightAK4Jet_drl2",TightAK4Jet_drl2)
    self.out.fillBranch("TightAK4Jet_pt",TightAK4Jet_pt)
    self.out.fillBranch("TightAK4Jet_eta",TightAK4Jet_eta)
    self.out.fillBranch("TightAK4Jet_phi",TightAK4Jet_phi)
    self.out.fillBranch("TightAK4Jet_mass",TightAK4Jet_mass)
    



    return True

HH12016apv = lambda: HHproducer1("2016apv")
HH12016 = lambda: HHproducer1("2016")
HH12017 = lambda: HHproducer1("2017")
HH12018 = lambda: HHproducer1("2018")
