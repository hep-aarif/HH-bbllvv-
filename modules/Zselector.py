import  ROOT
from  ROOT  import  TLorentzVector
ROOT.PyConfig.IgnoreCommandLineOptions  =  True

from  PhysicsTools.NanoAODTools.postprocessing.framework.datamodel  import  Collection
from  PhysicsTools.NanoAODTools.postprocessing.framework.eventloop  import  Module

class  Zselector(Module):

    def  __init__(self):
        pass

    def  beginJob(self):
        pass
    def  endJob(self):
        pass
        
    def  beginFile(self,  inputFile,  outputFile,  inputTree,  wrappedOutputTree):
        self.out  =  wrappedOutputTree

        self.out.branch("mu_channel",  "B")
        self.out.branch("l1_pt",  "F")
        self.out.branch("l1_rawpt",  "F")
        self.out.branch("l1_eta",  "F")
        self.out.branch("l1_phi",  "F")
        self.out.branch("l1_mass",  "F")
        self.out.branch("l1_rawmass",  "F")
        self.out.branch("l1_relPtl2",  "F")
        self.out.branch("l2_pt",  "F")
        self.out.branch("l2_rawpt",  "F")
        self.out.branch("l2_eta",  "F")
        self.out.branch("l2_phi",  "F")
        self.out.branch("l2_mass",  "F")
        self.out.branch("l2_rawmass",  "F")
        self.out.branch("l2_relPtl1",  "F")
        self.out.branch("drll",  "F")
        self.out.branch("detall",  "F")
        self.out.branch("dphill",  "F")
        self.out.branch("zlep_pt",  "F")
        self.out.branch("zlep_eta",  "F")
        self.out.branch("zlep_phi",  "F")
        self.out.branch("zlep_mass",  "F")

        

    def  endFile(self,  inputFile,  outputFile,  inputTree,  wrappedOutputTree):
        pass

    def  analyze(self,  event):

        muons  =  Collection(event,  'Muon')
        goodmuons  =  Collection(event,  'GoodMuon')
        goodeles  =  Collection(event,  'GoodElectron')
    

        mu_channel=False
        l1_pt=-99
        l1_rawpt=-99
        l1_eta=-99
        l1_phi=-99
        l1_mass=-99
        l1_rawmass=-99
        l1_relPtl2=-99
        l2_pt=-99
        l2_rawpt=-99
        l2_eta=-99
        l2_phi=-99
        l2_mass=-99
        l2_rawmass=-99
        l2_relPtl1=-99
        drll=-99
        detall=-99
        dphill=-99
        zlep_pt=-99
        zlep_eta=-99
        zlep_phi=-99
        zlep_mass=-99

        #  only  two  tight  leptons
        if  not  (event.nGoodMuon+event.nGoodElectron)==2:return  False
        if  not  (event.nGoodMuon==2  or  event.nGoodElectron==2):return  False

        l1v4_tmp=TLorentzVector()
        l2v4_tmp=TLorentzVector()
        l1v4raw_tmp=TLorentzVector()
        l2v4raw_tmp=TLorentzVector()

        if  event.nGoodMuon==2:
            mu_channel=True
            l1_pt=goodmuons[0].pt
            l1_rawpt=goodmuons[0].rawpt
            l1_eta=goodmuons[0].eta
            l1_phi=goodmuons[0].phi
            l1_mass=goodmuons[0].mass
            l1_rawmass=goodmuons[0].rawmass
            l2_pt=goodmuons[1].pt
            l2_rawpt=goodmuons[1].rawpt
            l2_eta=goodmuons[1].eta
            l2_phi=goodmuons[1].phi
            l2_mass=goodmuons[1].mass
            l2_rawmass=goodmuons[1].rawmass
        else:
            mu_channel=False
            l1_pt=goodeles[0].pt
            l1_rawpt=goodeles[0].rawpt
            l1_eta=goodeles[0].eta
            l1_phi=goodeles[0].phi
            l1_mass=goodeles[0].mass
            l1_rawmass=goodeles[0].rawmass
            l2_pt=goodeles[1].pt
            l2_rawpt=goodeles[1].rawpt
            l2_eta=goodeles[1].eta
            l2_phi=goodeles[1].phi
            l2_mass=goodeles[1].mass
            l2_rawmass=goodeles[1].rawmass

        l1v4_tmp.SetPtEtaPhiM(l1_pt,l1_eta,l1_phi,l1_mass)
        l2v4_tmp.SetPtEtaPhiM(l2_pt,l2_eta,l2_phi,l2_mass)
        l1v4raw_tmp.SetPtEtaPhiM(l1_rawpt,l1_eta,l1_phi,l1_rawmass)
        l2v4raw_tmp.SetPtEtaPhiM(l2_rawpt,l2_eta,l2_phi,l2_rawmass)
        l1_relPtl2=l1v4_tmp.Vect().Cross(l2v4_tmp.Vect().Unit()).Mag()
        l2_relPtl1=l2v4_tmp.Vect().Cross(l1v4_tmp.Vect().Unit()).Mag()
        drll=l1v4_tmp.DeltaR(l2v4_tmp)
        detall=abs(l1_eta  -  l2_eta)
        dphill=l1v4_tmp.DeltaPhi(l2v4_tmp)
        zlep_pt=(l1v4_tmp  +  l2v4_tmp).Pt()
        zlep_eta=(l1v4_tmp  +  l2v4_tmp).Eta()
        zlep_phi=(l1v4_tmp  +  l2v4_tmp).Phi()
        zlep_mass=(l1v4_tmp  +  l2v4_tmp).M()

       
        self.out.fillBranch("mu_channel",  mu_channel)
        self.out.fillBranch("l1_pt",  l1_pt)
        self.out.fillBranch("l1_rawpt",  l1_rawpt)
        self.out.fillBranch("l1_eta",  l1_eta)
        self.out.fillBranch("l1_phi",  l1_phi)
        self.out.fillBranch("l1_mass",  l1_mass)
        self.out.fillBranch("l1_rawmass",  l1_rawmass)
        self.out.fillBranch("l1_relPtl2",  l1_relPtl2)
        self.out.fillBranch("l2_pt",  l2_pt)
        self.out.fillBranch("l2_rawpt",  l2_rawpt)
        self.out.fillBranch("l2_eta",  l2_eta)
        self.out.fillBranch("l2_phi",  l2_phi)
        self.out.fillBranch("l2_mass",  l2_mass)
        self.out.fillBranch("l2_rawmass",  l2_rawmass)
        self.out.fillBranch("l2_relPtl1",  l2_relPtl1)
        self.out.fillBranch("drll",  drll)
        self.out.fillBranch("detall",  detall)
        self.out.fillBranch("dphill",  dphill)
        self.out.fillBranch("zlep_pt",  zlep_pt)
        self.out.fillBranch("zlep_eta",  zlep_eta)
        self.out.fillBranch("zlep_phi",  zlep_phi)
        self.out.fillBranch("zlep_mass",  zlep_mass)

        return  True

ZProducer  =  lambda:  Zselector()



