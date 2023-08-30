import ROOT

# Open the ROOT file and get the tree
fin1 = ROOT.TFile.Open('tree.root')
treein = fin1.Get('Events')

# Create a TH1D histograms
hist_z_mass = ROOT.TH1D('z_mass', 'z_mass', 30, 0, 120)
hist_z_pt = ROOT.TH1D('z_pt', 'z_pt', 20, 0, 1000)
hist_dr_l1l2 = ROOT.TH1D('dr_l1l2', 'dr_l1l2', 40, 0, 4)
hist_dr_l1j1 = ROOT.TH1D('dr_l1j1', 'dr_l1j1', 40, 0, 4)
hist_dr_l1j2 = ROOT.TH1D('dr_l1j2', 'dr_l1j2', 40, 0, 4)
hist_dr_l2j1 = ROOT.TH1D('dr_l2j1', 'dr_l2j1', 40, 0, 4)
hist_dr_l2j2 = ROOT.TH1D('dr_l2j2', 'dr_l2j2', 40, 0, 4)
hist_dr_j1j2 = ROOT.TH1D('dr_j1j2', 'dr_j1j2', 40, 0, 4)
hist_l1_pt = ROOT.TH1D('l1_pt', 'l1_pt', 20, 0, 1000)
hist_l2_pt = ROOT.TH1D('l2_pt', 'l2_pt', 20, 0, 1000)
hist_j1_pt = ROOT.TH1D('j1_pt', 'j1_pt', 20, 0, 1000)
hist_j2_pt = ROOT.TH1D('j2_pt', 'j2_pt', 20, 0, 1000)


z = ROOT.TLorentzVector()
l1 = ROOT.TLorentzVector()
l2 = ROOT.TLorentzVector()
j1 = ROOT.TLorentzVector()
j2 = ROOT.TLorentzVector()


# Loop through the events and fill the histogram
for ie in range(0, treein.GetEntriesFast()):

 treein.GetEntry(ie)
 z.SetPtEtaPhiM(treein.zlep_pt, treein.zlep_eta, treein.zlep_phi, treein.zlep_mass)
 l1.SetPtEtaPhiM(treein.l1_pt, treein.l1_eta, treein.l1_phi, treein.l1_mass)
 l2.SetPtEtaPhiM(treein.l2_pt, treein.l2_eta, treein.l2_phi, treein.l2_mass)
 j1.SetPtEtaPhiM(treein.h_j1_pt, treein.h_j1_eta, treein.h_j1_phi, treein.h_j1_mass)
 j2.SetPtEtaPhiM(treein.h_j2_pt, treein.h_j2_eta, treein.h_j2_phi, treein.h_j2_mass)
 hist_z_mass.Fill(z.M())
 hist_z_pt.Fill(z.Pt())
 hist_dr_l1l2.Fill(l1.DeltaR(l2))
 hist_dr_l1j1.Fill(l1.DeltaR(j1))
 hist_dr_l1j2.Fill(l1.DeltaR(j2))
 hist_dr_l2j1.Fill(l2.DeltaR(j1))
 hist_dr_l2j2.Fill(l2.DeltaR(j2))
 hist_dr_j1j2.Fill(j1.DeltaR(j2))
 hist_l1_pt.Fill(l1.Pt())
 hist_l2_pt.Fill(l2.Pt())
 hist_j1_pt.Fill(j1.Pt())
 hist_j2_pt.Fill(j2.Pt())

# Create a canvas
output_name = 'Output_histograms.root'
fileout = ROOT.TFile.Open(output_name, 'recreate')
fileout.cd()

hist_z_pt.Write()
hist_z_mass.Write()
hist_dr_l1l2.Write()
hist_dr_l1j1.Write()
hist_dr_l1j2.Write()
hist_dr_l2j1.Write()
hist_dr_l2j2.Write()
hist_dr_j1j2.Write()
hist_l1_pt.Write()
hist_l2_pt.Write()
hist_j1_pt.Write()
hist_j2_pt.Write()

fileout.Close()
