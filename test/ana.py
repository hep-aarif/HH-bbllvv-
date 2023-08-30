import ROOT

fin1=ROOT.TFile.Open('Mu_powheg.root')
treein1=fin1.Get('Events')

fin2=ROOT.TFile.Open('all_mg.root')
treein2=fin2.Get('Events')

h1=ROOT.TH1D('powheg_mll','powheg_mll',100,0,200)
h2=ROOT.TH1D('mg_mll','mg_mll',100,0,200)
l1=ROOT.TLorentzVector()
l2=ROOT.TLorentzVector()

nEvt_powheg_posi=470932
nEvt_powheg_nega=19568
powheg_weight=2.2155068e-06

for ie in range(0,treein1.GetEntriesFast()):
  if ie%1000==0:print('processing ',ie)
  treein1.GetEntry(ie)
  mu_arr=[]
  if treein1.nMuon<1:continue
  for im in range(0,treein1.nMuon):
    if not treein1.Muon_mediumId[im]:continue
    if treein1.Muon_pt[im]<20:continue
    mu_arr.append(im)
  if not len(mu_arr)==2:continue
  l1.SetPtEtaPhiM(treein1.Muon_pt[mu_arr[0]],treein1.Muon_eta[mu_arr[0]],treein1.Muon_phi[mu_arr[0]],treein1.Muon_mass[mu_arr[0]])
  l2.SetPtEtaPhiM(treein1.Muon_pt[mu_arr[1]],treein1.Muon_eta[mu_arr[1]],treein1.Muon_phi[mu_arr[1]],treein1.Muon_mass[mu_arr[1]])
  h1.Fill((l1+l2).M(),treein1.Generator_weight/(abs(treein1.Generator_weight)))

mg_weight=(2./treein2.GetEntriesFast())
for ie in range(0,treein2.GetEntriesFast()):
  if ie%1000==0:print('processing ',ie)
  treein2.GetEntry(ie)
  mu_arr=[]
  if treein2.nMuon<1:continue
  for im in range(0,treein2.nMuon):
    if not treein2.Muon_mediumId[im]:continue
    if treein2.Muon_pt[im]<20:continue
    mu_arr.append(im)
  if not len(mu_arr)==2:continue
  l1.SetPtEtaPhiM(treein2.Muon_pt[mu_arr[0]],treein2.Muon_eta[mu_arr[0]],treein2.Muon_phi[mu_arr[0]],treein2.Muon_mass[mu_arr[0]])
  l2.SetPtEtaPhiM(treein2.Muon_pt[mu_arr[1]],treein2.Muon_eta[mu_arr[1]],treein2.Muon_phi[mu_arr[1]],treein2.Muon_mass[mu_arr[1]])
  h2.Fill((l1+l2).M())

h1.Scale(powheg_weight)
h2.Scale(mg_weight)

c1=ROOT.TCanvas()
c1.cd()
h1.SetLineColor(ROOT.kRed)
h2.SetLineColor(ROOT.kBlue)
h1.Draw()
h2.Draw('same')
c1.SaveAs('aa.png')
c1.SaveAs('aa.root')


