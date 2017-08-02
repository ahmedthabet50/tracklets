void plotTracklets(TString data = "data.root", TString mc = "mc.root")
{

  TFile *file1 = new TFile(data,"read"); // Data
  TFile *file2 = new TFile(mc,"read"); //MC

  TList *MyOutputContainer1 = (TList *)file1->Get("Hist;1");
  TList *MyOutputContainer2 = (TList *)file2->Get("Hist;1");

  /***********************************************************************************/
  //  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  //=========================


  TH1F *fMultData_reco = (TH1F *)MyOutputContainer1->FindObject("fHistMult_reco");
  TH1F *fDataRecoZvtx = (TH1F *)MyOutputContainer1->FindObject("fHistRecoZvtx");
  TH1F *fMult_reco = (TH1F *)MyOutputContainer2->FindObject("fHistMult_reco");
  TH1F *fMult_prim = (TH1F *)MyOutputContainer2->FindObject("fHistMult_prim");
  TH1F *fMult_seco = (TH1F *)MyOutputContainer2->FindObject("fHistMult_seco");
  TH1F *fMult_fake = (TH1F *)MyOutputContainer2->FindObject("fHistMult_fake");
  TH1F *fRecoZvtx = (TH1F *)MyOutputContainer2->FindObject("fHistRecoZvtx");
  TH1F *fMCZvtx = (TH1F *)MyOutputContainer2->FindObject("fHistGenZvtx");


  TTree *trd = (TTree*)file1->Get("fTree_reco");
  TTree *tr = (TTree*)file2->Get("fTree_reco");
  TTree *tp = (TTree*)file2->Get("fTree_prim");
  TTree *ts = (TTree*)file2->Get("fTree_seco");
  TTree *tf = (TTree*)file2->Get("fTree_fake");


  if (!trd) return;

  TCanvas* cetaphi = new TCanvas("cetaphi", "cetaphi", 1100, 1100);
  cetaphi->Divide(2,2,0.01,0.01,0);

  cetaphi->cd(1);

  trd->Draw("fPhi_reco>>hphi");
  TH1F *hphi = (TH1F*)gPad->GetPrimitive("hphi");
  hphi->GetXaxis()->SetTitle("#phi (Reconstructed)");
  hphi->GetYaxis()->SetTitle("tracklets");

  cetaphi->cd(2);

  trd->Draw("fEta_reco>>heta(80, -2.0, 2.0)");
  TH1F *heta = (TH1F*)gPad->GetPrimitive("heta");
  heta->GetXaxis()->SetTitle("#eta (Reconstructed)");
  heta->GetYaxis()->SetTitle("tracklets");


  cetaphi->cd(3);

  fMultData_reco->Draw("p");

  cetaphi->cd(4);

  fDataRecoZvtx->Draw("p");

  cetaphi->SaveAs("reco_etaphimult.pdf");


  TCanvas* cpt = new TCanvas("cpt", "cpt", 1100, 1100);
  cpt->Divide(2,2,0.01,0.01,0);
  cpt_1->SetLogy();
  cpt_2->SetLogy();
  cpt_3->SetLogy();


  cpt->cd(1);


  tp->Draw("fPt_prim>>hptp");
  TH1F *hptp = (TH1F*)gPad->GetPrimitive("hptp");
  hptp->GetXaxis()->SetTitle("pT (primary tracklets)");
  hptp->GetYaxis()->SetTitle("tracklets");

  cpt->cd(2);

  ts->Draw("fPt_seco>>hpts");
  TH1F *hpts = (TH1F*)gPad->GetPrimitive("hpts");
  hpts->GetXaxis()->SetTitle("pT (Secondary  tracklets)");
  hpts->GetYaxis()->SetTitle("tracklets");

  cpt->cd(3);

  tf->Draw("fPt_fake>>hptf");
  TH1F *hptf = (TH1F*)gPad->GetPrimitive("hptf");
  hptf->GetXaxis()->SetTitle("pT (Combinotorial(Fake) tracklets)");
  hptf->GetYaxis()->SetTitle("tracklets");

  cpt->SaveAs("pT.pdf");


  TCanvas* corr = new TCanvas("corr", "corr", 1100, 1100);
  corr->Divide(2,2,0.01,0.01,0);
  corr_1->SetLogy();
  corr_2->SetLogy();
  corr_3->SetLogy();
  corr_4->SetLogy();


  corr->cd(1);

  tr->Draw("fDPhi_reco>>hdphir");
  TH1F *hdphi = (TH1F*)gPad->GetPrimitive("hdphir");
  hdphir->GetXaxis()->SetTitle("#Delta#phi (Reconstructed)");
  hdphir->GetYaxis()->SetTitle("Reconstructed tracklets");

  corr->cd(2);

  tr->Draw("fDTheta_reco>>hdthetar");
  TH1F *hdtheta = (TH1F*)gPad->GetPrimitive("hdthetar");
  hdthetar->GetXaxis()->SetTitle("#Delta#theta (Reconstructed)");
  hdthetar->GetYaxis()->SetTitle("Reconstructed tracklets");

  corr->cd(3);


  tp->Draw("fDPhi_prim>>hdphip");
  TH1F *hdphip = (TH1F*)gPad->GetPrimitive("hdphip");
  hdphip->GetXaxis()->SetTitle("#Delta#phi (primary)");
  hdphip->GetYaxis()->SetTitle("Reconstructed tracklets");


  corr->cd(4);

  tp->Draw("fDTheta_prim>>hthetap");
  TH1F *hthetap = (TH1F*)gPad->GetPrimitive("hthetap");
  hthetap->GetXaxis()->SetTitle("#Delta#theta (primary)");
  hthetap->GetYaxis()->SetTitle("Reconstructed tracklets");


  corr->SaveAs("DeltaEtaDeltaPhi.pdf");


  TCanvas* corr1 = new TCanvas("corr1", "corr1", 1100, 1100);
  corr1->Divide(2,2,0.01,0.01,0);
  corr1_4->SetLogy();


  corr1->cd(1);

  tp->Draw("fDTheta_prim:fDPhi_prim>>hetaphip(200,-0.1,0.1,200,-0.03,0.03)","","colz");
  TH2D *hetaphip = (TH2D*)gPad->GetPrimitive("hetaphip");
  hetaphip->GetXaxis()->SetTitle("#Delta#theta (primary)");
  hetaphip->GetYaxis()->SetTitle("#Delta#phi (primary)");

  corr1->cd(2);

  tr->Draw("fDTheta_reco:fDPhi_reco>>hetaphir(200,-0.1,0.1,200,-0.03,0.03)","","colz");
  TH2D *hetaphip = (TH2D*)gPad->GetPrimitive("hetaphir");
  hetaphir->GetXaxis()->SetTitle("#Delta#theta (Reconstructed)");
  hetaphir->GetYaxis()->SetTitle("#Delta#phi (Reconstructed)");

  corr1->cd(3);

  tp->Draw("fPt_prim:1/(1000*TMath::Abs(fDPhi_prim))>>hetaphi(100,0,5,100,0,5)","","colz");
  TH2D *hetaphi = (TH2D*)gPad->GetPrimitive("hetaphi");
  hetaphi->GetXaxis()->SetTitle("pT (primary)");
  hetaphi->GetYaxis()->SetTitle("1/|#phi|(mrad) (primary)");

  corr1->cd(4);

  tp->Draw("1000*TMath::Abs(fDPhi_prim)>>hdphi1");
  TH1F *hdphi1 = (TH1F*)gPad->GetPrimitive("hdphi1");
  hdphi1->GetXaxis()->SetTitle("|#Delta#phi| (mrad) (primary)");
  hdphi1->GetYaxis()->SetTitle("Reconstructed tracklets");


  corr1->SaveAs("correlation.pdf");


  TCanvas* cut1 = new TCanvas("cut1", "cut1", 900, 800);
  cut1->Divide(1,1,0.01,0.01,0);

  tp->Draw("fPt_prim>>hptcut1(200,0,5)","(1000*TMath::Abs(fDPhi_prim))<100");
  //TH1F *hptcut1 = (TH1F*)gPad->GetPrimitive("hptcut1");
  hptcut1->GetXaxis()->SetTitle("pT (primary)");
  hptcut1->GetYaxis()->SetTitle("Reconstructed tracklets");

  TH1F *hframe = (TH1F*)hptcut1->Clone("hframe");
  //  hptcut1->Scale(1/500);

  hframe->Draw();
  ((TH1F*)(gPad->GetListOfPrimitives()->At(0)))->SetLineColor(1); //

  tp->Draw("fPt_prim>>hptcut2","(1000*TMath::Abs(fDPhi_prim))<10","same");
  ((TH1F*)(gPad->GetListOfPrimitives()->At(1)))->SetLineColor(2); //

  tp->Draw("fPt_prim>>hptcut3","(1000*TMath::Abs(fDPhi_prim))<5","same");
  ((TH1F*)(gPad->GetListOfPrimitives()->At(2)))->SetLineColor(8); //

  tp->Draw("fPt_prim>>hptcut4","(1000*TMath::Abs(fDPhi_prim))<2","same");
  ((TH1F*)(gPad->GetListOfPrimitives()->At(3)))->SetLineColor(4); //

  tp->Draw("fPt_prim>>hptcut5","(1000*TMath::Abs(fDPhi_prim))<1","same");
  ((TH1F*)(gPad->GetListOfPrimitives()->At(4)))->SetLineColor(6); //

  tp->Draw("fPt_prim>>hptcut6","(1000*TMath::Abs(fDPhi_prim))<0.1","same");
  ((TH1F*)(gPad->GetListOfPrimitives()->At(5)))->SetLineColor(9); //

  TLegend *legend = new TLegend(0.42, 0.1+0.63, 0.7, 0.27+0.63);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(0.03);
  legend->AddEntry(hptcut1, "|#Delta#phi| < 100 , <pT> ~ 0.50", "l");
  legend->AddEntry(hptcut2, "|#Delta#phi| < 10, <pT> ~ 0.56", "l");
  legend->AddEntry(hptcut3, "|#Delta#phi| < 5, <pT> ~ 0.71", "l");
  legend->AddEntry(hptcut4, "|#Delta#phi| < 2, <pT> ~ 0.90", "l");
  legend->AddEntry(hptcut5, "|#Delta#phi| < 1, <pT> ~ 0.94", "l");
  legend->AddEntry(hptcut6, "|#Delta#phi| < 0.1, <pT> ~ 1.1", "l");
  legend->Draw("same");


  gPad->Modified(); gPad->Update(); // make sure the pad is redrawn

  cut1->SaveAs("pT_deltaphicut.pdf");


}
