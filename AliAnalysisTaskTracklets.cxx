/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
* Author: Prabi and Paolo Bartalini, July CCNU 2017               *
**************************************************************************/

/* AliAnalysisTaskTracklets source code
*
* simple task which can serve as a starting point for building an SPD tracklets analysis
*
*/

class TTree;
class AliPPVsMultUtils;
class AliESDtrackCuts;


#include <Riostream.h>
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"


#include "TVector3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TList.h"

#include "AliLog.h"
#include "AliVEvent.h"
#include "AliVVertex.h"
#include "AliVTrack.h"
//#include "AliV0vertexer.h"


#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDVZERO.h"

#include "AliESDInputHandler.h"
#include "AliMultiplicity.h"
#include "AliPPVsMultUtils.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskESDfilter.h"

#include "AliMCEventHandler.h"
#include "AliInputEventHandler.h"

#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include <AliHeader.h>

#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliMultInput.h"
#include "AliMultVariable.h"
#include "AliCentrality.h"
#include "AliMultiplicity.h"

#include "AliESDUtils.h"

#include "AliGenEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenPythiaEventHeader.h"


#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliPPVsMultUtils.h"
#include <AliESDVertex.h>
#include <AliMultiplicity.h>
#include <TTree.h>
#include <TMath.h>
#include <TDirectory.h>
#include <TBits.h>
#include <AliAnalysisFilter.h>

using std::cout;
using std::endl;

#include "AliAnalysisTaskTracklets.h"

class AliAnalysisTaskTracklets;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout' etc

ClassImp(AliAnalysisTaskTracklets) // classimp: necessary for root

AliAnalysisTaskTracklets::AliAnalysisTaskTracklets() : AliAnalysisTaskSE(),
fESD(0), fStack(0),fMCEvent(0), fOutputList(0),fWeight(1.),fEventCounter(0),fCuts (0), fUseMC(kFALSE),
fHistMultiplicityStd(0),fHistGenZvtx(0),fHistRecoZvtx(0),fPt_prim(0),fPt_seco(0),fPt_fake(0),fEta_prim(0),fPhi_prim(0),fDTheta_prim(0),fDPhi_prim(0),fPhi_seco(0),fPhi_fake(0),fDist_reco(0),fEta_reco(0),fDTheta_reco(0),fPhi_reco(0),fDPhi_reco(0),fHistMult_reco(0),fHistMult_prim(0),fHistMult_seco(0),fHistMult_fake(0)
{
  // default constructor, don't allocate memory here!  this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskTracklets::AliAnalysisTaskTracklets(const char* name) : AliAnalysisTaskSE(name),
fESD(0), fStack(0),fMCEvent(0), fOutputList(0),fWeight(1.),fEventCounter(0),fCuts (0), fUseMC(kFALSE),
fHistMultiplicityStd(0), fHistGenZvtx(0),fHistRecoZvtx(0), fPt_prim(0),fPt_seco(0),fPt_fake(0),fEta_prim(0),fPhi_prim(0),fDTheta_prim(0),fDPhi_prim(0),fPhi_seco(0),fPhi_fake(0),fDist_reco(0),fEta_reco(0),fDTheta_reco(0),fPhi_reco(0),fDPhi_reco(0),fHistMult_reco(0),fHistMult_prim(0),fHistMult_seco(0),fHistMult_fake(0)
{

  // constructor
  DefineInput(0, TChain::Class());    // define the input of the analysis: in this case you take a 'chain' of events
  // this chain is created by the analysis manager, so no need to worry about it, does its work automatically
  DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms
  DefineOutput(2, TTree::Class());    // define the ouptut of the analysis: in this case it's a list of histograms
  DefineOutput(3, TTree::Class());    // define the ouptut of the analysis: in this case it's a list of histograms
  DefineOutput(4, TTree::Class());    // define the ouptut of the analysis: in this case it's a list of histograms
  DefineOutput(5, TTree::Class());    // define the ouptut of the analysis: in this case it's a list of histograms

}
//_____________________________________________________________________________
AliAnalysisTaskTracklets::~AliAnalysisTaskTracklets()
{
  // destructor
  if(fOutputList) {
    delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    fOutputList = 0x0;
  }
  if (fTree_reco) {
    delete fTree_reco;
    fTree_reco = 0x0;
  }
  if (fTree_prim) {
    delete fTree_prim;
    fTree_prim = 0x0;
  }
  if (fTree_seco) {
    delete fTree_seco;
    fTree_seco = 0x0;
  }
  if (fTree_fake) {
    delete fTree_fake;
    fTree_fake = 0x0;
  }
  if(fCuts) {
    delete fCuts;
    fCuts = 0x0;
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskTracklets::UserCreateOutputObjects()
{

  AliAnalysisManager *mgr=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* handler = (AliInputEventHandler*) (mgr->GetInputEventHandler());
  handler->SetNeedField();

  fCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();

  fCuts->SetMinNClustersTPC(70);
  fCuts->SetMaxChi2PerClusterTPC(4);
  fCuts->SetMaxDCAToVertexXY(2.4);   // Alternatively use Gaussian fit to DCAxy distribution and use 7*sigma
  //fCuts->SetMaxDCAToVertexZ(3.2);
  fCuts->SetMaxDCAToVertexZ(2);
  fCuts->SetDCAToVertex2D(kTRUE);
  fCuts->SetMaxChi2TPCConstrainedGlobal(36);
  fCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny); // reduces secondaries by having atleast 1 spd hit

  fCuts->SetMaxFractionSharedTPCClusters(0.4);

  //fCuts->SetPtRange(0.15);
  fCuts->SetPtRange(0.5);
  fCuts->SetEtaRange(-0.8, 0.8);

  // create output objects

  OpenFile(1);
  fOutputList = new TList();          // this is a list which will contain all of your histograms
  // at the end of the analysis, the contents of this list are written  to the output file
  fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects and will delete them if requested

  if(! fEventCounter ){
    fEventCounter = new TH1I( "fEventCounter", ";Evt. Sel. Step;Count",16,0,16);
    fEventCounter->GetXaxis()->SetBinLabel(1, "Total processed Events");
    fEventCounter->GetXaxis()->SetBinLabel(2, "Events selected for Analysis");
    fEventCounter->GetXaxis()->SetBinLabel(3, "Pileup events in SPDInMultBins");
    fEventCounter->GetXaxis()->SetBinLabel(4, "Events not with INEL>0");
    fEventCounter->GetXaxis()->SetBinLabel(5, "Events not in Vertex cut");
    fEventCounter->GetXaxis()->SetBinLabel(6, "Events inconsistent with SPD&TrackVertx");
    fEventCounter->GetXaxis()->SetBinLabel(7, "Events without MinimumBias (kMB)");
    fEventCounter->GetXaxis()->SetBinLabel(8, "ClusterVsTracklet Cut rejected events");
    fEventCounter->GetXaxis()->SetBinLabel(9, "Out of Bunch events in 11 IR");
    fEventCounter->GetXaxis()->SetBinLabel(10, "Incomplete DAQ events");
    fEventCounter->GetXaxis()->SetBinLabel(11, "Events Rejected by events quality criteria");
    fEventCounter->GetXaxis()->SetBinLabel(12, "Events Rejected by V0 decision");
    fEventCounter->GetXaxis()->SetBinLabel(13, "Events Rejected by V0 Assymetry");
    fEventCounter->GetXaxis()->SetBinLabel(14, "Events Rejected by Fastor Online cut");
    fEventCounter->GetXaxis()->SetBinLabel(15, "Events Without BCMod4 == 2 ");
    fEventCounter->GetXaxis()->SetBinLabel(16, "Events with No Good Tracks");

    fOutputList->Add(fEventCounter);
  }

  // your histogram and TTree branches in the output file, add it to the list!

  fHistMultiplicityStd= new TH1I("fHistMultiplicityStd", "Standard Reference Multiplicity;Std Ref. Multiplicity  ; P(Nch) ",300,-0.5,299.5);       // create Reference Multiplicity histogram
  fHistGenZvtx = new TH1F("fHistGenZvtx", "GenZvtx; GenZvtx ; Events", 160, -20, 20);      // Generated Vertex Z
  fHistRecoZvtx = new TH1F("fHistRecoZvtx", "RecoZvtx; RecoZvtx ; Events", 160, -20, 20);      // Reconstructed Vertex Z
  fHistMult_reco= new TH1I("fHistMult_reco", "Reconstructed Tracklets Multiplicity;Multiplicity (Reconstructed)  ; P(Nch) ",300,-0.5,299.5);       // create Reference Multiplicity histogram
  fHistMult_prim= new TH1I("fHistMult_prim", "Primary Tracklets Multiplicity;Multiplicity (Primary) ; P(Nch) ",300,-0.5,299.5);       // create Reference Multiplicity histogram
  fHistMult_seco= new TH1I("fHistMult_seco", "Secondary Tracklets Multiplicity;Multiplicity  (Secondary); P(Nch) ",300,-0.5,299.5);       // create Reference Multiplicity histogram
  fHistMult_fake= new TH1I("fHistMult_fake", "Fakes Tracklets Multiplicity;Multiplicity  (Fakes); P(Nch) ",300,-0.5,299.5);       // create Reference Multiplicity histogram



  fTree_reco = new TTree("fTree_reco","fTree_reco");
  fTree_reco->Branch("fDist_reco",&fDist_reco,"fDist_reco/F");
  fTree_reco->Branch("fEta_reco",&fEta_reco,"fEta_reco/F");
  fTree_reco->Branch("fDTheta_reco",&fDTheta_reco,"fDTheta_reco/F");
  fTree_reco->Branch("fPhi_reco",&fPhi_reco,"fPhi_reco/F");
  fTree_reco->Branch("fDPhi_reco",&fDPhi_reco,"fDPhi_reco/F");


  fTree_prim = new TTree("fTree_prim","fTree_prim");
  fTree_prim->Branch("fPt_prim",&fPt_prim,"fPt_prim/F");
  fTree_prim->Branch("fEta_prim",&fEta_prim,"fEta_prim/F");
  fTree_prim->Branch("fPhi_prim",&fPhi_prim,"fPhi_prim/F");
  fTree_prim->Branch("fDTheta_prim",&fDTheta_prim,"fDTheta_prim/F");
  fTree_prim->Branch("fDPhi_prim",&fDPhi_prim,"fDPhi_prim/F");



  fTree_seco = new TTree("fTree_seco","fTree_seco");
  fTree_seco->Branch("fPt_seco",&fPt_seco,"fPt_seco/F");
  fTree_seco->Branch("fPhi_seco",&fPhi_seco,"fPhi_seco/F");


  fTree_fake = new TTree("fTree_fake","fTree_fake");
  fTree_fake->Branch("fPt_fake",&fPt_fake,"fPt_fake/F");
  fTree_fake->Branch("fPhi_fake",&fPhi_fake,"fPhi_fake/F");


  // your histogram in the output file, add it to the list!


  fHistMultiplicityStd->Sumw2();
  fHistGenZvtx->Sumw2();
  fHistRecoZvtx->Sumw2();
  fHistMult_reco->Sumw2();
  fHistMult_prim->Sumw2();
  fHistMult_seco->Sumw2();
  fHistMult_fake->Sumw2();


  fOutputList->Add(fHistMultiplicityStd);
  fOutputList->Add(fHistMult_reco);
  fOutputList->Add(fHistMult_prim);
  fOutputList->Add(fHistMult_seco);
  fOutputList->Add(fHistMult_fake);

  fOutputList->Add(fHistGenZvtx);
  fOutputList->Add(fHistRecoZvtx);



  PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the
  // fOutputList object. the manager will in the end take care of writing your output to file so it needs to know what's in the output

  PostData(2, fTree_prim);   //
  PostData(3, fTree_seco);   //
  PostData(4, fTree_fake);   //
  PostData(5, fTree_reco);   //

}
//_____________________________________________________________________________
void AliAnalysisTaskTracklets::UserExec(Option_t *)
{

  {
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if(!fESD) {
      AliWarning("ERROR: ESD event is not available \n");
      return;                                                 // if the pointer to the event is empty (getting it failed) skip this event
    }

    fEventCounter->Fill(0);

    ///< processed events counter
    if(HasNoPileupSPDInMultBins(fESD) == kFALSE){fEventCounter->Fill(2);}                  ///<  Rejects events tagged with IsPileupFromSPDInMultBins()
    if(IsINELgtZERO( fESD ) == kFALSE){fEventCounter->Fill(3);}                            ///<  INEL > 0 (with SPD tracklets atleast 1)
    if(HasAcceptedVertexPosition( fESD ) == kFALSE){fEventCounter->Fill(4);}               ///<  Checks for accepted vertex position (|eta|<10cm)
    if(HasNoInconsistentSPDandTrackVertices( fESD ) == kFALSE){fEventCounter->Fill(5);}    ///<  Checks for consistent SPD and track vertex (if track vertex exists)
    if(IsMinimumBias( fESD ) == kFALSE){fEventCounter->Fill(6);}                           ///<  kMB minimum bias trigger selection
    if(IsSPDClusterVsTrackletBG( fESD ) == kFALSE){fEventCounter->Fill(7);}                ///<  counts rejected events
    if(IsOutOfBunchPileup( fESD ) == kTRUE){fEventCounter->Fill(8);}                       ///<  counts rejected events
    if(IsInCompleteEvent( fESD ) == kTRUE){fEventCounter->Fill(9);}                        ///<  counts rejected events
    if(IsEventSelected( fESD ) == kFALSE){fEventCounter->Fill(10);}                        ///<  counts rejected events
    if(V0Decision( fESD ) == kFALSE){fEventCounter->Fill(11);}                      ///<  counts rejected events
    if(V0Asymmetry( fESD ) == kFALSE){fEventCounter->Fill(12);}                     ///<  counts rejected events
    if(IsOutOfBunchPileupFastor( fESD ) == kFALSE){fEventCounter->Fill(13);}        ///<  counts rejected events
    if(HasBCMod4( fESD ) == kFALSE){fEventCounter->Fill(14);}                       ///<  counts WITHOUT BCMOD4==2
    if(fCuts->CountAcceptedTracks(fESD)>0 == kFALSE){fEventCounter->Fill(15);}                       ///<  counts WITHOUT good any tracks

    if(!IsEventSelected(fESD) ) {
      PostData(1, fOutputList);   /// Event isn't selected, post output data, done here
      PostData(2, fTree_prim);   //
      PostData(3, fTree_seco);   //
      PostData(4, fTree_fake);   //
      PostData(5, fTree_reco);   //
      return;
    }

    fEventCounter->Fill(1);

    fHistMultiplicityStd->Fill(AliPPVsMultUtils::GetStandardReferenceMultiplicity(fESD,kFALSE),fWeight); // this will give you standard reference multiplicity

    /*
    if (fCuts->CountAcceptedTracks(fESD)>0) {
    TrackletsLoop(fESD);
  } else {
  //    AliWarning("<<<< 0 Tracks After All Cuts >>>>>>>");
  return;
}
*/
if (fUseMC){
  FillMCGen();  // Fill all the Generator level information
}

TrackletsLoop(fESD);

PostData(1, fOutputList); // stream the result of this event to the output manager which will write it to a file
PostData(2, fTree_prim);   //
PostData(3, fTree_seco);   //
PostData(4, fTree_fake);   //
PostData(5, fTree_reco);   //

}
}
//_____________________________________________________________________________
void AliAnalysisTaskTracklets::FillMCGen() {

  AliAnalysisManager *anmgr=AliAnalysisManager::GetAnalysisManager();
  AliMCEventHandler* eventHandler = 0;
  fMCEvent = 0;
  fStack = 0;
  eventHandler = (AliMCEventHandler*)anmgr->GetMCtruthEventHandler();
  if (!eventHandler) { printf("ERROR: Could not retrieve MC event handler\n"); return; }
  fMCEvent = eventHandler->MCEvent();
  if (!fMCEvent) { printf("ERROR: Could not retrieve MC event\n"); return; }
  fStack = fMCEvent->Stack();
  if (!fStack) { printf("Stack not available\n"); return; }


  AliGenEventHeader* mcGenH = 0;
  AliGenPythiaEventHeader* hPythia=0;

  mcGenH = fMCEvent->GenEventHeader();
  Float_t event_weight = fMCEvent->GenEventHeader()->EventWeight();

  if (mcGenH->InheritsFrom(AliGenCocktailEventHeader::Class())) {
    TList* headers = ((AliGenCocktailEventHeader*)mcGenH)->GetHeaders();
    TIter next(headers);
    AliGenEventHeader* mcGenH1 = 0;
    while ( (mcGenH1=(AliGenEventHeader*)next()) ) {
      if (mcGenH1->InheritsFrom(AliGenPythiaEventHeader::Class())) {
        hPythia = (AliGenPythiaEventHeader*)mcGenH1;
        break;
      }
    }
  }

  else {} // unknown generator

  TArrayF vtmc(3);
  mcGenH->PrimaryVertex(vtmc);
  for (int i=3;i--;) fVertexMC[i] = vtmc[i];
  fHistGenZvtx->Fill(fVertexMC[2],event_weight);

  const AliMultiplicity* mlt = fESD->GetMultiplicity();
  if (!mlt) {
    AliDebug(AliLog::kError, "AliMultiplicity not available");
    return;
  }

  int ntr = mlt->GetNumberOfTracklets();

  Float_t prim = 0;
  Float_t seco = 0;
  Float_t fake = 0;



  for (int itr=ntr;itr--;) {

    Float_t theta  = mlt->GetTheta(itr);
    Float_t phi  = mlt->GetPhi(itr);
    Float_t eta    = -TMath::Log(TMath::Tan(theta/2));
    //
    Float_t dtheta = mlt->GetDeltaTheta(itr);
    Float_t dphi   = mlt->GetDeltaPhi(itr);
    Float_t dist   = mlt->CalcDist(itr);

    if (fMCEvent) {

      int label0 = mlt->GetLabel(itr,0);
      int label1 = mlt->GetLabel(itr,1);
      int typeMC = 2; // comb.bg.
      if (label0 == label1) typeMC = fMCEvent->IsPhysicalPrimary(label0) ? 0:1; // prim or sec

    //  TParticle* part = fMCEvent->Particle(label0);
    AliMCParticle* part = (AliMCParticle*)fMCEvent->GetTrack(label0);

      if (part->Charge() == 0) continue;

      Float_t pt_tracklets = part->Pt();
      Float_t eta_tracklets = part->Eta();
      Float_t phi_tracklets = part->Phi();


      if (typeMC==0) {
        fPt_prim=pt_tracklets;
        fPhi_prim=phi_tracklets;
        fEta_prim=eta_tracklets;
        fDTheta_prim=dtheta;
        fDPhi_prim=dphi;
        prim++;
        fTree_prim->Fill();
        //  printf(" pt of primary tracklets   =======================     %f \n", pt_tracklets);
      }
      else if (typeMC==1) {
        fPt_seco=pt_tracklets;
        fPhi_seco=dphi;
        seco++;
        fTree_seco->Fill();
      }
      else {
        fPt_fake=pt_tracklets;
        fPhi_fake=dphi;
        fake++;
        fTree_fake->Fill();

      } //
    }
  }
if (prim > 0 )   fHistMult_prim->Fill(prim,fWeight);
if (seco > 0 )   fHistMult_seco->Fill(seco,fWeight);
if (fake > 0 )   fHistMult_fake->Fill(fake,fWeight);

}
//_____________________________________________________________________________________
void AliAnalysisTaskTracklets::Terminate(Option_t *)
{
  // terminate
  // called at the END of the analysis (when all events are processed)
}

//_____________________________________________________________________________
void AliAnalysisTaskTracklets::TrackletsLoop(AliESDEvent *fESD)
{

  const AliMultiplicity* mltr = fESD->GetMultiplicity();
  if (!mltr) {
    AliDebug(AliLog::kError, "AliMultiplicity not available");
    return;
  }


  for (int i=3;i--;) fVertexReco[i] = 0;

  const AliESDVertex* vtxSPD = fESD->GetPrimaryVertexSPD();
  const AliESDVertex* vtxTracks = fESD->GetPrimaryVertex();
  if ( !vtxSPD->IsFromVertexerZ() || (vtxSPD->GetDispersion()<0.02)) {
    fVertexReco[0] = vtxTracks->GetX(); // X from global tracks
    fVertexReco[1] = vtxTracks->GetY(); // Y
    fVertexReco[2] = vtxTracks->GetZ(); // Z
  }

  fHistRecoZvtx->Fill(fVertexReco[2],fWeight);



  int ntrk = mltr->GetNumberOfTracklets();
  Float_t reco = 0;


  for (int itrk=ntrk;itrk--;) {

    Float_t thetar  = mltr->GetTheta(itrk);
    Float_t phir  = mltr->GetPhi(itrk);
    Float_t etar    = -TMath::Log(TMath::Tan(thetar/2));
    //
    Float_t dthetar = mltr->GetDeltaTheta(itrk);
    Float_t dphir   = mltr->GetDeltaPhi(itrk);
    Float_t distr   = mltr->CalcDist(itrk);

    //if (dphi > TMath::Pi()) dphi = 2. * TMath::Pi() - dphi;

    fDist_reco=distr;
    fEta_reco=etar;
    fDTheta_reco=dthetar;
    fPhi_reco=phir;
    fDPhi_reco=dphir;
    reco++;
    fTree_reco->Fill();
  }

  fHistMult_reco->Fill(reco,fWeight);

}
//--------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskTracklets::IsMinimumBias(AliVEvent *fESD)
// Function to check for minimum-bias trigger (AliVEvent::kMB)
{
  // to reject events that aren't kMB/INT7/INT1
  UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  Bool_t IsTriggered = kFALSE;

  //    IsTriggered = (maskIsSelected & AliVEvent::kHighMultV0) == AliVEvent::kHighMultV0;   // CINT7 trigger
  //        IsTriggered = ((maskIsSelected & AliVEvent::kINT7) || (maskIsSelected & AliVEvent::kHighMultV0));   // CINT7 trigger
  IsTriggered = (maskIsSelected & AliVEvent::kINT7) == AliVEvent::kINT7;  // CINT7 trigger
  //  IsTriggered = (maskIsSelected & AliVEvent::kINT1) == AliVEvent::kINT1;  // CINT1 trigger
  //    IsTriggered = (maskIsSelected & AliVEvent::kINT10) == AliVEvent::kINT10;  // CINT10 trigger


  return IsTriggered;
}
//--------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskTracklets::IsINELgtZERO(AliVEvent *fESD)
// Function to check for INEL > 0 condition, need atleast one SPD tracklet
{
  Bool_t returnValue = kFALSE;
  AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(fESD);
  if (!esdevent) return kFALSE;
  if ( AliESDtrackCuts::GetReferenceMultiplicity(esdevent, AliESDtrackCuts::kTracklets, 1.0) >= 1 ) returnValue = kTRUE;

  return returnValue;
}

//--------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskTracklets::HasNoPileupSPDInMultBins(AliVEvent *fESD)
// Checks if No pileup from SPD (via IsPileupFromSPDInMultBins)
{
  Bool_t returnValue = kTRUE;

  AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(fESD);
  if (!esdevent) return kFALSE;
  if ( esdevent->IsPileupFromSPDInMultBins() == kTRUE ) returnValue = kFALSE;
  return returnValue;
}
//--------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskTracklets::HasNoInconsistentSPDandTrackVertices(AliVEvent *fESD)
{
  //Accepts events which has consistent SPD and Global tracks Z Vertices agreementent within 5 mm
  Bool_t returnValue = kTRUE;

  AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(fESD);
  if (!esdevent) return kFALSE;
  const AliESDVertex *lPrimaryVtxSPD    = NULL;
  const AliESDVertex *lPrimaryVtxTracks = NULL;

  lPrimaryVtxSPD    = esdevent->GetPrimaryVertexSPD   ();
  lPrimaryVtxTracks = esdevent->GetPrimaryVertexTracks();

  //Only continue if track vertex defined
  if( lPrimaryVtxTracks->GetStatus() && lPrimaryVtxSPD->GetStatus() ){
    //Copy-paste from refmult estimator
    // TODO value of displacement to be studied
    const Float_t maxDisplacement = 0.5;
    //check for displaced vertices
    Float_t displacement = TMath::Abs(lPrimaryVtxSPD->GetZ() - lPrimaryVtxTracks->GetZ());
    if (displacement > maxDisplacement) returnValue = kFALSE;
  }

  return returnValue;
}
//--------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskTracklets::HasAcceptedVertexPosition(AliVEvent *fESD)
// Will accept events only if best primary vertex Z position awailable with |z| < 10cm
{
  Bool_t returnValue = kFALSE;
  const AliVVertex *lPrimaryVtx = NULL;

  AliESDEvent *esdevent = dynamic_cast<AliESDEvent *>(fESD);
  if (!esdevent) return kFALSE;
  lPrimaryVtx = esdevent->GetPrimaryVertex();

  if ( TMath::Abs( lPrimaryVtx->GetZ() ) <= 10.0 ) returnValue = kTRUE;
  return returnValue;
}
//----------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskTracklets::IsSPDClusterVsTrackletBG(AliVEvent* fESD)
{
  // rejects BG background based on the cluster vs tracklet correlation
  // returns true if the event is BG slope is 4 and intercept is 65
  const AliVMultiplicity* mult = fESD->GetMultiplicity();
  // if (!mult) { AliFatal("No multiplicity object"); return 0;}
  Int_t ntracklet   = mult->GetNumberOfTracklets();
  Int_t spdClusters = fESD->GetNumberOfITSClusters(0) + fESD->GetNumberOfITSClusters(1);
  //return spdClusters < Float_t(fASPDCvsTCut) + Float_t(ntracklet)*fBSPDCvsTCut;
  return spdClusters < (65.0 + (ntracklet)*4.0);

}
//----------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskTracklets::IsOutOfBunchPileupFastor(AliVEvent* fESD)
{
  Float_t SPD_Online = (fESD->GetMultiplicity()->GetFastOrFiredChipMap().CountBits(400)+ fESD->GetMultiplicity()->GetFastOrFiredChipMap().CountBits(800));
  Float_t SPD_Offline = (fESD->GetMultiplicity()->GetFiredChipMap().CountBits(400) + fESD->GetMultiplicity()->GetFiredChipMap().CountBits(800));

  return SPD_Online >= (-20.589 + 0.73664*(SPD_Offline));
}
//----------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskTracklets::IsOutOfBunchPileup(AliVEvent* fESD)
{
  // rejects Outofbunch Pileup events based on the Interaction logic

  TBits fIR11 = fESD->GetHeader()->GetIRInt1InteractionMap();
  //TBits fIR22 = fESD->GetHeader()->GetIRInt2InteractionMap();

  Bool_t isOutOfBunchPileup = 0;
  for (Int_t i=1;i<=3;i++) { isOutOfBunchPileup|=fIR11.TestBitNumber(90-i);}
  for (Int_t i=8;i<=11;i++) { isOutOfBunchPileup|=fIR11.TestBitNumber(90+i);}

  //if(isOutOfBunchPileup == 1) return;

  //for (Int_t i=0;i<=7;i++) isOutOfBunchPileup|=fIR22.TestBitNumber(90-i);
  //for (Int_t i=4;i<=7;i++) isOutOfBunchPileup|=fIR22.TestBitNumber(90+i);

  return isOutOfBunchPileup;

}
//-----------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskTracklets::IsInCompleteEvent(AliVEvent* fESD)
{
  // rejects Incomplete events in the run
  Bool_t isIncomplete = kFALSE;
  if (fESD->IsIncompleteDAQ()) isIncomplete = kTRUE;

  return isIncomplete;

}
//-----------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskTracklets::V0Decision(AliVEvent* fESD)
{
  // rejects Incomplete events in the run
  AliVVZERO* vzero_decision = fESD->GetVZEROData();

  Bool_t isV0Decision = kFALSE;
  isV0Decision = ((vzero_decision->GetV0ADecision()==1) && (vzero_decision->GetV0CDecision()==1));

  return isV0Decision;
}
//-----------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskTracklets::V0Asymmetry(AliVEvent* fESD)
{
  // rejects Incomplete events in the run
  Bool_t isEventSelected_V0Asym = kTRUE;

  AliVVZERO* vzero_asym = fESD->GetVZEROData();

  Float_t v0c012 = vzero_asym->GetMRingV0C(0) + vzero_asym->GetMRingV0C(1) + vzero_asym->GetMRingV0C(2);
  Float_t v0c3   = vzero_asym->GetMRingV0C(3);

  isEventSelected_V0Asym &= vzero_asym->GetMTotV0C() < (330. + 100. * TMath::Power(vzero_asym->GetMTotV0A(), .2));
  isEventSelected_V0Asym &= (v0c012 < 160.) || (v0c3 > 12.*TMath::Power(.01*(v0c012 - 160.), 1.7));

  return isEventSelected_V0Asym;

}
//-----------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskTracklets::HasBCMod4(AliVEvent* fESD)
{
  // Bunch Crossing mod == 4

  Bool_t isBCMod4 = kFALSE;
  Int_t fBC = fESD->GetBunchCrossNumber();
  Int_t bcmod4 = fBC%4;

  if (bcmod4 == 2) isBCMod4 = kTRUE;

  //printf("fBC = %d, bcmod4 = %d, isBCMod4 = %d \n",fBC,bcmod4,isBCMod4);
  return isBCMod4;

}
//-----------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskTracklets::IsEventSelected(AliVEvent *fESD)

{
  Bool_t returnValue = kFALSE;

  if (
    HasNoPileupSPDInMultBins                 ( fESD ) == kTRUE &&
    //IsINELgtZERO                             ( fESD ) == kTRUE &&
    HasAcceptedVertexPosition                ( fESD ) == kTRUE  &&
    //HasNoInconsistentSPDandTrackVertices     ( fESD ) == kTRUE &&
    //IsSPDClusterVsTrackletBG                 ( fESD ) == kTRUE &&
    //IsOutOfBunchPileup                       ( fESD ) == kFALSE &&
    //IsOutOfBunchPileupFastor                 ( fESD ) == kTRUE &&
    //IsInCompleteEvent                        ( fESD ) == kFALSE &&
    //V0Decision                               ( fESD ) == kTRUE &&
    //V0Asymmetry                              ( fESD ) == kTRUE &&
    //    HasBCMod4                                 ( fESD ) == kTRUE &&
    IsMinimumBias                             ( fESD ) == kTRUE
  ) returnValue = kTRUE;
  return returnValue;
}
//-----------------------------------------------------------------------------------------

Float_t AliAnalysisTaskTracklets::DeltaPhi(Float_t phi1, Float_t phi2)
{
  Float_t dphi = TMath::Abs(phi1 - phi2);
  if (dphi > TMath::Pi()) dphi = 2. * TMath::Pi() - dphi;
  return dphi;
}

//-----------------------------------------------------------------------------------------
Float_t AliAnalysisTaskTracklets::DeltaR(Float_t phi1, Float_t eta1, Float_t phi2, Float_t eta2)
{
  Float_t dphi = DeltaPhi(phi1, phi2);
  Float_t deta = eta1 - eta2;
  return (TMath::Sqrt(dphi * dphi + deta * deta));

}
