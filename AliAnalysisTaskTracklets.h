/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* Add a description of your Tracklets analysis */
/* Author: Prabi and Paolo Bartalini, July CCNU 2017    */

#ifndef AliAnalysisTaskTracklets_H
#define AliAnalysisTaskTracklets_H

class AliESDtrackCuts;
class AliESDEvent;
class TList;
class TH1F;
class TH2F;
class TH1I;
class TTree;
class AliMultiplicity;

#include "AliAnalysisTaskSE.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"


class AliAnalysisTaskTracklets : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskTracklets();
  AliAnalysisTaskTracklets(const char *name);
  virtual                 ~AliAnalysisTaskTracklets();

  virtual void            UserCreateOutputObjects();
  virtual void            UserExec(Option_t* option);
  virtual void            Terminate(Option_t* option);

  void   TrackletsLoop(AliESDEvent *);
  static Bool_t IsMinimumBias(AliVEvent *);
  static Bool_t IsINELgtZERO(AliVEvent *);
  static Bool_t HasAcceptedVertexPosition(AliVEvent *);
  static Bool_t HasNoPileupSPDInMultBins(AliVEvent *);
  static Bool_t HasNoInconsistentSPDandTrackVertices(AliVEvent *);
  static Bool_t IsEventSelected(AliVEvent *);
  static Bool_t IsInCompleteEvent(AliVEvent *);
  static Bool_t IsSPDClusterVsTrackletBG(AliVEvent *);
  static Bool_t IsOutOfBunchPileup(AliVEvent *);
  static Bool_t IsOutOfBunchPileupFastor(AliVEvent *);
  static Bool_t HasBCMod4(AliVEvent *);

  static Bool_t V0Decision(AliVEvent *);
  static Bool_t V0Asymmetry(AliVEvent *);

  void       SetUseMC(Bool_t mc = kFALSE)              {fUseMC = mc;}   // use to analyse MC data

  virtual Float_t DeltaPhi(Float_t phi1, Float_t phi2);                               // difference in the phi of lead and track
  virtual Float_t DeltaR(Float_t phi1, Float_t eta1, Float_t phi2, Float_t eta2);   // distance between lead and track

protected:

  void       FillMCGen();
  Float_t fWeight;

private:
  AliESDtrackCuts*        fCuts;                                                //!
  AliESDEvent*            fESD;                                                 //! input ESD event
  AliStack*    fStack;                                                          //! MC stack
  AliMCEvent*  fMCEvent;                                                        //! MC Event
  Bool_t       fUseMC;                                                          // analyze MC events


  TList*                  fOutputList;                                          //! output list in the root file

  TH1I*                   fEventCounter;                                        //! Event Summary with cuts
  TH1I*                   fHistMultiplicityStd;                                 //! Std. Reference Multiplicity distribution
  TH1F*                   fHistGenZvtx;                                        //! Gen vertex Z
  TH1F*                   fHistRecoZvtx;                                        //! Reco vertex Z
  TH1I*                   fHistMult_reco;                                 //! reconstructed Multiplicity distribution
  TH1I*                   fHistMult_prim;                                 //! primary Multiplicity distribution
  TH1I*                   fHistMult_seco;                                 //! sec Reference Multiplicity distribution
  TH1I*                   fHistMult_fake;                                 //! fake. Reference Multiplicity distribution





  Float_t                 fVertexMC[3];                                           //! Gen vertex X,Y,Z
  Float_t                 fVertexReco[3];                                           //! Reconstructed vertex X,Y,Z


  TTree*   fTree_reco;                                                              //!
  TTree*   fTree_prim;                                                              //!
  TTree*   fTree_seco;                                                              //!
  TTree*   fTree_fake;                                                              //!

  Float_t fDist_reco;                                                            //! pT
  Float_t fPhi_reco;                                                            //!Phi
  Float_t fDPhi_reco;                                                            //!Phi
  Float_t fEta_reco;                                                            //!Eta
  Float_t fDTheta_reco;                                                            //Eta


  Float_t fPt_prim;                                                            //! pT
  Float_t fEta_prim;                                                           //!Eta
  Float_t fPhi_prim;                                                            //!Phi
  Float_t fDPhi_prim;                                                            //!Phi
  Float_t fDTheta_prim;                                                            //!Phi


  Float_t fPt_seco;                                                            //! pT
  Float_t fPhi_seco;                                                            //!Phi  Float_t fPt_tracklet;

  Float_t fPt_fake;                                                            //! pT
  Float_t fPhi_fake;                                                            //!Phi


  AliAnalysisTaskTracklets(const AliAnalysisTaskTracklets&);                    // not implemented
  AliAnalysisTaskTracklets& operator=(const AliAnalysisTaskTracklets&);         // not implemented

  ClassDef(AliAnalysisTaskTracklets, 1);
};

#endif
