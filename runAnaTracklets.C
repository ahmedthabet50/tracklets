#if !defined (__CINT__) || (defined(__MAKECINT__))
#include <iostream>
#include "AliAnalysisGrid.h"
#include "TSystem.h"
#include "TROOT.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisGrid.h"
#include "AliVEventHandler.h"
#include "AliPhysicsSelection.h"
#include "AliAnalysisAlien.h"
#include "AliMCEventHandler.h"
#include "AliPhysicsSelectionTask.h"
#include "TRegexp.h"
#include "TProof.h"
#include "AliESDInputHandler.h"
#include "AliOADBPhysicsSelection.h"
#include "TGrid.h"
#include "AliMultSelectionTask.h"

//#include "AliAnalysisTaskTracklets.h"

#endif

void runAnaTracklets(Bool_t local = kFALSE, Bool_t gridTest = kTRUE, Bool_t mcData = kTRUE, TString output = "mc" )
{

  // for arguments in the macro
  // 1.  set local = kTRUE if you want to run the analysis locally, or on grid (kFALSE) -> root -l 'runAnaTracklets.C(1,0,0,"data_local")'
  // 2. if you run on grid, specify test mode set gridTest = kTRUE  or full grid model (kFALSE) -> root -l 'runAnaTracklets.C(0,1,0,"data_test")'
  // 3. set mcData = 0 (Data) and mcData = 1 (MC)   -> root -l 'runAnaTracklets.C(0,0,0,"data_full")'
  // 4. Set desired name to your task and outputfile


  // Inorder to compile a class, tell root where to look for headers
  gROOT->ProcessLine(".include $ROOTSYS/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");

  //gSystem->Setenv("alien_CLOSE_SE","ALICE::WUHAN::SE");
  //if (!TGrid::Connect("alien://")) return;

  // create the analysis manager and ESD handler
  AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskTracklets");
  AliESDInputHandler *esdHandler = new AliESDInputHandler();
  esdHandler->SetNeedField();

  mgr->SetInputEventHandler(esdHandler);

  if (mcData){
    AliMCEventHandler* mcHandler = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mcHandler);
    //mgr->SetReadTR(kFALSE);
  }


  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask* physicsSelTask = AddTaskPhysicsSelection(mcData,kTRUE);

  //gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  //AliMultSelectionTask *MultSelection = AddTaskMultSelection();

  // compile the class (locally) before testing on the grid or for full grid job
  gROOT->LoadMacro("AliAnalysisTaskTracklets.cxx++g");
  // load the addtask macro
  gROOT->LoadMacro("AddTaskTracklets.C");
  // create an instance of your analysis task
  AliAnalysisTaskTracklets *task = AddTaskTracklets(output,mcData);

  if(!mgr->InitAnalysis()) return;
  mgr->SetDebugLevel(9);
  mgr->PrintStatus();
  mgr->SetUseProgressBar(1, 100);


  if(local) {
    // if you want to run locally, you need to define some input depending upon type of data

    if(mcData) {
      //    TString test_data= "2016l_mc.txt";
      TString test_data= "/Users/prabi/work/data/2015f_mc.txt";  // local path to local root_archive.zip files
    } else {
      //  TString test_data= "2016l_data.txt";
      TString test_data= "/Users/prabi/work/data/2015f.txt";   // local path to ESD files


    }

    TChain* chain = new TChain("esdTree");
    // add a few files to the chain (change this so that your local files are added)

    //chain->Add("AliESDs.root"); // start the analysis locally, reading the events from the tchain from list of files
    TFileCollection* file= new TFileCollection("data","",test_data);
    chain->AddFileInfoList(file->GetList());
    //    mgr->StartAnalysis("local", chain);
    mgr->StartAnalysis("local", chain,100,0); //2100

  } else {
    // if you want to run on grid,  you create and configure the plugin
    AliAnalysisAlien *alienHandler = new AliAnalysisAlien();

    alienHandler->SetUser("ppalni"); // username
    //  alienHandler->SetUser("xuz"); // username


    // also specify the include (header) paths on grid
    alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
    // make sure your source files get copied to grid
    alienHandler->SetAdditionalLibs("AliAnalysisTaskTracklets.cxx AliAnalysisTaskTracklets.h");
    alienHandler->SetAnalysisSource("AliAnalysisTaskTracklets.cxx");
    // select the aliphysics version. all other packages  are LOADED AUTOMATICALLY! you do not have to now specify aliroot version etc
    alienHandler->SetAliPhysicsVersion("vAN-20170725-1");
    // select the input data

    if(mcData){
      alienHandler->SetGridDataDir("/alice/sim/2017/LHC17e2");    // pp 5.02 TeV Pythia8
      //    alienHandler->SetGridDataDir("/alice/sim/2017/LHC17f2a_cent");  // pPb 5 TeV  EPOS

      alienHandler->SetDataPattern("/*/AliESDs.root"); // for testing over a single file segment
    }
    else{
      alienHandler->SetGridDataDir("/alice/data/2015/LHC15n"); // pp 5.02 TeV real data
      alienHandler->SetDataPattern("pass4/*/AliESDs.root"); // for testing over a single file segment

      //    alienHandler->SetGridDataDir("/alice/data/2016/LHC16q");  // p-Pb @ 5.02 TeV real data
      //    alienHandler->SetDataPattern("pass1_CENT_wSDD/*/AliESDs.root"); // for testing over a single file segment
      alienHandler->SetRunPrefix("000");
    }

    // runnumber

    alienHandler->AddRunNumber(244421);   // pp 5 TeV

    // alienHandler->AddRunNumber(265334);   // pPb  5 TeV

    // number of files per sub-job
    alienHandler->SetSplitMaxInputFileNumber(200); // previously it was 50 .. decrease this number to split master job into more sub-jobs
    alienHandler->SetExecutable("Tracklets.sh");
    // specify how many seconds your job may take
    alienHandler->SetMasterResubmitThreshold(99); // default (?) resubmit threshold
    alienHandler->SetTTL(80000);
    alienHandler->SetJDLName("Tracklets.jdl");

    alienHandler->SetOutputToRunNo(kTRUE);
    alienHandler->SetKeepLogs(kTRUE);
    // merging: run with kTRUE to merge on grid all the subjobs after re-running the jobs in SetRunMode("terminate")
    alienHandler->SetMaxMergeStages(1);
    alienHandler->SetMergeViaJDL(kTRUE);   // merge results on grid after finishing jobs (use terminate option when run this)

    alienHandler->SetOutputToRunNo(1); // run number wise output
    alienHandler->SetNrunsPerMaster(1); // run number wise masterjobs

  //  alienHandler->SetDefaultOutputs(kFALSE);
  //  alienHandler->SetOutputFiles("Tracklets.root");
//    alienHandler->SetMergeExcludes("EventStat_temp.root");

    // define the output folders

  //      alienHandler->SetGridWorkingDir("Data_July27");  // Data
    alienHandler->SetGridWorkingDir("MC_July27");  //MC


    // connect the alien plugin to the manager
    mgr->SetGridHandler(alienHandler);

    if(gridTest) {
      // speficy on how many files you want to run
      alienHandler->SetNtestFiles(2);
      // and launch the analysis
      alienHandler->SetRunMode("test");
      mgr->StartAnalysis("grid",10000,0);
    } else {
      // else launch the full grid analysis
      alienHandler->SetRunMode("full");
    //      alienHandler->SetRunMode("terminate");
      mgr->StartAnalysis("grid");
    }
  }
}
