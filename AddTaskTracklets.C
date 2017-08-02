///////////////////////////////////////////////////////////////////////
//                                                                    //
//            AddTaskTracklets Macro to run on grids                   //
//  Author: Prabi and Paolo Bartalini, July CCNU 2017                 //
//                                                                   //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

AliAnalysisTaskTracklets* AddTaskTracklets(TString name="Tracklets.root", Bool_t  useMC  = kFALSE)
{
  // get the manager via the static access member. since it's static, you don't need
  // an instance of the class to call the function

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    return 0x0;
  }
  // get the input event handler this handler is part of the managing system and feeds events to your task
  if (!mgr->GetInputEventHandler()) {
    return 0x0;
  }

  // now you create an instance of your task
  AliAnalysisTaskTracklets* tasktracklets = new AliAnalysisTaskTracklets("taskTracklets");
  if(!tasktracklets) return 0x0;

  // add your task to the manager
  mgr->AddTask(tasktracklets);
  tasktracklets->SetUseMC(useMC);

  // your task needs input: here you connect the manager to your task
  mgr->ConnectInput(tasktracklets,0,mgr->GetCommonInputContainer());
  // same for the output
  mgr->ConnectOutput(tasktracklets,1,mgr->CreateContainer("Hist", TList::Class(), AliAnalysisManager::kOutputContainer, name+".root"));
  mgr->ConnectOutput(tasktracklets,2,mgr->CreateContainer("Tree1", TTree::Class(), AliAnalysisManager::kOutputContainer, name+".root"));
  mgr->ConnectOutput(tasktracklets,3,mgr->CreateContainer("Tree2", TTree::Class(), AliAnalysisManager::kOutputContainer, name+".root"));
  mgr->ConnectOutput(tasktracklets,4,mgr->CreateContainer("Tree3", TTree::Class(), AliAnalysisManager::kOutputContainer, name+".root"));
  mgr->ConnectOutput(tasktracklets,5,mgr->CreateContainer("Tree4", TTree::Class(), AliAnalysisManager::kOutputContainer, name+".root"));



  // this macro returns a pointer to the task, this will be convenient when you run the analysis in an analysis train on grids
  return tasktracklets;
}
