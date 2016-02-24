#ifndef EVENT_INFO_INC
#define EVENT_INFO_INC

namespace ExoDiPhotons
{

  struct eventInfo_t {
    float pthat;
    float alphaqcd;
    float alphaqed;
    float qscale;
    float weight;
    int run;
    int LS;
    int evnum;
    int processid;
  };

  std::string eventInfoBranchDefString("pthat/F:alphaqcd:alphaqed:qscale:weight:run/I:LS:evnum:processid");
  
  void FillEventInfo(eventInfo_t &eventInfo, const edm::Event& iEvent) {
    eventInfo.run = iEvent.id().run();
    eventInfo.LS = iEvent.id().luminosityBlock();
    eventInfo.evnum = iEvent.id().event();
  }
  
  void InitEventInfo(eventInfo_t &eventInfo) {
    eventInfo.processid = (int) -99999.99;
    eventInfo.pthat = -99999.99;
    eventInfo.alphaqcd = -99999.99;
    eventInfo.alphaqed = -99999.99;
    eventInfo.qscale = -99999.99;
    eventInfo.weight = -99999.99;
    eventInfo.run = (int) -99999.99;
    eventInfo.LS = (int) -99999.99;
    eventInfo.evnum = (int) -99999.99;
  }
  
  
  
}

#endif
