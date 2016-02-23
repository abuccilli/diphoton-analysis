// -*- C++ -*-
//
// Package:    diphoton-analysis/ExoDiPhotonAnalyzer
// Class:      ExoDiPhotonAnalyzer
// 
/**\class ExoDiPhotonAnalyzer ExoDiPhotonAnalyzer.cc diphoton-analysis/ExoDiPhotonAnalyzer/plugins/ExoDiPhotonAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andrew Buccilli
//         Created:  Mon, 22 Feb 2016 12:27:26 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// from our CommomClasses
#include "diphoton-analysis/CommonClasses/interface/PhotonID.h"
#include "diphoton-analysis/CommonClasses/interface/PhotonInfo.h"

// for TFileService, trees
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TTree.h"

// for photons
#include "DataFormats/PatCandidates/interface/Photon.h"

// for genParticles
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class ExoDiPhotonAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ExoDiPhotonAnalyzer(const edm::ParameterSet&);
      ~ExoDiPhotonAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------

  // MiniAOD case data members
  edm::EDGetToken photonsMiniAODToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesMiniAODToken_;
  
  // main tree
  TTree *fTree;

  // photons
  ExoDiPhotons::photonInfo_t fPhotonInfo1; // leading
  ExoDiPhotons::photonInfo_t fPhotonInfo2; // sub-leading
  
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ExoDiPhotonAnalyzer::ExoDiPhotonAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   usesResource("TFileService");

   edm::Service<TFileService> fs;

   fTree = fs->make<TTree>("fTree","PhotonTree");
   fTree->Branch("Photon1",&fPhotonInfo1,ExoDiPhotons::photonBranchDefString.c_str());
   fTree->Branch("Photon2",&fPhotonInfo2,ExoDiPhotons::photonBranchDefString.c_str());

   // MiniAOD tokens
   photonsMiniAODToken_ = mayConsume<edm::View<pat::Photon> >
     (iConfig.getParameter<edm::InputTag>
      ("photonsMiniAOD"));
   
   genParticlesMiniAODToken_ = mayConsume<edm::View<reco::GenParticle> >
     (iConfig.getParameter<edm::InputTag>
      ("genParticlesMiniAOD"));
}


ExoDiPhotonAnalyzer::~ExoDiPhotonAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ExoDiPhotonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   using namespace pat;
   
   cout <<  "Run: " << iEvent.id().run() << ", LS: " <<  iEvent.id().luminosityBlock() << ", Event: " << iEvent.id().event() << endl;

   ExoDiPhotons::InitPhotonInfo(fPhotonInfo1);
   ExoDiPhotons::InitPhotonInfo(fPhotonInfo2);

   // Get pat::Photon
   edm::Handle<edm::View<pat::Photon> > photons;
   iEvent.getByToken(photonsMiniAODToken_,photons);

   // Get reco::GenParticle
   Handle<edm::View<reco::GenParticle> > genParticles;
   iEvent.getByToken(genParticlesMiniAODToken_,genParticles);


   //for (edm::View<pat::Photon>::const_iterator pho = photons->begin(); pho != photons->end(); ++pho) {
   for (size_t i = 0; i < photons->size(); ++i){
     const auto pho = photons->ptrAt(i);
     cout << "Photon: " << "pt = " << pho->pt() << "; eta = " << pho->eta() << "; phi = " << pho->phi() << endl;
   }
   
   cout << endl;
   
   //for (edm::View<reco::GenParticle>::const_iterator gen = genParticles->begin(); gen != genParticles->end(); ++gen) {
   for (size_t i = 0; i < genParticles->size(); ++i) {
     const auto gen = genParticles->ptrAt(i);
     if (gen->pt() > 50)
       cout << "GenParticle: " << "pt = " << gen->pt() << "; eta = " << gen->eta() << "; phi = " << gen->phi() << endl;
   }
   
   // fill our tree
   fTree->Fill();

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
ExoDiPhotonAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ExoDiPhotonAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ExoDiPhotonAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ExoDiPhotonAnalyzer);
