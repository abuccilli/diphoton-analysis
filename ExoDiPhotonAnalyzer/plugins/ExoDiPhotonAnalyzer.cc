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
#include "diphoton-analysis/CommonClasses/interface/EventInfo.h"
#include "diphoton-analysis/CommonClasses/interface/DiPhotonInfo.h"

// for TFileService, trees
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TTree.h"

// for ECAL topology
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

// for photons
#include "DataFormats/PatCandidates/interface/Photon.h"

// for diphotons

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

  // ECAL recHits
  edm::InputTag recHitsEBTag_;
  edm::InputTag recHitsEETag_;
  edm::EDGetTokenT<EcalRecHitCollection> recHitsEBToken;
  edm::EDGetTokenT<EcalRecHitCollection> recHitsEEToken;

  // rho token
  edm::EDGetTokenT<double> rhoToken_;

  // rho variable
  float rho_;
  
  // main tree
  TTree *fTree;

  // photons
  ExoDiPhotons::photonInfo_t fPhotonInfo1; // leading
  ExoDiPhotons::photonInfo_t fPhotonInfo2; // sub-leading

  // diphotons
  ExoDiPhotons::diphotonInfo_t fDiphotonInfo;
  
  // event
  ExoDiPhotons::eventInfo_t fEventInfo;
  
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
  : rhoToken_(consumes<double> (iConfig.getParameter<edm::InputTag>("rho")))

{
   //now do what ever initialization is needed
   usesResource("TFileService");

   edm::Service<TFileService> fs;

   fTree = fs->make<TTree>("fTree","PhotonTree");
   fTree->Branch("Event",&fEventInfo,ExoDiPhotons::eventBranchDefString.c_str());
   fTree->Branch("Photon1",&fPhotonInfo1,ExoDiPhotons::photonBranchDefString.c_str());
   fTree->Branch("Photon2",&fPhotonInfo2,ExoDiPhotons::photonBranchDefString.c_str());
   fTree->Branch("Diphoton",&fDiphotonInfo,ExoDiPhotons::diphotonBranchDefString.c_str());

   // MiniAOD tokens
   photonsMiniAODToken_ = mayConsume<edm::View<pat::Photon> >
     (iConfig.getParameter<edm::InputTag>
      ("photonsMiniAOD"));
   
   genParticlesMiniAODToken_ = mayConsume<edm::View<reco::GenParticle> >
     (iConfig.getParameter<edm::InputTag>
      ("genParticlesMiniAOD"));

   // ECAL RecHits
   recHitsEBTag_ = iConfig.getUntrackedParameter<edm::InputTag>("RecHitsEBTag",edm::InputTag("reducedEgamma:reducedEBRecHits"));
   recHitsEETag_ = iConfig.getUntrackedParameter<edm::InputTag>("RecHitsEETag",edm::InputTag("reducedEgamma:reducedEERecHits"));
   recHitsEBToken = consumes <edm::SortedCollection<EcalRecHit> > (recHitsEBTag_);
   recHitsEEToken = consumes <edm::SortedCollection<EcalRecHit> > (recHitsEETag_);
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

   ExoDiPhotons::InitEventInfo(fEventInfo);
   ExoDiPhotons::FillEventInfo(fEventInfo,iEvent);
   
   cout <<  "Run: " << iEvent.id().run() << ", LS: " <<  iEvent.id().luminosityBlock() << ", Event: " << iEvent.id().event() << endl;

   // Get rho
   edm::Handle< double > rhoH;
   iEvent.getByToken(rhoToken_,rhoH);
   rho_ = *rhoH;

   cout << "rho: " << rho_ << endl;
   
   // ECAL RecHits
   edm::Handle<EcalRecHitCollection> recHitsEB;
   iEvent.getByToken(recHitsEBToken,recHitsEB);   
   edm::Handle<EcalRecHitCollection> recHitsEE;
   iEvent.getByToken(recHitsEEToken,recHitsEE);

   if (!recHitsEB.isValid()) {
     return;
   }
   if (!recHitsEE.isValid()) {
     return;
   }
   
   
   // ECAL Topology
   const CaloSubdetectorTopology* subDetTopologyEB_;
   const CaloSubdetectorTopology* subDetTopologyEE_;
   edm::ESHandle<CaloTopology> caloTopology;
   iSetup.get<CaloTopologyRecord>().get(caloTopology);
   if (!caloTopology.isValid()) {
     return;
   }
   //const CaloTopology *topology = caloTopology.product();
   subDetTopologyEB_ = caloTopology->getSubdetectorTopology(DetId::Ecal,EcalBarrel);
   subDetTopologyEE_ = caloTopology->getSubdetectorTopology(DetId::Ecal,EcalEndcap);

   ExoDiPhotons::InitPhotonInfo(fPhotonInfo1);
   ExoDiPhotons::InitPhotonInfo(fPhotonInfo2);

   // Get pat::Photon
   edm::Handle<edm::View<pat::Photon> > photons;
   iEvent.getByToken(photonsMiniAODToken_,photons);

   if (!photons.isValid()) {
     return;
   }
   
   // Get reco::GenParticle
   Handle<edm::View<reco::GenParticle> > genParticles;
   iEvent.getByToken(genParticlesMiniAODToken_,genParticles);


   //for (edm::View<pat::Photon>::const_iterator pho = photons->begin(); pho != photons->end(); ++pho) {
   for (size_t i = 0; i < photons->size(); ++i){
     const auto pho = photons->ptrAt(i);
     cout << "Photon: " << "pt = " << pho->pt() << "; eta = " << pho->eta() << "; phi = " << pho->phi() << endl;
     cout << "isSat: " <<
       ExoDiPhotons::isSaturated(&(*pho), &(*recHitsEB), &(*recHitsEE), &(*subDetTopologyEB_), &(*subDetTopologyEE_))
	  << endl;
   }
   
   cout << endl;
   
   //for (edm::View<reco::GenParticle>::const_iterator gen = genParticles->begin(); gen != genParticles->end(); ++gen) {
   for (size_t i = 0; i < genParticles->size(); ++i) {
     const auto gen = genParticles->ptrAt(i);
     if (gen->pt() > 50)
       cout << "GenParticle: " << "pt = " << gen->pt() << "; eta = " << gen->eta() << "; phi = " << gen->phi() << endl;
   }

   ExoDiPhotons::InitDiphotonInfo(fDiphotonInfo);
   
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
