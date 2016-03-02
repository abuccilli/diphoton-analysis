#ifndef PHOTON_ID_INC
#define PHOTON_ID_INC

// for photons
#include "DataFormats/PatCandidates/interface/Photon.h"

// for saturation
#include "RecoCaloTools/Navigation/interface/CaloNavigator.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

namespace ExoDiPhotons
{
  // checking for saturated photons in 5x5 around seed crystal
  // considered saturated if any crystal is marked as saturated
  bool isSaturated(const pat::Photon *photon, const EcalRecHitCollection *recHitsEB, const EcalRecHitCollection *recHitsEE,
		   const CaloSubdetectorTopology* subDetTopologyEB_, const CaloSubdetectorTopology* subDetTopologyEE_) {
    using namespace std;
    
    bool isSat = false;
    DetId seedDetId = ((photon->superCluster())->seed())->seed();
    
    // check EB
    if (seedDetId.subdetId()==EcalBarrel) {
      CaloNavigator<DetId> cursor = CaloNavigator<DetId>(seedDetId,subDetTopologyEB_);
      for (int i = -2; i <= 2; ++i) {
	for (int j = -2; j <= 2; ++j) {
	  cursor.home();
	  cursor.offsetBy(i,j);
	  EcalRecHitCollection::const_iterator it = recHitsEB->find(*cursor);
	  if(it != recHitsEB->end()) {
	    cout << "Energy of (" << i << ", " << j << "): " << it-> energy()
		 << ", kSaturated: " << it->checkFlag(EcalRecHit::kSaturated)
		 << ", kDead: " << it->checkFlag(EcalRecHit::kDead)
		 << ", kKilled: " << it->checkFlag(EcalRecHit::kKilled)
		 << endl;
	    if (it->checkFlag(EcalRecHit::kSaturated) && !it->checkFlag(EcalRecHit::kDead) && !it->checkFlag(EcalRecHit::kKilled)) {
	      isSat = true;
	    }
	  }	  
	}
      }
    }
    // check EE
    else if (seedDetId.subdetId()==EcalEndcap) {
      CaloNavigator<DetId> cursor = CaloNavigator<DetId>(seedDetId,subDetTopologyEE_);
      for (int i = -2; i <= 2; ++i) {
	for (int j = -2; j <= 2; ++j) {
	  cursor.home();
	  cursor.offsetBy(i,j);
	  EcalRecHitCollection::const_iterator it = recHitsEE->find(*cursor);
	  if(it != recHitsEB->end()) {
	    cout << "Energy of (" << i << ", " << j << "): " << it->energy()
		 << ", kSaturated: " << it->checkFlag(EcalRecHit::kSaturated)
		 << ", kDead: " << it->checkFlag(EcalRecHit::kDead)
		 << ", kKilled: " << it->checkFlag(EcalRecHit::kKilled)
		 << endl;
	    if (it->checkFlag(EcalRecHit::kSaturated) && !it->checkFlag(EcalRecHit::kDead) && !it->checkFlag(EcalRecHit::kKilled)) {
	      isSat = true;
	    }
	  }
	}
      }
    }
    return isSat;
  }
  
  // Using High pT ID V2 cuts.
  double corPhoIso(const pat::Photon *photon, double phoIso, double rho)
  {
    double phoPt = photon->pt();
    double phoEta = fabs(photon->superCluster()->eta());
    
    if (phoEta < 1.4442) {
      if (phoEta < 0.9) {
	return (2.5 + phoIso - rho*0.17 - 0.0045*phoPt);
      }
      else {
	return (2.5 + phoIso - rho*0.14 - 0.0045*phoPt);
      }
    } // end EB
    
    else if (1.566 < phoEta && phoEta < 2.5) {
      if (phoEta < 2.0) {
	return (2.5 + phoIso - rho*0.11 - 0.0045*phoPt);
      }
      else if (phoEta < 2.2) {
	return (2.5 + phoIso - rho*0.14 - 0.003*phoPt);
      }
      else {
	return (2.5 + phoIso - rho*0.22 - 0.003*phoPt);
      }
    } // end EE

    else {
      return 99999.99;
    }
  }
  
  // Using High pT ID V2 cuts.
  bool passHighPtID(const pat::Photon *photon, double cHIso, double phoIso, double rho, double sigmaIeIe, bool isPhoSat)
  {
    bool passHadTowerOverEmCut = false;
    bool passCHIsoCut = false;
    bool passCorPhoIso = false;
    bool passSigmaIeIe = false;
    bool passCSEV = photon->passElectronVeto();
    double phoEta = fabs(photon->superCluster()->eta());
    
    if (phoEta < 1.4442) {
      passHadTowerOverEmCut = photon->hadTowOverEm() < 0.05;
      passCHIsoCut = cHIso < 5;
      passCorPhoIso = corPhoIso(photon,phoIso,rho) < 2.75;
      if (isPhoSat) {
	passSigmaIeIe = sigmaIeIe < 0.0112;
      }
      else {
	passSigmaIeIe = sigmaIeIe < 0.0105;
      }
    } // end EB

    else if (1.566 < phoEta && phoEta < 2.5) {
      passHadTowerOverEmCut = photon->hadTowOverEm() < 0.05;
      passCHIsoCut = cHIso < 5;
      passCorPhoIso = corPhoIso(photon,phoIso,rho) < 2.0;
      if (isPhoSat) {
	passSigmaIeIe = sigmaIeIe < 0.03;
      }
      else {
	passSigmaIeIe = sigmaIeIe < 0.028;
      }
    } // end EE    
    return passHadTowerOverEmCut && passCHIsoCut && passCorPhoIso && passSigmaIeIe && passCSEV;
  }
  
}

#endif
