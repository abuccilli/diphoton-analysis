#ifndef PHOTON_ID_INC
#define PHOTON_ID_INC

//for photons
//#include "DataFormats/EgammaCandidates/interface/Photon.h"
//#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

namespace ExoDiPhotons
{
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
  bool passHighPtID(const pat::Photon *photon, double cHIso, double phoIso, double rho, double sigmaIeIe, bool isSaturated, bool passCSEV)
  {
    bool passHadTowerOverEmCut = false;
    bool passCHIsoCut = false;
    bool passCorPhoIso = false;
    bool passSigmaIeIe = false;
    double phoEta = fabs(photon->superCluster()->eta());
    
    if (phoEta < 1.4442) {
      passHadTowerOverEmCut = photon->hadTowOverEm() < 0.05;
      passCHIsoCut = cHIso < 5;
      passCorPhoIso = corPhoIso(photon,phoIso,rho) < 2.75;
      if (isSaturated) {
	passSigmaIeIe = sigmaIeIe < 0.0112;
      }
      else {
	passSigmaIeIe = sigmaIeIe < 0.0105;
      }
    } // end EB

    else if (1.566 < phoEta && phoEta < 2.5) {
      passHadTowerOverEmCut = photon->hadTowOverEm() < 0.05;
      passCHIsoCut = cHIso < 5;
      passCorPhoIso = corPhoIso(photon,phoIso,rho) < 2.75;
      if (isSaturated) {
	passSigmaIeIe = sigmaIeIe < 0.03;
      }
      else {
	passSigmaIeIe = sigmaIeIe < 0.0112;
      }
    } // end EE    
    return passHadTowerOverEmCut && passCHIsoCut && passCorPhoIso && passSigmaIeIe && passCSEV;
  }
  
}

#endif
