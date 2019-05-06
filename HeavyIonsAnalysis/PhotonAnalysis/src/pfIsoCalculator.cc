#include "DataFormats/Math/interface/deltaPhi.h"

#include "HeavyIonsAnalysis/PhotonAnalysis/src/pfIsoCalculator.h"

pfIsoCalculator::pfIsoCalculator(const edm::Event &iEvent,
                                 const edm::EDGetTokenT<edm::View<reco::PFCandidate> > pfCandidates_,
                                 const math::XYZPoint& pv)
{
  iEvent.getByToken(pfCandidates_, candidatesView);

  vtx_ = pv;
}

double pfIsoCalculator::getPfIso(const reco::Photon& photon, int pfId,
                                 double r1, double r2, double threshold,
                                 double jWidth)
{
  double photonEta = photon.eta();
  double photonPhi = photon.phi();
  double totalEt = 0;

  for (const auto& pf : *candidatesView) {
    if ( pf.particleId() != pfId )   continue;
    double pfEta = pf.eta();
    double pfPhi = pf.phi();

    double dEta = std::abs(photonEta - pfEta);
    double dPhi = reco::deltaPhi(pfPhi, photonPhi);
    double dR2 = dEta*dEta + dPhi*dPhi;
    double pfPt = pf.pt();

    // remove the photon itself
    if ( pf.superClusterRef() == photon.superCluster() ) continue;

    if(pf.particleId() == reco::PFCandidate::h){
      float dz = std::abs(pf.vz() - vtx_.z());
      if (dz > 0.2) continue;
      double dxy = ((vtx_.x() - pf.vx())*pf.py() + (pf.vy() - vtx_.y())*pf.px()) / pf.pt();
      if (std::abs(dxy) > 0.1) continue;
    }

    // Jurassic Cone /////
    if (dR2 > r1 *r1) continue;
    if (dR2 < r2 * r2) continue;
    if (std::abs(dEta) < jWidth)  continue;
    if (pfPt < threshold) continue;
    totalEt += pfPt;
  }

  return totalEt;
}

double pfIsoCalculator::getPfIsoSubUE(const reco::Photon& photon, int pfId,
                                      double r1, double r2, double threshold,
                                      double jWidth)
{
  double photonEta = photon.eta();
  double photonPhi = photon.phi();
  double totalEt = 0;

  for (const auto& pf : *candidatesView) {
    if ( pf.particleId() != pfId )   continue;
    double pfEta = pf.eta();
    double pfPhi = pf.phi();

    double dEta = std::abs(photonEta - pfEta);
    if (dEta > r1) continue;
    if (dEta < jWidth)  continue;

    // remove the photon itself
    if ( pf.superClusterRef() == photon.superCluster() ) continue;

    if(pf.particleId() == reco::PFCandidate::h){
      float dz = std::abs(pf.vz() - vtx_.z());
      if (dz > 0.2) continue;
      double dxy = ((vtx_.x() - pf.vx())*pf.py() + (pf.vy() - vtx_.y())*pf.px()) / pf.pt();
      if (std::abs(dxy) > 0.1) continue;
    }

    double dPhi = reco::deltaPhi(pfPhi, photonPhi);
    double dR2 = dEta*dEta + dPhi*dPhi;
    double pfPt = pf.pt();

    // Jurassic Cone /////
    if (dR2 < r2 * r2) continue;
    if (pfPt < threshold) continue;
    totalEt += pfPt;
  }

  double areaStrip = 4*M_PI*(r1-jWidth);   // strip area over which UE is estimated
  double areaCone = M_PI*r1*r1;       // area inside which isolation is to be applied
  if (jWidth > 0) {
    // calculate overlap area between disk with radius r2 and rectangle of width jwidth
    double angTmp = std::acos(jWidth / r1);
    double lenTmp = std::sqrt(r1*r1 - jWidth*jWidth) * 2;
    double areaTwoTriangles = lenTmp * jWidth;
    double areaTwoArcs = (M_PI - 2*angTmp) * r1 * r1;
    areaCone -= (areaTwoTriangles + areaTwoArcs);
  }

  double areaInnerCone = 0;
  if (r2 > jWidth) {
    areaInnerCone = M_PI*r2*r2;
    if (jWidth > 0) {
      double angTmp = std::acos(jWidth / r2);
      double lenTmp = std::sqrt(r2*r2 - jWidth*jWidth) * 2;
      double areaTwoTriangles = lenTmp * jWidth;
      double areaTwoArcs = (M_PI - 2*angTmp) * r2 * r2;
      areaInnerCone -= (areaTwoTriangles + areaTwoArcs);
    }
  }
  areaStrip -= areaInnerCone;
  areaCone -= areaInnerCone;

  double subEt = getPfIso(photon, pfId, r1, r2, threshold, jWidth) - totalEt * (areaCone / areaStrip);
  return subEt;
}

double pfIsoCalculator::getPfIso(const reco::GsfElectron& ele, int pfId,
                                 double r1, double r2, double threshold)
{
  double eleEta = ele.eta();
  double elePhi = ele.phi();
  double totalEt = 0.;

  for (const auto& pf : *candidatesView) {
    if ( pf.particleId() != pfId )   continue;
    double pfEta = pf.eta();
    double pfPhi = pf.phi();

    double dEta = std::abs(eleEta - pfEta);
    double dPhi = reco::deltaPhi(pfPhi, elePhi);
    double dR2 = dEta*dEta + dPhi*dPhi;
    double pfPt = pf.pt();

    // remove electron itself
    if (pf.particleId() == reco::PFCandidate::e) continue;

    if (pf.particleId() == reco::PFCandidate::h) {
      float dz = std::abs(pf.vz() - vtx_.z());
      if (dz > 0.2) continue;
      double dxy = ((vtx_.x() - pf.vx())*pf.py() + (pf.vy() - vtx_.y())*pf.px()) / pf.pt();
      if (std::abs(dxy) > 0.1) continue;
    }

    // inside the cone size
    if (dR2 > r1 * r1) continue;
    if (dR2 < r2 * r2) continue;
    if (pfPt < threshold) continue;
    totalEt += pfPt;

  }

  return totalEt;
}
