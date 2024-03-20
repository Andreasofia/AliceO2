// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "ITSStudies/Efficiency.h"
#include "ITSStudies/ITSStudiesConfigParam.h"
#include "CommonUtils/TreeStreamRedirector.h"
#include "DataFormatsITS/TrackITS.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/TopologyDictionary.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/Task.h"
#include "ITSBase/GeometryTGeo.h"
#include "ITStracking/IOUtils.h"
#include "ReconstructionDataFormats/DCA.h"
#include "SimulationDataFormat/MCTrack.h"
#include "Steer/MCKinematicsReader.h"
#include "ReconstructionDataFormats/TrackParametrization.h"

#include <TH1.h>
#include <TH1D.h>
#include <TH1I.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TObjArray.h>
#include <THStack.h>
#include <TString.h>

namespace o2::its::study
{
using namespace o2::framework;
using namespace o2::globaltracking;

using GTrackID = o2::dataformats::GlobalTrackID;

class EfficiencyStudy : public Task
{
 public:
  EfficiencyStudy(std::shared_ptr<DataRequest> dr,
                  mask_t src,
                  bool useMC,
                  std::shared_ptr<o2::steer::MCKinematicsReader> kineReader,
                  std::shared_ptr<o2::base::GRPGeomRequest> gr) : mDataRequest(dr), mTracksSrc(src), mUseMC(useMC), mKineReader(kineReader), mGGCCDBRequest(gr){};

  ~EfficiencyStudy() final = default;
  void init(InitContext&) final;
  void run(ProcessingContext&) final;
  void endOfStream(EndOfStreamContext&) final;
  void finaliseCCDB(ConcreteDataMatcher&, void*) final;
  void initialiseRun(o2::globaltracking::RecoContainer&);
  int getDCAClusterTrackMC();
  void studyDCAcutsMC();
  void studyClusterSelectionMC();
  void process(o2::globaltracking::RecoContainer&);
  void setClusterDictionary(const o2::itsmft::TopologyDictionary* d) { mDict = d; }

 private:
  void updateTimeDependentParams(ProcessingContext& pc);
  bool mVerboseOutput = false;
  bool mUseMC;
  std::string mOutFileName;
  double b;
  std::shared_ptr<o2::steer::MCKinematicsReader> mKineReader;
  GeometryTGeo* mGeometry;
  const o2::itsmft::TopologyDictionary* mDict = nullptr;
  int mNLayers = 3;

  // Spans
  gsl::span<const o2::itsmft::ROFRecord> mTracksROFRecords;
  gsl::span<const o2::itsmft::ROFRecord> mClustersROFRecords;
  gsl::span<const o2::its::TrackITS> mTracks;
  gsl::span<const o2::MCCompLabel> mTracksMCLabels;
  gsl::span<const o2::itsmft::CompClusterExt> mClusters;
  gsl::span<const unsigned char> mClusPatterns;
  gsl::span<const int> mInputITSidxs;
  const o2::dataformats::MCLabelContainer* mClustersMCLCont;
  std::vector<o2::BaseCluster<float>> mITSClustersArray;

  // Data
  GTrackID::mask_t mTracksSrc{};
  std::shared_ptr<DataRequest> mDataRequest;
  unsigned short mMask = 0x7f;

  // Utils
  std::shared_ptr<o2::base::GRPGeomRequest> mGGCCDBRequest;
  std::unique_ptr<TFile> mOutFile;

  //// Histos
  // Distance betweeen original and duplicated clusters
  std::unique_ptr<TH1D> mDistanceClustersX;
  std::unique_ptr<TH1D> mDistanceClustersY;
  std::unique_ptr<TH1D> mDistanceClustersZ;
  std::unique_ptr<TH1D> mDistanceClusters;
  // DCA betweeen track and original cluster
  std::unique_ptr<TH1D> mDCAxyOriginal;
  std::unique_ptr<TH1D> mDCAzOriginal;
  // DCA betweeen track and duplicated cluster
  std::unique_ptr<TH1D> mDCAxyDuplicated;
  std::unique_ptr<TH1D> mDCAzDuplicated;
  // phi of the cluster
  std::unique_ptr<TH1D> mPhiOriginal;
  std::unique_ptr<TH1D> mPhiDuplicated;
  // position of the clusters
  std::unique_ptr<TH3D> m3DClusterPositions;
  std::unique_ptr<TH2D> m2DClusterOriginalPositions;
  std::unique_ptr<TH2D> m2DClusterDuplicatedPositions;

  std::unique_ptr<TH1D> mEfficiencyGoodMatch;
  std::unique_ptr<TH1D> mEfficiencyFakeMatch;
  std::unique_ptr<TH1D> mEfficiencyTotal;
};

void EfficiencyStudy::init(InitContext& ic)
{
  LOGP(info, "--------------- init");

  o2::base::GRPGeomHelper::instance().setRequest(mGGCCDBRequest);

  auto& pars = o2::its::study::ITSEfficiencyParamConfig::Instance();
  mOutFileName = pars.outFileName;
  b = pars.b;

  mDistanceClustersX = std::make_unique<TH1D>("distanceClustersX", ";Distance x (cm); ", 100, 0, 1);
  mDistanceClustersY = std::make_unique<TH1D>("distanceClustersY", ";Distance y (cm); ", 100, 0, 1);
  mDistanceClustersZ = std::make_unique<TH1D>("distanceClustersZ", ";Distance z (cm); ", 100, 0, 1);
  mDistanceClusters = std::make_unique<TH1D>("distanceClusters", ";Distance (cm); ", 100, 0, 1);

  mDCAxyOriginal = std::make_unique<TH1D>("dcaXYOriginal", "Distance between track and original cluster ;DCA xy (cm); ", 100, -0.05, 0.05);
  mDCAzOriginal = std::make_unique<TH1D>("dcaZOriginal", "Distance between track and original cluster ;DCA z (cm); ", 100, -0.05, 0.05);
  mDCAxyDuplicated = std::make_unique<TH1D>("dcaXYDuplicated", "Distance between track and duplicated cluster  ;DCA xy (cm); ", 100, -0.05, 0.05);
  mDCAzDuplicated = std::make_unique<TH1D>("dcaZDuplicated", "Distance between track and duplicated cluster  ;DCA z (cm); ", 100, -0.05, 0.05);

  mPhiOriginal = std::make_unique<TH1D>("phiOriginal", ";phi (deg); ", 120, -180, 180);
  mPhiDuplicated = std::make_unique<TH1D>("phiDuplicated", ";phi (deg); ", 120, -180, 180);

  m3DClusterPositions = std::make_unique<TH3D>("3DClusterPositions", ";x (cm);y (cm);z (cm)", 200, -10, 10, 200, -10, 10, 400, -20, 20);
  m2DClusterOriginalPositions = std::make_unique<TH2D>("m2DClusterOriginalPositions", ";x (cm);y (cm)", 400, -10, 10, 400, -6, 6);
  m2DClusterDuplicatedPositions = std::make_unique<TH2D>("m2DClusterDuplicatedPositions", ";x (cm);y (cm)", 400, -10, 10, 400, -6, 6);

  mEfficiencyGoodMatch = std::make_unique<TH1D>("mEfficiencyGoodMatch", ";#sigma(DCA) cut;Efficiency;", 20, 0.5, 20.5);
  mEfficiencyFakeMatch = std::make_unique<TH1D>("mEfficiencyFakeMatch", ";#sigma(DCA) cut;Efficiency;", 20, 0.5, 20.5);
  mEfficiencyTotal = std::make_unique<TH1D>("mEfficiencyTotal", ";#sigma(DCA) cut;Efficiency;", 20, 0.5, 20.5);
  gStyle->SetPalette(55);
}

void EfficiencyStudy::run(ProcessingContext& pc)
{
  LOGP(info, "--------------- run");
  o2::globaltracking::RecoContainer recoData;
  recoData.collectData(pc, *mDataRequest.get());

  updateTimeDependentParams(pc); // Make sure this is called after recoData.collectData, which may load some conditions
  initialiseRun(recoData);
  process(recoData);
}

void EfficiencyStudy::initialiseRun(o2::globaltracking::RecoContainer& recoData)
{
  LOGP(info, "--------------- initialiseRun");
  if (mUseMC) {
    mTracksMCLabels = recoData.getITSTracksMCLabels();
    mClustersMCLCont = recoData.getITSClustersMCLabels();
  }

  mTracksROFRecords = recoData.getITSTracksROFRecords();
  mTracks = recoData.getITSTracks();
  mClusters = recoData.getITSClusters();
  mClustersROFRecords = recoData.getITSClustersROFRecords();
  mClusPatterns = recoData.getITSClustersPatterns();
  mInputITSidxs = recoData.getITSTracksClusterRefs();
  mITSClustersArray.reserve(mClusters.size());
  auto pattIt = mClusPatterns.begin();
  o2::its::ioutils::convertCompactClusters(mClusters, pattIt, mITSClustersArray, mDict); // clusters converted to 3D spacepoints
}

int EfficiencyStudy::getDCAClusterTrackMC()
{
  // get the DCA between the clusters and the track from MC and fill histograms: distance between original and duplicated cluster, DCA, phi, clusters
  LOGP(info, "--------------- getDCAClusterTrackMC");

  mOutFile->mkdir("DistanceClusters/");
  mOutFile->mkdir("DCA/");

  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  o2::gpu::gpustd::array<float, 2> clusOriginalDCA, clusDuplicatedDCA;
  auto propagator = o2::base::Propagator::Instance();

  unsigned int rofIndexTrack = 0;
  unsigned int rofNEntriesTrack = 0;
  unsigned int rofIndexClus = 0;
  unsigned int rofNEntriesClus = 0;
  int nLabels = 0;
  unsigned int totClus = 0;

  int duplicated = 0;

  std::unordered_map<o2::MCCompLabel, std::vector<int>> label_vecClus[mClustersROFRecords.size()][mNLayers]; // array of maps nRofs x Nlayers -> {label, vec(iClus)} where vec(iClus) are the clusters that share the same label

  for (unsigned int iROF = 0; iROF < mTracksROFRecords.size(); iROF++) { // loop on ROFRecords array
    rofIndexTrack = mTracksROFRecords[iROF].getFirstEntry();
    rofNEntriesTrack = mTracksROFRecords[iROF].getNEntries();

    rofIndexClus = mClustersROFRecords[iROF].getFirstEntry();
    rofNEntriesClus = mClustersROFRecords[iROF].getNEntries();

    for (unsigned int iTrack = rofIndexTrack; iTrack < rofIndexTrack + rofNEntriesTrack; iTrack++) { // loop on tracks per ROF
      auto track = mTracks[iTrack];
      o2::track::TrackParCov trackParCov = mTracks[iTrack];
      int firstClus = track.getFirstClusterEntry(); // get the first cluster of the track
      int ncl = track.getNumberOfClusters();        // get the number of clusters of the track

      auto& tracklab = mTracksMCLabels[iTrack];
      if (tracklab.isFake())
        continue;

      for (int iclTrack = firstClus; iclTrack < firstClus + ncl; iclTrack++) { // loop on clusters associated to the track
        auto& clusOriginal = mClusters[mInputITSidxs[iclTrack]];
        auto clusOriginalPoint = mITSClustersArray[mInputITSidxs[iclTrack]]; // cluster spacepoint in the tracking system
        auto staveOriginal = mGeometry->getStave(clusOriginal.getSensorID());
        auto chipOriginal = mGeometry->getChipIdInStave(clusOriginal.getSensorID());

        auto layer = mGeometry->getLayer(clusOriginal.getSensorID());
        if (layer >= mNLayers)
          continue;                                                            // checking only selected layers
        auto labsTrack = mClustersMCLCont->getLabels(mInputITSidxs[iclTrack]); // get labels of the cluster associated to the track

        o2::math_utils::Point3D<float> clusOriginalPointTrack = {clusOriginalPoint.getX(), clusOriginalPoint.getY(), clusOriginalPoint.getZ()};
        o2::math_utils::Point3D<float> clusOriginalPointGlob = mGeometry->getMatrixT2G(clusOriginal.getSensorID()) * clusOriginalPointTrack;

        auto phiOriginal = std::atan2(clusOriginalPointGlob.y(), clusOriginalPointGlob.x()) * 180 / M_PI;
        mPhiOriginal->Fill(phiOriginal);
        m3DClusterPositions->Fill(clusOriginalPointGlob.x(), clusOriginalPointGlob.y(), clusOriginalPointGlob.z());
        m2DClusterOriginalPositions->Fill(clusOriginalPointGlob.x(), clusOriginalPointGlob.y());

        for (auto& labT : labsTrack) { // for each valid label iterate over ALL the clusters in the ROF to see if there are duplicates
          if (labT != tracklab)
            continue;
          nLabels++;
          if (labT.isValid()) {
            for (unsigned int iClus = rofIndexClus; iClus < rofIndexClus + rofNEntriesClus; iClus++) { // iteration over ALL the clusters in the ROF
              auto clusDuplicated = mClusters[iClus];
              auto clusDuplicatedPoint = mITSClustersArray[iClus];

              auto layerClus = mGeometry->getLayer(clusDuplicated.getSensorID());
              if (layerClus != layer)
                continue;

              o2::math_utils::Point3D<float> clusDuplicatedPointTrack = {clusDuplicatedPoint.getX(), clusDuplicatedPoint.getY(), clusDuplicatedPoint.getZ()};
              o2::math_utils::Point3D<float> clusDuplicatedPointGlob = mGeometry->getMatrixT2G(clusDuplicated.getSensorID()) * clusDuplicatedPointTrack;

              auto phiDuplicated = std::atan2(clusDuplicatedPointGlob.y(), clusDuplicatedPointGlob.x()) * 180 / M_PI;

              auto labsClus = mClustersMCLCont->getLabels(iClus); // ideally I can have more than one label per cluster
              for (auto labC : labsClus) {
                if (labC == labT) {
                  label_vecClus[iROF][layerClus][labT].push_back(iClus); // same cluster: label from the track = label from the cluster
                                                                         // if a duplicate cluster is found, propagate the track to the duplicate cluster and compute the distance from the original cluster
                                                                         // if (clusOriginalPointGlob != clusDuplicatedPointGlob) { /// check that the duplicated cluster is not the original one just counted twice
                                                                         // if (clusDuplicated.getSensorID() != clusOriginal.getSensorID()) { /// check that the duplicated cluster is not the original one just counted twice

                  // applying constraints: the cluster should be on the same layer, should be on an adjacent stave and on the same or adjacent chip position
                  if (clusDuplicated.getSensorID() == clusOriginal.getSensorID())
                    continue;
                  auto layerDuplicated = mGeometry->getLayer(clusDuplicated.getSensorID());
                  if (layerDuplicated != layerClus)
                    continue;
                  auto staveDuplicated = mGeometry->getStave(clusDuplicated.getSensorID());
                  if (abs(staveDuplicated - staveOriginal) != 1)
                    continue;
                  auto chipDuplicated = mGeometry->getChipIdInStave(clusDuplicated.getSensorID());
                  if (abs(chipDuplicated - chipOriginal) > 1)
                    continue;

                  duplicated++; 
                  mPhiDuplicated->Fill(phiDuplicated);
                  m3DClusterPositions->Fill(clusDuplicatedPointGlob.x(), clusDuplicatedPointGlob.y(), clusDuplicatedPointGlob.z());
                  m2DClusterDuplicatedPositions->Fill(clusDuplicatedPointGlob.x(), clusDuplicatedPointGlob.y());

                  /// compute the distance between original and dubplicated cluster
                  mDistanceClustersX->Fill(abs(clusOriginalPointGlob.x() - clusDuplicatedPointGlob.x()));
                  mDistanceClustersY->Fill(abs(clusOriginalPointGlob.y() - clusDuplicatedPointGlob.y()));
                  mDistanceClustersZ->Fill(abs(clusOriginalPointGlob.z() - clusDuplicatedPointGlob.z()));
                  mDistanceClusters->Fill(std::hypot(clusOriginalPointGlob.x() - clusDuplicatedPointGlob.x(), clusOriginalPointGlob.y() - clusDuplicatedPointGlob.y(), clusOriginalPointGlob.z() - clusDuplicatedPointGlob.z()));

                  /// Compute the DCA between the cluster location and the track

                  /// first propagate to the original cluster
                  trackParCov.rotate(mGeometry->getSensorRefAlpha(clusOriginal.getSensorID()));
                  if (propagator->propagateToDCA(clusOriginalPointGlob, trackParCov, b, 2.f, matCorr, &clusOriginalDCA)) {
                    mDCAxyOriginal->Fill(clusOriginalDCA[0]);
                    mDCAzOriginal->Fill(clusOriginalDCA[1]);
                  }
                  /// then propagate to the duplicated cluster
                  trackParCov.rotate(mGeometry->getSensorRefAlpha(clusDuplicated.getSensorID()));
                  if (propagator->propagateToDCA(clusDuplicatedPointGlob, trackParCov, b, 2.f, matCorr, &clusDuplicatedDCA)) {
                    mDCAxyDuplicated->Fill(clusDuplicatedDCA[0]);
                    mDCAzDuplicated->Fill(clusDuplicatedDCA[1]);
                  }
                  ///////////////////////////////////////////////////////
                }
              }
            }
          }
        }
      } // end loop on clusters
      totClus += ncl;
    } // end loop on tracks per ROF
  }   // end loop on ROFRecords array
  LOGP(info, "Total number of clusters: {} ", totClus);
  LOGP(info, "total nLabels: {}", nLabels);
  LOGP(info, "Number of duplicated clusters: {}", duplicated);

  mPhiOriginal->Write();
  mPhiDuplicated->Write();
  m3DClusterPositions->Write();
  m2DClusterOriginalPositions->Write();
  m2DClusterDuplicatedPositions->Write();

  mOutFile->cd("DistanceClusters");
  mDistanceClustersX->Write();
  mDistanceClustersY->Write();
  mDistanceClustersZ->Write();
  mDistanceClusters->Write();

  mOutFile->cd("DCA");
  mDCAxyOriginal->Write();
  mDCAzOriginal->Write();
  mDCAxyDuplicated->Write();
  mDCAzDuplicated->Write();

  if (mVerboseOutput && mUseMC) {
    // printing the duplicates
    for (unsigned int iROF = 0; iROF < mClustersROFRecords.size(); iROF++) {
      LOGP(info, "°°°°°°°°°°°°°°°°°°°°°°°° ROF {} °°°°°°°°°°°°°°°°°°°°°°°°", iROF);
      for (unsigned int lay = 0; lay < mNLayers; lay++) {
        LOGP(info, "°°°°°°°°°°°°°°°°°°°°°°°° LAYER {} °°°°°°°°°°°°°°°°°°°°°°°°", lay);
        for (auto& it : label_vecClus[iROF][lay]) {
          if (it.second.size() <= 1)
            continue; // printing only duplicates
          std::cout << " \n++++++++++++ Label: ";
          auto label = it.first;
          it.first.print();
          for (auto iClus : it.second) {
            auto name = mGeometry->getSymbolicName(mClusters[iClus].getSensorID());
            auto chipid = mClusters[iClus].getChipID();
            auto clus = mClusters[iClus];
            auto clusPoint = mITSClustersArray[iClus];

            o2::math_utils::Point3D<float> clusPointTrack = {clusPoint.getX(), clusPoint.getY(), clusPoint.getZ()};
            o2::math_utils::Point3D<float> clusPointGlob = mGeometry->getMatrixT2G(clus.getSensorID()) * clusPointTrack;
            std::cout << "ROF: " << iROF << ", iClus: " << iClus << " -> chip: " << chipid << " = " << name << std::endl;
            LOGP(info, "LOCtrack: {} {} {}", clusPointTrack.x(), clusPointTrack.y(), clusPointTrack.z());
            LOGP(info, "LOCglob {} {} {}", clusPointGlob.x(), clusPointGlob.y(), clusPointGlob.z());
          }
        }
      }
    }
  }
  return duplicated;
}

void EfficiencyStudy::studyDCAcutsMC()
{

  LOGP(info, "--------------- studyDCAcutsMC");

  int duplicated = getDCAClusterTrackMC();

  double meanDCAxyDuplicated = mDCAxyDuplicated->GetMean();
  double meanDCAzDuplicated = mDCAzDuplicated->GetMean();
  double sigmaDCAxyDuplicated = mDCAxyDuplicated->GetRMS();
  double sigmaDCAzDuplicated = mDCAzDuplicated->GetRMS();

  LOGP(info, "meanDCAxyDuplicated: {}, meanDCAzDuplicated: {}, sigmaDCAxyDuplicated: {}, sigmaDCAzDuplicated: {}", meanDCAxyDuplicated, meanDCAzDuplicated, sigmaDCAxyDuplicated, sigmaDCAzDuplicated);

  // now we have the DCA distribution:
  //  ->iterate again over tracks and over duplicated clusters and find the matching ones basing on DCA cuts (1 sigma, 2 sigma,...)
  //  then control if the matching ones according to the DCA matches have the same label
  //  if yes, then we have a good match -> increase the good match counter
  //  if not, keep it as a fake match -> increase the fake match counter
  //  the efficiency of each one will be match counter / total of the duplicated clusters
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  o2::gpu::gpustd::array<float, 2> clusOriginalDCA, clusDuplicatedDCA;
  auto propagator = o2::base::Propagator::Instance();

  unsigned int rofIndexTrack = 0;
  unsigned int rofNEntriesTrack = 0;
  unsigned int rofIndexClus = 0;
  unsigned int rofNEntriesClus = 0;
  int nLabels = 0;
  unsigned int totClus = 0;

  unsigned int nDCAMatches[20] = {0};
  unsigned int nGoodMatches[20] = {0};
  unsigned int nFakeMatches[20] = {0};
  unsigned int nOriginalMatches = 0;

  for (unsigned int iROF = 0; iROF < mTracksROFRecords.size(); iROF++) { // loop on ROFRecords array
    rofIndexTrack = mTracksROFRecords[iROF].getFirstEntry();
    rofNEntriesTrack = mTracksROFRecords[iROF].getNEntries();

    rofIndexClus = mClustersROFRecords[iROF].getFirstEntry();
    rofNEntriesClus = mClustersROFRecords[iROF].getNEntries();

    for (unsigned int iTrack = rofIndexTrack; iTrack < rofIndexTrack + rofNEntriesTrack; iTrack++) { // loop on tracks per ROF
      auto track = mTracks[iTrack];
      o2::track::TrackParCov trackParCov = mTracks[iTrack];
      int firstClus = track.getFirstClusterEntry(); // get the first cluster of the track
      int ncl = track.getNumberOfClusters();        // get the number of clusters of the track

      auto& tracklab = mTracksMCLabels[iTrack];
      if (tracklab.isFake())
        continue;

      if (mVerboseOutput) {
        LOGP(info, "--------- track Label: ");
        tracklab.print();
      }

      for (int iclTrack = firstClus; iclTrack < firstClus + ncl; iclTrack++) { // loop on clusters associated to the track to extract layer, stave and chip to restrict the possible matches to be searched with the DCA cut

        auto& clusOriginal = mClusters[mInputITSidxs[iclTrack]];
        auto layerOriginal = mGeometry->getLayer(clusOriginal.getSensorID());
        if (layerOriginal >= mNLayers)
          continue;
        nOriginalMatches++;                                                       // checking only selected layers
        auto labsOriginal = mClustersMCLCont->getLabels(mInputITSidxs[iclTrack]); // get labels of the cluster associated to the track (original)
        auto staveOriginal = mGeometry->getStave(clusOriginal.getSensorID());
        auto chipOriginal = mGeometry->getChipIdInStave(clusOriginal.getSensorID());

        /// for each oroginal cluster iterate over all the possible "adjacten" clusters (stave +-1, chip =,+-1) and calculate the DCA with the track. Then compare the cluster label with the track label to see if it is a true or fake match
        for (unsigned int iClus = rofIndexClus; iClus < rofIndexClus + rofNEntriesClus; iClus++) { // iteration over ALL the clusters in the ROF
          auto clusDuplicated = mClusters[iClus];
          //// applying constraints: the cluster should be on the same layer, should be on an adjacent stave and on the same or adjacent chip position
          if (clusDuplicated.getSensorID() == clusOriginal.getSensorID())
            continue;
          auto layerDuplicated = mGeometry->getLayer(clusDuplicated.getSensorID());
          if (layerDuplicated != layerOriginal)
            continue;
          auto staveDuplicated = mGeometry->getStave(clusDuplicated.getSensorID());
          if (abs(staveDuplicated - staveOriginal) != 1)
            continue;
          auto chipDuplicated = mGeometry->getChipIdInStave(clusDuplicated.getSensorID());
          if (abs(chipDuplicated - chipOriginal) > 1)
            continue;

          auto labsDuplicated = mClustersMCLCont->getLabels(iClus);

          /// if the cheks are passed, then calculate the DCA
          auto clusDuplicatedPoint = mITSClustersArray[iClus];

          o2::math_utils::Point3D<float> clusDuplicatedPointTrack = {clusDuplicatedPoint.getX(), clusDuplicatedPoint.getY(), clusDuplicatedPoint.getZ()};
          o2::math_utils::Point3D<float> clusDuplicatedPointGlob = mGeometry->getMatrixT2G(clusDuplicated.getSensorID()) * clusDuplicatedPointTrack;

          /// Compute the DCA between the duplicated cluster location and the track
          trackParCov.rotate(mGeometry->getSensorRefAlpha(clusDuplicated.getSensorID()));
          if (propagator->propagateToDCA(clusDuplicatedPointGlob, trackParCov, b, 2.f, matCorr, &clusDuplicatedDCA)) { // check if the propagation fails
            if (mVerboseOutput)
              LOGP(info, "Propagation ok");
            /// checking the DCA for 15 different sigma ranges
            for (int i = 0; i < 20; i++) {

              // if (abs(meanDCAxyDuplicated - clusDuplicatedDCA[0]) < (i+1)*sigmaDCAxyDuplicated){ // check if the DCA is within the cut i*sigma
              if (abs(meanDCAxyDuplicated - clusDuplicatedDCA[0]) < (i + 1) * sigmaDCAxyDuplicated && abs(meanDCAzDuplicated - clusDuplicatedDCA[1]) < (i + 1) * sigmaDCAzDuplicated) { // check if the DCA is within the cut i*sigma
                if (mVerboseOutput)
                  LOGP(info, "Check DCA ok: {} < {}; {} < {}", abs(meanDCAxyDuplicated - clusDuplicatedDCA[0]), (i + 1) * sigmaDCAxyDuplicated, abs(meanDCAzDuplicated - clusDuplicatedDCA[1]), (i + 1) * sigmaDCAzDuplicated);
                nDCAMatches[i]++;
                bool isGoodMatch = false;

                for (auto labD : labsDuplicated) { // at this point the track has been matched with a duplicated cluster based on the DCA cut. Now we check if the matching is good ore not based on the label
                  if (mVerboseOutput) {
                    LOGP(info, "tracklab, labD:");
                    tracklab.print();
                    labD.print();
                  }
                  if (labD == tracklab) { //// check if the label of the origial cluster is equal to the label of the duplicated cluster among all the labels for a cluster
                    isGoodMatch = true;
                    continue;
                  }
                }
                if (isGoodMatch)
                  nGoodMatches[i]++;
                else
                  nFakeMatches[i]++;
              } else if (mVerboseOutput)
                LOGP(info, "Check DCA failed");
            }
          } else if (mVerboseOutput)
            LOGP(info, "Propagation failed");
        } // end loop on all the clusters in the rof
      }   // end loop on clusters associated to the track
    }     // end loop on tracks per ROF
  }       // end loop on ROFRecords array

  for (int i = 0; i < 20; i++) {
    LOGP(info, "Cut: {} sigma -> number of duplicated clusters: {} nDCAMatches: {} nGoodMatches: {} nFakeMatches: {}", i + 1, duplicated, nDCAMatches[i], nGoodMatches[i], nFakeMatches[i]);
    mEfficiencyGoodMatch->SetBinContent(i + 1, double(nGoodMatches[i]) / double(duplicated));
    mEfficiencyFakeMatch->SetBinContent(i + 1, double(nFakeMatches[i]) / double(duplicated));
    mEfficiencyTotal->SetBinContent(i + 1, double(nGoodMatches[i] + nFakeMatches[i]) / double(duplicated));
  }

  mOutFile->mkdir("Efficiency/");
  mOutFile->cd("Efficiency/");
  mEfficiencyGoodMatch->Write();
  mEfficiencyFakeMatch->Write();
  mEfficiencyTotal->Write();

  mEfficiencyGoodMatch->GetYaxis()->SetRangeUser(-0.1, 1.1);
  mEfficiencyFakeMatch->GetYaxis()->SetRangeUser(-0.1, 1.1);
  mEfficiencyTotal->GetYaxis()->SetRangeUser(-0.1, 1.1);

  TCanvas c;
  auto leg = std::make_unique<TLegend>(0.75, 0.45, 0.89, 0.65);
  leg->AddEntry(mEfficiencyGoodMatch.get(), "#frac{# good matches}{# tot duplicated clusters}", "p");
  leg->AddEntry(mEfficiencyFakeMatch.get(), "#frac{# fake matches}{# tot duplicated clusters}", "p");
  leg->AddEntry(mEfficiencyTotal.get(), "#frac{# tot matches}{# tot duplicated clusters}", "p");

  mEfficiencyGoodMatch->Draw("lp PLC PMC");
  mEfficiencyFakeMatch->Draw("same lp PLC PMC");
  mEfficiencyTotal->Draw("same lp PLC PMC");
  leg->Draw("same");
  c.Write();
  c.SaveAs("prova.png");
}

void EfficiencyStudy::studyClusterSelectionMC()
{
  // study to find a good selection method for the duplicated cluster, to be used for non MC data 
  // iterate over tracks an associated clusters, and find the closer cluster that is not the original one
  // control that the DCA is within 1,2,...,15 sigma -> compute the efficiency for each cut

  LOGP(info, "--------------- studyClusterSelection");

  int duplicated = getDCAClusterTrackMC();

  double meanDCAxyDuplicated = mDCAxyDuplicated->GetMean();
  double meanDCAzDuplicated = mDCAzDuplicated->GetMean();
  double sigmaDCAxyDuplicated = mDCAxyDuplicated->GetRMS();
  double sigmaDCAzDuplicated = mDCAzDuplicated->GetRMS();

  LOGP(info, "meanDCAxyDuplicated: {}, meanDCAzDuplicated: {}, sigmaDCAxyDuplicated: {}, sigmaDCAzDuplicated: {}", meanDCAxyDuplicated, meanDCAzDuplicated, sigmaDCAxyDuplicated, sigmaDCAzDuplicated);


  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  o2::gpu::gpustd::array<float, 2> clusOriginalDCA, clusDuplicatedDCA;
  auto propagator = o2::base::Propagator::Instance();

  short unsigned int mNLayers = 3;
  unsigned int rofIndexTrack = 0;
  unsigned int rofNEntriesTrack = 0;
  unsigned int rofIndexClus = 0;
  unsigned int rofNEntriesClus = 0;
  int nLabels = 0;
  unsigned int totClus = 0;

  unsigned int nDCAMatches[15] = {0};
  unsigned int nGoodMatches[15] = {0};
  unsigned int nFakeMatches[15] = {0};
  unsigned int nOriginalMatches = 0;

  std::map<std::tuple<int, double, o2::MCCompLabel>, bool> clusterMatches = {};

  for (unsigned int iROF = 0; iROF < mTracksROFRecords.size(); iROF++) { // loop on ROFRecords array
    rofIndexTrack = mTracksROFRecords[iROF].getFirstEntry();
    rofNEntriesTrack = mTracksROFRecords[iROF].getNEntries();

    rofIndexClus = mClustersROFRecords[iROF].getFirstEntry();
    rofNEntriesClus = mClustersROFRecords[iROF].getNEntries();

    for (unsigned int iTrack = rofIndexTrack; iTrack < rofIndexTrack + rofNEntriesTrack; iTrack++) { // loop on tracks per ROF
      auto track = mTracks[iTrack];
      o2::track::TrackParCov trackParCov = mTracks[iTrack];
      int firstClus = track.getFirstClusterEntry(); // get the first cluster of the track
      int ncl = track.getNumberOfClusters();        // get the number of clusters of the track

      auto& tracklab = mTracksMCLabels[iTrack];
      if (tracklab.isFake())
        continue;

      if (mVerboseOutput) {
        LOGP(info, "--------- track Label: ");
        tracklab.print();
      }

      for (int iclTrack = firstClus; iclTrack < firstClus + ncl; iclTrack++) { // loop on clusters associated to the track to extract layer, stave and chip to restrict the possible matches to be searched with the DCA cut
        LOGP(info, "New cluster");
        auto& clusOriginal = mClusters[mInputITSidxs[iclTrack]];
        auto layerOriginal = mGeometry->getLayer(clusOriginal.getSensorID());
        if (layerOriginal >= mNLayers)
          continue;
        nOriginalMatches++;                                                       // checking only selected layers
        auto labsOriginal = mClustersMCLCont->getLabels(mInputITSidxs[iclTrack]); // get labels of the cluster associated to the track (original)
        auto staveOriginal = mGeometry->getStave(clusOriginal.getSensorID());
        auto chipOriginal = mGeometry->getChipIdInStave(clusOriginal.getSensorID());

        // std::vector<std::tuple<int,double,gsl::span<const o2::MCCompLabel>>> vec_clusID_rDCA_label;
        std::tuple<int, double, gsl::span<const o2::MCCompLabel>> clusID_rDCA_label = {0, 999., gsl::span<const o2::MCCompLabel>()}; // inizializing tuple with dummy values

        // double minDCA = 999.;

        /// for each oroginal cluster iterate over all the possible "adjacten" clusters (stave +-1, chip =,+-1) and calculate the DCA with the track. Then choose the closest one.
        for (unsigned int iClus = rofIndexClus; iClus < rofIndexClus + rofNEntriesClus; iClus++) { // iteration over ALL the clusters in the ROF
          auto clusDuplicated = mClusters[iClus];

          if (clusDuplicated.getSensorID() == clusOriginal.getSensorID())
            continue;

          //// applying constraints: the cluster should be on the same layer, should be on an adjacent stave and on the same or adjacent chip position
          auto layerDuplicated = mGeometry->getLayer(clusDuplicated.getSensorID());
          if (layerDuplicated != layerOriginal)
            continue;
          auto staveDuplicated = mGeometry->getStave(clusDuplicated.getSensorID());
          if (abs(staveDuplicated - staveOriginal) != 1)
            continue;
          auto chipDuplicated = mGeometry->getChipIdInStave(clusDuplicated.getSensorID());
          if (abs(chipDuplicated - chipOriginal) > 1)
            continue;

          auto labsDuplicated = mClustersMCLCont->getLabels(iClus);

          /// if the cheks are passed, then calculate the DCA
          auto clusDuplicatedPoint = mITSClustersArray[iClus];

          o2::math_utils::Point3D<float> clusDuplicatedPointTrack = {clusDuplicatedPoint.getX(), clusDuplicatedPoint.getY(), clusDuplicatedPoint.getZ()};
          o2::math_utils::Point3D<float> clusDuplicatedPointGlob = mGeometry->getMatrixT2G(clusDuplicated.getSensorID()) * clusDuplicatedPointTrack;

          /// Compute the DCA between the duplicated cluster location and the track
          trackParCov.rotate(mGeometry->getSensorRefAlpha(clusDuplicated.getSensorID()));
          if (!propagator->propagateToDCA(clusDuplicatedPointGlob, trackParCov, b, 2.f, matCorr, &clusDuplicatedDCA)) { // check if the propagation fails
            continue;
          }

          if (!(abs(meanDCAxyDuplicated - clusDuplicatedDCA[0]) < 10 * sigmaDCAxyDuplicated && abs(meanDCAzDuplicated - clusDuplicatedDCA[1]) < 10 * sigmaDCAzDuplicated)) {
            continue;
          }

          // if (mVerboseOutput)
          //   LOGP(info, "Propagation ok");
          double rDCA = std::hypot(clusDuplicatedDCA[0], clusDuplicatedDCA[1]);

          // for (auto lab: std::get<2>(clusID_rDCA_label)){
          //   lab.print();
          // }

          if (rDCA < std::get<1>(clusID_rDCA_label)) {
            // minDCA = rDCA;
            // LOGP(info, "+++++++++++++++++ NEW min found: {} {} <---- {}", iClus, rDCA, std::get<1>(clusID_rDCA_label));
            clusID_rDCA_label = {iClus, rDCA, labsDuplicated};
            // for (auto lab: labsDuplicated){
            //   lab.print();
            // }
          }
          // else
          // LOGP(info, " new value bigger, skipping (old, new)= {} {} ", std::get<1>(clusID_rDCA_label), rDCA);
          // std::tuple<int, double, gsl::span<const o2::MCCompLabel>> tuple_clusID_rDCA_label = {iClus, rDCA, labsDuplicated};
          // vec_clusID_rDCA_label.push_back(tuple_clusID_rDCA_label);

          // here clusID_rDCA_label is updated with the closest cluster to the track other than the original one
          // checking if it is a good or fake match looking at the labels
          bool isGood = false;
          for (auto lab : std::get<2>(clusID_rDCA_label)) {
            if (lab == tracklab) {
              std::tuple<int, double, o2::MCCompLabel> tuple_clusID_rDCA_label = {std::get<0>(clusID_rDCA_label), std::get<1>(clusID_rDCA_label), lab};
              clusterMatches[tuple_clusID_rDCA_label] = true;
              LOGP(info, " Duplicated cluster found: clus, rdca, label: {} {} {} TRUE", std::get<0>(clusID_rDCA_label), std::get<1>(clusID_rDCA_label), lab.getTrackID());
              isGood = true;
              continue;
            }
          }
          if (!isGood) {
            std::tuple<int, double, o2::MCCompLabel> tuple_clusID_rDCA_label = {std::get<0>(clusID_rDCA_label), std::get<1>(clusID_rDCA_label), std::get<2>(clusID_rDCA_label)[0]};
            LOGP(info, " Duplicated cluster found: clus, rdca, label: {} {} {} FALSE", std::get<0>(clusID_rDCA_label), std::get<1>(clusID_rDCA_label), std::get<2>(clusID_rDCA_label)[0].getTrackID());
            clusterMatches[tuple_clusID_rDCA_label] = false;
          }

        } // end loop on all the clusters in the rof
      }   // end loop on clusters associated to the track
    }     // end loop on tracks per ROF
  }       // end loop on ROFRecords array

  // for (auto el:clusterMatches ){
  //   LOGP(info, " clus, rdca, good, label: {} {} {}", std::get<0>(el.first),  std::get<1>(el.first), el.second  );
  //   std::get<2>(el.first).print();
  // }

  for (int i = 0; i < 15; i++) {
    LOGP(info, "Cut: {} sigma -> number of duplicated clusters: {} nDCAMatches: {} nGoodMatches: {} nFakeMatches: {}", i + 1, duplicated, nDCAMatches[i], nGoodMatches[i], nFakeMatches[i]);
    mEfficiencyGoodMatch->SetBinContent(i + 1, double(nGoodMatches[i]) / double(duplicated));
    mEfficiencyFakeMatch->SetBinContent(i + 1, double(nFakeMatches[i]) / double(duplicated));
    mEfficiencyTotal->SetBinContent(i + 1, double(nGoodMatches[i] + nFakeMatches[i]) / double(duplicated));
  }

  mOutFile->mkdir("Efficiency/");
  mOutFile->cd("Efficiency/");
  mEfficiencyGoodMatch->Write();
  mEfficiencyFakeMatch->Write();
  mEfficiencyTotal->Write();

  mEfficiencyGoodMatch->GetYaxis()->SetRangeUser(-0.1, 1.1);
  mEfficiencyFakeMatch->GetYaxis()->SetRangeUser(-0.1, 1.1);
  mEfficiencyTotal->GetYaxis()->SetRangeUser(-0.1, 1.1);

  TCanvas c;
  auto leg = std::make_unique<TLegend>(0.75, 0.45, 0.89, 0.65);
  leg->AddEntry(mEfficiencyGoodMatch.get(), "#frac{# good matches}{# tot duplicated clusters}", "p");
  leg->AddEntry(mEfficiencyFakeMatch.get(), "#frac{# fake matches}{# tot duplicated clusters}", "p");
  leg->AddEntry(mEfficiencyTotal.get(), "#frac{# tot matches}{# tot duplicated clusters}", "p");

  mEfficiencyGoodMatch->Draw("lp PLC PMC");
  mEfficiencyFakeMatch->Draw("same lp PLC PMC");
  mEfficiencyTotal->Draw("same lp PLC PMC");
  leg->Draw("same");
  c.Write();
  c.SaveAs("prova.png");
}

void EfficiencyStudy::process(o2::globaltracking::RecoContainer& recoData)
{
  LOGP(info, "--------------- process");

  mOutFile = std::make_unique<TFile>(mOutFileName.c_str(), "recreate");

  if (mUseMC) {
    studyDCAcutsMC();
    // studyClusterSelectionMC();
  }

  LOGP(info, "** Found in {} rofs:\n\t- {} clusters\n\t",
       mClustersROFRecords.size(), mClusters.size());

  if (mUseMC) {
    LOGP(info, "mClusters size: {}, mClustersROFRecords size: {}, mClustersMCLCont size: {}, mClustersconverted size: {} ", mClusters.size(), mClustersROFRecords.size(), mClustersMCLCont->getNElements(), mITSClustersArray.size());
    LOGP(info, "mTracks size: {}, mTracksROFRecords size: {}, mTracksMCLabels size: {}", mTracks.size(), mTracksROFRecords.size(), mTracksMCLabels.size());
  } else {
    LOGP(info, "mClusters size: {}, mClustersROFRecords size: {}, mClustersconverted size: {} ", mClusters.size(), mClustersROFRecords.size(), mITSClustersArray.size());
    LOGP(info, "mTracks size: {}, mTracksROFRecords size: {}", mTracks.size(), mTracksROFRecords.size());
  }
}

void EfficiencyStudy::updateTimeDependentParams(ProcessingContext& pc)
{
  static bool initOnceDone = false;
  o2::base::GRPGeomHelper::instance().checkUpdates(pc);
  if (!initOnceDone) { // this params need to be queried only once
    initOnceDone = true;
    mGeometry = GeometryTGeo::Instance();
    mGeometry->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::T2GRot, o2::math_utils::TransformType::T2G, o2::math_utils::TransformType::L2G));
  }
}

void EfficiencyStudy::endOfStream(EndOfStreamContext& ec)
{
}

void EfficiencyStudy::finaliseCCDB(ConcreteDataMatcher& matcher, void* obj)
{
  if (o2::base::GRPGeomHelper::instance().finaliseCCDB(matcher, obj)) {
    return;
  }
  if (matcher == ConcreteDataMatcher("ITS", "CLUSDICT", 0)) {
    setClusterDictionary((const o2::itsmft::TopologyDictionary*)obj);
    return;
  }
}

DataProcessorSpec getEfficiencyStudy(mask_t srcTracksMask, mask_t srcClustersMask, bool useMC, std::shared_ptr<o2::steer::MCKinematicsReader> kineReader)
{
  std::vector<OutputSpec> outputs;
  auto dataRequest = std::make_shared<DataRequest>();
  dataRequest->requestTracks(srcTracksMask, useMC);
  dataRequest->requestClusters(srcClustersMask, useMC);

  auto ggRequest = std::make_shared<o2::base::GRPGeomRequest>(false,                             // orbitResetTime
                                                              true,                              // GRPECS=true
                                                              false,                             // GRPLHCIF
                                                              true,                              // GRPMagField
                                                              true,                              // askMatLUT
                                                              o2::base::GRPGeomRequest::Aligned, // geometry
                                                              dataRequest->inputs,
                                                              true);
  return DataProcessorSpec{
    "its-efficiency-study",
    dataRequest->inputs,
    outputs,
    AlgorithmSpec{adaptFromTask<EfficiencyStudy>(dataRequest, srcTracksMask, useMC, kineReader, ggRequest)},
    Options{}};
}

} // namespace o2::its::study