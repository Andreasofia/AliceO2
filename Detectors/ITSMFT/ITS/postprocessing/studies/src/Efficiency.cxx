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
                  std::shared_ptr<o2::base::GRPGeomRequest> gr) : mDataRequest(dr), mTracksSrc(src), mKineReader(kineReader), mGGCCDBRequest(gr){};

  ~EfficiencyStudy() final = default;
  void init(InitContext&) final;
  void run(ProcessingContext&) final;
  void endOfStream(EndOfStreamContext&) final;
  void finaliseCCDB(ConcreteDataMatcher&, void*) final;
  void initialiseRun(o2::globaltracking::RecoContainer&);
  void process(o2::globaltracking::RecoContainer&);
  void setClusterDictionary(const o2::itsmft::TopologyDictionary* d) { mDict = d; }

 private:
  void updateTimeDependentParams(ProcessingContext& pc);
  std::string mOutFileName;
  double b;
  std::shared_ptr<o2::steer::MCKinematicsReader> mKineReader;
  GeometryTGeo* mGeometry;
  const o2::itsmft::TopologyDictionary* mDict = nullptr;

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
  // Distance betweeen track and original cluster, method 1
  std::unique_ptr<TH1D> mDCAxOriginalMethod1;
  std::unique_ptr<TH1D> mDCAyOriginalMethod1;
  std::unique_ptr<TH1D> mDCAzOriginalMethod1;
  std::unique_ptr<TH1D> mDCAOriginalMethod1;
  // Distance betweeen track and duplicated cluster, method 1
  std::unique_ptr<TH1D> mDCAxDuplicatedMethod1;
  std::unique_ptr<TH1D> mDCAyDuplicatedMethod1;
  std::unique_ptr<TH1D> mDCAzDuplicatedMethod1;
  std::unique_ptr<TH1D> mDCADuplicatedMethod1;
  // DCA betweeen track and original cluster, method 2
  std::unique_ptr<TH1D> mDCAxyOriginalMethod2;
  std::unique_ptr<TH1D> mDCAzOriginalMethod2;
  std::unique_ptr<TH1D> mDCAOriginalMethod2;
  // DCA betweeen track and duplicated cluster, method 2
  std::unique_ptr<TH1D> mDCAxyDuplicatedMethod2;
  std::unique_ptr<TH1D> mDCAzDuplicatedMethod2;
  std::unique_ptr<TH1D> mDCADuplicatedMethod2;
  // phi of the cluster
  std::unique_ptr<TH1D> mPhi;
};

void EfficiencyStudy::init(InitContext& ic)
{
  LOGP(info, "--------------- init");

  o2::base::GRPGeomHelper::instance().setRequest(mGGCCDBRequest);

  auto& pars = o2::its::study::ITSEfficiencyParamConfig::Instance();
  mOutFileName = pars.outFileName;
  b = pars.b;

  mDistanceClustersX = std::make_unique<TH1D>("distanceClustersX", ";Distance x (cm); ", 100, 0, 5);
  mDistanceClustersY = std::make_unique<TH1D>("distanceClustersY", ";Distance y (cm); ", 100, 0, 5);
  mDistanceClustersZ = std::make_unique<TH1D>("distanceClustersZ", ";Distance z (cm); ", 100, 0, 5);
  mDistanceClusters = std::make_unique<TH1D>("distanceClusters", ";Distance (cm); ", 100, 0, 5);

  mDCAxOriginalMethod1 = std::make_unique<TH1D>("dcaXOriginalMethod1", "Method 1: distance between track and original cluster ;DCA x (cm); ", 100, 0, 0.5);
  mDCAyOriginalMethod1 = std::make_unique<TH1D>("dcaYOriginalMethod1", "Method 1: distance between track and original cluster ;DCA y (cm); ", 100, 0, 5);
  mDCAzOriginalMethod1 = std::make_unique<TH1D>("dcaZOriginalMethod1", "Method 1: distance between track and original cluster ;DCA z (cm); ", 100, 0, 5);
  mDCAOriginalMethod1 = std::make_unique<TH1D>("dcaOriginalMethod1", "Method 1: distance between track and original cluster ;DCA (cm); ", 100, 0, 5);
  mDCAxDuplicatedMethod1 = std::make_unique<TH1D>("dcaXDuplicatedMethod1", "Method 1: distance between track and duplicated cluster  ;DCA x (cm); ", 100, 0, 5);
  mDCAyDuplicatedMethod1 = std::make_unique<TH1D>("dcaYDuplicatedMethod1", "Method 1: distance between track and duplicated cluster  ;DCA y (cm); ", 100, 0, 5);
  mDCAzDuplicatedMethod1 = std::make_unique<TH1D>("dcaZDuplicatedMethod1", "Method 1: distance between track and duplicated cluster  ;DCA z (cm); ", 100, 0, 5);
  mDCADuplicatedMethod1 = std::make_unique<TH1D>("dcaDuplicatedMethod1", "Method 1: distance between track and duplicated cluster  ;DCA (cm); ", 100, 0, 5);

  mDCAxyOriginalMethod2 = std::make_unique<TH1D>("dcaXYOriginalMethod2", "Method 1: distance between track and original cluster ;DCA xy (cm); ", 100, -5, 5);
  mDCAzOriginalMethod2 = std::make_unique<TH1D>("dcaZOriginalMethod2", "Method 1: distance between track and original cluster ;DCA z (cm); ", 100, -5, 5);
  mDCAOriginalMethod2 = std::make_unique<TH1D>("dcaOriginalMethod2", "Method 1: distance between track and original cluster ;DCA (cm); ", 100, -5, 5);
  mDCAxyDuplicatedMethod2 = std::make_unique<TH1D>("dcaXYDuplicatedMethod2", "Method 1: distance between track and duplicated cluster  ;DCA xy (cm); ", 100, -5, 5);
  mDCAzDuplicatedMethod2 = std::make_unique<TH1D>("dcaZDuplicatedMethod2", "Method 1: distance between track and duplicated cluster  ;DCA z (cm); ", 100, -5, 5);
  mDCADuplicatedMethod2 = std::make_unique<TH1D>("dcaDuplicatedMethod2", "Method 1: distance between track and duplicated cluster  ;DCA (cm); ", 100, 0, 10);

  mPhi = std::make_unique<TH1D>("phi", ";phi (deg); ", 360, 0, 360);
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
  mTracksROFRecords = recoData.getITSTracksROFRecords();
  mTracks = recoData.getITSTracks();
  mTracksMCLabels = recoData.getITSTracksMCLabels();
  mClusters = recoData.getITSClusters();
  mClustersROFRecords = recoData.getITSClustersROFRecords();
  mClustersMCLCont = recoData.getITSClustersMCLabels();
  mClusPatterns = recoData.getITSClustersPatterns();
  mInputITSidxs = recoData.getITSTracksClusterRefs();
  mITSClustersArray.reserve(mClusters.size());
  auto pattIt = mClusPatterns.begin();
  o2::its::ioutils::convertCompactClusters(mClusters, pattIt, mITSClustersArray, mDict);
}

void EfficiencyStudy::process(o2::globaltracking::RecoContainer& recoData)
{
  mOutFile = std::make_unique<TFile>(mOutFileName.c_str(), "recreate");
  mOutFile->mkdir("DistanceClusters/");
  mOutFile->mkdir("Method 1/");
  mOutFile->mkdir("Method 2/");

  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  o2::itsmft::ClusterPattern patt;
  o2::itsmft::ClusterPattern pattoriginal;
  o2::itsmft::ClusterPattern pattduplicated;
  auto pattIt = mClusPatterns.begin();
  auto pattItoriginal = mClusPatterns.begin();
  auto pattItduplicated = mClusPatterns.begin();
  o2::gpu::gpustd::array<float, 2> clusOriginalDCA, clusDuplicatedDCA;

  ////starting from tracks
  LOGP(info, "--------------- process");

  LOGP(info, "** Found in {} rofs:\n\t- {} clusters\n\t",
       mClustersROFRecords.size(), mClusters.size());

  LOGP(info, "mClusters size: {}, mClustersROFRecords size: {}, mClustersMCLCont size: {}, mClustersconverted size: {} ", mClusters.size(), mClustersROFRecords.size(), mClustersMCLCont->getNElements(), mITSClustersArray.size());
  LOGP(info, "mTracks size: {}, mTracksROFRecords size: {}, mTracksMCLabels size: {}", mTracks.size(), mTracksROFRecords.size(), mTracksMCLabels.size());

  short unsigned int nLayers = 3;
  unsigned int rofIndexTrack = 0;
  unsigned int rofNEntriesTrack = 0;
  unsigned int rofIndexClus = 0;
  unsigned int rofNEntriesClus = 0;
  unsigned int totalRofEntries = 0;
  int nLabels = 0;
  unsigned int totPrecClus = 0;

  std::unordered_map<o2::MCCompLabel, std::vector<int>> label_vecClus[mClustersROFRecords.size()][nLayers]; // array of maps nRofs x Nlayers -> {label, vec(iClus)} where vec(iClus) are the clusters that share the same label

  for (unsigned int iROF = 0; iROF < mTracksROFRecords.size(); iROF++) { // loop on ROFRecords array
    rofIndexTrack = mTracksROFRecords[iROF].getFirstEntry();
    rofNEntriesTrack = mTracksROFRecords[iROF].getNEntries();

    rofIndexClus = mClustersROFRecords[iROF].getFirstEntry();
    rofNEntriesClus = mClustersROFRecords[iROF].getNEntries();

    for (unsigned int iTrack = rofIndexTrack; iTrack < rofIndexTrack + rofNEntriesTrack; iTrack++) { // loop on tracks per ROF
      auto track = mTracks[iTrack];
      int firstClus = track.getFirstClusterEntry(); // get the first cluster of the track
      int ncl = track.getNumberOfClusters();        // get the number of clusters of the track

      auto& tracklab = mTracksMCLabels[iTrack];
      if (tracklab.isFake())
        continue;

      for (int iclTrack = firstClus; iclTrack < firstClus + ncl; iclTrack++) { // loop on clusters associated to the track
        auto& clus = mClusters[mInputITSidxs[iclTrack]];

        auto layer = mGeometry->getLayer(clus.getSensorID());
        if (layer >= nLayers)
          continue;                                                            // checking only selected layers
        auto labsTrack = mClustersMCLCont->getLabels(mInputITSidxs[iclTrack]); // get labels of clusters associated to the track

        /////////////////
        auto pattID = clus.getPatternID();
        if (pattID == o2::itsmft::CompCluster::InvalidPatternID || mDict->isGroup(pattID)) {
          patt.acquirePattern(pattIt);
        } else {
          patt = mDict->getPattern(pattID);
        }

        o2::math_utils::Point3D<float> locC;
        locC = mDict->getClusterCoordinates(clus, patt, false);
        auto clID = clus.getSensorID();

        LOGP(info, "geom {}", (void*)mGeometry);
        LOGP(info, "cliID: {}, ", clID);
        LOGP(info, "X, Y, Z: {} {} {}", locC.x(), locC.y(), locC.z());
        auto matrixLToG = mGeometry->getMatrixL2G(clID);

        // auto trkXYZ = mGeometry->getMatrixT2L(clID) ^ locC;
        // o2::math_utils::Point3D<float>

        matrixLToG.print();
        LOGP(info, "X, Y, Z: {} {} {}", locC.x(), locC.y(), locC.z());

        auto gloXYZ = mGeometry->getMatrixL2G(clID) * locC;

        LOGP(info, "CHIAO");
        // LOGP(info, "X, Y: {} {} ", gloXYZ.x(), gloXYZ.y());  //////////// !!!!!!!!!!!!!!!! se uso gloXYZ crasha
        // if (gloXYZ.y()!=0 && gloXYZ.x()!=0){
        //   auto phi = o2::gpu::CAMath::ATan2(gloXYZ.y(), gloXYZ.x());
        //   mPhi->Fill(phi);
        // }
        //////////////////

        for (auto& labT : labsTrack) { // for each valid label iterate over ALL the clusters in the ROF to see if there are duplicates
          nLabels++;
          if (labT.isValid()) {
            for (unsigned int iClus = rofIndexClus; iClus < rofIndexClus + rofNEntriesClus; iClus++) { // iteration over ALL the clusters in the ROF
              auto cluster = mClusters[iClus];
              // auto clusPoint = mITSClustersArray[iClus];
              // LOGP(info, "Cluster x, y, z: {} {} {} ", clusPoint.getX(), clusPoint.getY(), clusPoint.getZ());
              // auto phi = std::atan2( clusPoint.getY(), clusPoint.getX()) * 180 / M_PI;
              // mPhi->Fill(phi);
              auto layerClus = mGeometry->getLayer(cluster.getSensorID());
              if (layerClus != layer)
                continue;                                         // we are interested only in duplicated clusters on the same layer
              auto labsClus = mClustersMCLCont->getLabels(iClus); // ideally I can have more than one label per cluster
              for (auto labC : labsClus) {
                if (labT == labC) {
                  label_vecClus[iROF][layerClus][labT].push_back(iClus); // same cluster: label from the track = label from the cluster
                  // if a duplicate cluster is found, propagate the track to the duplicate cluster and compute the distance from the original cluster
                  auto pattIDoriginal = clus.getPatternID();
                  if (pattIDoriginal == o2::itsmft::CompCluster::InvalidPatternID || mDict->isGroup(pattIDoriginal)) {
                    pattoriginal.acquirePattern(pattItoriginal);
                  } else {
                    pattoriginal = mDict->getPattern(pattIDoriginal);
                  }

                  auto pattIDduplicated = cluster.getPatternID();
                  if (pattIDduplicated == o2::itsmft::CompCluster::InvalidPatternID || mDict->isGroup(pattIDduplicated)) {
                    pattduplicated.acquirePattern(pattItduplicated);
                  } else {
                    pattduplicated = mDict->getPattern(pattIDduplicated);
                  }

                  o2::math_utils::Point3D<float> locCoriginal, locCduplicate;
                  locCoriginal = mDict->getClusterCoordinates(clus, pattoriginal, false);
                  o2::math_utils::Point3D<float> trkXYZoriginal = mGeometry->getMatrixT2L(clus.getSensorID()) ^ locCoriginal;

                  locCduplicate = mDict->getClusterCoordinates(cluster, pattduplicated, false);
                  o2::math_utils::Point3D<float> trkXYZduplicate = mGeometry->getMatrixT2L(cluster.getSensorID()) ^ locCduplicate;

                  LOGP(info, "Posizione del cluster: {} {} {}", locCoriginal.x(), locCoriginal.y(), locCoriginal.z());

                  if (locCoriginal != locCduplicate) { /// check that the duplicated cluster is not the original one just counted twice

                    // LOGP(info, "cluster: {} {} {} -> {} {} {}", locCoriginal.x(), locCoriginal.y(), locCoriginal.z(), locCduplicate.x(), locCduplicate.y(), locCduplicate.z());
                    // LOGP(info, "track: {} {} {} -> {} {} {}", trkXYZoriginal.x(), trkXYZoriginal.y(), trkXYZoriginal.z(), trkXYZduplicate.x(), trkXYZduplicate.y(), trkXYZduplicate.z());

                    LOGP(info, "Distanza tra i cluster: {} {} {}", locCoriginal.x() - locCduplicate.x(), locCoriginal.y() - locCduplicate.y(), locCoriginal.z() - locCduplicate.z());
                    mDistanceClustersX->Fill(abs(locCoriginal.x() - locCduplicate.x()));
                    mDistanceClustersY->Fill(abs(locCoriginal.y() - locCduplicate.y()));
                    mDistanceClustersZ->Fill(abs(locCoriginal.z() - locCduplicate.z()));
                    mDistanceClusters->Fill(std::hypot(locCoriginal.x() - locCduplicate.x(), locCoriginal.y() - locCduplicate.y(), locCoriginal.z() - locCduplicate.z()));

                    /// method 1 //////////////   -> propagate the track to the cluster and manually compute the distance between the track anda the cluster
                    LOGP(info, "track prima: {} {} {} ----- clsXYZ: {}, {}, {} ---> DCA: {} {} {}", track.getX(), track.getY(), track.getZ(), trkXYZoriginal.x(), trkXYZoriginal.y(), trkXYZoriginal.z(), trkXYZoriginal.x() - track.getX(), trkXYZoriginal.y() - track.getY(), trkXYZoriginal.z() - track.getZ());

                    track.propagateTo(trkXYZoriginal.x(), b);
                    LOGP(info, "track dopo: {} {} {} ----- clsXYZ: {}, {}, {} ---> -> DCA: {} {} {} ", track.getX(), track.getY(), track.getZ(), trkXYZoriginal.x(), trkXYZoriginal.y(), trkXYZoriginal.z(), trkXYZoriginal.x() - track.getX(), trkXYZoriginal.y() - track.getY(), trkXYZoriginal.z() - track.getZ());
                    mDCAxOriginalMethod1->Fill(abs(trkXYZoriginal.x() - track.getX()));
                    mDCAyOriginalMethod1->Fill(abs(trkXYZoriginal.y() - track.getY()));
                    mDCAzOriginalMethod1->Fill(abs(trkXYZoriginal.z() - track.getZ()));
                    mDCAOriginalMethod1->Fill(std::hypot(trkXYZoriginal.x() - track.getX(), trkXYZoriginal.y() - track.getY(), trkXYZoriginal.z() - track.getZ()));

                    track.propagateTo(trkXYZduplicate.x(), b);
                    LOGP(info, "track dopo ancora: {} {} {} ----- clsXYZ: {}, {}, {} ---> -> DCA: {} {} {}", track.getX(), track.getY(), track.getZ(), trkXYZduplicate.x(), trkXYZduplicate.y(), trkXYZduplicate.z(), trkXYZduplicate.x() - track.getX(), trkXYZduplicate.y() - track.getY(), trkXYZduplicate.z() - track.getZ());
                    mDCAxDuplicatedMethod1->Fill(abs(trkXYZduplicate.x() - track.getX()));
                    mDCAyDuplicatedMethod1->Fill(abs(trkXYZduplicate.y() - track.getY()));
                    mDCAzDuplicatedMethod1->Fill(abs(trkXYZduplicate.z() - track.getZ()));
                    mDCADuplicatedMethod1->Fill(std::hypot(trkXYZduplicate.x() - track.getX(), trkXYZduplicate.y() - track.getY(), trkXYZduplicate.z() - track.getZ()));
                    ///////////////////////////////

                    /// method 2 ///////////////// -> returns directly the DCA between the cluster location in the track reference system and the track
                    o2::track::TrackParCov trackParCov = mTracks[iTrack];
                    auto propagator = o2::base::Propagator::Instance();
                    LOGP(info, "track prima: {} {} {}", trackParCov.getX(), trackParCov.getY(), trackParCov.getZ());
                    if (propagator->propagateToDCA(trkXYZoriginal, trackParCov, b, 2.f, matCorr, &clusOriginalDCA)) {
                      LOGP(info, "track dopo: {} {} {} --> DCA: {} {}", trackParCov.getX(), trackParCov.getY(), trackParCov.getZ(), clusOriginalDCA[0], clusOriginalDCA[1]);
                      mDCAxyOriginalMethod2->Fill(clusOriginalDCA[0]);
                      mDCAzOriginalMethod2->Fill(clusOriginalDCA[1]);
                      mDCAOriginalMethod2->Fill(std::hypot(clusOriginalDCA[0], clusOriginalDCA[1]));
                    }
                    if (propagator->propagateToDCA(trkXYZduplicate, trackParCov, b, 2.f, matCorr, &clusDuplicatedDCA)) {
                      LOGP(info, "track dopo ancora: {} {} {} --> DCA: {} {}", trackParCov.getX(), trackParCov.getY(), trackParCov.getZ(), clusDuplicatedDCA[0], clusDuplicatedDCA[1]);
                      mDCAxyDuplicatedMethod2->Fill(clusDuplicatedDCA[0]);
                      mDCAzDuplicatedMethod2->Fill(clusDuplicatedDCA[1]);
                      mDCADuplicatedMethod2->Fill(std::hypot(clusDuplicatedDCA[0], clusDuplicatedDCA[1]));
                    }
                    //////////////////
                  }
                }
              }
            }
          }
        }
      } // end loop on clusters
      totPrecClus += ncl;
    } // end loop on tracks per ROF
  }   // end loop on ROFRecords array

  mPhi->Write();

  mOutFile->cd("DistanceClusters");
  mDistanceClustersX->Write();
  mDistanceClustersY->Write();
  mDistanceClustersZ->Write();
  mDistanceClusters->Write();

  mOutFile->cd("Method 1");
  mDCAxOriginalMethod1->Write();
  mDCAyOriginalMethod1->Write();
  mDCAzOriginalMethod1->Write();
  mDCAOriginalMethod1->Write();
  mDCAxDuplicatedMethod1->Write();
  mDCAyDuplicatedMethod1->Write();
  mDCAzDuplicatedMethod1->Write();
  mDCADuplicatedMethod1->Write();

  mOutFile->cd("Method 2");
  mDCAxyOriginalMethod2->Write();
  mDCAzOriginalMethod2->Write();
  mDCAOriginalMethod2->Write();
  mDCAxyDuplicatedMethod2->Write();
  mDCAzDuplicatedMethod2->Write();
  mDCADuplicatedMethod2->Write();

  LOGP(info, "Total number of clusters: {} ", totPrecClus);
  LOGP(info, "total rof entries: {}", totalRofEntries);
  LOGP(info, "total nLabels: {}", nLabels);

  // printing the duplicates
  // for (unsigned int iROF = 0; iROF < mClustersROFRecords.size(); iROF++) {
  //   LOGP(info, "°°°°°°°°°°°°°°°°°°°°°°°° ROF {} °°°°°°°°°°°°°°°°°°°°°°°°", iROF);
  //   for (unsigned int lay=0; lay<nLayers; lay++){
  //     LOGP(info, "°°°°°°°°°°°°°°°°°°°°°°°° LAYER {} °°°°°°°°°°°°°°°°°°°°°°°°", lay);
  //     for (auto & it : label_vecClus[iROF][lay]) {
  //       if (it.second.size()<=1) continue;  // printing only duplicates
  //       std::cout<<" \n++++++++++++ Label: ";
  //       auto label = it.first;
  //       it.first.print();
  //       for (auto iClus : it.second) {
  //         auto name = mGeometry->getSymbolicName(mClusters[iClus].getSensorID());
  //         auto chipid = mClusters[iClus].getChipID();

  //         auto pattID = mClusters[iClus].getPatternID();
  //         if (pattID == o2::itsmft::CompCluster::InvalidPatternID || mDict->isGroup(pattID)) {
  //           patt.acquirePattern(pattIt);
  //         } else {
  //           patt = mDict->getPattern(pattID);
  //         }

  //         o2::math_utils::Point3D<float> locC;
  //         locC= mDict->getClusterCoordinates(mClusters[iClus], patt, false);
  //         std::cout<< "ROF: "<<iROF<<", iClus: "<<iClus<<" -> chip: "<< chipid <<" = "<<name<<std::endl;
  //         LOGP(info, "LOCC: {} {} {}", locCoriginal.x(), locCoriginal.y(), locCoriginal.z());

  //       }
  //     }
  //   }
  // }
}

void EfficiencyStudy::updateTimeDependentParams(ProcessingContext& pc)
{
  static bool initOnceDone = false;
  o2::base::GRPGeomHelper::instance().checkUpdates(pc);
  if (!initOnceDone) { // this params need to be queried only once
    initOnceDone = true;
    mGeometry = GeometryTGeo::Instance();
    mGeometry->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::T2GRot, o2::math_utils::TransformType::T2G));
  }
}

void EfficiencyStudy::endOfStream(EndOfStreamContext& ec)
{
  // saveOutput();
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