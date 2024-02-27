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
  bool mVerboseOutput = false;
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
  // std::unique_ptr<TH1D> mDCAxOriginalMethod1;
  // std::unique_ptr<TH1D> mDCAyOriginalMethod1;
  // std::unique_ptr<TH1D> mDCAzOriginalMethod1;
  // std::unique_ptr<TH1D> mDCAOriginalMethod1;
  // Distance betweeen track and duplicated cluster, method 1
  // std::unique_ptr<TH1D> mDCAxDuplicatedMethod1;
  // std::unique_ptr<TH1D> mDCAyDuplicatedMethod1;
  // std::unique_ptr<TH1D> mDCAzDuplicatedMethod1;
  // std::unique_ptr<TH1D> mDCADuplicatedMethod1;
  // DCA betweeen track and original cluster, method 2
  std::unique_ptr<TH1D> mDCAxyOriginalMethod2;
  std::unique_ptr<TH1D> mDCAzOriginalMethod2;
  // DCA betweeen track and duplicated cluster, method 2
  std::unique_ptr<TH1D> mDCAxyDuplicatedMethod2;
  std::unique_ptr<TH1D> mDCAzDuplicatedMethod2;
  // phi of the cluster
  std::unique_ptr<TH1D> mPhi;
  // position of the clusters
  std::unique_ptr<TH3D> m3DClusterPositions;
  std::unique_ptr<TH2D> m2DClusterPositions;
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

  // mDCAxOriginalMethod1 = std::make_unique<TH1D>("dcaXOriginalMethod1", "Method 1: distance between track and original cluster ;DCA x (cm); ", 100, -0.05, 0.05);
  // mDCAyOriginalMethod1 = std::make_unique<TH1D>("dcaYOriginalMethod1", "Method 1: distance between track and original cluster ;DCA y (cm); ", 100, -0.05, 0.05);
  // mDCAzOriginalMethod1 = std::make_unique<TH1D>("dcaZOriginalMethod1", "Method 1: distance between track and original cluster ;DCA z (cm); ", 100, -0.05, 0.05);
  // mDCAOriginalMethod1 = std::make_unique<TH1D>("dcaOriginalMethod1", "Method 1: distance between track and original cluster ;DCA (cm); ",  100, -0.05, 0.05);
  // mDCAxDuplicatedMethod1 = std::make_unique<TH1D>("dcaXDuplicatedMethod1", "Method 1: distance between track and duplicated cluster  ;DCA x (cm); ",  100, -0.05, 0.05);
  // mDCAyDuplicatedMethod1 = std::make_unique<TH1D>("dcaYDuplicatedMethod1", "Method 1: distance between track and duplicated cluster  ;DCA y (cm); ",  100, -0.05, 0.05);
  // mDCAzDuplicatedMethod1 = std::make_unique<TH1D>("dcaZDuplicatedMethod1", "Method 1: distance between track and duplicated cluster  ;DCA z (cm); ",  100, -0.05, 0.05);
  // mDCADuplicatedMethod1 = std::make_unique<TH1D>("dcaDuplicatedMethod1", "Method 1: distance between track and duplicated cluster  ;DCA (cm); ",  100, -0.05, 0.05);

  mDCAxyOriginalMethod2 = std::make_unique<TH1D>("dcaXYOriginalMethod2", "Method 1: distance between track and original cluster ;DCA xy (cm); ", 100, -0.05, 0.05);
  mDCAzOriginalMethod2 = std::make_unique<TH1D>("dcaZOriginalMethod2", "Method 1: distance between track and original cluster ;DCA z (cm); ", 100, -0.05, 0.05);
  mDCAxyDuplicatedMethod2 = std::make_unique<TH1D>("dcaXYDuplicatedMethod2", "Method 1: distance between track and duplicated cluster  ;DCA xy (cm); ", 100, -0.05, 0.05);
  mDCAzDuplicatedMethod2 = std::make_unique<TH1D>("dcaZDuplicatedMethod2", "Method 1: distance between track and duplicated cluster  ;DCA z (cm); ", 100, -0.05, 0.05);

  mPhi = std::make_unique<TH1D>("phi", ";phi (deg); ", 120, -180, 180);

  m3DClusterPositions = std::make_unique<TH3D>("3DClusterPositions", ";x (cm);y (cm);z (cm)", 200, -10, 10, 200, -10, 10, 400, -20, 20);
  m2DClusterPositions = std::make_unique<TH2D>("3DClusterPositions", ";x (cm);y (cm)", 200, -10, 10, 200, -6, 6);
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
  o2::its::ioutils::convertCompactClusters(mClusters, pattIt, mITSClustersArray, mDict);  // clusters converted to 3D spacepoints
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
      o2::track::TrackParCov trackParCov = mTracks[iTrack];
      int firstClus = track.getFirstClusterEntry(); // get the first cluster of the track
      int ncl = track.getNumberOfClusters();        // get the number of clusters of the track

      auto& tracklab = mTracksMCLabels[iTrack];
      if (tracklab.isFake())
        continue;

      for (int iclTrack = firstClus; iclTrack < firstClus + ncl; iclTrack++) { // loop on clusters associated to the track
        auto& clusOriginal = mClusters[mInputITSidxs[iclTrack]];
        auto clusOriginalPoint = mITSClustersArray[mInputITSidxs[iclTrack]]; // cluster spacepoint in the tracking system

        auto layer = mGeometry->getLayer(clusOriginal.getSensorID());
        if (layer >= nLayers)
          continue;                                                            // checking only selected layers
        auto labsTrack = mClustersMCLCont->getLabels(mInputITSidxs[iclTrack]); // get labels of clusters associated to the track

        o2::math_utils::Point3D<float> clusOriginalPointTrack = {clusOriginalPoint.getX(), clusOriginalPoint.getY(), clusOriginalPoint.getZ()};
        o2::math_utils::Point3D<float> clusOriginalPointGlob = mGeometry->getMatrixT2G(clusOriginal.getSensorID()) * clusOriginalPointTrack;

        auto phiOriginal = std::atan2( clusOriginalPointGlob.y(), clusOriginalPointGlob.x()) * 180 / M_PI;
        mPhi->Fill(phiOriginal);
        m3DClusterPositions->Fill(clusOriginalPointGlob.x(), clusOriginalPointGlob.y(), clusOriginalPointGlob.z());
        m2DClusterPositions->Fill(clusOriginalPointGlob.x(), clusOriginalPointGlob.y());

        for (auto& labT : labsTrack) { // for each valid label iterate over ALL the clusters in the ROF to see if there are duplicates
          nLabels++;
          if (labT.isValid()) {
            for (unsigned int iClus = rofIndexClus; iClus < rofIndexClus + rofNEntriesClus; iClus++) { // iteration over ALL the clusters in the ROF
              auto clusDuplicated = mClusters[iClus];
              auto clusDuplicatedPoint = mITSClustersArray[iClus];
              
              auto layerClus = mGeometry->getLayer(clusDuplicated.getSensorID());
              if (layerClus != layer)
                continue;                                         // we are interested only in duplicated clusters on the same layer
              
              o2::math_utils::Point3D<float> clusDuplicatedPointTrack = {clusDuplicatedPoint.getX(), clusDuplicatedPoint.getY(), clusDuplicatedPoint.getZ()};
              o2::math_utils::Point3D<float> clusDuplicatedPointGlob = mGeometry->getMatrixT2G(clusDuplicated.getSensorID()) * clusDuplicatedPointTrack;

              auto phiDuplicated = std::atan2( clusDuplicatedPointGlob.y(), clusDuplicatedPointGlob.x()) * 180 / M_PI;
              mPhi->Fill(phiDuplicated);
              m3DClusterPositions->Fill(clusDuplicatedPointGlob.x(), clusDuplicatedPointGlob.y(), clusDuplicatedPointGlob.z());
              m2DClusterPositions->Fill(clusDuplicatedPointGlob.x(), clusDuplicatedPointGlob.y());
              

              auto labsClus = mClustersMCLCont->getLabels(iClus); // ideally I can have more than one label per cluster
              for (auto labC : labsClus) {
                if (labT == labC) {
                  label_vecClus[iROF][layerClus][labT].push_back(iClus); // same cluster: label from the track = label from the cluster
                  // if a duplicate cluster is found, propagate the track to the duplicate cluster and compute the distance from the original cluster
                  if (clusOriginalPointGlob != clusDuplicatedPointGlob) { /// check that the duplicated cluster is not the original one just counted twice

                    /// compute the distance between original and dubplicated cluster
                    mDistanceClustersX->Fill(abs(clusOriginalPointGlob.x() - clusDuplicatedPointGlob.x()));
                    mDistanceClustersY->Fill(abs(clusOriginalPointGlob.y() - clusDuplicatedPointGlob.y()));
                    mDistanceClustersZ->Fill(abs(clusOriginalPointGlob.z() - clusDuplicatedPointGlob.z()));
                    mDistanceClusters->Fill(std::hypot(clusOriginalPointGlob.x() - clusDuplicatedPointGlob.x(), clusOriginalPointGlob.y() - clusDuplicatedPointGlob.y(), clusOriginalPointGlob.z() - clusDuplicatedPointGlob.z()));

                    /// method 1 //////////////   -> propagate the track to the cluster and manually compute the distance between the track and the cluster
                    /// first propagate to the original cluster
                    // track.rotate(mGeometry->getSensorRefAlpha(clusOriginal.getSensorID()));
                    // if (track.propagateTo(clusOriginalPointTrack.x(), b)){
                    //   mDCAxOriginalMethod1->Fill(abs(clusOriginalPointTrack.x() - track.getX()));
                    //   mDCAyOriginalMethod1->Fill(abs(clusOriginalPointTrack.y() - track.getY()));
                    //   mDCAzOriginalMethod1->Fill(abs(clusOriginalPointTrack.z() - track.getZ()));
                    //   mDCAOriginalMethod1->Fill(std::hypot(clusOriginalPointTrack.x() - track.getX(), clusOriginalPointTrack.y() - track.getY(), clusOriginalPointTrack.z() - track.getZ()));
                    // }
                    // // then propagate to the duplicated cluster
                    // track.rotate(mGeometry->getSensorRefAlpha(clusDuplicated.getSensorID()));
                    // if (track.propagateTo(clusDuplicatedPointTrack.x(), b)){
                    //   mDCAxDuplicatedMethod1->Fill(abs(clusDuplicatedPointTrack.x() - track.getX()));
                    //   mDCAyDuplicatedMethod1->Fill(abs(clusDuplicatedPointTrack.y() - track.getY()));
                    //   mDCAzDuplicatedMethod1->Fill(abs(clusDuplicatedPointTrack.z() - track.getZ()));
                    //   mDCADuplicatedMethod1->Fill(std::hypot(clusDuplicatedPointTrack.x() - track.getX(), clusDuplicatedPointTrack.y() - track.getY(), clusDuplicatedPointTrack.z() - track.getZ()));
                    // }
                    ///////////////////////////////////////////////////////

                    /// method 2 ///////////////// -> returns directly the DCA between the cluster location and the track
                    auto propagator = o2::base::Propagator::Instance();
                    /// first propagate to the original cluster
                    trackParCov.rotate(mGeometry->getSensorRefAlpha(clusOriginal.getSensorID()));
                    if (propagator->propagateToDCA(clusOriginalPointGlob, trackParCov, b, 2.f, matCorr, &clusOriginalDCA)) {
                      mDCAxyOriginalMethod2->Fill(clusOriginalDCA[0]);
                      mDCAzOriginalMethod2->Fill(clusOriginalDCA[1]);
                    }
                    /// then propagate to the duplicated cluster
                    trackParCov.rotate(mGeometry->getSensorRefAlpha(clusDuplicated.getSensorID()));
                    if (propagator->propagateToDCA(clusDuplicatedPointGlob, trackParCov, b, 2.f, matCorr, &clusDuplicatedDCA)) {
                      mDCAxyDuplicatedMethod2->Fill(clusDuplicatedDCA[0]);
                      mDCAzDuplicatedMethod2->Fill(clusDuplicatedDCA[1]);
                    }
                    ///////////////////////////////////////////////////////
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
  m3DClusterPositions->Write();
  m2DClusterPositions->Write();

  mOutFile->cd("DistanceClusters");
  mDistanceClustersX->Write();
  mDistanceClustersY->Write();
  mDistanceClustersZ->Write();
  mDistanceClusters->Write();

  // mOutFile->cd("Method 1");
  // mDCAxOriginalMethod1->Write();
  // mDCAyOriginalMethod1->Write();
  // mDCAzOriginalMethod1->Write();
  // mDCAOriginalMethod1->Write();
  // mDCAxDuplicatedMethod1->Write();
  // mDCAyDuplicatedMethod1->Write();
  // mDCAzDuplicatedMethod1->Write();
  // mDCADuplicatedMethod1->Write();

  mOutFile->cd("Method 2");
  mDCAxyOriginalMethod2->Write();
  mDCAzOriginalMethod2->Write();
  mDCAxyDuplicatedMethod2->Write();
  mDCAzDuplicatedMethod2->Write();

  LOGP(info, "Total number of clusters: {} ", totPrecClus);
  LOGP(info, "total rof entries: {}", totalRofEntries);
  LOGP(info, "total nLabels: {}", nLabels);

  if (mVerboseOutput){
  // printing the duplicates
    for (unsigned int iROF = 0; iROF < mClustersROFRecords.size(); iROF++) {
      LOGP(info, "°°°°°°°°°°°°°°°°°°°°°°°° ROF {} °°°°°°°°°°°°°°°°°°°°°°°°", iROF);
      for (unsigned int lay=0; lay<nLayers; lay++){
        LOGP(info, "°°°°°°°°°°°°°°°°°°°°°°°° LAYER {} °°°°°°°°°°°°°°°°°°°°°°°°", lay);
        for (auto & it : label_vecClus[iROF][lay]) {
          if (it.second.size()<=1) continue;  // printing only duplicates
          std::cout<<" \n++++++++++++ Label: ";
          auto label = it.first;
          it.first.print();
          for (auto iClus : it.second) {
            auto name = mGeometry->getSymbolicName(mClusters[iClus].getSensorID());
            auto chipid = mClusters[iClus].getChipID();
            auto clus = mClusters[iClus];
            auto clusPoint = mITSClustersArray[iClus];

            o2::math_utils::Point3D<float> clusPointTrack = {clusPoint.getX(), clusPoint.getY(), clusPoint.getZ()};
            o2::math_utils::Point3D<float> clusPointGlob = mGeometry->getMatrixT2G(clus.getSensorID()) * clusPointTrack;
            std::cout<< "ROF: "<<iROF<<", iClus: "<<iClus<<" -> chip: "<< chipid <<" = "<<name<<std::endl;
            LOGP(info, "LOCtrack: {} {} {}", clusPointTrack.x(), clusPointTrack.y(), clusPointTrack.z());
            LOGP(info, "LOCglob {} {} {}", clusPointGlob.x(), clusPointGlob.y(), clusPointGlob.z());
          }
        }
      }
    }
  }
}

void EfficiencyStudy::updateTimeDependentParams(ProcessingContext& pc)
{
  static bool initOnceDone = false;
  o2::base::GRPGeomHelper::instance().checkUpdates(pc);
  if (!initOnceDone) { // this params need to be queried only once
    initOnceDone = true;
    mGeometry = GeometryTGeo::Instance();
    mGeometry->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::T2GRot, o2::math_utils::TransformType::T2G,o2::math_utils::TransformType::L2G ));
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