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
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "Framework/Task.h"
#include "ITSBase/GeometryTGeo.h"
#include "SimulationDataFormat/MCTrack.h"
#include "Steer/MCKinematicsReader.h"

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
  struct ParticleInfo {
    int event;
    int pdg;
    float pt;
    float eta;
    float phi;
    int mother;
    int first;
    float vx;
    float vy;
    float vz;
    int unsigned short clusters = 0u;
    unsigned char isReco = 0u;
    unsigned char isFake = 0u;
    bool isPrimary = false;
    unsigned char storedStatus = 2; /// not stored = 2, fake = 1, good = 0
    const char* prodProcessName;
    int prodProcess;
    o2::its::TrackITS track;
  };

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

 private:
  void updateTimeDependentParams(ProcessingContext& pc);
  std::string mOutFileName = "ITS_efficiencyStudy.root";
  std::shared_ptr<o2::steer::MCKinematicsReader> mKineReader;
  GeometryTGeo* mGeometry;

  // Spans
  gsl::span<const o2::itsmft::ROFRecord> mTracksROFRecords;
  gsl::span<const o2::itsmft::ROFRecord> mClustersROFRecords;
  gsl::span<const o2::its::TrackITS> mTracks;
  gsl::span<const o2::MCCompLabel> mTracksMCLabels;
  gsl::span<const o2::itsmft::CompClusterExt> mClusters;
  gsl::span<const int> mInputITSidxs;
  const o2::dataformats::MCLabelContainer* mClustersMCLCont;

  // Data
  GTrackID::mask_t mTracksSrc{};
  std::shared_ptr<DataRequest> mDataRequest;
  std::vector<std::vector<std::vector<ParticleInfo>>> mParticleInfo; // src/event/track
  unsigned short mMask = 0x7f;

  // Utils
  std::shared_ptr<o2::base::GRPGeomRequest> mGGCCDBRequest;
};

void EfficiencyStudy::init(InitContext& ic)
{
  LOGP(info, "--------------- init");

  o2::base::GRPGeomHelper::instance().setRequest(mGGCCDBRequest);

  auto& pars = o2::its::study::ITSCheckTracksParamConfig::Instance();
  mOutFileName = pars.outFileName;
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
  mInputITSidxs = recoData.getITSTracksClusterRefs();

  // LOGP(info, "** Found in {} rofs:\n\t- {} clusters\n\t",
  //       mClustersROFRecords.size(),  mClusters.size());
}

void EfficiencyStudy::process(o2::globaltracking::RecoContainer& recoData)
{
  LOGP(info, "--------------- process");

  // std::unordered_set<o2::MCCompLabel> clsLabels;
  // std::unordered_map<o2::MCCompLabel, std::tuple<int,std::vector<o2::itsmft::CompClusterExt>>> mapLab_countClus;
  std::unordered_map<o2::MCCompLabel, std::vector<o2::itsmft::CompClusterExt>> mapLab_vecClus;
  std::unordered_map<o2::MCCompLabel, std::vector<std::tuple<o2::itsmft::CompClusterExt,unsigned int >>> mapLab_vecClusROF;
  std::unordered_map<o2::MCCompLabel, std::unordered_set<int>> mapLab_setClus;
  std::unordered_map<o2::MCCompLabel, std::unordered_set<int>> mapLab_setROF;

  std::unordered_map<o2::MCCompLabel, std::unordered_set<std::string>> mapLab_setClusROF;
  // std::unordered_map<o2::MCCompLabel, std::unordered_set<std::tuple<int,unsigned int >>> mapLab_setClusROF;

  std::unordered_map<int, std::vector<std::tuple<o2::MCCompLabel, o2::itsmft::CompClusterExt, unsigned int>>> mapLay_labClusRof;  // {layer, vec(label, cluster, ROF)}

  unsigned int rofIndex = 0;
  unsigned int rofNEntries = 0;

  for (unsigned int iROF = 0; iROF < mTracksROFRecords.size(); iROF++) {
    rofIndex = mTracksROFRecords[iROF].getFirstEntry();
    rofNEntries = mTracksROFRecords[iROF].getNEntries();
    LOGP(info, "\n");
    LOGP(info, "Track ROF: {}, entries: {}", iROF, rofNEntries);

    for (unsigned int iTrack=rofIndex; iTrack < rofIndex+rofNEntries; iTrack++){
      auto ITStrack = recoData.getITSTrack(iTrack);

      int firstClus = ITStrack.getFirstClusterEntry(); // get the first cluster of the track
      int ncl = ITStrack.getNumberOfClusters();        // get the number of clusters of the track


      auto& tracklab = mTracksMCLabels[iTrack];
      if (tracklab.isFake())
        continue;

      LOGP(info, "\n");
      LOGP(info, "Track number: {}, track ID: {}, event ID: {}, source ID: {}, isFake: {}, N Clusters: {}", iTrack, tracklab.getTrackID(), tracklab.getEventID(), tracklab.getSourceID(), tracklab.isFake(), ncl);
      LOGP(info, "From cluster labels: {}, {} ", iTrack, ncl);

      for (int icl = 0; icl < ncl; icl++) { // loop on clusters associated to the track
        auto& clus = mClusters[mInputITSidxs[firstClus + icl]];
        auto layer = mGeometry->getLayer(clus.getSensorID());

        std::tuple<o2::MCCompLabel, o2::itsmft::CompClusterExt, unsigned int> lab_clus;
        auto labs = mClustersMCLCont->getLabels(mInputITSidxs[firstClus + icl]); // get labels of clusters associated to the track
        for (auto lab : labs) {
          lab_clus = std::make_tuple(lab, clus, iROF);
          mapLay_labClusRof[layer].push_back(lab_clus);
          LOGP(info, "Cluster # {}: Track ID from clus: {}, Track ID from track: {}, Event ID:{}, Source ID:{}, chip: {}, index {}, rof: {}", icl, lab.getTrackID(), tracklab.getTrackID(), lab.getEventID(), lab.getSourceID(), clus.getChipID(), mInputITSidxs[firstClus + icl], iROF);
        }

      } // end loop on clusters



    } // end loop on tracks

  } // end loop on ROFRecords array
  



////////////////////////////

// for (unsigned int iTrack{0}; iTrack < mTracks.size(); ++iTrack) { // loop over the tracks
//   auto ITStrack = recoData.getITSTrack(iTrack);                   // get the ITS track
//   // auto ITStrackROFRecords = recoData.getITSTracksROFRecords();
//   // LOGP(info, "size: {}", ITStrackROFRecords.size());

//   int firstClus = ITStrack.getFirstClusterEntry(); // get the first cluster of the track
//   int ncl = ITStrack.getNumberOfClusters();        // get the number of clusters of the track

//   auto& tracklab = mTracksMCLabels[iTrack];
//   if (tracklab.isFake())
//     continue;

//   LOGP(info, "\n");
//   LOGP(info, "Track number: {}, track ID: {}, event ID: {}, source ID: {}, isFake: {}, N Clusters: {}", iTrack, tracklab.getTrackID(), tracklab.getEventID(), tracklab.getSourceID(), tracklab.isFake(), ncl);
//   LOGP(info, "From cluster labels: {}, {} ", iTrack, ncl);

//   for (int icl = 0; icl < ncl; icl++) { // loop on clusters associated to the track
//     auto& clus = mClusters[mInputITSidxs[firstClus + icl]];
//     auto layer = mGeometry->getLayer(clus.getSensorID());

//     std::tuple<o2::MCCompLabel, o2::itsmft::CompClusterExt> lab_clus;
//     auto labs = mClustersMCLCont->getLabels(mInputITSidxs[firstClus + icl]); // get labels of clusters associated to the track
//     for (auto lab : labs) {
//       lab_clus = std::make_tuple(lab, clus);
//       mapLay_labClus[layer].push_back(lab_clus);
//       LOGP(info, "Cluster # {}: Track ID from clus: {}, Track ID from track: {}, Event ID:{}, Source ID:{}, chip: {}, index {}", icl, lab.getTrackID(), tracklab.getTrackID(), lab.getEventID(), lab.getSourceID(), clus.getChipID(), mInputITSidxs[firstClus + icl]);
//     }

//   } // loop on clusters
// } // loop on tracks

// Going layer by layer
for (auto labclsROFVector : mapLay_labClusRof) {                                       // loop on layers -> each value is a vector of tuples (label, std::vector<cluster>)
  mapLab_vecClus.clear();                                                        // reusing the same map for each layer
  mapLab_vecClusROF.clear();                                                        // reusing the same map for each layer
  mapLab_setClus.clear();                                                        // reusing the same map for each layer
  mapLab_setROF.clear();                                                        // reusing the same map for each layer
  mapLab_setClusROF.clear();                                                        // reusing the same map for each layer
  // mapLab_setClusROF.clear();                                                        // reusing the same map for each layer
  for (auto labcls : labclsROFVector.second) {                                      // loop on vector of tuples (label, cluster) -> each element is a tuple (label, cluster)
    mapLab_vecClus[std::get<0>(labcls)].push_back(std::get<1>(labcls));          // mapLab_vecClus[label].push_back(cluster)
    mapLab_vecClusROF[std::get<0>(labcls)].push_back(std::make_tuple(std::get<1>(labcls), std::get<2>(labcls)));          // mapLab_vecClus[label].push_back(cluster)
    mapLab_setClus[std::get<0>(labcls)].insert(std::get<1>(labcls).getChipID()); // mapLab_setClus[label].insert(chipid)
    mapLab_setROF[std::get<0>(labcls)].insert(std::get<2>(labcls)); // mapLab_setClus[label].insert(chipid)

    auto layer = mGeometry->getLayer(std::get<1>(labcls).getSensorID());
    auto stave = mGeometry->getStave(std::get<1>(labcls).getSensorID());
    auto chipInStave1 = mGeometry->getChipIdInStave(std::get<1>(labcls).getSensorID());
    auto mod = mGeometry->getModule(std::get<1>(labcls).getSensorID());
    auto hs = mGeometry->getHalfStave(std::get<1>(labcls).getSensorID());
    auto name = mGeometry->getSymbolicName(std::get<1>(labcls).getSensorID());
    auto chipInStave2 = mGeometry->getChipIdInModule(std::get<1>(labcls).getSensorID());
    // std::string chipROF ="chip: "+std::to_string(std::get<1>(labcls).getChipID()) + "---> L" + std::to_string(layer) + "_S" + std::to_string(stave) + "_C" + std::to_string(chipInStave1)+"/"+std::to_string(chipInStave2) + " ROF:"+std::to_string(std::get<2>(labcls));
    std::string chipROF ="chip: "+std::to_string(std::get<1>(labcls).getChipID()) + "---> L" + name + " ROF:"+std::to_string(std::get<2>(labcls));
    // std::string clusROF = std::to_string(std::get<1>(labcls).getChipID()) +"_"+std::to_string(std::get<2>(labcls));
    // mapLab_setClusROF[std::get<0>(labcls)].insert(clusROF); // mapLab_setClus[label].insert(chipid)
    mapLab_setClusROF[std::get<0>(labcls)].insert(chipROF); // mapLab_setClus[label].insert(chipid)
    
    // mapLab_setClusROF[std::get<0>(labcls)].insert(std::make_tuple(std::get<1>(labcls).getChipID(), std::get<2>(labcls))); // mapLab_setClus[label].insert(chipid)
  }

  // if (labclsVector.first!=0) continue;

  // // Printing the number of clusters per label
  LOGP(info, "Layer {}:", labclsROFVector.first);
  // for (auto it : mapLab_vecClusROF) {
  //   LOGP(info, "Label: Track ID: {}, Event ID: {}, Source ID: {}, Count: {}",
  //        it.first.getTrackID(),
  //        it.first.getEventID(),
  //        it.first.getSourceID(),
  //        it.second.size());
  //   for (auto clus : it.second) {
  //     LOGP(info, "chipID: {}, ROF: {}", std::get<0>(clus).getChipID(), std::get<1>(clus));
  //   }
  // }

  
  //// print with the rof
  for (auto it: mapLab_setClusROF){
    if (it.second.size()==1) continue;
    LOGP(info, "Label: Track ID: {}, Event ID: {}, Source ID: {}, Count: {}",
        it.first.getTrackID(),
        it.first.getEventID(),
        it.first.getSourceID(),
        it.second.size());
    for (auto clusrof: it.second){
      LOGP(info, "Chip_ROF: {}", clusrof);
    }
  }

  
  //// print without the rof
  LOGP(info, "Only cluster repeated in different chips:\n");
  for (auto it: mapLab_setClus){
    if (it.second.size()==1) continue;
    LOGP(info, "Label: Track ID: {}, Event ID: {}, Source ID: {}, Count: {}",
        it.first.getTrackID(),
        it.first.getEventID(),
        it.first.getSourceID(),
        it.second.size());

    for (auto id: it.second){
      LOGP(info, "chipID: {}", id);
    }
  }


} // end loop on layers

// LOGP(info, "** Found {} mTracksMCLabels  :",mTracksMCLabels.size());

// LOGP(info, "** Found {} mClustersMCLCont  :",mClustersMCLCont->getIndexedSize());
// LOGP(info, "** Found {} mInputITSidxs  :",mInputITSidxs.size());

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
