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
}

void EfficiencyStudy::process(o2::globaltracking::RecoContainer& recoData)
{

  ////starting from tracks
  LOGP(info, "--------------- process");

  LOGP(info, "** Found in {} rofs:\n\t- {} clusters\n\t",
       mClustersROFRecords.size(), mClusters.size());

  LOGP(info, "mClusters size: {}, mClustersROFRecords size: {}, mClustersMCLCont size: {}", mClusters.size(), mClustersROFRecords.size(), mClustersMCLCont->getNElements());
  LOGP(info, "mTracks size: {}, mTracksROFRecords size: {}, mTracksMCLabels size: {}", mTracks.size(), mTracksROFRecords.size(), mTracksMCLabels.size());

  short unsigned int nLayers = 3;
  unsigned int rofIndexTrack = 0;
  unsigned int rofNEntriesTrack = 0;  
  unsigned int rofIndexClus = 0;
  unsigned int rofNEntriesClus = 0;
  unsigned int totalRofEntries = 0;
  int nLabels = 0;
  unsigned int totPrecClus = 0;


  std::unordered_map<o2::MCCompLabel, std::vector<int>> label_vecClus[mClustersROFRecords.size()][nLayers];  // array of maps nRofs x Nlayers -> {label, vec(iClus)} where vec(iClus) are the clusters that share the same label

  for (unsigned int iROF = 0; iROF < mTracksROFRecords.size(); iROF++) {  // loop on ROFRecords array
    rofIndexTrack = mTracksROFRecords[iROF].getFirstEntry();
    rofNEntriesTrack = mTracksROFRecords[iROF].getNEntries();

    rofIndexClus = mClustersROFRecords[iROF].getFirstEntry();
    rofNEntriesClus = mClustersROFRecords[iROF].getNEntries();


    for (unsigned int iTrack=rofIndexTrack; iTrack < rofIndexTrack+rofNEntriesTrack; iTrack++){  // loop on tracks per ROF
      auto track = mTracks[iTrack];
      int firstClus = track.getFirstClusterEntry(); // get the first cluster of the track
      int ncl = track.getNumberOfClusters();        // get the number of clusters of the track

      auto& tracklab = mTracksMCLabels[iTrack]; 
      if (tracklab.isFake())
        continue;
      
      for (int iclTrack =firstClus; iclTrack < firstClus+ ncl; iclTrack++) { // loop on clusters associated to the track
        auto& clus = mClusters[mInputITSidxs[iclTrack]];
        auto layer = mGeometry->getLayer(clus.getSensorID());
        if (layer>=nLayers) continue; // checking only selected layers
        auto labsTrack = mClustersMCLCont->getLabels(mInputITSidxs[iclTrack]); // get labels of clusters associated to the track
      
        for (auto& labT : labsTrack) {    // for each valid label iterate over ALL the clusters in the ROF to see if there are duplicates
          nLabels++;
          if (labT.isValid()) {
            for (unsigned int iClus=rofIndexClus; iClus < rofIndexClus + rofNEntriesClus; iClus++){  // iteration over ALL the clusters in thre ROF
              auto cluster = mClusters[iClus];
              auto layerClus = mGeometry->getLayer(cluster.getSensorID());
              if (layerClus!=layer) continue;  // we are interested only in duplicated clusters on the same layer
              auto labsClus = mClustersMCLCont->getLabels(iClus); // ideally I can have more than one label per cluster
              for (auto labC: labsClus){
                if (labT == labC) {
                  label_vecClus[iROF][layerClus][labT].push_back(iClus);   //same cluster: label from the track = label from the cluster
                }
              }
            }
          }
        }
      } // end loop on clusters
      totPrecClus+=ncl;
    } // end loop on tracks per ROF
  } // end loop on ROFRecords array
  
  LOGP(info, "Total number of clusters: {} ", totPrecClus);
  LOGP(info, "total rof entries: {}", totalRofEntries);
  LOGP(info, "total nLabels: {}", nLabels);

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
          std::cout<< "ROF: "<<iROF<<", iClus: "<<iClus<<" -> chip: "<< chipid <<" = "<<name<<std::endl;
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