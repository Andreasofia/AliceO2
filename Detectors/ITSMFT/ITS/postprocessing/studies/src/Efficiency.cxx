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

#include <TEfficiency.h>
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
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TObjArray.h>
#include <THStack.h>
#include <TString.h>
#include <numeric>

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
  void stileEfficiencyGraph(std::unique_ptr<TEfficiency>& eff, const char* name, const char* title, const int markerStyle , const double markersize, const int markercolor , const int linercolor);
  int getDCAClusterTrackMC(int countDuplicated);
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
  float mrangesPt[3][2] = {{0,0.5}, {0.5,2}, {2,7.5}};


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
  int mDuplicated_layer[3] = {0};

  //// Histos
  // Distance betweeen original and duplicated clusters
  std::unique_ptr<TH1D> mDistanceClustersX[3];
  std::unique_ptr<TH1D> mDistanceClustersY[3];
  std::unique_ptr<TH1D> mDistanceClustersZ[3];
  std::unique_ptr<TH1D> mDistanceClusters[3];
  // DCA betweeen track and original cluster
  std::unique_ptr<TH1D> mDCAxyOriginal[3];
  std::unique_ptr<TH1D> mDCAzOriginal[3];
  // DCA betweeen track and duplicated cluster
  std::unique_ptr<TH1D> mDCAxyDuplicated;
  std::unique_ptr<TH1D> mDCAzDuplicated;

  // DCA betweeen track and duplicated cluster per layer
  std::unique_ptr<TH1D> mDCAxyDuplicated_layer[3];
  std::unique_ptr<TH1D> mDCAzDuplicated_layer[3];

  // phi, eta, pt of the cluster
  std::unique_ptr<TH1D> mPhiOriginal[3];
  std::unique_ptr<TH1D> mEtaOriginal[3];
  std::unique_ptr<TH1D> mPtOriginal[3];
  std::unique_ptr<TH1D> mPtDuplicated[3];

  // position of the clusters
  std::unique_ptr<TH3D> m3DClusterPositions;
  std::unique_ptr<TH2D> m2DClusterOriginalPositions;
  std::unique_ptr<TH2D> m2DClusterDuplicatedPositions;

  // Efficiency histos
  std::unique_ptr<TH1D> mEfficiencyGoodMatch;
  std::unique_ptr<TH1D> mEfficiencyFakeMatch;
  std::unique_ptr<TH1D> mEfficiencyTotal;
  std::unique_ptr<TH1D> mEfficiencyGoodMatch_layer[3];
  std::unique_ptr<TH1D> mEfficiencyFakeMatch_layer[3];
  std::unique_ptr<TH1D> mEfficiencyTotal_layer[3];

  // phi, eta, pt of the duplicated cluster per layer
  TH2D * mPt_EtaDupl[3];

  // duplicated per layer and per cut
  std::unique_ptr<TH1D> mDuplicatedEtaAllPt[3];
  std::unique_ptr<TH1D> mDuplicatedEta[3][3];
  std::unique_ptr<TH1D> mDuplicatedPhi[3];
  TH1D* mDuplicatedPt[3];

  // matches per layer and per cut
  std::unique_ptr<TH1D> mNGoodMatchesEtaAllPt[3];
  std::unique_ptr<TH1D> mNGoodMatchesEta[3][3];
  std::unique_ptr<TH1D> mNGoodMatchesPhi[3];

  std::unique_ptr<TH1D> mNFakeMatchesEtaAllPt[3];
  std::unique_ptr<TH1D> mNFakeMatchesEta[3][3];
  std::unique_ptr<TH1D> mNFakeMatchesPhi[3];
  
  TH1D* mNGoodMatchesPt[3];
  TH1D* mNFakeMatchesPt[3];

  // calculating the efficiency with TEfficiency class
  std::unique_ptr<TEfficiency> mEffPtGood[3];
  std::unique_ptr<TEfficiency> mEffPtFake[3];

  std::unique_ptr<TEfficiency> mEffEtaGoodAllPt[3];
  std::unique_ptr<TEfficiency> mEffEtaGood[3][3];
  std::unique_ptr<TEfficiency> mEffEtaFakeAllPt[3];
  std::unique_ptr<TEfficiency> mEffEtaFake[3][3];

  std::unique_ptr<TEfficiency> mEffPhiGood[3];
  std::unique_ptr<TEfficiency> mEffPhiFake[3];
};

void EfficiencyStudy::init(InitContext& ic)
{
  LOGP(info, "--------------- init");

  o2::base::GRPGeomHelper::instance().setRequest(mGGCCDBRequest);

  auto& pars = o2::its::study::ITSEfficiencyParamConfig::Instance();
  mOutFileName = pars.outFileName;
  b = pars.b;

  int nbPt = 75;
  double xbins[nbPt + 1], ptcutl = 0.0001, ptcuth = 7.5;
  double a = std::log(ptcuth / ptcutl) / nbPt;
  for (int i = 0; i <= nbPt; i++)
    xbins[i] = ptcutl * std::exp(i * a);

  mDCAxyDuplicated = std::make_unique<TH1D>("dcaXYDuplicated", "Distance between track and duplicated cluster  ;DCA xy (cm); ", 400, -0.2, 0.2);
  mDCAzDuplicated = std::make_unique<TH1D>("dcaZDuplicated", "Distance between track and duplicated cluster  ;DCA z (cm); ", 400, -0.2, 0.2);

  m3DClusterPositions = std::make_unique<TH3D>("3DClusterPositions", ";x (cm);y (cm);z (cm)", 200, -10, 10, 200, -10, 10, 400, -20, 20);
  m2DClusterOriginalPositions = std::make_unique<TH2D>("m2DClusterOriginalPositions", ";x (cm);y (cm)", 400, -10, 10, 400, -6, 6);
  m2DClusterDuplicatedPositions = std::make_unique<TH2D>("m2DClusterDuplicatedPositions", ";x (cm);y (cm)", 400, -10, 10, 400, -6, 6);

  mEfficiencyGoodMatch = std::make_unique<TH1D>("mEfficiencyGoodMatch", ";#sigma(DCA) cut;Efficiency;", 20, 0.5, 20.5);
  mEfficiencyFakeMatch = std::make_unique<TH1D>("mEfficiencyFakeMatch", ";#sigma(DCA) cut;Efficiency;", 20, 0.5, 20.5);
  mEfficiencyTotal = std::make_unique<TH1D>("mEfficiencyTotal", ";#sigma(DCA) cut;Efficiency;", 20, 0.5, 20.5);
  
  for (int i=0; i<mNLayers; i++) {

    mDistanceClustersX[i] = std::make_unique<TH1D>(Form("distanceClustersX_L%d",i), ";Distance x (cm); ", 100, 0, 1);
    mDistanceClustersY[i] = std::make_unique<TH1D>(Form("distanceClustersY_L%d",i), ";Distance y (cm); ", 100, 0, 1);
    mDistanceClustersZ[i] = std::make_unique<TH1D>(Form("distanceClustersZ_L%d",i), ";Distance z (cm); ", 100, 0, 1);
    mDistanceClusters[i] = std::make_unique<TH1D>(Form("distanceClusters_L%d",i), ";Distance (cm); ", 100, 0, 1);

    mDCAxyOriginal[i] = std::make_unique<TH1D>(Form("dcaXYOriginal_L%d",i), "Distance between track and original cluster ;DCA xy (cm); ", 400, -0.2, 0.2);
    mDCAzOriginal[i] = std::make_unique<TH1D>(Form("dcaZOriginal_L%d",i), "Distance between track and original cluster ;DCA z (cm); ", 400, -0.2, 0.2);


    mPhiOriginal[i] = std::make_unique<TH1D>(Form("phiOriginal_L%d",i), ";phi (deg); ", 120, 0, 360);
    mEtaOriginal[i] = std::make_unique<TH1D>(Form("etaOriginal_L%d",i), ";eta (deg); ", 100, -2, 2);
    mPtOriginal[i] = std::make_unique<TH1D>(Form("ptOriginal_L%d",i), ";pt (GeV); ", 100, 0, 10);
    
    mPtDuplicated[i] = std::make_unique<TH1D>(Form("ptDuplicated_L%d",i), ";pt (GeV); ", 100, 0, 10);

    mDCAxyDuplicated_layer[i] = std::make_unique<TH1D>(Form("dcaXYDuplicated_layer_L%d",i), "Distance between track and duplicated cluster  ;DCA xy (cm); ", 400, -0.2, 0.2);
    mDCAzDuplicated_layer[i] = std::make_unique<TH1D>(Form("dcaZDuplicated_layer_L%d",i), "Distance between track and duplicated cluster  ;DCA z (cm); ", 400, -0.2, 0.2);

    mEfficiencyGoodMatch_layer[i] = std::make_unique<TH1D>(Form("mEfficiencyGoodMatch_layer_L%d",i), ";#sigma(DCA) cut;Efficiency;", 20, 0.5, 20.5);
    mEfficiencyFakeMatch_layer[i] = std::make_unique<TH1D>(Form("mEfficiencyFakeMatch_layer_L%d",i), ";#sigma(DCA) cut;Efficiency;", 20, 0.5, 20.5);
    mEfficiencyTotal_layer[i] = std::make_unique<TH1D>(Form("mEfficiencyTotal_layer_L%d",i), ";#sigma(DCA) cut;Efficiency;", 20, 0.5, 20.5);
    
    mPt_EtaDupl[i] = new TH2D(Form("mPt_EtaDupl_L%d",i), ";#it{p}_{T} (GeV/c);#eta; ", 100, 0, 10, 100, -2, 2); 

    mDuplicatedPt[i] = new TH1D(Form("mDuplicatedPt_log_L%d",i), Form("; #it{p}_{T} (GeV/c); Number of duplciated clusters L%d",i), nbPt, xbins);
    mDuplicatedPt[i]->Sumw2();
    mNGoodMatchesPt[i] = new TH1D(Form("mNGoodMatchesPt_L%d",i), Form("; #it{p}_{T} (GeV/c); Number of good matches L%d",i), nbPt, xbins);
    mNGoodMatchesPt[i]->Sumw2();
    mNFakeMatchesPt[i] = new TH1D(Form("mNFakeMatchesPt_L%d",i), Form("; #it{p}_{T} (GeV/c); Number of fake matches L%d",i), nbPt, xbins);
    mNFakeMatchesPt[i]->Sumw2();

    mDuplicatedEtaAllPt[i] = std::make_unique<TH1D>(Form("mDuplicatedEtaAllPt_L%d",i), Form("; #eta; Number of duplicated clusters L%d",i), 40, -2, 2);
    mNGoodMatchesEtaAllPt[i] = std::make_unique<TH1D>(Form("mNGoodMatchesEtaAllPt_L%d",i), Form("; #eta; Number of good matches L%d",i), 40, -2, 2);
    mNFakeMatchesEtaAllPt[i] = std::make_unique<TH1D>(Form("mNFakeMatchesEtaAllPt_L%d",i), Form("; #eta; Number of fake matches L%d",i), 40, -2, 2);

    mDuplicatedPhi[i] = std::make_unique<TH1D>(Form("mDuplicatedPhi_L%d",i), Form("; #phi (deg); Number of duplicated clusters L%d",i), 120, 0, 360);
    mNGoodMatchesPhi[i] = std::make_unique<TH1D>(Form("mNGoodMatchesPhi_L%d",i), Form("; #phi (deg); Number of good matches L%d",i), 120, 0, 360);
    mNFakeMatchesPhi[i] = std::make_unique<TH1D>(Form("mNFakeMatchesPhi_L%d",i), Form("; #phi (deg); Number of fake matches L%d",i), 120, 0, 360);
  

    for (int j=0; j<3; j++){
      mDuplicatedEta[i][j] = std::make_unique<TH1D>(Form("mDuplicatedEta_L%d_pt%d",i,j), Form("%f < #it{p}_{T} < %f GeV/c; #eta; Number of duplicated clusters L%d", mrangesPt[j][0], mrangesPt[j][1],i), 40, -2, 2);
      mNGoodMatchesEta[i][j] = std::make_unique<TH1D>(Form("mNGoodMatchesEta_L%d_pt%d",i,j), Form("%f < #it{p}_{T} < %f GeV/c; #eta; Number of good matches L%d",i, mrangesPt[j][0], mrangesPt[j][1],i), 40, -2, 2);
      mNFakeMatchesEta[i][j] = std::make_unique<TH1D>(Form("mNFakeMatchesEta_L%d_pt%d",i,j), Form("%f < #it{p}_{T} < %f GeV/c; #eta; Number of fake matches L%d",i, mrangesPt[j][0], mrangesPt[j][1],i), 40, -2, 2);
    }
  }
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

void EfficiencyStudy::stileEfficiencyGraph(std::unique_ptr<TEfficiency>& eff, const char* name, const char* title, const int markerStyle = kFullCircle, const double markersize = 1, const int markercolor = kBlack, const int linecolor = kBlack)
{
  eff->SetName(name);
  eff->SetTitle(title);
  eff->SetMarkerStyle(markerStyle);
  eff->SetMarkerSize(markersize);
  eff->SetMarkerColor(markercolor);
  eff->SetLineColor(linecolor);
}
int EfficiencyStudy::getDCAClusterTrackMC(int countDuplicated = 0)
{
  // get the DCA between the clusters and the track from MC and fill histograms: distance between original and duplicated cluster, DCA, phi, clusters
  LOGP(info, "--------------- getDCAClusterTrackMC");

  mOutFile->mkdir("DistanceClusters/");
  mOutFile->mkdir("DCA/");
  mOutFile->mkdir("Pt_Eta_Phi/");
  
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

      auto pt = trackParCov.getPt();
      auto eta = trackParCov.getEta();
      auto phi = trackParCov.getPhi()*180/M_PI;

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

        mPhiOriginal[layer]->Fill(phi);
        mPtOriginal[layer]->Fill(pt);
        mEtaOriginal[layer]->Fill(eta);
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
                  
                  if (countDuplicated == 0){
                    mDuplicated_layer[layerDuplicated]++; // This has to be incremented at the first call
                  }

                  if (countDuplicated == 1){            
                    mDuplicatedPt[layerDuplicated]->Fill(pt);
                    mDuplicatedEtaAllPt[layerDuplicated]->Fill(eta);
                    for (int ipt = 0; ipt < 3; ipt++){
                      if (pt>=mrangesPt[ipt][0] && pt<mrangesPt[ipt][1]){
                        mDuplicatedEta[layerDuplicated][ipt]->Fill(eta);
                      }
                    }
                    
                    mDuplicatedPhi[layerDuplicated]->Fill(phi);
                    mPtDuplicated[layerClus]->Fill(pt);
                    mPt_EtaDupl[layerClus]->Fill(pt,eta);
                  }

                  m3DClusterPositions->Fill(clusDuplicatedPointGlob.x(), clusDuplicatedPointGlob.y(), clusDuplicatedPointGlob.z());
                  m2DClusterDuplicatedPositions->Fill(clusDuplicatedPointGlob.x(), clusDuplicatedPointGlob.y());

                  /// compute the distance between original and dubplicated cluster
                  mDistanceClustersX[layerClus]->Fill(abs(clusOriginalPointGlob.x() - clusDuplicatedPointGlob.x()));
                  mDistanceClustersY[layerClus]->Fill(abs(clusOriginalPointGlob.y() - clusDuplicatedPointGlob.y()));
                  mDistanceClustersZ[layerClus]->Fill(abs(clusOriginalPointGlob.z() - clusDuplicatedPointGlob.z()));
                  mDistanceClusters[layerClus]->Fill(std::hypot(clusOriginalPointGlob.x() - clusDuplicatedPointGlob.x(), clusOriginalPointGlob.y() - clusDuplicatedPointGlob.y(), clusOriginalPointGlob.z() - clusDuplicatedPointGlob.z()));

                  /// Compute the DCA between the cluster location and the track

                  /// first propagate to the original cluster
                  trackParCov.rotate(mGeometry->getSensorRefAlpha(clusOriginal.getSensorID()));
                  if (propagator->propagateToDCA(clusOriginalPointGlob, trackParCov, b, 2.f, matCorr, &clusOriginalDCA)) {
                    mDCAxyOriginal[layerClus]->Fill(clusOriginalDCA[0]);
                    mDCAzOriginal[layerClus]->Fill(clusOriginalDCA[1]);
                  }
                  /// then propagate to the duplicated cluster
                  trackParCov.rotate(mGeometry->getSensorRefAlpha(clusDuplicated.getSensorID()));
                  if (propagator->propagateToDCA(clusDuplicatedPointGlob, trackParCov, b, 2.f, matCorr, &clusDuplicatedDCA)) {
                    mDCAxyDuplicated->Fill(clusDuplicatedDCA[0]);
                    mDCAzDuplicated->Fill(clusDuplicatedDCA[1]);
                    mDCAxyDuplicated_layer[layerDuplicated]->Fill(clusDuplicatedDCA[0]);
                    mDCAzDuplicated_layer[layerDuplicated]->Fill(clusDuplicatedDCA[1]);             
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
  if (countDuplicated == 1){
    m3DClusterPositions->Write();
    m2DClusterOriginalPositions->Write();
    m2DClusterDuplicatedPositions->Write();

    mOutFile->cd("DistanceClusters");
    for (int i=0; i<mNLayers; i++) {
      mDistanceClustersX[i]->Write();
      mDistanceClustersY[i]->Write();
      mDistanceClustersZ[i]->Write();
      mDistanceClusters[i]->Write();
    }

    mOutFile->cd("DCA");
    mDCAxyDuplicated->Write();
    mDCAzDuplicated->Write();
    for (int i=0; i<mNLayers; i++) {
      mDCAxyDuplicated_layer[i]->Write();
      mDCAzDuplicated_layer[i]->Write();

      mDCAxyOriginal[i]->Write();
      mDCAzOriginal[i]->Write();
    }

    mOutFile->cd("Pt_Eta_Phi/");
    for (int i=0; i<mNLayers; i++) {
      mPhiOriginal[i]->Write();
      mDuplicatedPhi[i]->Write();
      mPtOriginal[i]->Write();
      mPtDuplicated[i]->Write();
      mDuplicatedPt[i]->Write();
      mEtaOriginal[i]->Write();
      mDuplicatedEtaAllPt[i]->Write(); 
      for (int p=0; p<3; p++){
        mDuplicatedEta[i][p]->Write();
      }
      mPt_EtaDupl[i]->Write();
    }
  }


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

  int duplicated = getDCAClusterTrackMC(0);

  double meanDCAxyDuplicated[mNLayers] = {0};
  double meanDCAzDuplicated[mNLayers] = {0};
  double sigmaDCAxyDuplicated[mNLayers] = {0};
  double sigmaDCAzDuplicated[mNLayers] = {0};

  for (int l=0; l<mNLayers; l++) {
    meanDCAxyDuplicated[l] = mDCAxyDuplicated_layer[l]->GetMean();
    meanDCAzDuplicated[l] = mDCAzDuplicated_layer[l]->GetMean();
    sigmaDCAxyDuplicated[l] = mDCAxyDuplicated_layer[l]->GetRMS();
    sigmaDCAzDuplicated[l] = mDCAzDuplicated_layer[l]->GetRMS();
    
  }

  for (int l=0; l<mNLayers; l++) {
    LOGP(info, "meanDCAxyDuplicated L{}: {}, meanDCAzDuplicated: {}, sigmaDCAxyDuplicated: {}, sigmaDCAzDuplicated: {}",l, meanDCAxyDuplicated[l], meanDCAzDuplicated[l], sigmaDCAxyDuplicated[l], sigmaDCAzDuplicated[l]);
  }
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

  unsigned int nGoodMatches_layer[mNLayers][20] = {0};
  unsigned int nFakeMatches_layer[mNLayers][20] = {0};

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
        auto labsOriginal = mClustersMCLCont->getLabels(mInputITSidxs[iclTrack]); // get labels of the cluster associated to the track (original)
        auto staveOriginal = mGeometry->getStave(clusOriginal.getSensorID());
        auto chipOriginal = mGeometry->getChipIdInStave(clusOriginal.getSensorID());



        for (auto& labT : labsOriginal) { // for each valid label iterate over ALL the clusters in the ROF to see if there are duplicates
          if (labT != tracklab)
            continue;
          if (!labT.isValid()) 
            continue;

          /// for each oroginal cluster iterate over all the possible "adjacent" clusters (stave +-1, chip =,+-1) and calculate the DCA with the track. Then compare the cluster label with the track label to see if it is a true or fake match
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
                if (abs(meanDCAxyDuplicated[layerDuplicated] - clusDuplicatedDCA[0]) < (i + 1) * sigmaDCAxyDuplicated[layerDuplicated] && abs(meanDCAzDuplicated[layerDuplicated] - clusDuplicatedDCA[1]) < (i + 1) * sigmaDCAzDuplicated[layerDuplicated]) { // check if the DCA is within the cut i*sigma
                  if (mVerboseOutput)
                    LOGP(info, "Check DCA ok: {} < {}; {} < {}", abs(meanDCAxyDuplicated[layerDuplicated] - clusDuplicatedDCA[0]), (i + 1) * sigmaDCAxyDuplicated[layerDuplicated], abs(meanDCAzDuplicated[layerDuplicated] - clusDuplicatedDCA[1]), (i + 1) * sigmaDCAzDuplicated[layerDuplicated]);
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
                  if (isGoodMatch){
                    nGoodMatches[i]++;
                    nGoodMatches_layer[layerDuplicated][i]++;
                  }
                  else {
                    nFakeMatches[i]++;
                    nFakeMatches_layer[layerDuplicated][i]++;
                  }
                } else if (mVerboseOutput)
                  LOGP(info, "Check DCA failed");
              }
            } else if (mVerboseOutput)
              LOGP(info, "Propagation failed");
          } // end loop on all the clusters in the rof
        }
      }   // end loop on clusters associated to the track
    }     // end loop on tracks per ROF
  }       // end loop on ROFRecords array

  for (int i = 0; i < 20; i++) {
    LOGP(info, "Cut: {} sigma -> number of duplicated clusters: {} nDCAMatches: {} nGoodMatches: {} nFakeMatches: {}", i + 1, duplicated, nDCAMatches[i], nGoodMatches[i], nFakeMatches[i]);
    mEfficiencyGoodMatch->SetBinContent(i+1, nGoodMatches[i]);
    mEfficiencyFakeMatch->SetBinContent(i+1, nFakeMatches[i]);
    mEfficiencyTotal->SetBinContent(i+1, double(nGoodMatches[i] + nFakeMatches[i]));

    for (int l=0; l<mNLayers; l++){
      mEfficiencyGoodMatch_layer[l]->SetBinContent(i+1, nGoodMatches_layer[l][i]);
      mEfficiencyFakeMatch_layer[l]->SetBinContent(i+1, nFakeMatches_layer[l][i]);
      mEfficiencyTotal_layer[l]->SetBinContent(i+1, double(nGoodMatches_layer[l][i] + nFakeMatches_layer[l][i]));
    }
  }

  for (int l=0; l<mNLayers; l++){
    mEfficiencyGoodMatch_layer[l]->Scale(1./double(mDuplicated_layer[l]), "b");
    mEfficiencyFakeMatch_layer[l]->Scale(1./double(mDuplicated_layer[l]), "b");
    mEfficiencyTotal_layer[l]->Scale(1./double(mDuplicated_layer[l]), "b");
  }

  mEfficiencyGoodMatch->Scale(1./double(duplicated), "b");
  mEfficiencyFakeMatch->Scale(1./double(duplicated), "b");
  mEfficiencyTotal->Scale(1./double(duplicated), "b");

  mOutFile->mkdir("Efficiency/");
  mOutFile->cd("Efficiency/");
  mEfficiencyGoodMatch->Write();
  mEfficiencyFakeMatch->Write();
  mEfficiencyTotal->Write();
  for (int l=0; l<mNLayers; l++){
    mEfficiencyGoodMatch_layer[l]->Write();
    mEfficiencyFakeMatch_layer[l]->Write();
    mEfficiencyTotal_layer[l]->Write();
    
    mEfficiencyGoodMatch_layer[l]->GetYaxis()->SetRangeUser(-0.1, 1.1);
    mEfficiencyFakeMatch_layer[l]->GetYaxis()->SetRangeUser(-0.1, 1.1);
    mEfficiencyTotal_layer[l]->GetYaxis()->SetRangeUser(-0.1, 1.1);
  }

  mEfficiencyGoodMatch->GetYaxis()->SetRangeUser(-0.1, 1.1);
  mEfficiencyFakeMatch->GetYaxis()->SetRangeUser(-0.1, 1.1);
  mEfficiencyTotal->GetYaxis()->SetRangeUser(-0.1, 1.1);

  TCanvas c;
  c.SetName("EffVsDCA_allLayers");
  auto leg = std::make_unique<TLegend>(0.75, 0.45, 0.89, 0.65);
  leg->AddEntry(mEfficiencyGoodMatch.get(), "#frac{# good matches}{# tot duplicated clusters}", "p");
  leg->AddEntry(mEfficiencyFakeMatch.get(), "#frac{# fake matches}{# tot duplicated clusters}", "p");
  leg->AddEntry(mEfficiencyTotal.get(), "#frac{# tot matches}{# tot duplicated clusters}", "p");

  mEfficiencyGoodMatch->Draw("P l E1_NOSTAT PLC PMC ");
  mEfficiencyFakeMatch->Draw("same P l E1_NOSTAT  PLC PMC");
  mEfficiencyTotal->Draw("same P l E1_NOSTAT  PLC PMC");
  leg->Draw("same");
  c.Write();
  c.SaveAs("prova.png");

  TCanvas cc[mNLayers];
  for (int l=0; l<mNLayers; l++){
    cc[l].cd();
    cc[l].SetName(Form("EffVsDCA_layer_L%d", l));

    auto leg = std::make_unique<TLegend>(0.75, 0.45, 0.89, 0.65);
    leg->AddEntry(mEfficiencyGoodMatch_layer[l].get(), "#frac{# good matches}{# tot duplicated clusters}", "p");
    leg->AddEntry(mEfficiencyFakeMatch_layer[l].get(), "#frac{# fake matches}{# tot duplicated clusters}", "p");
    leg->AddEntry(mEfficiencyTotal_layer[l].get(), "#frac{# tot matches}{# tot duplicated clusters}", "p");
    
    mEfficiencyGoodMatch_layer[l]->SetLineColor(kBlue+3);
    mEfficiencyGoodMatch_layer[l]->SetMarkerColor(kBlue+3);
    mEfficiencyGoodMatch_layer[l]->Draw("P l E1_NOSTAT");
    mEfficiencyFakeMatch_layer[l]->SetLineColor(kAzure+7);
    mEfficiencyFakeMatch_layer[l]->SetMarkerColor(kAzure+7);
    mEfficiencyFakeMatch_layer[l]->Draw("same P l E1_NOSTAT");
    mEfficiencyTotal_layer[l]->SetLineColor(kGreen+1);
    mEfficiencyTotal_layer[l]->SetMarkerColor(kGreen+1);
    mEfficiencyTotal_layer[l]->Draw("same P l E1_NOSTAT");
    leg->Draw("same");
    cc[l].Write();
    cc[l].SaveAs(Form("provaLayer%d.png", l));
  }


}

void EfficiencyStudy::studyClusterSelectionMC()
{
  // study to find a good selection method for the duplicated cluster, to be used for non MC data 
  // iterate over tracks an associated clusters, and find the closer cluster that is not the original one applying cuts on staveID and chipID
  // fix the DCA < 10 sigma, then compute the efficiency for each bin of pt, eta and phi

  LOGP(info, "--------------- studyClusterSelection");

  int duplicated = getDCAClusterTrackMC(1);

  std::cout<<"duplicated: "<<duplicated<<std::endl;

  double meanDCAxyDuplicated[mNLayers] = {0};
  double meanDCAzDuplicated[mNLayers] = {0};
  double sigmaDCAxyDuplicated[mNLayers] = {0};
  double sigmaDCAzDuplicated[mNLayers] = {0};
  
  for (int l=0; l<mNLayers; l++) {
    meanDCAxyDuplicated[l] = mDCAxyDuplicated_layer[l]->GetMean();
    meanDCAzDuplicated[l] = mDCAzDuplicated_layer[l]->GetMean();
    sigmaDCAxyDuplicated[l] = mDCAxyDuplicated_layer[l]->GetRMS();
    sigmaDCAzDuplicated[l] = mDCAzDuplicated_layer[l]->GetRMS();
    
  }

  for (int l=0; l<mNLayers; l++) {
    LOGP(info, "meanDCAxyDuplicated L{}: {}, meanDCAzDuplicated: {}, sigmaDCAxyDuplicated: {}, sigmaDCAzDuplicated: {}",l, meanDCAxyDuplicated[l], meanDCAzDuplicated[l], sigmaDCAxyDuplicated[l], sigmaDCAzDuplicated[l]);
  }

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

  std::map<std::tuple<int, double, o2::MCCompLabel>, bool> clusterMatchesPtEta[100][100] = {};


  for (unsigned int iROF = 0; iROF < mTracksROFRecords.size(); iROF++) { // loop on ROFRecords array
    rofIndexTrack = mTracksROFRecords[iROF].getFirstEntry();
    rofNEntriesTrack = mTracksROFRecords[iROF].getNEntries();

    rofIndexClus = mClustersROFRecords[iROF].getFirstEntry();
    rofNEntriesClus = mClustersROFRecords[iROF].getNEntries();

    //////calculcating efficiency vs pt, eta, phi 
    for (unsigned int iTrack = rofIndexTrack; iTrack < rofIndexTrack + rofNEntriesTrack; iTrack++) { // loop on tracks per ROF
      auto track = mTracks[iTrack];
      o2::track::TrackParCov trackParCov = mTracks[iTrack];
      
      int firstClus = track.getFirstClusterEntry(); // get the first cluster of the track
      int ncl = track.getNumberOfClusters();        // get the number of clusters of the track

      auto& tracklab = mTracksMCLabels[iTrack];
      if (tracklab.isFake())
        continue;
      
      auto pt = trackParCov.getPt();
      auto eta = trackParCov.getEta();
      auto phi = trackParCov.getPhi()*180/M_PI;


      if (mVerboseOutput) {
        LOGP(info, "--------- track Label: ");
        tracklab.print();
      }
      for (int iclTrack = firstClus; iclTrack < firstClus + ncl; iclTrack++) { // loop on clusters associated to the track to extract layer, stave and chip to restrict the possible matches to be searched with the DCA cut
        // LOGP(info, "New cluster");
        auto& clusOriginal = mClusters[mInputITSidxs[iclTrack]];
        auto layerOriginal = mGeometry->getLayer(clusOriginal.getSensorID());
        if (layerOriginal >= mNLayers)
          continue;

        auto labsOriginal = mClustersMCLCont->getLabels(mInputITSidxs[iclTrack]); // get labels of the cluster associated to the track (original)
        auto staveOriginal = mGeometry->getStave(clusOriginal.getSensorID());
        auto chipOriginal = mGeometry->getChipIdInStave(clusOriginal.getSensorID());

        std::tuple<int, double, gsl::span<const o2::MCCompLabel>> clusID_rDCA_label = {0, 999., gsl::span<const o2::MCCompLabel>()}; // inizializing tuple with dummy values

        bool adjacentFound = 0;
          /// for each oroginal cluster iterate over all the possible "adjacten" clusters (stave +-1, chip =,+-1) and calculate the DCA with the track. Then choose the closest one.
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
          if (!propagator->propagateToDCA(clusDuplicatedPointGlob, trackParCov, b, 2.f, matCorr, &clusDuplicatedDCA)) { // check if the propagation fails
            continue;
          }

          /// Imposing that the distance between the original cluster and the duplicated one is less than x sigma
          if (!(abs(meanDCAxyDuplicated[layerDuplicated] - clusDuplicatedDCA[0]) < 5 * sigmaDCAxyDuplicated[layerDuplicated] && abs(meanDCAzDuplicated[layerDuplicated] - clusDuplicatedDCA[1]) < 5 * sigmaDCAzDuplicated[layerDuplicated])) {
            continue;
          }

          if (mVerboseOutput)
            LOGP(info, "Propagation ok");
          double rDCA = std::hypot(clusDuplicatedDCA[0], clusDuplicatedDCA[1]);

          // taking the closest cluster within x sigma
          if (rDCA < std::get<1>(clusID_rDCA_label)) {   // updating the closest cluster 
            clusID_rDCA_label = {iClus, rDCA, labsDuplicated};
          }
          adjacentFound = 1;
        } // end loop on all the clusters in the rof 

        // here clusID_rDCA_label is updated with the closest cluster to the track other than the original one
        // checking if it is a good or fake match looking at the labels

        if (!adjacentFound)
          continue;

        bool isGood = false;
        for (auto lab : std::get<2>(clusID_rDCA_label)) {
          if (lab == tracklab) {
            isGood = true;
            mNGoodMatchesPt[layerOriginal]->Fill(pt);
            mNGoodMatchesEtaAllPt[layerOriginal]->Fill(eta);
            for (int ipt = 0; ipt < 3; ipt++){
              if (pt>=mrangesPt[ipt][0] && pt<mrangesPt[ipt][1]){
                mNGoodMatchesEta[layerOriginal][ipt]->Fill(eta);
              }
            }
            mNGoodMatchesPhi[layerOriginal]->Fill(phi);
            continue;
          }
        }
        if (!isGood) {

          mNFakeMatchesPt[layerOriginal]->Fill(pt);
          mNFakeMatchesEtaAllPt[layerOriginal]->Fill(eta);
          for (int ipt = 0; ipt < 3; ipt++){
            if (pt>=mrangesPt[ipt][0] && pt<mrangesPt[ipt][1]){
              mNFakeMatchesEta[layerOriginal][ipt]->Fill(eta);
            }
          }
          mNFakeMatchesPhi[layerOriginal]->Fill(phi);
        }

      }   // end loop on clusters associated to the track
    }     // end loop on tracks per ROF
  }       // end loop on ROFRecords array


  mOutFile->mkdir("EfficiencyCuts/");
  mOutFile->cd("EfficiencyCuts/");

  TH1D* axpt = new TH1D("axpt","",1,0,7.5);
  TH1D* axetaAllPt = new TH1D("axetaAllPt","",1,-2,2);
  TH1D* axeta[3];
  for (int ipt=0; ipt<3; ipt++){
    axeta[ipt] = new TH1D(Form("axeta%d",ipt),Form("axeta%d",ipt),1,-2,2);
  }
  TH1D* axphi = new TH1D("axphi","",1,0,360);

  TCanvas *effPt[3];
  TCanvas *effEtaAllPt[3];
  TCanvas *effEta[3][3];
  TCanvas *effPhi[3];
  for (int l=0; l<3; l++){
    if (mVerboseOutput) std::cout<<"Pt L"<<l<<"\n\n";
  
    effPt[l]= new TCanvas(Form("effPt_L%d",l));

    mEffPtGood[l] = std::make_unique<TEfficiency>(*mNGoodMatchesPt[l], *mDuplicatedPt[l]);
    stileEfficiencyGraph(mEffPtGood[l], Form("mEffPtGood_L%d",l), Form("L%d;#it{p}_{T} (GeV/#it{c});Efficiency",l ), kFullDiamond, 1,kGreen+3, kGreen+3 );

    for(int ibin=1; ibin<=mNFakeMatchesPt[l]->GetNbinsX(); ibin++){
      if (mNFakeMatchesPt[l]->GetBinContent(ibin) > mDuplicatedPt[l]->GetBinContent(ibin)){
        std::cout<<"--- Pt: Npass = "<<mNFakeMatchesPt[l]->GetBinContent(ibin)<<",  Nall = "<<mDuplicatedPt[l]->GetBinContent(ibin)<<" for ibin = "<<ibin<<std::endl;
        mNFakeMatchesPt[l]->SetBinContent(ibin, mDuplicatedPt[l]->GetBinContent(ibin));
      }
    }
    mEffPtFake[l] = std::make_unique<TEfficiency>(*mNFakeMatchesPt[l], *mDuplicatedPt[l]);
    stileEfficiencyGraph(mEffPtFake[l], Form("mEffPtFake_L%d",l), Form("L%d;#it{p}_{T} (GeV/#it{c});Efficiency",l ), kFullDiamond, 1, kRed+1, kRed+1);

    axpt->SetTitle(Form("L%d;#it{p}_{T} (GeV/#it{c});Efficiency",l));
    axpt->GetYaxis()->SetRangeUser(-0.1,1.1);
    axpt->GetXaxis()->SetRangeUser(0,7.5);
    axpt->Draw();
    mEffPtGood[l]->Draw("same p");
    mEffPtFake[l]->Draw("same p");
  
    auto legpt = std::make_unique<TLegend>(0.70, 0.15, 0.89, 0.35);
    legpt->AddEntry(mEffPtGood[l].get(), "#frac{# good matches}{# tot duplicated clusters}", "pl");
    legpt->AddEntry(mEffPtFake[l].get(), "#frac{# fake matches}{# tot duplicated clusters}", "pl");
    legpt->Draw("same");
    effPt[l]->Write();

    if (mVerboseOutput) std::cout<<"Eta L"<<l<<"\n\n";
  
    effEtaAllPt[l]= new TCanvas(Form("effEtaAllPt_L%d",l));

    mEffEtaGoodAllPt[l] = std::make_unique<TEfficiency>(*mNGoodMatchesEtaAllPt[l], *mDuplicatedEtaAllPt[l]);
    stileEfficiencyGraph(mEffEtaGoodAllPt[l], Form("mEffEtaGoodAllPt_L%d",l), Form("L%d;#eta;Efficiency",l ), kFullDiamond, 1, kGreen+3, kGreen+3);
    
    for(int ibin=1; ibin<=mNFakeMatchesEtaAllPt[l]->GetNbinsX(); ibin++){
      if (mNFakeMatchesEtaAllPt[l]->GetBinContent(ibin) > mDuplicatedEtaAllPt[l]->GetBinContent(ibin)){
        std::cout<<"--- EtaAllPt: Npass = "<<mNFakeMatchesEtaAllPt[l]->GetBinContent(ibin)<<",  Nall = "<<mDuplicatedEtaAllPt[l]->GetBinContent(ibin)<<" for ibin = "<<ibin<<std::endl;
        mNFakeMatchesEtaAllPt[l]->SetBinContent(ibin, mDuplicatedEtaAllPt[l]->GetBinContent(ibin));
      }
    }
    mEffEtaFakeAllPt[l] = std::make_unique<TEfficiency>(*mNFakeMatchesEtaAllPt[l], *mDuplicatedEtaAllPt[l]);
    stileEfficiencyGraph(mEffEtaFakeAllPt[l], Form("mEffEtaFakeAllPt_L%d",l), Form("L%d;#eta;Efficiency",l ), kFullDiamond, 1, kRed+1, kRed+1);

    axetaAllPt->SetTitle(Form("L%d;#eta;Efficiency",l));
    axetaAllPt->GetYaxis()->SetRangeUser(-0.1,1.1);

    axetaAllPt->Draw();
    mEffEtaGoodAllPt[l]->Draw("same p");
    mEffEtaFakeAllPt[l]->Draw("same p");

    auto legEta = std::make_unique<TLegend>(0.70, 0.15, 0.89, 0.35);
    legEta->AddEntry(mEffEtaGoodAllPt[l].get(), "#frac{# good matches}{# tot duplicated clusters}", "pl");
    legEta->AddEntry(mEffEtaFakeAllPt[l].get(), "#frac{# fake matches}{# tot duplicated clusters}", "pl");
    legEta->Draw("same");
    effEtaAllPt[l]->Write();

    /// eta in different pt ranges 
    for (int ipt=0; ipt<3; ipt++){
      effEta[l][ipt]= new TCanvas(Form("effEta_L%d_pt%d",l,ipt));

      mEffEtaGood[l][ipt] = std::make_unique<TEfficiency>(*mNGoodMatchesEta[l][ipt], *mDuplicatedEta[l][ipt]);
      stileEfficiencyGraph(mEffEtaGood[l][ipt], Form("mEffEtaGood_L%d_pt%d",l,ipt), Form("L%d     %.1f #leq #it{p}_{T} < %.1f GeV/#it{c};#eta;Efficiency",l, mrangesPt[ipt][0], mrangesPt[ipt][1] ), kFullDiamond, 1, kGreen+3, kGreen+3);
      
      for(int ibin=1; ibin<=mNFakeMatchesEta[l][ipt]->GetNbinsX(); ibin++){
        if (mNFakeMatchesEta[l][ipt]->GetBinContent(ibin) > mDuplicatedEta[l][ipt]->GetBinContent(ibin)){
          std::cout<<"--- Eta : Npass = "<<mNFakeMatchesEta[l][ipt]->GetBinContent(ibin)<<",  Nall = "<<mDuplicatedEta[l][ipt]->GetBinContent(ibin)<<" for ibin = "<<ibin<<std::endl;
          mNFakeMatchesEta[l][ipt]->SetBinContent(ibin, mDuplicatedEta[l][ipt]->GetBinContent(ibin));
        }
      }

      mEffEtaFake[l][ipt] = std::make_unique<TEfficiency>(*mNFakeMatchesEta[l][ipt], *mDuplicatedEta[l][ipt]);
      stileEfficiencyGraph(mEffEtaFake[l][ipt], Form("mEffEtaFake_L%d_pt%d",l,ipt), Form("L%d    %.1f #leq #it{p}_{T} < %.1f GeV/#it{c};#eta;Efficiency",l, mrangesPt[ipt][0], mrangesPt[ipt][1] ), kFullDiamond, 1, kRed+1, kRed+1);

      axeta[ipt]->SetTitle(Form("L%d     %.1f #leq #it{p}_{T} < %.1f GeV/#it{c};#eta;Efficiency",l, mrangesPt[ipt][0], mrangesPt[ipt][1] ));
      axeta[ipt]->GetYaxis()->SetRangeUser(-0.1,1.1);

      axeta[ipt]->Draw();
      mEffEtaGood[l][ipt]->Draw("same p");
      mEffEtaFake[l][ipt]->Draw("same p");

      auto legEta = std::make_unique<TLegend>(0.70, 0.15, 0.89, 0.35);
      legEta->AddEntry(mEffEtaGood[l][ipt].get(), "#frac{# good matches}{# tot duplicated clusters}", "pl");
      legEta->AddEntry(mEffEtaFake[l][ipt].get(), "#frac{# fake matches}{# tot duplicated clusters}", "pl");
      legEta->Draw("same");
      effEta[l][ipt]->Write();      
    }
    
    if (mVerboseOutput) std::cout<<"Phi L"<<l<<"\n\n";
  
    effPhi[l]= new TCanvas(Form("effPhi_L%d",l));

    mEffPhiGood[l] = std::make_unique<TEfficiency>(*mNGoodMatchesPhi[l], *mDuplicatedPhi[l]);
    stileEfficiencyGraph(mEffPhiGood[l], Form("mEffPhiGood_L%d",l), Form("L%d;#phi;Efficiency",l ), kFullDiamond, 1, kGreen+3, kGreen+3);

    for(int ibin=1; ibin<=mNFakeMatchesPhi[l]->GetNbinsX(); ibin++){
      if (mNFakeMatchesPhi[l]->GetBinContent(ibin) > mDuplicatedPhi[l]->GetBinContent(ibin)){
        std::cout<<"--- phi Npass = "<<mNFakeMatchesPhi[l]->GetBinContent(ibin)<<",  Nall = "<<mDuplicatedPhi[l]->GetBinContent(ibin)<<" for ibin = "<<ibin<<std::endl;
        mNFakeMatchesPhi[l]->SetBinContent(ibin, mDuplicatedPhi[l]->GetBinContent(ibin));
      }
    }
    mEffPhiFake[l] = std::make_unique<TEfficiency>(*mNFakeMatchesPhi[l], *mDuplicatedPhi[l]);
    stileEfficiencyGraph(mEffPhiFake[l], Form("mEffPhiFake_L%d",l), Form("L%d;#phi;Efficiency",l ), kFullDiamond, 1, kRed+1, kRed+1);

    axphi->SetTitle(Form("L%d;#phi;Efficiency",l));
    axphi->GetYaxis()->SetRangeUser(-0.1,1.1);
    axphi->Draw();
    mEffPhiGood[l]->Draw("same p");
    mEffPhiFake[l]->Draw("same p");

    auto legPhi = std::make_unique<TLegend>(0.70, 0.15, 0.89, 0.35);
    legPhi->AddEntry(mEffPhiGood[l].get(), "#frac{# good matches}{# tot duplicated clusters}", "pl");
    legPhi->AddEntry(mEffPhiFake[l].get(), "#frac{# fake matches}{# tot duplicated clusters}", "pl");
    legPhi->Draw("same");
    effPhi[l]->Write();

  }
}

void EfficiencyStudy::process(o2::globaltracking::RecoContainer& recoData)
{
  LOGP(info, "--------------- process");

  mOutFile = std::make_unique<TFile>(mOutFileName.c_str(), "recreate");

  if (mUseMC) {
    studyDCAcutsMC();
    studyClusterSelectionMC();
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