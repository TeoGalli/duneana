////////////////////////////////////////////////////////////////////////
// Class:       MichelElectron
// Plugin Type: analyzer (Unknown Unknown)
// File:        MichelElectron_module.cc
//
// Generated at Tue Apr  2 04:47:23 2024 by Matteo Galli using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// Additional framework includes
#include "art_root_io/TFileService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RecoBase/Hit.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"

// ROOT includes
#include <TTree.h>

namespace dune {
  class MichelElectron;
}


class dune::MichelElectron : public art::EDAnalyzer {
public:
  explicit MichelElectron(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MichelElectron(MichelElectron const&) = delete;
  MichelElectron(MichelElectron&&) = delete;
  MichelElectron& operator=(MichelElectron const&) = delete;
  MichelElectron& operator=(MichelElectron&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  // Create output Tree
  TTree *fTree;

  // Tree variabiles
  unsigned int fEventID;
  unsigned int Nhits;
  unsigned int Nmctruth;
  // unsigned int Nmctrack;
  unsigned int Ntrack;
  unsigned int NTrackers;
  unsigned int Npfp;
  std::vector<float> nuPDG_truth;
  std::vector<float> ccnc_truth;
  std::vector<float> PDG_pfp;
  std::vector<float> Daughters_pfp;

  std::string fMCParticleModuleLabel;
  std::string fSimChannelModuleLabel;
  std::string fHitsModuleLabel;
  std::string fGenieGenModuleLabel;
  // std::string fMCTrackModuleLabel;
  std::string fTrackModuleLabel;
  std::string fPFPModuleLabel;
};


dune::MichelElectron::MichelElectron(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fMCParticleModuleLabel(p.get<std::string>("MCParticleModuleLabel")),
  fSimChannelModuleLabel(p.get<std::string>("SimChannelModuleLabel")),
  fHitsModuleLabel(p.get<std::string>("HitsModuleLabel")),
  fGenieGenModuleLabel(p.get<std::string>("GenieGenModuleLabel")),
  // fMCTrackModuleLabel(p.get<std::string>("MCTrackModuleLabel")),
  fTrackModuleLabel(p.get<std::string>("TrackModuleLabel")),
  fPFPModuleLabel(p.get<std::string>("PFPModuleLabel"))  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void dune::MichelElectron::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  // Set the event id
  fEventID = e.id().event();

  // Set all counter to 0 for current event
  Nhits = 0;
  Nmctruth = 0;
  // Nmctrack = 0;
  Ntrack = 0;
  NTrackers = 0;
  Npfp = 0;

  // Truth information
  std::vector<art::Ptr<simb::MCParticle>> mcparticlelist;
  if (fMCParticleModuleLabel.size()) {
    art::ValidHandle<std::vector<simb::MCParticle>> mcparticleHandle = e.getValidHandle<std::vector<simb::MCParticle>>(fMCParticleModuleLabel);
    art::fill_ptr_vector(mcparticlelist, mcparticleHandle);
  }

  std::vector<art::Ptr<sim::SimChannel>> simchannellist;
  if (fSimChannelModuleLabel.size()) {
    art::ValidHandle<std::vector<sim::SimChannel>> simchannelHandle = e.getValidHandle<std::vector<sim::SimChannel>>(fSimChannelModuleLabel);
    art::fill_ptr_vector(simchannellist, simchannelHandle);
  }

  //Get Hit information
  art::ValidHandle<std::vector<recob::Hit>> hitHandle = e.getValidHandle<std::vector<recob::Hit>>(fHitsModuleLabel);
  std::vector<art::Ptr<recob::Hit>> hitlist;

  if (hitHandle.isValid())
   art::fill_ptr_vector(hitlist, hitHandle);

  Nhits = hitlist.size(); 

  //Get MC truth information
  art::ValidHandle<std::vector<simb::MCTruth>> mctruthHandle = e.getValidHandle<std::vector<simb::MCTruth>>(fGenieGenModuleLabel);
  std::vector<art::Ptr<simb::MCTruth>> mclist;

  if (mctruthHandle.isValid())
    art::fill_ptr_vector(mclist, mctruthHandle);

  Nmctruth = mclist.size();
  nuPDG_truth.clear();
  ccnc_truth.clear();

  nuPDG_truth.resize(Nmctruth);
  ccnc_truth.resize(Nmctruth);

  for (size_t i =0; i<Nmctruth; i++)
  {
    nuPDG_truth[i] = mclist[i]->GetNeutrino().Nu().PdgCode();
    ccnc_truth[i] = mclist[i]->GetNeutrino().CCNC();
  }
  
  //Get Track information

  NTrackers = fTrackModuleLabel.size();

  // std::vector< art::ValidHandle< std::vector<recob::Track> > > trackHandle(NTrackers);
  // std::vector< std::vector< art::Ptr<recob::Track> > > tracklist(NTrackers);
  // for (unsigned int i=0; i<NTrackers; ++i){
  //   trackHandle[i] = e.getValidHandle<std::vector<recob::Track>>(fTrackModuleLabel);
  //   if (trackHandle[i])
  //     art::fill_ptr_vector(tracklist[i], trackHandle[i]);
  // }

  art::ValidHandle<std::vector<recob::Track>> trackHandle = e.getValidHandle<std::vector<recob::Track>>(fTrackModuleLabel);
  std::vector<art::Ptr<recob::Track>> tracklist;

  art::fill_ptr_vector(tracklist, trackHandle);

  Ntrack = tracklist.size();

  //Get PFParticle information
  art::ValidHandle<std::vector<recob::PFParticle>> pfpHandle = e.getValidHandle<std::vector<recob::PFParticle>>(fPFPModuleLabel);
  std::vector<art::Ptr<recob::PFParticle>> pfplist;

  art::fill_ptr_vector(pfplist, pfpHandle);

  Npfp = pfplist.size();
  PDG_pfp.clear();
  Daughters_pfp.clear();

  PDG_pfp.resize(Npfp);
  Daughters_pfp.resize(Npfp);

  for (size_t i = 0; i<Npfp; i++)
  {
    PDG_pfp[i] = pfplist[i]->PdgCode();
    Daughters_pfp[i] = pfplist[i]->NumDaughters();
  }
  // Fill tree
  fTree->Fill();
}

void dune::MichelElectron::beginJob()
{
  // Implementation of optional member function here.
  // Get the TFileService to create the output TTree for us
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("tree", "Output TTree");

  // Add branches to TTree
  fTree->Branch("eventID", &fEventID);
  fTree->Branch("no_hits", &Nhits);
  // fTree->Branch("no_mctracks", &Nmctrack);
  fTree->Branch("no_trackers", &NTrackers);
  fTree->Branch("no_tracks", &Ntrack);
  fTree->Branch("nuPDG_truth", &nuPDG_truth);
  fTree->Branch("ccnc_truth", &ccnc_truth);
  fTree->Branch("PDG_pfp", &PDG_pfp);
  fTree->Branch("Daughters_pfp", &Daughters_pfp);
}

void dune::MichelElectron::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(dune::MichelElectron)
