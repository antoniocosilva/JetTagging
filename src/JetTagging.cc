#include "JetTagging.h"

/// Cluster/Calorimeter includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calotrigger/CaloTriggerInfo.h>

#include <phool/phool.h>

/// Jet includes
#include <g4jets/Jet.h>
#include <g4jets/JetMap.h>

/// Tracking includes
#include <g4vertex/GlobalVertex.h>
#include <g4vertex/GlobalVertexMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>

/// Truth evaluation includes
#include <g4eval/JetEvalStack.h>
#include <g4eval/SvtxEvalStack.h>

/// HEPMC truth includes
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

/// Fun4All includes
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <KFParticle.h>
#include <kfparticle_sphenix/KFParticle_Container.h>

#include <g4jets/FastJetAlgo.h>
#include <g4jets/Jet.h>
#include <g4jets/Jetv1.h>
//#include <jetbackground/FastJetAlgoSub.h>

#include <fastjet/ClusterSequence.hh>
#include <fastjet/FunctionOfPseudoJet.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

#include <g4main/PHG4Particle.h>            // for PHG4Particle
#include <g4main/PHG4TruthInfoContainer.h>  // for PHG4TruthInfoContainer
#include <g4main/PHG4VtxPoint.h>            // for PHG4VtxPoint
#include <trackbase_historic/SvtxPHG4ParticleMap_v1.h>
#include <kfparticle_sphenix/KFParticle_truthAndDetTools.h>

/// ROOT includes
#include <TFile.h>
#include <TH1.h>
#include <TH1I.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>

/// C++ includes
#include <cassert>
#include <sstream>
#include <string>

using namespace std;

/**
 * JetTagging is a class developed to reconstruct jets containing a D-meson
 * The class can be adapted to tag jets using any kind of particle
 * Author: Antonio Silva (antonio.sphenix@gmail.com)
 */

/**
 * Constructor of module
 */
JetTagging::JetTagging(const std::string &name, const std::string &filename)
  : SubsysReco(name)
  , m_outfilename(filename)
  , m_hm(nullptr)
  , m_track_minpt(0.)
  , m_track_maxpt(9999.)
  , m_track_mineta(-1.1)
  , m_track_maxeta(1.1)
  , m_EMCal_cluster_minpt(0.)
  , m_EMCal_cluster_maxpt(9999.)
  , m_EMCal_cluster_mineta(-1.1)
  , m_EMCal_cluster_maxeta(1.1)
  , m_HCal_cluster_minpt(0.)
  , m_HCal_cluster_maxpt(9999.)
  , m_HCal_cluster_mineta(-1.1)
  , m_HCal_cluster_maxeta(1.1)
  , m_add_tracks(true)
  , m_add_EMCal_clusters(false)
  , m_add_HCal_clusters(false)
  , m_jetr(0.4)
  , m_jetalgo(fastjet::antikt_algorithm)
  , m_recomb_scheme(fastjet::pt_scheme)
  , m_tag_pdg(421)
  , m_qualy_plots(false)
{
  /// Initialize variables and trees so we don't accidentally access
  /// memory that was never allocated
  initializeVariables();
  initializeTrees();
}

/**
 * Destructor of module
 */
JetTagging::~JetTagging()
{
  delete m_hm;
  delete m_taggedjettree;
}

/**
 * Initialize the module and prepare looping over events
 */
int JetTagging::Init(PHCompositeNode *topNode)
{
  if (Verbosity() > 5)
  {
    cout << "Beginning Init in JetTagging" << endl;
  }

  m_outfile = new TFile(m_outfilename.c_str(), "RECREATE");

  m_eventcount_h = new TH1I("eventcount_h", "Event Count", 2, -0.5, 1.5);
  m_eventcount_h->GetXaxis()->SetBinLabel(1,"N ev anal");
  m_eventcount_h->GetXaxis()->SetBinLabel(2,"N D candidates");

  m_rec_tagpart_pt = new TH1D("rec_tagpart_pt", ";#it{p}_{T,D} (GeV/#it{c});Entries", 20, 0., 20.);
  m_gen_withrec_tagpart_pt = new TH1D("gen_withrec_tagpart_pt", ";#it{p}_{T,D} (GeV/#it{c});Entries", 20, 0., 20.);
  m_gen_tagpart_pt = new TH1D("gen_tagpart_pt", ";#it{p}_{T,D} (GeV/#it{c});Entries", 20, 0., 20.);

  m_rec_tracks_pt = new TH1D("rec_tracks_pt", ";#it{p}_{T} (GeV/#it{c});Entries", 200, 0., 20.);
  m_gen_tracks_pt = new TH1D("gen_tracks_pt", ";#it{p}_{T} (GeV/#it{c});Entries", 200, 0., 20.);

  m_rec_emcal_clusters_pt = new TH1D("rec_emcal_clusters_pt", ";#it{p}_{T} (GeV/#it{c});Entries", 200, 0., 20.);
  m_rec_hcalin_clusters_pt = new TH1D("rec_hcalin_clusters_pt", ";#it{p}_{T} (GeV/#it{c});Entries", 200, 0., 20.);
  m_rec_hcalout_clusters_pt = new TH1D("rec_hcalout_clusters_pt", ";#it{p}_{T} (GeV/#it{c});Entries", 200, 0., 20.);

  return 0;
}

/**
 * Main workhorse function where each event is looped over and
 * data from each event is collected from the node tree for analysis
 */
int JetTagging::process_event(PHCompositeNode *topNode)
{
  m_eventcount_h->Fill(0);


  KFParticle_Container* kfContainer = new KFParticle_Container();

  PHNodeIterator nodeIter(topNode);
  PHNode *findNode = dynamic_cast<PHNode *>(nodeIter.findFirst("DST"));
  if (findNode)
  {
    findNode = dynamic_cast<PHNode *>(nodeIter.findFirst("reconstructedParticles_KFParticle_Container"));
    if (findNode)
    {
      kfContainer = findNode::getClass<KFParticle_Container>(topNode, "reconstructedParticles_KFParticle_Container");
      //cout << "****************** CONTAINER FOUND with size: " << kfContainer->size() << endl;
    }
  }

  KFParticle *D0cand;

  KFParticle *D0dau[2];

  int nDaughters = 2;//D0cand->NDaughters() is returning 0, bug?

  for(unsigned int i = 0; i < kfContainer->size(); i++)
  {
    D0cand = kfContainer->get(i);
    if(TMath::Abs(D0cand->GetPDG()) == m_tag_pdg)
    {
      //cout << "Rec D0 px: " << D0cand->Px() << " py: " << D0cand->Py() << " pz: " << D0cand->Pz() << endl;
      m_eventcount_h->Fill(1);
      D0dau[0] = kfContainer->get(i+1);
      D0dau[1] = kfContainer->get(i+2);

      //cout << "D0: " << D0cand->GetPDG();
      cout << "D0dau0: " << D0dau[0]->GetPDG() << endl;
      cout << "D0dau0 ID: " << D0dau[0]->Id() << endl;
      cout << "D0dau1: " << D0dau[1]->GetPDG() << endl;
      cout << "D0dau1 ID: " << D0dau[1]->Id() << endl;

      recTaggedJets(topNode, D0cand, D0dau, nDaughters);
      recMCTaggedJets(topNode, D0dau, nDaughters);
      m_rec_tagpart_pt->Fill(D0cand->GetPt());
      m_taggedjettree->Fill();
      i += nDaughters; //Go to the next D meson
    }
  }

  if((kfContainer->size() > 0) && m_qualy_plots) doMCLoop(topNode);



  return Fun4AllReturnCodes::EVENT_OK;
}

/**
 * End the module and finish any data collection. Clean up any remaining
 * loose ends
 */
int JetTagging::End(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    cout << "Ending JetTagging analysis package" << endl;
  }

  /// Change to the outfile
  m_outfile->cd();


  m_taggedjettree->Write();

  m_eventcount_h->Write();

  m_rec_tagpart_pt->Write();
  m_gen_withrec_tagpart_pt->Write();
  m_gen_tagpart_pt->Write();

  m_rec_tracks_pt->Write();
  m_gen_tracks_pt->Write();

  m_rec_emcal_clusters_pt->Write();
  m_rec_hcalin_clusters_pt->Write();
  m_rec_hcalout_clusters_pt->Write();

  /// Write and close the outfile
  m_outfile->Write();
  m_outfile->Close();

  /// Clean up pointers and associated histos/trees in TFile
  delete m_outfile;

  if (Verbosity() > 1)
  {
    cout << "Finished JetTagging analysis package" << endl;
  }

  return 0;
}

void JetTagging::recTaggedJets(PHCompositeNode *topNode, KFParticle *Tag, KFParticle *TagDecays[], int nDecays)
{
  fastjet::JetDefinition* jetdef = new fastjet::JetDefinition(m_jetalgo, m_jetr, m_recomb_scheme, fastjet::Best);

  std::vector<fastjet::PseudoJet> particles;

  Jet* taggedJet = new Jetv1();

  fastjet::PseudoJet fjTag(Tag->Px(), Tag->Py(), Tag->Pz(), Tag->E());
  fjTag.set_user_index(0); //index 0 is the tag particle
  particles.push_back(fjTag);
  //int idpart = 1; //index 0 is the tag particle

  if(m_add_tracks) addTracks(topNode, particles, TagDecays, nDecays);
  if(m_add_EMCal_clusters || m_add_HCal_clusters) addClusters(topNode, particles);

  fastjet::ClusterSequence jetFinder(particles, *jetdef);
  std::vector<fastjet::PseudoJet> fastjets = jetFinder.inclusive_jets();

  delete jetdef;

  for (unsigned int ijet = 0; ijet < fastjets.size(); ++ijet)
  {
    std::vector<fastjet::PseudoJet> comps = fastjets[ijet].constituents();
    for (unsigned int icomp = 0; icomp < comps.size(); ++icomp)
    {
      if(comps[icomp].user_index() == 0)
      {
        //cout << "D0 found in jet #" << ijet << endl;
        taggedJet->set_px(fastjets[ijet].px());
        taggedJet->set_py(fastjets[ijet].py());
        taggedJet->set_pz(fastjets[ijet].pz());
        taggedJet->set_e(fastjets[ijet].e());
        break;
      }
    }
  }

  //cout << "D0 pt: " << Tag->GetPt() << endl;

  m_tagpartpx = Tag->Px();
  m_tagpartpy = Tag->Py();
  m_tagpartpz = Tag->Pz();
  m_tagpartpt = Tag->GetPt();
  m_tagparteta = Tag->GetEta();
  m_tagpartphi = Tag->GetPhi();
  m_tagpartm = Tag->GetMass();

  //cout << "D0-Jet pt: " << taggedJet->get_pt() << endl;

  m_tagjetpx = taggedJet->get_px();
  m_tagjetpy = taggedJet->get_py();
  m_tagjetpz = taggedJet->get_pz();
  m_tagjetpt = taggedJet->get_pt();
  m_tagjeteta = taggedJet->get_eta();
  m_tagjetphi = taggedJet->get_phi();



  return;

}

void JetTagging::addTracks(PHCompositeNode *topNode, std::vector<fastjet::PseudoJet> &particles, KFParticle *TagDecays[], int nDecays)
{
  SvtxTrackMap *trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");

  int idpart = particles.size();
  //cout << "Tracks idpart = " << idpart << endl;

  if (!trackmap)
  {
    cout << PHWHERE
         << "SvtxTrackMap node is missing, can't collect tracks"
         << endl;
    return;
  }

  for (SvtxTrackMap::Iter iter = trackmap->begin();
       iter != trackmap->end();
       ++iter)
  {
    SvtxTrack *track = iter->second;

    if(!isAcceptableTrack(track)) continue;

    if(m_qualy_plots) m_rec_tracks_pt->Fill(track->get_pt());

    if(isDecay(track, TagDecays, nDecays)) continue;

    fastjet::PseudoJet fjTrack(track->get_px(), track->get_py(), track->get_pz(), 0.);

    fjTrack.set_user_index(idpart);
    particles.push_back(fjTrack);
    idpart++;
  }


}

bool JetTagging::isAcceptableTrack(SvtxTrack *track)
{
  if((track->get_pt() < m_track_minpt) || (track->get_pt() > m_track_maxpt)) return false;
  if((track->get_eta() < m_track_mineta) || (track->get_eta() > m_track_maxeta)) return false;

  return true;
}


void JetTagging::addClusters(PHCompositeNode *topNode, std::vector<fastjet::PseudoJet> &particles)
{

  GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!vertexmap)
  {
    cout << "JetTagging::getEmcalClusters - Fatal Error - GlobalVertexMap node is missing. Please turn on the do_global flag in the main macro in order to reconstruct the global vertex." << endl;
    assert(vertexmap);  // force quit

    return;
  }

  if (vertexmap->empty())
  {
    cout << "JetTagging::getEmcalClusters - Fatal Error - GlobalVertexMap node is empty. Please turn on the do_global flag in the main macro in order to reconstruct the global vertex." << endl;
    return;
  }

  GlobalVertex *vtx = vertexmap->begin()->second;
  if (vtx == nullptr)
    return;

  // Get the current position in the particles vector
  int idpart = particles.size();

  //EMCAL Clusters
  if(m_add_EMCal_clusters)
  {
    RawClusterContainer *clustersEMC = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_CEMC");

    if (!clustersEMC)
    {
      cout << PHWHERE
           << "EMCal cluster node is missing, can't collect EMCal clusters"
           << endl;
      return;
    }

    RawClusterContainer::ConstRange begin_end_EMC = clustersEMC->getClusters();
    RawClusterContainer::ConstIterator clusIter_EMC;


    /// Loop over the EMCal clusters
    for (clusIter_EMC = begin_end_EMC.first;
         clusIter_EMC != begin_end_EMC.second;
         ++clusIter_EMC)
    {

      const RawCluster *cluster = clusIter_EMC->second;

      CLHEP::Hep3Vector vertex(vtx->get_x(), vtx->get_y(), vtx->get_z());
      CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*cluster, vertex);

      if(m_qualy_plots)
      {
        //No pt cut, but within eta acceptance
        if((E_vec_cluster.pseudoRapidity() > m_EMCal_cluster_mineta) && (E_vec_cluster.pseudoRapidity() < m_EMCal_cluster_maxeta)) m_rec_emcal_clusters_pt->Fill(E_vec_cluster.perp());
      }

      if(!isAcceptableEMCalCluster(E_vec_cluster)) continue;

      double cluster_e = E_vec_cluster.mag();

      double cluster_pt = E_vec_cluster.perp();

      double cluster_phi = E_vec_cluster.getPhi();

      double cluster_px = cluster_pt * cos(cluster_phi);
      double cluster_py = cluster_pt * sin(cluster_phi);
      double cluster_pz = sqrt(cluster_e * cluster_e - cluster_px * cluster_px - cluster_py * cluster_py);

      fastjet::PseudoJet fjCluster(cluster_px, cluster_py, cluster_pz, cluster_e);

      fjCluster.set_user_index(idpart);
      particles.push_back(fjCluster);
      idpart++;
    }
  }

  if(m_add_HCal_clusters)
  {
    //HCAL Clusters
    RawClusterContainer *clustersHCALIN = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_HCALIN");

    if (!clustersHCALIN)
    {
      cout << PHWHERE
           << "EMCal cluster node is missing, can't collect EMCal clusters"
           << endl;
      return;
    }

    RawClusterContainer::ConstRange begin_end_HCALIN = clustersHCALIN->getClusters();
    RawClusterContainer::ConstIterator clusIter_HCALIN;

    //cout << "idpart before HCALIN: " << idpart << endl;

    /// Loop over the EMCal clusters
    for (clusIter_HCALIN = begin_end_HCALIN.first;
         clusIter_HCALIN != begin_end_HCALIN.second;
         ++clusIter_HCALIN)
    {
      /// Get this cluster
      const RawCluster *cluster = clusIter_HCALIN->second;

      CLHEP::Hep3Vector vertex(vtx->get_x(), vtx->get_y(), vtx->get_z());
      CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*cluster, vertex);

      if(m_qualy_plots)
      {
        //No pt cut, but within eta acceptance
        if((E_vec_cluster.pseudoRapidity() > m_EMCal_cluster_mineta) && (E_vec_cluster.pseudoRapidity() < m_EMCal_cluster_maxeta)) m_rec_hcalin_clusters_pt->Fill(E_vec_cluster.perp());
      }

      if(!isAcceptableHCalCluster(E_vec_cluster)) continue;

      double cluster_e = E_vec_cluster.mag();
      double cluster_pt = E_vec_cluster.perp();
      double cluster_phi = E_vec_cluster.getPhi();

      double cluster_px = cluster_pt * cos(cluster_phi);
      double cluster_py = cluster_pt * sin(cluster_phi);
      double cluster_pz = sqrt(cluster_e * cluster_e - cluster_px * cluster_px - cluster_py * cluster_py);

      fastjet::PseudoJet fjCluster(cluster_px, cluster_py, cluster_pz, cluster_e);

      fjCluster.set_user_index(idpart);
      particles.push_back(fjCluster);
      idpart++;


    }

    RawClusterContainer *clustersHCALOUT = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_HCALOUT");

    if (!clustersHCALOUT)
    {
      cout << PHWHERE
           << "EMCal cluster node is missing, can't collect EMCal clusters"
           << endl;
      return;
    }

    RawClusterContainer::ConstRange begin_end_HCALOUT = clustersHCALOUT->getClusters();
    RawClusterContainer::ConstIterator clusIter_HCALOUT;

    //cout << "idpart before HCALOUT: " << idpart << endl;

    /// Loop over the EMCal clusters
    for (clusIter_HCALOUT = begin_end_HCALOUT.first;
         clusIter_HCALOUT != begin_end_HCALOUT.second;
         ++clusIter_HCALOUT)
    {
      /// Get this cluster
      const RawCluster *cluster = clusIter_HCALOUT->second;

      /// Get cluster characteristics
      /// This helper class determines the photon characteristics
      /// depending on the vertex position
      /// This is important for e.g. eta determination and E_T determination
      CLHEP::Hep3Vector vertex(vtx->get_x(), vtx->get_y(), vtx->get_z());
      CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*cluster, vertex);

      if(m_qualy_plots)
      {
        //No pt cut, but within eta acceptance
        if((E_vec_cluster.pseudoRapidity() > m_EMCal_cluster_mineta) && (E_vec_cluster.pseudoRapidity() < m_EMCal_cluster_maxeta)) m_rec_hcalout_clusters_pt->Fill(E_vec_cluster.perp());
      }

      if(!isAcceptableHCalCluster(E_vec_cluster)) continue;

      double cluster_e = E_vec_cluster.mag();
      //double cluster_eta = E_vec_cluster.pseudoRapidity();
      //double cluster_theta = E_vec_cluster.getTheta();
      double cluster_pt = E_vec_cluster.perp();
      double cluster_phi = E_vec_cluster.getPhi();

      double cluster_px = cluster_pt * cos(cluster_phi);
      double cluster_py = cluster_pt * sin(cluster_phi);
      double cluster_pz = sqrt(cluster_e * cluster_e - cluster_px * cluster_px - cluster_py * cluster_py);

      fastjet::PseudoJet fjCluster(cluster_px, cluster_py, cluster_pz, cluster_e);

      fjCluster.set_user_index(idpart);
      particles.push_back(fjCluster);
      idpart++;


    }
  }



  //cout << "idpart after HCALOUT: " << idpart << endl;


}

bool JetTagging::isAcceptableEMCalCluster(CLHEP::Hep3Vector &E_vec_cluster)
{
  if((E_vec_cluster.perp() < m_EMCal_cluster_minpt) || (E_vec_cluster.perp() > m_EMCal_cluster_maxpt)) return false;
  if((E_vec_cluster.pseudoRapidity() < m_EMCal_cluster_mineta) || (E_vec_cluster.pseudoRapidity() > m_EMCal_cluster_maxeta)) return false;

  return true;
}

bool JetTagging::isAcceptableHCalCluster(CLHEP::Hep3Vector &E_vec_cluster)
{
  if((E_vec_cluster.perp() < m_HCal_cluster_minpt) || (E_vec_cluster.perp() > m_HCal_cluster_maxpt)) return false;
  if((E_vec_cluster.pseudoRapidity() < m_HCal_cluster_mineta) || (E_vec_cluster.pseudoRapidity() > m_HCal_cluster_maxeta)) return false;

  return true;
}

bool JetTagging::isSameParticle(SvtxTrack *track, KFParticle *particle)
{
  if(TMath::Abs(track->get_px()-particle->Px()) < 0.00001
     && TMath::Abs(track->get_py()-particle->Py()) < 0.00001
     && TMath::Abs(track->get_pz()-particle->Pz()) < 0.00001) return true;

  return false;
}

bool JetTagging::isDecay(SvtxTrack *track, KFParticle *decays[], int nDecays)
{
  for(int idecay = 0; idecay < nDecays; idecay++)
  {
    //if(track->get_id() == decays[idecay]->Id()) return true; //KFParticle is not storing SvtxTrack IDs
    if(isSameParticle(track, decays[idecay])) return true;
  }
  return false;
}

void JetTagging::recMCTaggedJets(PHCompositeNode *topNode, KFParticle *decays[], int nDecays)
{
  fastjet::JetDefinition* jetdef = new fastjet::JetDefinition(m_jetalgo, m_jetr, m_recomb_scheme, fastjet::Best);

  std::vector<fastjet::PseudoJet> particles;

  Jet* mcTaggedJet = new Jetv1();

  PHG4Particle *mcDaughters[nDecays];

  HepMC::GenParticle *mcTag = findMCTag(topNode, decays, nDecays, mcDaughters);

  if(!mcTag) return;

  fastjet::PseudoJet fjMCTag(mcTag->momentum().px(), mcTag->momentum().py(), mcTag->momentum().pz(), mcTag->momentum().e());
  fjMCTag.set_user_index(0); //index 0 is the tag particle
  particles.push_back(fjMCTag);

  //HepMC particle loop
  PHHepMCGenEventMap *hepmceventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");

  /// If the node was not properly put on the tree, return
  if (!hepmceventmap)
  {
    cout << PHWHERE
         << "HEPMC event map node is missing, can't collected HEPMC truth particles"
         << endl;
    return;
  }

  PHHepMCGenEvent *hepmcevent = hepmceventmap->get(1);

  if(!hepmcevent) return;

  HepMC::GenEvent* hepMCevent = hepmcevent->getEvent();

  if(!hepMCevent) return;

  int idpart = particles.size();

  TDatabasePDG *database = TDatabasePDG::Instance();
  TParticlePDG *partPDG = 0;

  for (HepMC::GenEvent::particle_const_iterator p = hepMCevent->particles_begin(); p != hepMCevent->particles_end(); ++p)
  {
    if(((*p)->momentum().eta() < m_track_mineta) || ((*p)->momentum().eta() > m_track_maxeta)) continue; //Maybe make it depend on charge? (Track or cluster)
    partPDG = database->GetParticle((*p)->pdg_id());
    double hepmcPartPt = TMath::Sqrt(((*p)->momentum().px()*(*p)->momentum().px()) + ((*p)->momentum().py()*(*p)->momentum().py()));
    if(partPDG->Charge() == 0)
    {

      if((hepmcPartPt < m_EMCal_cluster_minpt) || (hepmcPartPt > m_EMCal_cluster_maxpt)) continue;
    }
    else
    {
      if((hepmcPartPt < m_track_minpt) || (hepmcPartPt > m_track_maxpt)) continue;
    }
    if((*p)->status() > 1) continue; //Get only final state particles

    if(isDecay((*p), mcDaughters, nDecays)) continue; //remove daughters

    if(mcTag->barcode() == (*p)->barcode()) continue; //if D0 is not stable, this is redundant

    fastjet::PseudoJet fjMCParticle((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e());
    fjMCParticle.set_user_index(idpart);
    particles.push_back(fjMCParticle);
    idpart++;
  }


  fastjet::ClusterSequence jetFinder(particles, *jetdef);
  std::vector<fastjet::PseudoJet> mcfastjets = jetFinder.inclusive_jets();

  delete jetdef;

  for (unsigned int ijet = 0; ijet < mcfastjets.size(); ++ijet)
  {
    std::vector<fastjet::PseudoJet> comps = mcfastjets[ijet].constituents();
    for (unsigned int icomp = 0; icomp < comps.size(); ++icomp)
    {
      if(comps[icomp].user_index() == 0)
      {
        //cout << "D0 found in jet #" << ijet << endl;
        mcTaggedJet->set_px(mcfastjets[ijet].px());
        mcTaggedJet->set_py(mcfastjets[ijet].py());
        mcTaggedJet->set_pz(mcfastjets[ijet].pz());
        mcTaggedJet->set_e(mcfastjets[ijet].e());
        break;
      }
    }
  }

  m_truth_tagpartpx = mcTag->momentum().px();
  m_truth_tagpartpy = mcTag->momentum().py();
  m_truth_tagpartpz = mcTag->momentum().pz();
  m_truth_tagpartpt = sqrt(m_truth_tagpartpx * m_truth_tagpartpx + m_truth_tagpartpy * m_truth_tagpartpy);
  m_truth_tagparteta = atanh(m_truth_tagpartpz / mcTag->momentum().e());
  m_truth_tagpartphi = atan(m_truth_tagpartpy / m_truth_tagpartpx);


  m_truth_tagjetpx = mcTaggedJet->get_px();
  m_truth_tagjetpy = mcTaggedJet->get_py();
  m_truth_tagjetpz = mcTaggedJet->get_pz();
  m_truth_tagjetpt = mcTaggedJet->get_pt();
  m_truth_tagjeteta = mcTaggedJet->get_eta();
  m_truth_tagjetphi = mcTaggedJet->get_phi();

  m_gen_withrec_tagpart_pt->Fill(m_truth_tagpartpt);

}

HepMC::GenParticle* JetTagging::findMCTag(PHCompositeNode *topNode, KFParticle *decays[], int nDecays, PHG4Particle *mcDaughters[])
{
  PHG4Particle *g4particle = nullptr;
  SvtxTrack *dauTrack = nullptr;
  HepMC::GenParticle *mcTag[nDecays];
  int mcDauID = 0;

  PHNodeIterator nodeIter(topNode);
  PHNode *findNode = dynamic_cast<PHNode *>(nodeIter.findFirst("SvtxPHG4ParticleMap"));
  PHG4TruthInfoContainer *truthinfo = nullptr;
  if (findNode)
  {
    findNode = dynamic_cast<PHNode *>(nodeIter.findFirst("G4TruthInfo"));
    if (findNode)
    {
      truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    }
    else
    {
      std::cout << "KFParticle truth matching: G4TruthInfo does not exist" << std::endl;
      return 0;
    }
  }
  //Truth map
  SvtxPHG4ParticleMap_v1 *dst_reco_truth_map = findNode::getClass<SvtxPHG4ParticleMap_v1>(topNode, "SvtxPHG4ParticleMap");

  if(!dst_reco_truth_map) return 0;


  SvtxTrackMap *trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");

  if (!trackmap)
  {
    cout << PHWHERE
         << "SvtxTrackMap node is missing, can't collect tracks"
         << endl;
    return 0;
  }


  for (SvtxTrackMap::Iter iter = trackmap->begin();
       iter != trackmap->end();
       ++iter)
  {
    SvtxTrack *track = iter->second;

    for(int idecay = 0; idecay < nDecays; idecay++)
    {
      if(isSameParticle(track, decays[idecay]))
      {
        dauTrack = track;
        std::map<float, std::set<int>> truth_set = dst_reco_truth_map->get(dauTrack->get_id());
        const auto& best_weight = truth_set.rbegin();
        int best_truth_id =  *best_weight->second.rbegin();
        g4particle = truthinfo->GetParticle(best_truth_id);


        mcDaughters[mcDauID] = g4particle;
        mcDauID++;

        mcTag[idecay] = getMother(topNode, g4particle);
        if(mcTag[idecay] == 0) return 0;

      }
    }
  }
  //check is decays are from the same mother, otherwise it is background
  for(int idecay = 1; idecay < nDecays; idecay++)
  {
    if(mcTag[idecay]->barcode() != mcTag[idecay-1]->barcode()) return 0;
  }

  return mcTag[0];

}

HepMC::GenParticle* JetTagging::getMother(PHCompositeNode *topNode, PHG4Particle *g4daughter)
{
  PHHepMCGenEventMap *hepmceventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");

  /// If the node was not properly put on the tree, return
  if (!hepmceventmap)
  {
    cout << PHWHERE
         << "HEPMC event map node is missing, can't collected HEPMC truth particles"
         << endl;
    return 0;
  }

  PHHepMCGenEvent *hepmcevent = hepmceventmap->get(1);
  if(!hepmcevent) return 0;

  HepMC::GenEvent* hepMCGenEvent = hepmcevent->getEvent();
  if(!hepMCGenEvent) return 0;

  HepMC::GenParticle *mcDaughter = 0;

  if(g4daughter->get_barcode() > 0) mcDaughter = hepMCGenEvent->barcode_to_particle(g4daughter->get_barcode());
  if(!mcDaughter) return 0;

  HepMC::GenVertex* TagVertex = mcDaughter->production_vertex();
  for (HepMC::GenVertex::particle_iterator it = TagVertex->particles_begin(HepMC::ancestors); it != TagVertex->particles_end(HepMC::ancestors); ++it)
  {
    if(TMath::Abs((*it)->pdg_id()) == m_tag_pdg) return (*it);
  }

  return 0;

}

bool JetTagging::isDecay(HepMC::GenParticle *particle, PHG4Particle *decays[], int nDecays)
{
  for(int idecay = 0; idecay < nDecays; idecay++)
  {
    //if(isSameParticle(particle, decays[idecay])) return true;
    if(particle->barcode() == decays[idecay]->get_barcode()) return true;
  }
  return false;
}

void JetTagging::doMCLoop(PHCompositeNode *topNode)
{
  PHHepMCGenEventMap *hepmceventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");

  /// If the node was not properly put on the tree, return
  if (!hepmceventmap)
  {
    cout << PHWHERE
         << "HEPMC event map node is missing, can't collected HEPMC truth particles"
         << endl;
    return;
  }

  cout << "hepmceventmap size: " << hepmceventmap->size() << endl;

  PHHepMCGenEvent *hepmcevent = hepmceventmap->get(1);

  if(!hepmcevent) return;

  HepMC::GenEvent* hepMCGenEvent = hepmcevent->getEvent();

  if(!hepMCGenEvent) return;

  TDatabasePDG *database = TDatabasePDG::Instance();
  TParticlePDG *partPDG = 0;

  for (HepMC::GenEvent::particle_const_iterator p = hepMCGenEvent->particles_begin(); p != hepMCGenEvent->particles_end(); ++p)
  {
    if(((*p)->momentum().eta() < m_track_mineta) || ((*p)->momentum().eta() > m_track_maxeta)) continue;

    partPDG = database->GetParticle((*p)->pdg_id());
    double hepmcPartPt = TMath::Sqrt(((*p)->momentum().px()*(*p)->momentum().px()) + ((*p)->momentum().py()*(*p)->momentum().py()));

    if(TMath::Abs((*p)->pdg_id()) == m_tag_pdg)
    {
      m_gen_tagpart_pt->Fill(hepmcPartPt); //For now, D0 is stable. This should change
      HepMC::GenVertex* TagVertex = (*p)->end_vertex();
      for (HepMC::GenVertex::particle_iterator it = TagVertex->particles_begin(HepMC::descendants); it != TagVertex->particles_end(HepMC::descendants); ++it)
         cout << "daughter PDG: " << (*it)->pdg_id() << "daughter barcode: " << (*it)->barcode() << " px: " << (*it)->momentum().px() << " py: " << (*it)->momentum().py() << " pz: " << (*it)->momentum().pz() << endl;
      continue; //Remove this stable D0
    }

    if((*p)->status() > 1) continue;



    if(partPDG->Charge() != 0)
    {
      if((hepmcPartPt < m_track_minpt) || (hepmcPartPt > m_track_maxpt)) continue;
      m_gen_tracks_pt->Fill(hepmcPartPt);
    }
  }

}




/**
 * This function puts all of the tree branch assignments in one place so as to not
 * clutter up the JetTagging::Init function.
 */
void JetTagging::initializeTrees()
{


  m_taggedjettree = new TTree("m_taggedjettree", "A tree with Tagged-Jet info");
  m_taggedjettree->Branch("m_tagpartpx", &m_tagpartpx, "m_tagpartpx/D");
  m_taggedjettree->Branch("m_tagpartpy", &m_tagpartpy, "m_tagpartpy/D");
  m_taggedjettree->Branch("m_tagpartpz", &m_tagpartpz, "m_tagpartpz/D");
  m_taggedjettree->Branch("m_tagpartpt", &m_tagpartpt, "m_tagpartpt/D");
  m_taggedjettree->Branch("m_tagparteta", &m_tagparteta, "m_tagparteta/D");
  m_taggedjettree->Branch("m_tagpartphi", &m_tagpartphi, "m_tagpartphi/D");
  m_taggedjettree->Branch("m_tagpartm", &m_tagpartm, "m_tagpartm/D");
  m_taggedjettree->Branch("m_tagjetpx", &m_tagjetpx, "m_tagjetpx/D");
  m_taggedjettree->Branch("m_tagjetpy", &m_tagjetpy, "m_tagjetpy/D");
  m_taggedjettree->Branch("m_tagjetpz", &m_tagjetpz, "m_tagjetpz/D");
  m_taggedjettree->Branch("m_tagjetpt", &m_tagjetpt, "m_tagjetpt/D");
  m_taggedjettree->Branch("m_tagjeteta", &m_tagjeteta, "m_tagjeteta/D");
  m_taggedjettree->Branch("m_tagjetphi", &m_tagjetphi, "m_tagjetphi/D");

  m_taggedjettree->Branch("m_truth_tagpartpx", &m_truth_tagpartpx, "m_truth_tagpartpx/D");
  m_taggedjettree->Branch("m_truth_tagpartpy", &m_truth_tagpartpy, "m_truth_tagpartpy/D");
  m_taggedjettree->Branch("m_truth_tagpartpz", &m_truth_tagpartpz, "m_truth_tagpartpz/D");
  m_taggedjettree->Branch("m_truth_tagpartpt", &m_truth_tagpartpt, "m_truth_tagpartpt/D");
  m_taggedjettree->Branch("m_truth_tagparteta", &m_truth_tagparteta, "m_truth_tagparteta/D");
  m_taggedjettree->Branch("m_truth_tagpartphi", &m_truth_tagpartphi, "m_truth_tagpartphi/D");
  m_taggedjettree->Branch("m_truth_tagjetpx", &m_truth_tagjetpx, "m_truth_tagjetpx/D");
  m_taggedjettree->Branch("m_truth_tagjetpy", &m_truth_tagjetpy, "m_truth_tagjetpy/D");
  m_taggedjettree->Branch("m_truth_tagjetpz", &m_truth_tagjetpz, "m_truth_tagjetpz/D");
  m_taggedjettree->Branch("m_truth_tagjetpt", &m_truth_tagjetpt, "m_truth_tagjetpt/D");
  m_taggedjettree->Branch("m_truth_tagjeteta", &m_truth_tagjeteta, "m_truth_tagjeteta/D");
  m_taggedjettree->Branch("m_truth_tagjetphi", &m_truth_tagjetphi, "m_truth_tagjetphi/D");

}

/**
 * This function initializes all of the member variables in this class so that there
 * are no variables that might not be set before e.g. writing them to the output
 * trees.
 */
void JetTagging::initializeVariables()
{
  m_outfile = new TFile();
}
