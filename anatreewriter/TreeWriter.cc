#include "TreeWriter.h"


#include <fun4all/PHTFileServer.h>
#include <fun4all/Fun4AllServer.h>

#include <ffaobjects/EventHeader.h>
#include <ffaobjects/EventHeaderv1.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/GlobalVertexMapv1.h>

#include <centrality/CentralityInfo.h>
#include <centrality/CentralityInfov1.h>
#include <centrality/CentralityInfov2.h>

#include <eventplaneinfo/Eventplaneinfo.h>
#include <eventplaneinfo/EventplaneinfoMap.h>

#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

#include <epd/EpdGeom.h>

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>


#include <mbd/MbdOut.h>
#include <mbd/MbdOutV1.h>
#include <mbd/MbdOutV2.h>
#include <mbd/MbdPmtHit.h>
#include <mbd/MbdPmtContainerV1.h>


#include <ffarawobjects/Gl1Packet.h>
#include <ffarawobjects/Gl1Packetv2.h>

#include <jetbase/Jet.h>
#include <jetbase/Jetv2.h>
#include <jetbase/JetContainer.h>
#include <jetbase/JetContainerv1.h>

#include <jetbackground/TowerRho.h>
#include <jetbackground/TowerRhov1.h>
#include <jetbackground/TowerBackground.h>
#include <jetbackground/TowerBackgroundv1.h>

#include <TTree.h>
#include <TLorentzVector.h>
#include <TRandom3.h>

// standard includes
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <utility>
#include <vector>
#include <cassert>


int TreeWriter::Init( PHCompositeNode * /*topNode*/ )
{
  
  // create output file
  PHTFileServer::get().open( m_output_filename, "RECREATE" );

  if ( Verbosity () > 0 ) 
  {
    std::cout << "TreeWriter::Init - opening file " << m_output_filename << std::endl;
  } 

  // Reset counters
  m_run_tree = new TTree( "RunTree", "RunTree" );
  m_run_tree -> Branch( "num_events", &m_nevents, "num_events/I" );
  m_run_tree -> Branch( "weight", &m_weight, "weight/F" );

  // create tree
  m_tree = new TTree( "EventTree", "EventTree" );
  m_tree -> Branch( "event_id", &m_event_id, "event_id/I" );
  m_event_id = -1;

  // gl1 info
  if ( !m_gl1_node.empty() ) 
  {
    m_tree -> Branch( "s_triggervec", &s_triggervec, "s_triggervec/l" );
    m_tree -> Branch( "l_triggervec", &l_triggervec, "l_triggervec/l" );
    if ( Verbosity() > 0 ) 
    {
      std::cout << "TreeWriter::Init - GL1 node: " << m_gl1_node << std::endl;
    }
  }

  // zvtx
  if ( !m_zvrtx_node.empty() ) 
  {
    m_tree -> Branch( "zvrtx", &m_zvtx, "zvrtx/F" );
    if ( Verbosity() > 0 ) 
    {
      std::cout << "TreeWriter::Init - Registered zvtx" << std::endl;
    }
  }

  // centrality
  if ( !m_cent_node.empty() ) 
  {
    m_tree->Branch("cent", &m_cent, "cent/I");
    if ( Verbosity() > 0 ) 
    {
      std::cout << "TreeWriter::Init - Registered centrality" << std::endl;
    }
  }

  // eventplane info
  if ( !m_eventplane_node.empty() )
  {
    m_tree -> Branch( "psi_shifted_NS", &m_psi_shifted_NS);
    m_tree -> Branch( "psi_shifted_S", &m_psi_shifted_S);
    m_tree -> Branch( "psi_shifted_N", &m_psi_shifted_N);
    m_tree -> Branch( "psi_NS", &m_psi_NS);
    m_tree -> Branch( "psi_S", &m_psi_S);
    m_tree -> Branch( "psi_N", &m_psi_N);
    if ( Verbosity() > 0 ) 
    {
      std::cout << "TreeWriter::Init - Registered Event Plane nodes: " << m_eventplane_node << std::endl;
    }
  }

  // sepd infoK
  if ( !m_sepd_node.empty() )
  {
    m_tree -> Branch( "sepd_energy", &m_sepd_energy);
    m_tree -> Branch( "sepd_isgood", &m_sepd_isgood);
    m_tree -> Branch( "sepd_arm", &m_sepd_arm);
    m_tree -> Branch( "sepd_radius", &m_sepd_radius);
    m_tree -> Branch( "sepd_phi", &m_sepd_phi);
    if ( Verbosity() > 0 ) 
    {
      std::cout << "TreeWriter::Init - Registered sEPD nodes: " << m_sepd_node << std::endl;
    }
  }
  
  // event header info
  if ( !m_eventhead_node.empty() ) 
  {
    m_tree -> Branch( "b", &m_b, "b/F" );
    m_tree -> Branch( "ep_angle", &m_ep_angle, "ep_angle/F" );
    m_tree -> Branch( "ecc", &m_ecc, "ecc/F" );
    m_tree -> Branch( "psi1", &m_psi1, "psi1/F" );
    m_tree -> Branch( "psi2", &m_psi2, "psi2/F" );
    m_tree -> Branch( "psi3", &m_psi3, "psi3/F" );
    m_tree -> Branch( "psi4", &m_psi4, "psi4/F" );
    m_tree -> Branch( "psi5", &m_psi5, "psi5/F" );
    m_tree -> Branch( "psi6", &m_psi6, "psi6/F" );
    m_tree -> Branch( "ncoll", &m_ncoll, "ncoll/F" );
    m_tree -> Branch( "npart", &m_npart, "npart/F" );
    if ( Verbosity() > 0 ) 
    {
      std::cout << "TreeWriter::Init - Event header node: " << m_eventhead_node << std::endl;
    }
  }

  // MBD
  if ( !m_mbd_node.empty() ) 
  {
    m_tree -> Branch( "mbd_q_N", &m_mbd_q_N, "mbd_q_N/F" );
    m_tree -> Branch( "mbd_q_S", &m_mbd_q_S, "mbd_q_S/F" );
    m_tree -> Branch( "mbd_t_N", &m_mbd_t_N, "mbd_t_N/F" );
    m_tree -> Branch( "mbd_t_S", &m_mbd_t_S, "mbd_t_S/F" );
    m_tree -> Branch( "mbd_n_N", &m_mbd_n_N, "mbd_n_N/F" );
    m_tree -> Branch( "mbd_n_S", &m_mbd_n_S, "mbd_n_S/F" );
    if ( Verbosity() > 0 ) 
    {
      std::cout << "TreeWriter::Init - MBD node: " << m_mbd_node << std::endl;
    }
  }

  // cemc
  if ( !m_cemc_node.empty() ) 
  {
    m_tree -> Branch( "cemc_E", &m_cemc_E, Form("cemc_E[%d][%d]/F", k_ieta, k_iphi) );
    m_tree -> Branch( "cemc_isgood", &m_cemc_isgood, Form("cemc_isgood[%d][%d]/F", k_ieta, k_iphi) );
    m_tree -> Branch( "cemc_time", &m_cemc_time, Form("cemc_time[%d][%d]/F", k_ieta, k_iphi) );
    m_tree -> Branch( "cemc_eta", &m_cemc_eta, Form("cemc_eta[%d][%d]/F", k_ieta, k_iphi) );
    m_tree -> Branch( "cemc_phi", &m_cemc_phi, Form("cemc_phi[%d][%d]/F", k_ieta, k_iphi) );
    if ( Verbosity() > 0 ) 
    {
      std::cout << "TreeWriter::Init - Registered CEMC" << std::endl;
    } 

  }

  // hcalin
  if ( !m_hcalin_node.empty() )
  {
    m_tree -> Branch( "hcalin_E", &m_hcalin_E, Form("hcalin_E[%d][%d]/F", k_ieta, k_iphi) );
    m_tree -> Branch( "hcalin_isgood", &m_hcalin_isgood, Form("hcalin_isgood[%d][%d]/F", k_ieta, k_iphi) );
    m_tree -> Branch( "hcalin_time", &m_hcalin_time, Form("hcalin_time[%d][%d]/F", k_ieta, k_iphi) );
    m_tree -> Branch( "hcalin_eta", &m_hcalin_eta, Form("hcalin_eta[%d][%d]/F", k_ieta, k_iphi) );
    m_tree -> Branch( "hcalin_phi", &m_hcalin_phi, Form("hcalin_phi[%d][%d]/F", k_ieta, k_iphi) );
    if ( Verbosity() > 0 ) 
    {
      std::cout << "TreeWriter::Init - Registered HCALIN" << std::endl;
    } 
  }

  // hcalout
  if ( !m_hcalout_node.empty() )
  {
    m_tree -> Branch( "hcalout_E", &m_hcalout_E, Form("hcalout_E[%d][%d]/F", k_ieta, k_iphi) );
    m_tree -> Branch( "hcalout_isgood", &m_hcalout_isgood, Form("hcalout_isgood[%d][%d]/F", k_ieta, k_iphi) );
    m_tree -> Branch( "hcalout_time", &m_hcalout_time, Form("hcalout_time[%d][%d]/F", k_ieta, k_iphi) );
    m_tree -> Branch( "hcalout_eta", &m_hcalout_eta, Form("hcalout_eta[%d][%d]/F", k_ieta, k_iphi) );
    m_tree -> Branch( "hcalout_phi", &m_hcalout_phi, Form("hcalout_phi[%d][%d]/F", k_ieta, k_iphi) );
    if ( Verbosity() >  0 ) 
    {
      std::cout << "TreeWriter::Init - Registered HCALOUT" << std::endl;
    }
  }

  // cemc sub1
  if ( !m_cemc_sub1_node.empty() ) 
  {
    m_tree -> Branch( "cemc_sub1_E", &m_cemc_sub1_E, Form("cemc_sub1_E[%d][%d]/F", k_ieta, k_iphi) );
    m_tree -> Branch( "cemc_sub1_isgood", &m_cemc_sub1_isgood, Form("cemc_sub1_isgood[%d][%d]/F", k_ieta, k_iphi) );
    m_tree -> Branch( "cemc_sub1_time", &m_cemc_sub1_time, Form("cemc_sub1_time[%d][%d]/F", k_ieta, k_iphi) );
    m_tree -> Branch( "cemc_sub1_eta", &m_cemc_sub1_eta, Form("cemc_sub1_eta[%d][%d]/F", k_ieta, k_iphi) );
    m_tree -> Branch( "cemc_sub1_phi", &m_cemc_sub1_phi, Form("cemc_sub1_phi[%d][%d]/F", k_ieta, k_iphi) );
    if ( Verbosity() > 0 ) 
    {
      std::cout << "TreeWriter::Init - Registered CEMC sub1" << std::endl;
    }
  }

  // hcalin sub1
  if ( !m_hcalin_sub1_node.empty() )
  {
    m_tree -> Branch( "hcalin_sub1_E", &m_hcalin_sub1_E, Form("hcalin_sub1_E[%d][%d]/F", k_ieta, k_iphi) );
    m_tree -> Branch( "hcalin_sub1_isgood", &m_hcalin_sub1_isgood, Form("hcalin_sub1_isgood[%d][%d]/F", k_ieta, k_iphi) );
    m_tree -> Branch( "hcalin_sub1_time", &m_hcalin_sub1_time, Form("hcalin_sub1_time[%d][%d]/F", k_ieta, k_iphi) );
    m_tree -> Branch( "hcalin_sub1_eta", &m_hcalin_sub1_eta, Form("hcalin_sub1_eta[%d][%d]/F", k_ieta, k_iphi) );
    m_tree -> Branch( "hcalin_sub1_phi", &m_hcalin_sub1_phi, Form("hcalin_sub1_phi[%d][%d]/F", k_ieta, k_iphi) );
    if ( Verbosity() > 0 ) 
    {
      std::cout << "TreeWriter::Init - Registered HCALIN sub1" << std::endl;
    } 
  } 

  // hcalout sub1
  if ( !m_hcalout_sub1_node.empty() )
  {
    m_tree -> Branch( "hcalout_sub1_E", &m_hcalout_sub1_E, Form("hcalout_sub1_E[%d][%d]/F", k_ieta, k_iphi) );
    m_tree -> Branch( "hcalout_sub1_isgood", &m_hcalout_sub1_isgood, Form("hcalout_sub1_isgood[%d][%d]/F", k_ieta, k_iphi) );   
    m_tree -> Branch( "hcalout_sub1_time", &m_hcalout_sub1_time, Form("hcalout_sub1_time[%d][%d]/F", k_ieta, k_iphi) );
    m_tree -> Branch( "hcalout_sub1_eta", &m_hcalout_sub1_eta, Form("hcalout_sub1_eta[%d][%d]/F", k_ieta, k_iphi) );
    m_tree -> Branch( "hcalout_sub1_phi", &m_hcalout_sub1_phi, Form("hcalout_sub1_phi[%d][%d]/F", k_ieta, k_iphi) );
    if ( Verbosity() > 0 ) 
    {
      std::cout << "TreeWriter::Init - Registered HCALOUT sub1" << std::endl;
    }
  }

  // rho
  for ( unsigned int i = 0; i < m_rho_nodes.size(); ++i ) 
  {
    m_tree -> Branch( m_rho_nodes[i].c_str(), &m_rho_val[i], (m_rho_nodes[i] + "/F").c_str() );
    m_tree -> Branch( (m_rho_nodes[i] + "_std").c_str(), &m_std_rho_val[i], (m_rho_nodes[i] + "_std/F").c_str() );
    if ( Verbosity() > 0 ) 
    {
      std::cout << "TreeWriter::Init - Registered rho nodes: " << m_rho_nodes[i] << std::endl;
    }
  }
  


  // towrer background
  if ( m_do_towerbkgd )
  {
    m_tree -> Branch( "sub1_towerbkgd_cemc", &m_sub1_towerbkgd_cemc, Form("sub1_towerbkgd_cemc[%d]/F", k_ieta) );
    m_tree -> Branch( "sub1_towerbkgd_hcalin", &m_sub1_towerbkgd_hcalin, Form("sub1_towerbkgd_hcalin[%d]/F", k_ieta) );
    m_tree -> Branch( "sub1_towerbkgd_hcalout", &m_sub1_towerbkgd_hcalout, Form("sub1_towerbkgd_hcalout[%d]/F", k_ieta) );
    m_tree -> Branch( "sub2_towerbkgd_cemc", &m_sub2_towerbkgd_cemc, Form("sub2_towerbkgd_cemc[%d]/F", k_ieta) );
    m_tree -> Branch( "sub2_towerbkgd_hcalin", &m_sub2_towerbkgd_hcalin, Form("sub2_towerbkgd_hcalin[%d]/F", k_ieta) );
    m_tree -> Branch( "sub2_towerbkgd_hcalout", &m_sub2_towerbkgd_hcalout, Form("sub2_towerbkgd_hcalout[%d]/F", k_ieta) );
    m_tree -> Branch( "sub1_v2", &m_sub1_v2, "sub1_v2/F" );
    m_tree -> Branch( "sub2_v2", &m_sub2_v2, "sub2_v2/F" );
    m_tree -> Branch( "sub1_flowfaliure", &m_sub1_flowfaliure, "sub1_flowfaliure/I" );
    m_tree -> Branch( "sub2_flowfaliure", &m_sub2_flowfaliure, "sub2_flowfaliure/I" );
    if ( Verbosity() > 0 ) 
    {
      std::cout << "TreeWriter::Init - Registered tower background info" << std::endl;
    }
  }
  // seeds
  if ( !m_rawseed_node.empty() ) 
  {
    m_tree -> Branch( "raw_seed_E", &m_raw_seed_E );
    m_tree -> Branch( "raw_seed_phi", &m_raw_seed_phi );
    m_tree -> Branch( "raw_seed_eta", &m_raw_seed_eta );
    m_tree -> Branch( "raw_seed_pT", &m_raw_seed_pT );
    m_tree -> Branch( "raw_seed_avg_super_E", &m_mean_super_tower_E );
    m_tree -> Branch( "raw_seed_max_super_E", &m_max_super_tower_E );
    m_tree -> Branch( "raw_seed_comp_ieta", &m_raw_seed_comp_ieta );
    m_tree -> Branch( "raw_seed_comp_iphi", &m_raw_seed_comp_iphi );
    m_tree -> Branch( "raw_seed_comp_caloid", &m_raw_seed_comp_caloid );
    m_tree -> Branch( "raw_seed_comp_status", &m_raw_seed_comp_status );
    m_tree -> Branch( "raw_seed_comp_E", &m_raw_seed_comp_E );
    m_tree -> Branch( "raw_seed_comp_eta", &m_raw_seed_comp_eta );
    m_tree -> Branch( "raw_seed_comp_phi", &m_raw_seed_comp_phi );
    m_tree -> Branch( "raw_seed_super_tower_E", &m_raw_seed_super_tower_E );

    if ( Verbosity() > 0 ) 
    {
      std::cout << "TreeWriter::Init - Registered seed node: " << m_rawseed_node << std::endl;
    }
  }

  // raw jets
  for ( unsigned int i = 0; i < m_rawjet_nodes.size(); ++i ) 
  {

    m_raw_jet_trees[i] = new TTree( m_rawjet_nodes[i].c_str(), m_rawjet_nodes[i].c_str() );
    m_raw_jet_trees[i] -> Branch( "event_id", &m_event_id, "event_id/I" );
    m_raw_jet_trees[i] -> Branch( "jet_E", &m_raw_jet_E );
    m_raw_jet_trees[i] -> Branch( "jet_phi", &m_raw_jet_phi );
    m_raw_jet_trees[i] -> Branch( "jet_eta", &m_raw_jet_eta );
    m_raw_jet_trees[i] -> Branch( "jet_pT", &m_raw_jet_pT );
    m_raw_jet_trees[i] -> Branch( "jet_comp_ieta", &m_raw_jet_comp_ieta );
    m_raw_jet_trees[i] -> Branch( "jet_comp_iphi", &m_raw_jet_comp_iphi );
    m_raw_jet_trees[i] -> Branch( "jet_comp_caloid", &m_raw_jet_comp_caloid );
    m_raw_jet_trees[i] -> Branch( "jet_comp_status", &m_raw_jet_comp_status );
    m_raw_jet_trees[i] -> Branch( "jet_comp_E", &m_raw_jet_comp_E );
    m_raw_jet_trees[i] -> Branch( "jet_comp_eta", &m_raw_jet_comp_eta );
    m_raw_jet_trees[i] -> Branch( "jet_comp_phi", &m_raw_jet_comp_phi );
    if ( Verbosity() > 0 ) 
    {
      std::cout << "AnaTreeWriter::Init - Registered Jet nodes: " << m_rawjet_nodes[i] << std::endl;
    }
  }

  // sub1 jets
  for ( unsigned int i = 0; i < m_sub1jet_nodes.size(); ++i ) 
  {
    m_sub1_jet_trees[i] = new TTree( m_sub1jet_nodes[i].c_str(), m_sub1jet_nodes[i].c_str() );
    m_sub1_jet_trees[i] -> Branch( "event_id", &m_event_id, "event_id/I" );
    m_sub1_jet_trees[i] -> Branch( "jet_E", &m_sub1_jet_E );
    m_sub1_jet_trees[i] -> Branch( "jet_phi", &m_sub1_jet_phi );
    m_sub1_jet_trees[i] -> Branch( "jet_eta", &m_sub1_jet_eta );
    m_sub1_jet_trees[i] -> Branch( "jet_pT", &m_sub1_jet_pT );
    m_sub1_jet_trees[i] -> Branch( "jet_unsub_pT", &m_sub1_jet_unsub_pT );
    m_sub1_jet_trees[i] -> Branch( "jet_unsub_E", &m_sub1_jet_unsub_E );
    m_sub1_jet_trees[i] -> Branch( "jet_comp_ieta", &m_sub1_jet_comp_ieta );
    m_sub1_jet_trees[i] -> Branch( "jet_comp_iphi", &m_sub1_jet_comp_iphi );
    m_sub1_jet_trees[i] -> Branch( "jet_comp_caloid", &m_sub1_jet_comp_caloid );
    m_sub1_jet_trees[i] -> Branch( "jet_comp_status", &m_sub1_jet_comp_status );
    m_sub1_jet_trees[i] -> Branch( "jet_comp_E", &m_sub1_jet_comp_E );
    m_sub1_jet_trees[i] -> Branch( "jet_comp_eta", &m_sub1_jet_comp_eta );
    m_sub1_jet_trees[i] -> Branch( "jet_comp_phi", &m_sub1_jet_comp_phi );
    if ( Verbosity() > 0 ) 
    {
      std::cout << "AnaTreeWriter::Init - Registered Jet nodes: " << m_sub1jet_nodes[i] << std::endl;
    }
  }

  // truth jets
  for ( unsigned int i = 0; i < m_truthjet_nodes.size(); ++i ) 
  {
    m_truth_jet_trees[i] = new TTree( m_truthjet_nodes[i].c_str(), m_truthjet_nodes[i].c_str() );
    m_truth_jet_trees[i] -> Branch( "event_id", &m_event_id, "event_id/I" );
    m_truth_jet_trees[i] -> Branch( "jet_E", &m_truth_jet_E );
    m_truth_jet_trees[i] -> Branch( "jet_phi", &m_truth_jet_phi );
    m_truth_jet_trees[i] -> Branch( "jet_eta", &m_truth_jet_eta );
    m_truth_jet_trees[i] -> Branch( "jet_pT", &m_truth_jet_pT );

    if ( Verbosity() > 0 ) 
    {
      std::cout << "AnaTreeWriter::Init - Registered Truth Jet nodes: " << m_truthjet_nodes[i] << std::endl;
    }
  }

  if ( !m_g4truth_node.empty() )
  {
    m_tree -> Branch( "g4truth_zvtx", &m_g4truth_zvtx, "g4truth_zvtx/F" );
    m_tree -> Branch( "g4truth_id", &m_g4truth_id );
    m_tree -> Branch( "g4truth_status", &m_g4truth_status );
    m_tree -> Branch( "g4truth_px", &m_g4truth_px );
    m_tree -> Branch( "g4truth_py", &m_g4truth_py );
    m_tree -> Branch( "g4truth_pz", &m_g4truth_pz );
    m_tree -> Branch( "g4truth_E", &m_g4truth_E );
    m_tree -> Branch( "g4truth_eta", &m_g4truth_eta );
    m_tree -> Branch( "g4truth_phi", &m_g4truth_phi );
    m_tree -> Branch( "g4truth_pT", &m_g4truth_pT );
    if ( m_do_flow )
    {
      m_tree -> Branch( "g4truth_v2", &m_g4truth_v2 );
      m_tree -> Branch( "g4truth_v3", &m_g4truth_v3 );
      m_tree -> Branch( "g4truth_v4", &m_g4truth_v4 );
      m_tree -> Branch( "g4truth_v5", &m_g4truth_v5 );
      m_tree -> Branch( "g4truth_v6", &m_g4truth_v6 );
      m_tree -> Branch( "m_g4truth_v2reco", &m_g4truth_v2reco );
      m_tree -> Branch( "m_g4truth_v3reco", &m_g4truth_v3reco );
      m_tree -> Branch( "m_g4truth_v4reco", &m_g4truth_v4reco );
      m_tree -> Branch( "m_g4truth_v5reco", &m_g4truth_v5reco );
      m_tree -> Branch( "m_g4truth_v6reco", &m_g4truth_v6reco );
    }
  
    if ( Verbosity() > 0 ) 
    {
      std::cout << "TreeWriter::Init - Registered g4 truth node: " << m_g4truth_node << std::endl;
    }
  }

  m_rand = new TRandom3(0);
  if ( Verbosity () > 0 )
  {
    std::cout << "TreeWriter::Init - done" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TreeWriter::process_event( PHCompositeNode *topNode )
{

  m_event_id++; 

  
  if(Verbosity() > 1) 
  {
    std::cout << "TreeWriter::process_event - Process event " << m_event_id << std::endl;
  }

  if ( !m_gl1_node.empty() )
  {  // get GL1
    auto res = GetGL1(topNode);
    if ( res != Fun4AllReturnCodes::EVENT_OK ) { return res; }
  }

  if ( !m_cent_node.empty() ) 
  { // get centrality
    auto res = GetCentInfo(topNode);
    if ( res != Fun4AllReturnCodes::EVENT_OK ) { return res; }
  }

  if ( !m_zvrtx_node.empty() ) 
  { // get zvtx
    auto res = GetZvtx(topNode);
    if ( res != Fun4AllReturnCodes::EVENT_OK ) { return res; }
  }

  if ( !m_eventhead_node.empty() ) 
  { // get event header info
    auto res = GetEventHeaderInfo(topNode);
    if ( res != Fun4AllReturnCodes::EVENT_OK ) { return res; }
  }
 
  if ( !m_eventplane_node.empty() ) 
  { // get eventplane info
    auto res = GetEventPlaneInfo(topNode);
    if ( res != Fun4AllReturnCodes::EVENT_OK ) { return res; }
  }

  if ( !m_sepd_node.empty() ) 
  { // get sepd info
    auto res = GetSepdInfo(topNode);
    if ( res != Fun4AllReturnCodes::EVENT_OK ) { return res; }
  }

  
  if ( !m_mbd_node.empty() ) 
  { // get MBD
    auto res = GetMbdInfo(topNode);
    if ( res != Fun4AllReturnCodes::EVENT_OK ) { return res; }
  }

  if ( !m_cemc_node.empty() || !m_hcalin_node.empty() || !m_hcalout_node.empty() ) 
  { // get calo info
    auto res = GetRawCaloInfo(topNode);
    if ( res != Fun4AllReturnCodes::EVENT_OK ) { return res; }
  }

  if ( !m_cemc_sub1_node.empty() || !m_hcalin_sub1_node.empty() || !m_hcalout_sub1_node.empty() ) 
  { // get calo sub1 info
    auto res = GetSubCaloInfo(topNode);
    if ( res != Fun4AllReturnCodes::EVENT_OK ) { return res; }
  }

  if ( m_do_towerbkgd ) 
  { // get tower background info
    auto res = GetTowerBkgdInfo(topNode);
    if ( res != Fun4AllReturnCodes::EVENT_OK ) { return res; }
  }

  if( m_rho_nodes.size() > 0 )
  {
    auto res = GetRhoInfo(topNode);
    if ( res != Fun4AllReturnCodes::EVENT_OK ) { return res; }
  }


  if ( !m_rawseed_node.empty() ) 
  { // get seed info
    auto res = GetSeedInfo(topNode);
    if ( res != Fun4AllReturnCodes::EVENT_OK ) { return res; }
  }

  for ( unsigned int i = 0; i < m_rawjet_nodes.size(); ++i ) 
  { // reset raw jet info
    auto res = GetRawJetInfo(topNode, i);
    if ( res != Fun4AllReturnCodes::EVENT_OK ) { return res; }
  }
  

  for ( unsigned int i = 0; i < m_sub1jet_nodes.size(); ++i ) 
  { // reset sub1 jet info
    auto res = GetSub1JetInfo(topNode, i);
    if ( res != Fun4AllReturnCodes::EVENT_OK ) { return res; }
  }

  for ( unsigned int i = 0; i < m_truthjet_nodes.size(); ++i ) 
  { // reset truth jet info
    auto res = GetTruthJetInfo(topNode, i);
    if ( res != Fun4AllReturnCodes::EVENT_OK ) { return res; }
  }

  if ( !m_g4truth_node.empty() ) 
  { // get g4 truth info
    auto res = GetG4TruthInfo(topNode);
    if ( res != Fun4AllReturnCodes::EVENT_OK ) { return res; }
  }

  // fill tree
  m_tree->Fill();
  
  return Fun4AllReturnCodes::EVENT_OK;

}

int TreeWriter::End( PHCompositeNode * /*topNode*/ )
{
  
  if( Verbosity() > 0 ) 
  {
    std::cout << "TreeWriter::EndRun - End run " << std::endl;
    std::cout << "TreeWriter::EndRun - Writing to " << m_output_filename << std::endl;
  }

  PHTFileServer::get().cd(m_output_filename); 

  m_tree->Write();
  for ( unsigned int i = 0; i < m_rawjet_nodes.size(); ++i ) 
  {
    m_raw_jet_trees[i]->Write();
  }
  for ( unsigned int i = 0; i < m_sub1jet_nodes.size(); ++i ) 
  {
    m_sub1_jet_trees[i]->Write();
  }
  for (unsigned int i = 0; i < m_truthjet_nodes.size(); ++i ) 
  {
    m_truth_jet_trees[i]->Write();
  }

  m_nevents = m_event_id+1;
  m_run_tree->Fill();
  m_run_tree->Write();

  if ( Verbosity() > 0 ) 
  {
    std::cout << "TreeWriter::EndRun - Number of events: " << m_nevents << std::endl;
    std::cout << "TreeWriter::EndRun - Writing run tree" << std::endl;
  }
  PHTFileServer::get().close();
 
  if ( Verbosity () > 0 ) 
  {
    std::cout << "TreeWriter::EndRun - done" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TreeWriter::GetGL1( PHCompositeNode *topNode )
{

  ResetGL1();
  auto gl1 = findNode::getClass< Gl1Packetv2 >( topNode, m_gl1_node );
  if( !gl1 ) 
  {
    std::cout << PHWHERE << " No GL1 packet found! Abort." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // get GL1 trigger vectors
  s_triggervec = gl1->getScaledVector();
  l_triggervec = gl1->getLiveVector();

  if ( Verbosity() > 1 ) 
  {
    std::cout << PHWHERE << " - s_triggervec = " << std::hex << s_triggervec << ", l_triggervec = " << l_triggervec << std::dec << std::endl;
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int TreeWriter::GetZvtx( PHCompositeNode *topNode )
{
    // get zvtx
    ResetZvtx();
    auto vertexmap = findNode::getClass<GlobalVertexMap>( topNode, m_zvrtx_node );
    if ( !vertexmap  || vertexmap->empty() ) 
    {
      std::cout << PHWHERE << "" << m_zvrtx_node << " node missing, skipping event." << std::endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }
    // for 
    // loop over all global vertices and identify
    // for ( size_t i = 0; i < vertexmap->size(); ++i ) 
    // {
    //   auto vtxpair = vertexmap->get(i);
    //   std::cout << "Vertex ID: " << vtxpair->get_id() << " Details (x,y,z): \n" ;
    //   vtxpair->identify();
    // }
    
    auto vtx = vertexmap->begin()->second;
    m_zvtx = NAN;
    if ( vtx ) 
    {
      m_zvtx = vtx->get_z();
    } 

    if ( std::isnan(m_zvtx) || m_zvtx > 1e3) 
    {
      static bool once = true;
      if (once) 
      {
        once = false;
        std::cout << PHWHERE << "vertex is " << m_zvtx << ". Drop all tower inputs (further vertex warning will be suppressed)." << std::endl;
      }
      return Fun4AllReturnCodes::ABORTEVENT;
    }


    if ( Verbosity() > 1 ) 
    {
      std::cout << PHWHERE << " - zvtx = " << m_zvtx << std::endl;
    }

    return Fun4AllReturnCodes::EVENT_OK;
}

int TreeWriter::GetCentInfo( PHCompositeNode *topNode )
{
  // get centrality
  ResetCent();
  auto cent_node = findNode::getClass< CentralityInfo >( topNode, m_cent_node );
  if ( !cent_node ) 
  {
    std::cout << PHWHERE << m_cent_node << " node missing, Abort!." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if ( cent_node ) 
  {
    m_cent = (int)(cent_node->get_centrality_bin(CentralityInfo::PROP::mbd_NS));
  } 
  else 
  {
    std::cout << PHWHERE << "CentralityInfo is missing." << std::endl;
  }
  

  if ( Verbosity() > 1 ) 
  {
    std::cout << PHWHERE << "- Centrality = " << m_cent << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TreeWriter::GetMbdInfo( PHCompositeNode *topNode )
{
  
  // get MBD info
  ResetMbd();
  auto mbd = findNode::getClass<MbdOutV2>( topNode, m_mbd_node );
  if ( !mbd ) 
  {
    static bool once = true;
    if ( once ) 
    {
      std::cout << PHWHERE << m_mbd_node << " node missing, supressing further warnings." << std::endl;
      once = false;
    }
    return Fun4AllReturnCodes::EVENT_OK; // skipping event
  }

  m_mbd_q_N = mbd->get_q(k_mbd_north_code);
  m_mbd_q_S = mbd->get_q(k_mbd_south_code);
  m_mbd_t_N = mbd->get_time(k_mbd_north_code);
  m_mbd_t_S = mbd->get_time(k_mbd_south_code);
  m_mbd_n_N = mbd->get_npmt(k_mbd_north_code);
  m_mbd_n_S = mbd->get_npmt(k_mbd_south_code);
  

  if ( Verbosity() > 1 )
  {
    std::cout << PHWHERE << "- MBD info N:(q,t,n), S:(q,t,n) = (" << m_mbd_q_N << ", " << m_mbd_t_N << ", " << m_mbd_n_N << "), (" << m_mbd_q_S << ", " << m_mbd_t_S << ", " << m_mbd_n_S << ")" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;

}

int TreeWriter::GetEventHeaderInfo( PHCompositeNode *topNode )
{
  // get event header info
  ResetEventHeader();
  auto eventhead = findNode::getClass<EventHeader>( topNode, m_eventhead_node );
  if ( !eventhead ) 
  {
    std::cout << PHWHERE << " Input node " << m_eventhead_node << " Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_b = eventhead->get_ImpactParameter();
  m_ep_angle = eventhead->get_EventPlaneAngle();
  m_ecc = eventhead->get_eccentricity();
  m_psi1 = eventhead->get_FlowPsiN(1);
  m_psi2 = eventhead->get_FlowPsiN(2);
  m_psi3 = eventhead->get_FlowPsiN(3);
  m_psi4 = eventhead->get_FlowPsiN(4);
  m_psi5 = eventhead->get_FlowPsiN(5);
  m_psi6 = eventhead->get_FlowPsiN(6);
  m_ncoll = eventhead->get_ncoll();
  m_npart = eventhead->get_npart();
  
  if ( Verbosity() > 1 ) 
  {
    std::cout << PHWHERE << " - b = " << m_b << ", ep_angle = " << m_ep_angle << ", ecc = " << m_ecc << ", psi2 = " << m_psi2 << ", ncoll = " << m_ncoll << ", npart = " << m_npart << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TreeWriter::GetRhoInfo( PHCompositeNode *topNode )
{

  ResetRhoArrays();
  // get rho nodes
  for ( unsigned int i = 0; i < m_rho_nodes.size(); ++i ) 
  {

    auto rho = findNode::getClass<TowerRhov1>(topNode, m_rho_nodes[i]);
    if ( !rho ) 
    {
      std::cout << PHWHERE << " Input node " << m_rho_nodes[i] << " Node missing, doing nothing." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

    m_rho_val[i] = rho->get_rho();
    m_std_rho_val[i] = rho->get_sigma();
  }

  if ( Verbosity() > 1 ) 
  {
    for ( unsigned int i = 0; i < m_rho_nodes.size(); ++i ) 
    {
      std::cout << PHWHERE << " - Rho from " << m_rho_nodes[i] << " = " << m_rho_val[i] << " +/- " << m_std_rho_val[i] << std::endl;
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}



int TreeWriter::GetRawCaloInfo( PHCompositeNode *topNode )
{

  ResetCaloArrays();

  // get tower info containers
  auto towerinfosEM3 = findNode::getClass< TowerInfoContainer >( topNode, m_cemc_node);
  auto towerinfosIH3 = findNode::getClass< TowerInfoContainer >( topNode, m_hcalin_node);
  auto towerinfosOH3 = findNode::getClass< TowerInfoContainer >( topNode, m_hcalout_node);
  if( !towerinfosIH3 && !m_cemc_node.empty() )
  {
    std::cout << PHWHERE << " TowerInfoContainer for CEMC is missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN; // fatal error
  }
  if( !towerinfosIH3 && !m_hcalin_node.empty() )
  {
    std::cout << PHWHERE << " TowerInfoContainer for HCALIN is missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN; // fatal error
  }
  if( !towerinfosOH3 && !m_hcalout_node.empty() )
  {
    std::cout << PHWHERE << " TowerInfoContainer for HCALOUT is missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN; // fatal error
  }

  // get geometry containers
  auto geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  auto geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT"); 
  auto geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  if ( !geomIH && !m_hcalin_node.empty() )
  {
    std::cout << PHWHERE << " RawTowerGeomContainer for HCALIN is missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN; // fatal error
  }
  if ( !geomOH && !m_hcalout_node.empty() )
  {
    std::cout << PHWHERE << " RawTowerGeomContainer for HCALOUT is missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN; // fatal error
  }
  if ( !geomEM && !m_cemc_node.empty() )
  {
    std::cout << PHWHERE << " RawTowerGeomContainer for CEMC is missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN; // fatal error
  }

  // fill EMCal
  if ( towerinfosEM3 )
  {
   
    unsigned int ntowers = towerinfosEM3->size();
    for ( unsigned int ichannel = 0; ichannel < ntowers; ichannel++ ) 
    {
      auto tower = towerinfosEM3->get_tower_at_channel(ichannel);
      assert(tower);
      
      unsigned int key = towerinfosEM3->encode_key(ichannel);
      int ieta = towerinfosEM3->getTowerEtaBin(key);
      int iphi = towerinfosEM3->getTowerPhiBin(key);
      float this_eta = geomIH->get_etacenter(ieta);
      float this_phi = geomIH->get_phicenter(iphi);
      int this_isgood = tower->get_isGood();
      float this_E = tower->get_energy();
      float this_time = tower->get_time();

      m_cemc_E[ieta][iphi] = this_E;
      m_cemc_isgood[ieta][iphi] = this_isgood;
      m_cemc_time[ieta][iphi] = this_time;
      m_cemc_eta[ieta][iphi] = this_eta;
      m_cemc_phi[ieta][iphi] = this_phi;
    }
  }
  // fill HCalIN
  if ( towerinfosIH3 )
  {
    unsigned int ntowers = towerinfosIH3->size();
    for ( unsigned int ichannel = 0; ichannel < ntowers; ichannel++ ) 
    {
      auto tower = towerinfosIH3->get_tower_at_channel(ichannel);
      assert(tower);
      
      unsigned int key = towerinfosIH3->encode_key(ichannel);
      int ieta = towerinfosIH3->getTowerEtaBin(key);
      int iphi = towerinfosIH3->getTowerPhiBin(key);
      float this_eta = geomIH->get_etacenter(ieta);
      float this_phi = geomIH->get_phicenter(iphi);
      int this_isgood = tower->get_isGood();
      float this_E = tower->get_energy();
      float this_time = tower->get_time();

      m_hcalin_E[ieta][iphi] = this_E;
      m_hcalin_isgood[ieta][iphi] = this_isgood;
      m_hcalin_time[ieta][iphi] = this_time;
      m_hcalin_eta[ieta][iphi] = this_eta;
      m_hcalin_phi[ieta][iphi] = this_phi;
    }
  }
  // fill HCalOUT
  if ( towerinfosOH3 )
  {
    unsigned int ntowers = towerinfosOH3->size();
    for ( unsigned int ichannel = 0; ichannel < ntowers; ichannel++ ) 
    {
      auto tower = towerinfosOH3->get_tower_at_channel(ichannel);
      assert(tower);  
      unsigned int key = towerinfosOH3->encode_key(ichannel);
      int ieta = towerinfosOH3->getTowerEtaBin(key);
      int iphi = towerinfosOH3->getTowerPhiBin(key);
      float this_eta = geomOH->get_etacenter(ieta);
      float this_phi = geomOH->get_phicenter(iphi);
      int this_isgood = tower->get_isGood();
      float this_E = tower->get_energy();
      float this_time = tower->get_time();

      m_hcalout_E[ieta][iphi] = this_E;
      m_hcalout_isgood[ieta][iphi] = this_isgood;
      m_hcalout_time[ieta][iphi] = this_time;
      m_hcalout_eta[ieta][iphi] = this_eta;
      m_hcalout_phi[ieta][iphi] = this_phi;
    }
  }
  
  // printout
  if ( Verbosity() > 1 ) 
  {
    std::cout << PHWHERE << " - Filled calo info: CEMC: " << towerinfosEM3->size() << " towers, HCALIN: " << towerinfosIH3->size() << " towers, HCALOUT: " << towerinfosOH3->size() << " towers." << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;

}

int TreeWriter::GetSubCaloInfo( PHCompositeNode *topNode )
{

  ResetCaloSub1Arrays();
    
  
  // get tower info containers
  auto towerinfosEM3 = findNode::getClass< TowerInfoContainer >( topNode, m_cemc_sub1_node);
  auto towerinfosIH3 = findNode::getClass< TowerInfoContainer >( topNode, m_hcalin_sub1_node);
  auto towerinfosOH3 = findNode::getClass< TowerInfoContainer >( topNode, m_hcalout_sub1_node);
  if( !towerinfosIH3 && !m_cemc_sub1_node.empty() )
  {
    std::cout << PHWHERE << " TowerInfoContainer for CEMC sub1 is missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN; // fatal error
  }
  if( !towerinfosIH3 && !m_hcalin_sub1_node.empty() )
  {
    std::cout << PHWHERE << " TowerInfoContainer for HCALIN sub1 is missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN; // fatal error
  }
  if( !towerinfosOH3 && !m_hcalout_sub1_node.empty()  )
  {
    std::cout << PHWHERE << " TowerInfoContainer for HCALOUT sub1 is missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN; // fatal error
  }
  
  // get geometry containers
  auto geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  auto geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT"); 
  auto geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  if ( !geomIH && !m_hcalin_sub1_node.empty() )
  {
    std::cout << PHWHERE << " RawTowerGeomContainer for HCALIN is missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN; // fatal error
  }
  if ( !geomOH && !m_hcalout_sub1_node.empty() )
  {
    std::cout << PHWHERE << " RawTowerGeomContainer for HCALOUT is missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN; // fatal error
  }
  if ( !geomEM && !m_cemc_sub1_node.empty() )
  {
    std::cout << PHWHERE << " RawTowerGeomContainer for CEMC is missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN; // fatal error
  }


  // fill EMCal
  if ( towerinfosEM3 )
  {
   
    unsigned int ntowers = towerinfosEM3->size();
    for ( unsigned int ichannel = 0; ichannel < ntowers; ichannel++ ) 
    {
      auto tower = towerinfosEM3->get_tower_at_channel(ichannel);
      assert(tower);
      
      unsigned int key = towerinfosEM3->encode_key(ichannel);
      int ieta = towerinfosEM3->getTowerEtaBin(key);
      int iphi = towerinfosEM3->getTowerPhiBin(key);
      float this_eta = geomIH->get_etacenter(ieta);
      float this_phi = geomIH->get_phicenter(iphi);
      int this_isgood = tower->get_isGood();
      float this_E = tower->get_energy();
      float this_time = tower->get_time();

      m_cemc_sub1_E[ieta][iphi] = this_E;
      m_cemc_sub1_isgood[ieta][iphi] = this_isgood;
      m_cemc_sub1_time[ieta][iphi] = this_time;
      m_cemc_sub1_eta[ieta][iphi] = this_eta;
      m_cemc_sub1_phi[ieta][iphi] = this_phi;
    }
  }
  // fill HCalIN
  if ( towerinfosIH3 )
  {
    unsigned int ntowers = towerinfosIH3->size();
    for ( unsigned int ichannel = 0; ichannel < ntowers; ichannel++ ) 
    {
      auto tower = towerinfosIH3->get_tower_at_channel(ichannel);
      assert(tower);
      
      unsigned int key = towerinfosIH3->encode_key(ichannel);
      int ieta = towerinfosIH3->getTowerEtaBin(key);
      int iphi = towerinfosIH3->getTowerPhiBin(key);
      float this_eta = geomIH->get_etacenter(ieta);
      float this_phi = geomIH->get_phicenter(iphi);
      int this_isgood = tower->get_isGood();
      float this_E = tower->get_energy();
      float this_time = tower->get_time();

      m_hcalin_sub1_E[ieta][iphi] = this_E;
      m_hcalin_sub1_isgood[ieta][iphi] = this_isgood;
      m_hcalin_sub1_time[ieta][iphi] = this_time;
      m_hcalin_sub1_eta[ieta][iphi] = this_eta;
      m_hcalin_sub1_phi[ieta][iphi] = this_phi;
    }
  }
  // fill HCalOUT
  if ( towerinfosOH3 )
  {
    unsigned int ntowers = towerinfosOH3->size();
    for ( unsigned int ichannel = 0; ichannel < ntowers; ichannel++ ) 
    {
      auto tower = towerinfosOH3->get_tower_at_channel(ichannel);
      assert(tower);  
      unsigned int key = towerinfosOH3->encode_key(ichannel);
      int ieta = towerinfosOH3->getTowerEtaBin(key);
      int iphi = towerinfosOH3->getTowerPhiBin(key);
      float this_eta = geomOH->get_etacenter(ieta);
      float this_phi = geomOH->get_phicenter(iphi);
      int this_isgood = tower->get_isGood();
      float this_E = tower->get_energy();
      float this_time = tower->get_time();

      m_hcalout_sub1_E[ieta][iphi] = this_E;
      m_hcalout_sub1_isgood[ieta][iphi] = this_isgood;
      m_hcalout_sub1_time[ieta][iphi] = this_time;
      m_hcalout_sub1_eta[ieta][iphi] = this_eta;
      m_hcalout_sub1_phi[ieta][iphi] = this_phi;
    }
  }
  
  // printout
  if ( Verbosity() > 1 ) 
  {
    std::cout << PHWHERE << " - Filled calo info: CEMC: " << towerinfosEM3->size() << " towers, HCALIN: " << towerinfosIH3->size() << " towers, HCALOUT: " << towerinfosOH3->size() << " towers." << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;

}

// int TreeWriter::GetTowerBkgdInfo( PHCompositeNode *topNode )
// {

//   ResetTowerBkgd();
  
//   auto tower_background_sub1 = findNode::getClass<TowerBackgroundv1>(topNode, m_towerbkgd_node_sub1);
//   if ( !tower_background_sub1 )
//   {
//     std::cout << PHWHERE << " TowerBackgroundv1 node is missing, skipping." << std::endl;
//     return Fun4AllReturnCodes::ABORTRUN; // fatal error
//   }
//   auto tower_background_sub2 = findNode::getClass<TowerBackgroundv1>(topNode, m_towerbkgd_node_sub2);
//   if ( !tower_background_sub2 )
//   {
//     std::cout << PHWHERE << " TowerBackgroundv1 node is missing, skipping." << std::endl;
//     return Fun4AllReturnCodes::ABORTRUN; // fatal error
//   }

//   for ( int ilayer = 0 ; ilayer < 3; ilayer++ )
//   {
//     std::vector<float> this_ue_sub1 = tower_background_sub1->get_UE(ilayer);
//     std::vector<float> this_ue_sub2 = tower_background_sub2->get_UE(ilayer);
//     for ( int ieta = 0; ieta < 24; ieta++ )
//     {
//       if ( ilayer == 0 ) 
//       {
//         m_sub1_towerbkgd_cemc[ieta] = this_ue_sub1[ieta];
//         m_sub2_towerbkgd_cemc[ieta] = this_ue_sub2[ieta];
//       }
//       else if ( ilayer == 1 ) 
//       {
//         m_sub1_towerbkgd_hcalin[ieta] = this_ue_sub1[ieta];
//         m_sub2_towerbkgd_hcalin[ieta] = this_ue_sub2[ieta];
//       }
//       else if ( ilayer == 2 ) 
//       {
//         m_sub1_towerbkgd_hcalout[ieta] = this_ue_sub1[ieta];
//         m_sub2_towerbkgd_hcalout[ieta] = this_ue_sub2[ieta];
//       }
//     }
//   } // end loop over layers

//   if ( Verbosity() > 1 ) 
//   {
//     std::cout << PHWHERE << " - Tower background from " << m_towerbkgd_node_sub1 << " and " << m_towerbkgd_node_sub2 << std::endl;
//   }


//   return Fun4AllReturnCodes::EVENT_OK;

// }

int TreeWriter::GetTowerBkgdInfo( PHCompositeNode *topNode )
{

  ResetTowerBkgd();
  
  auto tower_background_sub1 = findNode::getClass<TowerBackgroundv1>(topNode, m_towerbkgd_node_sub1);
  if ( !tower_background_sub1 )
  {
    std::cout << PHWHERE << " TowerBackgroundv1 node is missing, skipping." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN; // fatal error
  }
  auto tower_background_sub2 = findNode::getClass<TowerBackgroundv1>(topNode, m_towerbkgd_node_sub2);
  if ( !tower_background_sub2 )
  {
    std::cout << PHWHERE << " TowerBackgroundv1 node is missing, skipping." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN; // fatal error
  }

  m_sub1_v2 = tower_background_sub1->get_v2();
  m_sub2_v2 = tower_background_sub2->get_v2();
  m_sub1_flowfaliure = tower_background_sub1->get_flow_failure_flag() == true ? 1 : 0;
  m_sub2_flowfaliure = tower_background_sub2->get_flow_failure_flag() == true ? 1 : 0;
  for ( int ilayer = 0 ; ilayer < 3; ilayer++ )
  {
    std::vector<float> this_ue_sub1 = tower_background_sub1->get_UE(ilayer);
    std::vector<float> this_ue_sub2 = tower_background_sub2->get_UE(ilayer);
    for ( int ieta = 0; ieta < 24; ieta++ )
    {
      if ( ilayer == 0 ) 
      {
        m_sub1_towerbkgd_cemc[ieta] = this_ue_sub1[ieta];
        m_sub2_towerbkgd_cemc[ieta] = this_ue_sub2[ieta];
      }
      else if ( ilayer == 1 ) 
      {
        m_sub1_towerbkgd_hcalin[ieta] = this_ue_sub1[ieta];
        m_sub2_towerbkgd_hcalin[ieta] = this_ue_sub2[ieta];
      }
      else if ( ilayer == 2 ) 
      {
        m_sub1_towerbkgd_hcalout[ieta] = this_ue_sub1[ieta];
        m_sub2_towerbkgd_hcalout[ieta] = this_ue_sub2[ieta];
      }
    }
  } // end loop over layers

  if ( Verbosity() > 1 ) 
  {
    std::cout << PHWHERE << " - Tower background from " << m_towerbkgd_node_sub1 << " and " << m_towerbkgd_node_sub2 << std::endl;
    std::cout << "   sub1 v2: " << m_sub1_v2 << ", flow failure: " << m_sub1_flowfaliure << std::endl;
  }


  return Fun4AllReturnCodes::EVENT_OK;

}

int TreeWriter::GetSeedInfo( PHCompositeNode *topNode )
{


  ResetRawSeed();

  // get tower info containers
  auto towerinfosEM3 = findNode::getClass< TowerInfoContainer >( topNode, m_cemc_node);
  auto towerinfosIH3 = findNode::getClass< TowerInfoContainer >( topNode, m_hcalin_node);
  auto towerinfosOH3 = findNode::getClass< TowerInfoContainer >( topNode, m_hcalout_node);
  if( !towerinfosIH3 || !towerinfosOH3 || !towerinfosEM3 )
  {
    std::cout
      << PHWHERE
      << " One of the following nodes is missing: "
      << m_cemc_node << ", "
      << m_hcalin_node << ", "
      << m_hcalout_node << "."
      << " Skipping calo filling."
      << std::endl;
    return Fun4AllReturnCodes::ABORTRUN; // fatal error
  }
  
  // get geometry containers
  auto geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  auto geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT"); 
  if ( !geomIH || !geomOH ) 
  {
    std::cout << PHWHERE << " RawTowerGeomContainer for HCALIN or HCALOUT is missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN; // fatal error
  }
  
  auto seeds = findNode::getClass<JetContainer>( topNode, m_rawseed_node );
  if ( !seeds )
  {
    std::cout << PHWHERE << " Input node " << m_rawseed_node << " Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT; 
  }


  for ( auto seed : *seeds )
  {

    float this_pt = seed->get_pt();
    float this_eta = seed->get_eta();
    float this_phi = seed->get_phi();
    float this_e = seed->get_e();

    std::map < int, double > constituent_ETsum {} ;
    std::vector<int> tower_ieta {};
    std::vector<int> tower_iphi {};
    std::vector<int> tower_caloid {};
    std::vector<float> tower_E {};
    std::vector<float> tower_eta {};
    std::vector<float> tower_phi {};
    std::vector<int> tower_status {};
    std::vector<float> super_tower_E {};

    for ( const auto &comp : seed->get_comp_vec() )
    {
      int this_comp_ieta = -999;
      int this_comp_iphi = -999;
      float this_comp_E = 0;
      float this_comp_eT = 0;
      float this_comp_eta = 0;
      float this_comp_phi = 0;
      int this_comp_status = -1;
      int this_comp_caloid = comp.first;

      bool is_hcalin = ( comp.first == 5 || comp.first == 26 );
      bool is_hcalout = ( comp.first == 7 || comp.first == 27 );
      bool is_cemc = ( comp.first == 13 || comp.first == 28 );
      if ( is_cemc )
      {
        auto tower = towerinfosEM3->get_tower_at_channel( comp.second );
        assert(tower);
        unsigned int towerkey = towerinfosEM3->encode_key( comp.second );
        this_comp_ieta = towerinfosEM3->getTowerEtaBin(towerkey);
        this_comp_iphi = towerinfosEM3->getTowerPhiBin(towerkey);
        const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid( RawTowerDefs::CalorimeterId::HCALIN,  this_comp_ieta,  this_comp_iphi);
        auto geom = geomIH->get_tower_geometry(key);
        if ( geom )
        {
          this_comp_eta = geom->get_eta();
          this_comp_phi = geom->get_phi();
        }
        this_comp_E = tower->get_energy();
        this_comp_eT = this_comp_E / cosh(this_comp_eta);
        this_comp_status = tower->get_isGood();
      }
      else if ( is_hcalin )
      {
        auto tower = towerinfosIH3->get_tower_at_channel( comp.second );
        assert(tower);
        unsigned int towerkey = towerinfosIH3->encode_key( comp.second );
        this_comp_ieta = towerinfosIH3->getTowerEtaBin(towerkey);
        this_comp_iphi = towerinfosIH3->getTowerPhiBin(towerkey);
        const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid( RawTowerDefs::CalorimeterId::HCALIN,  this_comp_ieta,  this_comp_iphi);
        auto geom = geomIH->get_tower_geometry(key);
        if ( geom )
        {
          this_comp_eta = geom->get_eta();
          this_comp_phi = geom->get_phi();
        }
        this_comp_E = tower->get_energy();
        this_comp_eT = this_comp_E / cosh(this_comp_eta);
        this_comp_status = tower->get_isGood();
      }
      else if ( is_hcalout )
      {
        auto tower = towerinfosOH3->get_tower_at_channel( comp.second );
        assert(tower);
        unsigned int towerkey = towerinfosOH3->encode_key( comp.second );
        this_comp_ieta = towerinfosOH3->getTowerEtaBin(towerkey);
        this_comp_iphi = towerinfosOH3->getTowerPhiBin(towerkey);
        const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid( RawTowerDefs::CalorimeterId::HCALOUT,  this_comp_ieta,  this_comp_iphi);
        auto geom = geomOH->get_tower_geometry(key);
        if ( geom )
        {
          this_comp_eta = geom->get_eta();
          this_comp_phi = geom->get_phi();
        }
        this_comp_E = tower->get_energy();
        this_comp_eT = this_comp_E / cosh(this_comp_eta);
        this_comp_status = tower->get_isGood();
      }
      else
      {
        std::cout << PHWHERE << " Warning: seed constituent caloid " << comp.first << " not recognized, skipping." << std::endl;
        continue;
      }

      if ( this_comp_status != 1 ) 
      {
        continue; // skip bad towers
      }


      int comp_ikey = (1000 * this_comp_ieta) + this_comp_iphi;
      constituent_ETsum[comp_ikey] += this_comp_eT;


      tower_ieta.push_back(this_comp_ieta);
      tower_iphi.push_back(this_comp_iphi);
      tower_caloid.push_back(this_comp_caloid);
      tower_E.push_back(this_comp_E);
      tower_eta.push_back(this_comp_eta);
      tower_phi.push_back(this_comp_phi);
      tower_status.push_back(this_comp_status);
      
    } // end loop over constituents

    float constituent_max_ET = 0;
    float constituent_sum_ET = 0;
    int nconstituents = 0;
    for (auto &map_iter : constituent_ETsum)
    {
      nconstituents++;
      constituent_sum_ET += map_iter.second;
      constituent_max_ET = std::max<double>(map_iter.second, constituent_max_ET);
      super_tower_E.push_back(map_iter.second);
    }
    
    float mean_constituent_ET = 0;
    if ( nconstituents != 0 ) 
    {
      mean_constituent_ET = constituent_sum_ET / nconstituents;
    }

    
  
    m_raw_seed_E.push_back(this_e);
    m_raw_seed_eta.push_back(this_eta);
    m_raw_seed_phi.push_back(this_phi);
    m_raw_seed_pT.push_back(this_pt);
    m_raw_seed_comp_ieta.push_back(tower_ieta);
    m_raw_seed_comp_iphi.push_back(tower_iphi);
    m_raw_seed_comp_caloid.push_back(tower_caloid);
    m_raw_seed_comp_status.push_back(tower_status);
    m_raw_seed_comp_E.push_back(tower_E);
    m_raw_seed_comp_eta.push_back(tower_eta);
    m_raw_seed_comp_phi.push_back(tower_phi);
    m_raw_seed_super_tower_E.push_back(super_tower_E);
    m_mean_super_tower_E.push_back(mean_constituent_ET);
    m_max_super_tower_E.push_back(constituent_max_ET);

  } // end loop over seeds

  if ( Verbosity() > 0 ) 
  {
    std::cout << PHWHERE << " - Found " << seeds->size() << " seeds in node " << m_rawseed_node << std::endl;
  }


  return Fun4AllReturnCodes::EVENT_OK;

}

int TreeWriter::GetRawJetInfo( PHCompositeNode *topNode , const int idx)
{

  const std::string & node_name = m_rawjet_nodes[idx];
 
  ResetRawJet();

  // get tower info containers
  auto towerinfosEM3 = findNode::getClass< TowerInfoContainer >( topNode, "TOWERINFO_CALIB_CEMC_RETOWER");
  auto towerinfosIH3 = findNode::getClass< TowerInfoContainer >( topNode, "TOWERINFO_CALIB_HCALIN" );
  auto towerinfosOH3 = findNode::getClass< TowerInfoContainer >( topNode,  "TOWERINFO_CALIB_HCALOUT" );
  if( !towerinfosIH3 || !towerinfosOH3 || !towerinfosEM3 )
  {
    std::cout
      << PHWHERE
      << " One of the following nodes is missing: "
      << "TOWERINFO_CALIB_CEMC_RETOWER, "
      << "TOWERINFO_CALIB_HCALIN, "
      << "TOWERINFO_CALIB_HCALOUT."
      << " Skipping calo filling."
      << std::endl;
    return Fun4AllReturnCodes::ABORTRUN; // fatal error
  }
  
  // get geometry containers
  auto geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  auto geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT"); 
  if ( !geomIH || !geomOH ) 
  {
    std::cout << PHWHERE << " RawTowerGeomContainer for HCALIN or HCALOUT is missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN; // fatal error
  }
  
  // get raw jets
  auto jets = findNode::getClass<JetContainer>( topNode, node_name );
  if ( !jets )
  {
    std::cout << PHWHERE << " Input node " << node_name << " Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT; 
  }

  for ( auto jet : *jets )
  {

    float this_pt = jet->get_pt();
    float this_eta = jet->get_eta();
    float this_phi = jet->get_phi();
    float this_e = jet->get_e();

    std::vector<int> tower_ieta {};
    std::vector<int> tower_iphi {};
    std::vector<int> tower_caloid {};
    std::vector<float> tower_E {};
    std::vector<float> tower_eta {};
    std::vector<float> tower_phi {};
    std::vector<int> tower_status {};

    for ( const auto &comp : jet->get_comp_vec() )
    {
      int this_comp_ieta = -999;
      int this_comp_iphi = -999;
      float this_comp_E = 0;
      float this_comp_eta = 0;
      float this_comp_phi = 0;
      int this_comp_status = -1;
      int this_comp_caloid = comp.first;

      bool is_hcalin = ( comp.first == Jet::SRC::HCALIN_TOWERINFO || comp.first == Jet::SRC::HCALIN_TOWERINFO_SUB1 );
      bool is_hcalout = ( comp.first == Jet::SRC::HCALOUT_TOWERINFO || comp.first == Jet::SRC::HCALOUT_TOWERINFO_SUB1 );
      bool is_cemc = ( comp.first == Jet::SRC::CEMC_TOWERINFO || comp.first == Jet::SRC::CEMC_TOWERINFO_RETOWER || comp.first == Jet::SRC::CEMC_TOWERINFO_SUB1);
      if ( is_cemc )
      {
        auto tower = towerinfosEM3->get_tower_at_channel( comp.second );
        assert(tower);
        unsigned int towerkey = towerinfosEM3->encode_key( comp.second );
        this_comp_ieta = towerinfosEM3->getTowerEtaBin(towerkey);
        this_comp_iphi = towerinfosEM3->getTowerPhiBin(towerkey);
        const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid( RawTowerDefs::CalorimeterId::HCALIN,  this_comp_ieta,  this_comp_iphi);
        auto geom = geomIH->get_tower_geometry(key);
        if ( geom )
        {
          this_comp_eta = geom->get_eta();
          this_comp_phi = geom->get_phi();
        }
        this_comp_E = tower->get_energy();
        this_comp_status = tower->get_isGood();
      }
      else if ( is_hcalin )
      {
        auto tower = towerinfosIH3->get_tower_at_channel( comp.second );
        assert(tower);
        unsigned int towerkey = towerinfosIH3->encode_key( comp.second );
        this_comp_ieta = towerinfosIH3->getTowerEtaBin(towerkey);
        this_comp_iphi = towerinfosIH3->getTowerPhiBin(towerkey);
        const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid( RawTowerDefs::CalorimeterId::HCALIN,  this_comp_ieta,  this_comp_iphi);
        auto geom = geomIH->get_tower_geometry(key);
        if ( geom )
        {
          this_comp_eta = geom->get_eta();
          this_comp_phi = geom->get_phi();
        }
        this_comp_E = tower->get_energy();
        this_comp_status = tower->get_isGood();
      }
      else if ( is_hcalout )
      {
        auto tower = towerinfosOH3->get_tower_at_channel( comp.second );
        assert(tower);
        unsigned int towerkey = towerinfosOH3->encode_key( comp.second );
        this_comp_ieta = towerinfosOH3->getTowerEtaBin(towerkey);
        this_comp_iphi = towerinfosOH3->getTowerPhiBin(towerkey);
        const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid( RawTowerDefs::CalorimeterId::HCALOUT,  this_comp_ieta,  this_comp_iphi);
        auto geom = geomOH->get_tower_geometry(key);
        if ( geom )
        {
          this_comp_eta = geom->get_eta();
          this_comp_phi = geom->get_phi();
        }
        this_comp_E = tower->get_energy();
        this_comp_status = tower->get_isGood();
      }
      else
      {
        std::cout << PHWHERE << " Warning: jets constituent caloid " << comp.first << " not recognized, skipping." << std::endl;
        continue;
      }

      if ( this_comp_status != 1 ) 
      {
        continue; // skip bad towers
      }


      tower_ieta.push_back(this_comp_ieta);
      tower_iphi.push_back(this_comp_iphi);
      tower_caloid.push_back(this_comp_caloid);
      tower_E.push_back(this_comp_E);
      tower_eta.push_back(this_comp_eta);
      tower_phi.push_back(this_comp_phi);
      tower_status.push_back(this_comp_status);
      
    } // end loop over constituents


    m_raw_jet_E.push_back(this_e);
    m_raw_jet_eta.push_back(this_eta);
    m_raw_jet_phi.push_back(this_phi);
    m_raw_jet_pT.push_back(this_pt);
    m_raw_jet_comp_ieta.push_back(tower_ieta);
    m_raw_jet_comp_iphi.push_back(tower_iphi);
    m_raw_jet_comp_caloid.push_back(tower_caloid);
    m_raw_jet_comp_status.push_back(tower_status);
    m_raw_jet_comp_E.push_back(tower_E);
    m_raw_jet_comp_eta.push_back(tower_eta);
    m_raw_jet_comp_phi.push_back(tower_phi);

  } // end loop over jets

  if ( Verbosity() > 0 ) 
  {
    std::cout << PHWHERE << " - Found " << jets->size() << " raw jets in node " << node_name << std::endl;
  }

  m_raw_jet_trees[idx]->Fill();

  
  return Fun4AllReturnCodes::EVENT_OK;

}

int TreeWriter::GetSub1JetInfo( PHCompositeNode *topNode , const int idx)
{

  const std::string & node_name = m_sub1jet_nodes[idx];
 
  ResetSub1Jet();

  // get tower info containers
  auto towerinfosEM3 = findNode::getClass< TowerInfoContainer >( topNode, "TOWERINFO_CALIB_CEMC_RETOWER_SUB1");
  auto towerinfosIH3 = findNode::getClass< TowerInfoContainer >( topNode, "TOWERINFO_CALIB_HCALIN_SUB1" );
  auto towerinfosOH3 = findNode::getClass< TowerInfoContainer >( topNode, "TOWERINFO_CALIB_HCALOUT_SUB1" );
  if( !towerinfosIH3 || !towerinfosOH3 || !towerinfosEM3 )
  {
    std::cout
      << PHWHERE
      << " One of the following nodes is missing: "
      << "TOWERINFO_SUB1_CEMC_RETOWER, "
      << "TOWERINFO_SUB1_HCALIN, "
      << "TOWERINFO_SUB1_HCALOUT."
      << " Skipping calo filling."
      << std::endl;
    return Fun4AllReturnCodes::ABORTRUN; // fatal error
  }
  
  // get geometry containers
  auto geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  auto geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT"); 
  if ( !geomIH || !geomOH ) 
  {
    std::cout << PHWHERE << " RawTowerGeomContainer for HCALIN or HCALOUT is missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN; // fatal error
  }
  
  // get raw jets
  auto jets = findNode::getClass<JetContainer>( topNode, node_name );
  if ( !jets )
  {
    std::cout << PHWHERE << " Input node " << node_name << " Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT; 
  }

  auto tower_background_sub1 = findNode::getClass<TowerBackgroundv1>(topNode, "TowerInfoBackground_Sub1");
  if ( !tower_background_sub1 )
  {
    std::cout << PHWHERE << " TowerBackgroundv1 node is missing, skipping." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN; // fatal error
  }

  auto tower_background_sub2 = findNode::getClass<TowerBackgroundv1>(topNode, "TowerInfoBackground_Sub2");
  if ( !tower_background_sub2 )
  {
    std::cout << PHWHERE << " TowerBackgroundv1 node is missing, skipping." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN; // fatal error
  }

  for ( int ilayer = 0 ; ilayer < 3; ilayer++ )
  {
    std::vector<float> this_ue_sub1 = tower_background_sub1->get_UE(ilayer);
    std::vector<float> this_ue_sub2 = tower_background_sub2->get_UE(ilayer);
    for ( int ieta = 0; ieta < 24; ieta++ )
    {
      if ( ilayer == 0 ) 
      {
        m_sub1_towerbkgd_cemc[ieta] = this_ue_sub1[ieta];
        m_sub2_towerbkgd_cemc[ieta] = this_ue_sub2[ieta];
      }
      else if ( ilayer == 1 ) 
      {
        m_sub1_towerbkgd_hcalin[ieta] = this_ue_sub1[ieta];
        m_sub2_towerbkgd_hcalin[ieta] = this_ue_sub2[ieta];
      }
      else if ( ilayer == 2 ) 
      {
        m_sub1_towerbkgd_hcalout[ieta] = this_ue_sub1[ieta];
        m_sub2_towerbkgd_hcalout[ieta] = this_ue_sub2[ieta];
      }
    }
  } // end loop over layers


  for ( auto jet : *jets )
  {

    float this_pt = jet->get_pt();
    float this_eta = jet->get_eta();
    float this_phi = jet->get_phi();
    float this_e = jet->get_e();
    float unsub_px = 0, unsub_py = 0;
    float unsub_E = 0;


    std::vector<int> tower_ieta {};
    std::vector<int> tower_iphi {};
    std::vector<int> tower_caloid {};
    std::vector<float> tower_E {};
    std::vector<float> tower_eta {};
    std::vector<float> tower_phi {};
    std::vector<int> tower_status {};

    for ( const auto &comp : jet->get_comp_vec() )
    {
      int this_comp_ieta = -999;
      int this_comp_iphi = -999;
      float this_comp_E = 0;
      float this_comp_eta = 0;
      float this_comp_phi = 0;
      int this_comp_status = -1;
      int this_comp_caloid = comp.first;

      bool is_hcalin = ( comp.first == Jet::SRC::HCALIN_TOWERINFO || comp.first == Jet::SRC::HCALIN_TOWERINFO_SUB1 );
      bool is_hcalout = ( comp.first == Jet::SRC::HCALOUT_TOWERINFO || comp.first == Jet::SRC::HCALOUT_TOWERINFO_SUB1 );
      bool is_cemc = ( comp.first == Jet::SRC::CEMC_TOWERINFO || comp.first == Jet::SRC::CEMC_TOWERINFO_RETOWER || comp.first == Jet::SRC::CEMC_TOWERINFO_SUB1);
      if ( is_cemc )
      {
        auto tower = towerinfosEM3->get_tower_at_channel( comp.second );
        assert(tower);
        unsigned int towerkey = towerinfosEM3->encode_key( comp.second );
        this_comp_ieta = towerinfosEM3->getTowerEtaBin(towerkey);
        this_comp_iphi = towerinfosEM3->getTowerPhiBin(towerkey);
        const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid( RawTowerDefs::CalorimeterId::HCALIN,  this_comp_ieta,  this_comp_iphi);
        auto geom = geomIH->get_tower_geometry(key);
        if ( geom )
        {
          this_comp_eta = geom->get_eta();
          this_comp_phi = geom->get_phi();
        }
        this_comp_E = tower->get_energy();
        this_comp_status = tower->get_isGood();

        float this_ue = m_sub2_towerbkgd_cemc[this_comp_ieta];
        float this_unsub_pt = ( this_comp_E + this_ue ) / cosh(this_comp_eta);
        unsub_px += this_unsub_pt * cos(this_comp_phi);
        unsub_py += this_unsub_pt * sin(this_comp_phi);
        unsub_E += ( this_comp_E + this_ue );
      }
      else if ( is_hcalin )
      {
        auto tower = towerinfosIH3->get_tower_at_channel( comp.second );
        assert(tower);
        unsigned int towerkey = towerinfosIH3->encode_key( comp.second );
        this_comp_ieta = towerinfosIH3->getTowerEtaBin(towerkey);
        this_comp_iphi = towerinfosIH3->getTowerPhiBin(towerkey);
        const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid( RawTowerDefs::CalorimeterId::HCALIN,  this_comp_ieta,  this_comp_iphi);
        auto geom = geomIH->get_tower_geometry(key);
        if ( geom )
        {
          this_comp_eta = geom->get_eta();
          this_comp_phi = geom->get_phi();
        }
        this_comp_E = tower->get_energy();
        this_comp_status = tower->get_isGood();

        float this_ue = m_sub2_towerbkgd_hcalin[this_comp_ieta];
        float this_unsub_pt = ( this_comp_E + this_ue ) / cosh(this_comp_eta);
        unsub_px += this_unsub_pt * cos(this_comp_phi);
        unsub_py += this_unsub_pt * sin(this_comp_phi);
        unsub_E += ( this_comp_E + this_ue );
      }
      else if ( is_hcalout )
      {
        auto tower = towerinfosOH3->get_tower_at_channel( comp.second );
        assert(tower);
        unsigned int towerkey = towerinfosOH3->encode_key( comp.second );
        this_comp_ieta = towerinfosOH3->getTowerEtaBin(towerkey);
        this_comp_iphi = towerinfosOH3->getTowerPhiBin(towerkey);
        const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid( RawTowerDefs::CalorimeterId::HCALOUT,  this_comp_ieta,  this_comp_iphi);
        auto geom = geomOH->get_tower_geometry(key);
        if ( geom )
        {
          this_comp_eta = geom->get_eta();
          this_comp_phi = geom->get_phi();
        }
        this_comp_E = tower->get_energy();
        this_comp_status = tower->get_isGood();
        float this_ue = m_sub2_towerbkgd_hcalout[this_comp_ieta];
        float this_unsub_pt = ( this_comp_E + this_ue ) / cosh(this_comp_eta);
        unsub_px += this_unsub_pt * cos(this_comp_phi);
        unsub_py += this_unsub_pt * sin(this_comp_phi);
        unsub_E += ( this_comp_E + this_ue );
      }
      else
      {
        std::cout << PHWHERE << " Warning: jets constituent caloid " << comp.first << " not recognized, skipping." << std::endl;
        continue;
      }

      if ( this_comp_status != 1 ) 
      {
        continue; // skip bad towers
      }


      tower_ieta.push_back(this_comp_ieta);
      tower_iphi.push_back(this_comp_iphi);
      tower_caloid.push_back(this_comp_caloid);
      tower_E.push_back(this_comp_E);
      tower_eta.push_back(this_comp_eta);
      tower_phi.push_back(this_comp_phi);
      tower_status.push_back(this_comp_status);
      
    } // end loop over constituents


    float unsub_pt = sqrt( (unsub_px * unsub_px) + (unsub_py * unsub_py) );
    m_sub1_jet_E.push_back(this_e);
    m_sub1_jet_eta.push_back(this_eta);
    m_sub1_jet_phi.push_back(this_phi);
    m_sub1_jet_pT.push_back(this_pt);
    m_sub1_jet_unsub_pT.push_back(unsub_pt);
    m_sub1_jet_unsub_E.push_back(unsub_E);
    m_sub1_jet_comp_ieta.push_back(tower_ieta);
    m_sub1_jet_comp_iphi.push_back(tower_iphi);
    m_sub1_jet_comp_caloid.push_back(tower_caloid);
    m_sub1_jet_comp_status.push_back(tower_status);
    m_sub1_jet_comp_E.push_back(tower_E);
    m_sub1_jet_comp_eta.push_back(tower_eta);
    m_sub1_jet_comp_phi.push_back(tower_phi); 


  } // end loop over jets

  if ( Verbosity() > 0 ) 
  {
    std::cout << PHWHERE << " - Found " << jets->size() << " sub1 jets in node " << node_name << std::endl;
  }

  m_sub1_jet_trees[idx]->Fill();

  
  return Fun4AllReturnCodes::EVENT_OK;

}

int TreeWriter::GetTruthJetInfo( PHCompositeNode *topNode , const int idx)
{

  const std::string & node_name = m_truthjet_nodes[idx];
 
  ResetTruthJet();

  // get truth jets
  auto jets = findNode::getClass<JetContainer>( topNode, node_name );
  if ( !jets )
  {
    std::cout << PHWHERE << " Input node " << node_name << " Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT; 
  }

  for ( auto jet : *jets )
  {

    float this_pt = jet->get_pt();
    float this_eta = jet->get_eta();
    float this_phi = jet->get_phi();
    float this_e = jet->get_e();

    m_truth_jet_E.push_back(this_e);
    m_truth_jet_eta.push_back(this_eta);
    m_truth_jet_phi.push_back(this_phi);
    m_truth_jet_pT.push_back(this_pt);

  } // end loop over jets

  if ( Verbosity() > 0 ) 
  {
    std::cout << PHWHERE << " - Found " << jets->size() << " truth jets in node " << node_name << std::endl;
  }

  m_truth_jet_trees[idx]->Fill();

  
  return Fun4AllReturnCodes::EVENT_OK;

}

int TreeWriter::GetG4TruthInfo( PHCompositeNode *topNode )
{
  ResetG4Truth();

  // get g4 truth info
  PHG4TruthInfoContainer * truthinfo = findNode::getClass<PHG4TruthInfoContainer>( topNode, m_g4truth_node );
  if ( !truthinfo ) 
  {
    std::cout << PHWHERE << " Input node "<< m_g4truth_node << " Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  float _ptsum_v2 = 0;
  float _ptsum_v3 = 0;
  float _ptsum_v4 = 0;
  float _ptsum_v5 = 0;
  float _ptsum_v6 = 0;
  float _ptsum = 0;
  PHG4TruthInfoContainer::ConstRange range = truthinfo->GetSPHENIXPrimaryParticleRange();
  PHG4VtxPoint *gvertex = truthinfo->GetPrimaryVtx(truthinfo->GetPrimaryVertexIndex());
  m_g4truth_zvtx = gvertex->get_z();
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter) 
  {
    auto g4particle = iter->second;
    
    TLorentzVector t;
    t.SetPxPyPzE(g4particle->get_px(), g4particle->get_py(), g4particle->get_pz(), g4particle->get_e());

    float truth_pt = t.Pt();
    float truth_eta = t.Eta();
    float truth_phi = t.Phi();
    m_g4truth_eta.push_back(truth_eta);
    m_g4truth_phi.push_back(truth_phi);
    m_g4truth_pT.push_back(truth_pt);
    m_g4truth_px.push_back(g4particle->get_px());
    m_g4truth_py.push_back(g4particle->get_py());
    m_g4truth_pz.push_back(g4particle->get_pz());
    m_g4truth_E.push_back(g4particle->get_e());
    m_g4truth_id.push_back(g4particle->get_pid());
    m_g4truth_status.push_back(truthinfo->isEmbeded(g4particle->get_track_id()));

    if (m_do_flow && !m_eventhead_node.empty() )
    {
      float _v2, _v3, _v4, _v5, _v6;
      CalcFlow( m_b, truth_eta, truth_phi, _v2, _v3, _v4, _v5, _v6, m_rand , m_do_fluc, m_flow_scale);
      m_g4truth_v2.push_back(_v2);
      m_g4truth_v3.push_back(_v3);
      m_g4truth_v4.push_back(_v4);
      m_g4truth_v5.push_back(_v5);
      m_g4truth_v6.push_back(_v6);
      _ptsum_v2 += truth_pt * TMath::Cos(2.0 * (truth_phi - m_ep_angle));
      _ptsum_v3 += truth_pt * TMath::Cos(3.0 * (truth_phi - m_psi3 ));
      _ptsum_v4 += truth_pt * TMath::Cos(4.0 * (truth_phi - m_psi4 ));
      _ptsum_v5 += truth_pt * TMath::Cos(5.0 * (truth_phi - m_psi5 ));
      _ptsum_v6 += truth_pt * TMath::Cos(6.0 * (truth_phi - m_psi6 ));
      _ptsum += truth_pt;
    }
  }

  if (m_do_flow && _ptsum > 0)
  {
    m_g4truth_v2reco = _ptsum_v2 / _ptsum;
    m_g4truth_v3reco = _ptsum_v3 / _ptsum;
    m_g4truth_v4reco = _ptsum_v4 / _ptsum;
    m_g4truth_v5reco = _ptsum_v5 / _ptsum;
    m_g4truth_v6reco = _ptsum_v6 / _ptsum;
  }
  else if (m_do_flow)
  {
    m_g4truth_v2reco = 0;
    m_g4truth_v3reco = 0;
    m_g4truth_v4reco = 0;
    m_g4truth_v5reco = 0;
    m_g4truth_v6reco = 0;
  }
  

  if ( Verbosity() > 0 ) 
  {
    std::cout << PHWHERE << " - Found " << m_g4truth_pT.size() << " truth particles in node " << m_g4truth_node << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TreeWriter::GetSepdInfo( PHCompositeNode *topNode )
{
  // get sEPD info
  auto sepd = findNode::getClass<TowerInfoContainer>(topNode, m_sepd_node );

  if ( !sepd ) {
    std::cout << PHWHERE << m_sepd_node << " node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // reset sEPD info
  m_sepd_energy.clear();
  m_sepd_isgood.clear();
  m_sepd_arm.clear();
  m_sepd_radius.clear();
  m_sepd_phi.clear();

  int num_sepd = sepd->size();
  auto epdgeom = findNode::getClass<EpdGeom>(topNode, "TOWERGEOM_EPD");
  for ( int i = 0; i < num_sepd; ++i ) {
    auto tower = sepd->get_tower_at_channel(i);

    if ( !tower ) {
      std::cout << PHWHERE << "sEPD tower at channel " << i << " is missing, skipping." << std::endl;
      continue;
    }
    unsigned int key = TowerInfoDefs::encode_epd(i);
    int epd_isgood = (tower->get_isGood()?1:0);
    float epd_e = tower->get_energy();
    int epd_arm = TowerInfoDefs::get_epd_arm(key);
	  float epd_phi = epdgeom->get_phi(key);
	  float epd_r = epdgeom->get_r(key);
	
	
    m_sepd_energy.push_back(epd_e);
    m_sepd_isgood.push_back(epd_isgood);
    m_sepd_arm.push_back(epd_arm);
    m_sepd_radius.push_back(epd_r);
    m_sepd_phi.push_back(epd_phi);

  }
  if ( Verbosity() > 1 ) {
    std::cout << "AnaTreeWriter::GetSEPDInfo - sEPD - Num Towers: " << m_sepd_energy.size() << std::endl;
  }


  return Fun4AllReturnCodes::EVENT_OK;
}

int TreeWriter::GetEventPlaneInfo( PHCompositeNode *topNode )
{
  m_psi_shifted_NS.clear();
  m_psi_shifted_S.clear();
  m_psi_shifted_N.clear();
  m_psi_NS.clear();
  m_psi_S.clear();
  m_psi_N.clear();
  auto m_evpmap  = findNode::getClass<EventplaneinfoMap>(topNode,  m_eventplane_node );
  if ( !m_evpmap  || m_evpmap->empty() ) {
    std::cout << PHWHERE << "EventplaneinfoMap node missing, doing nothing." << std::endl;
  } else {

    auto EPDNS = m_evpmap->get(EventplaneinfoMap::sEPDNS);
    if ( EPDNS ) {
      if (EPDNS->get_shifted_psi(2)){
        m_psi_shifted_NS.push_back(EPDNS->get_shifted_psi(2));
        m_psi_NS.push_back(EPDNS->get_psi(2));
      } else {
        m_psi_shifted_NS.push_back(-999);
        m_psi_NS.push_back(-999);
      }

      if ( EPDNS->get_shifted_psi(3) ) {
        m_psi_shifted_NS.push_back(EPDNS->get_shifted_psi(3));
        m_psi_NS.push_back(EPDNS->get_psi(3));
      } else {
        m_psi_shifted_NS.push_back(-999);
        m_psi_NS.push_back(-999);
      }
    }
    auto EPDS = m_evpmap->get(EventplaneinfoMap::sEPDS);
    if ( EPDS ) {
      if (EPDS->get_shifted_psi(2)){
        m_psi_shifted_S.push_back(EPDS->get_shifted_psi(2));
        m_psi_S.push_back(EPDS->get_psi(2));
      } else {
        m_psi_shifted_S.push_back(-999);
        m_psi_S.push_back(-999);
      }
      if (EPDS->get_shifted_psi(3)){
        m_psi_shifted_S.push_back(EPDS->get_shifted_psi(3));
        m_psi_S.push_back(EPDS->get_psi(3));
      } else {
        m_psi_shifted_S.push_back(-999);
        m_psi_S.push_back(-999);
      }
    }
    auto EPDN = m_evpmap->get(EventplaneinfoMap::sEPDN);
    if ( EPDN ) {
      if (EPDN->get_shifted_psi(2)){
        m_psi_shifted_N.push_back(EPDN->get_shifted_psi(2));
        m_psi_N.push_back(EPDN->get_psi(2));
      } else {
        m_psi_shifted_N.push_back(-999);
        m_psi_N.push_back(-999);
      }
      if (EPDN->get_shifted_psi(3)){
        m_psi_shifted_N.push_back(EPDN->get_shifted_psi(3));
        m_psi_N.push_back(EPDN->get_psi(3));
      } else {
        m_psi_shifted_N.push_back(-999);
        m_psi_N.push_back(-999);
      }

  }

  }


  if ( Verbosity() > 1 ) {
    std::cout << "AnaTreeWriter::GetEventPlaneInfo - Event Plane info: "
              << "NS: (" << m_psi_shifted_NS[0] << ", " << m_psi_NS[0] << "), "
              << "S: (" << m_psi_shifted_S[0] << ", " << m_psi_S[0] << "), "
              << "N: (" << m_psi_shifted_N[0] << ", " << m_psi_N[0] << ")" 
              << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}


void TreeWriter::CalcFlow( float b , float eta, float pt, 
  float &v2 , float &v3, float &v4, float &v5, float &v6 , 
  TRandom3 * engine, bool do_fluc, float scale)
{
    
    float a1 = 0.4397 * std::exp(-(b - 4.526) * (b - 4.526) / 72.0) + 0.636;
    float a2 = 1.916 / (b + 2) + 0.1;
    float a3 = 4.79 * 0.0001 * (b - 0.621) * (b - 10.172) * (b - 23) + 1.2;  // this is >0 for b>0
    float a4 = 0.135 * std::exp(-0.5 * (b - 10.855) * (b - 10.855) / 4.607 / 4.607) + 0.0120;
    float temp1 = std::pow(pt, a1) / (1 + std::exp((pt - 3.0) / a3));
    float temp2 = std::pow(pt + 0.1, -a2) / (1 + std::exp(-(pt - 4.5) / a3));
    float temp3 = 0.01 / (1 + std::exp(-(pt - 4.5) / a3));

    v2 = (a4 * (temp1 + temp2) + temp3) * exp(-0.5 * eta * eta / 3.43 / 3.43);
    
    float v2_sqrt = std::sqrt(v2);
    float fb = (
        0.97 
        + (1.06 * std::exp(-0.5 * b * b / 3.2 / 3.2))
        ) * v2_sqrt;
    float gb = (
        1.096 
        + (1.36 * std::exp(-0.5 * b * b / 3.0 / 3.0))
        ) * v2_sqrt;

    v3 = std::pow( fb, 3); 
    v4 = std::pow( gb, 4); 
    v5 = std::pow( gb, 5); 
    v6 = std::pow( gb, 6); 

    if ( do_fluc )
    {
        const std::array<float, 8> coeffs = { -7.00411e-09, 4.24567e-07, -9.87748e-06, 0.000112689, -0.000694686, 0.002413930, -0.00324709, 0.0107906};
        float sigma = 0.0;
        for ( float coeff : coeffs ) { sigma = sigma * b + coeff; }
        if (sigma < 0.0) { sigma = 0.0; }

        const float s2 = sigma;
        const float s3 = 1.5*std::pow(fb,3) * std::sqrt(v2) * sigma;
        const float s4 = 2.0*std::pow(gb,4) * v2 * sigma;
        const float s5 = 2.5*std::pow(gb,5) * std::sqrt(v2*v2*v2) * sigma;
        const float s6 = 3.0*std::pow(gb,6) * v2*v2 * sigma;

        float u;
        float v;
        float s;
        do 
        {
            u = 2.0*engine->Uniform() - 1.0;
            v = 2.0*engine->Uniform() - 1.0;
            s = u*u + v*v;
        } while (s >= 1.0 || s == 0.0);

        const float mul = std::sqrt(-2.0*std::log(s)/s);
        const float z0  = u*mul;
        const float z1  = v*mul;
        if ( v2 != 0.0 )
        {
            float _v2 = std::hypot(v2 + s2*z0,  s2*z1);
            v2 = std::clamp(_v2, 0.0F, 1.0F);
        }
        if ( v3 != 0.0 )
        {
            float _v3 = std::hypot(v3 + s3*z0,  s3*z1);
            v3 = std::clamp(_v3, 0.0F, 1.0F);
        }
        if ( v4 != 0.0 )
        {
            float _v4 = std::hypot(v4 + s4*z0,  s4*z1);
            v4 = std::clamp(_v4, 0.0F, 1.0F);
        }
        if ( v5 != 0.0 )
        {
            float _v5 = std::hypot(v5 + s5*z0,  s5*z1);
            v5 = std::clamp(_v5, 0.0F, 1.0F);
        }
        if ( v6 != 0.0 )
        {
            float _v6 = std::hypot(v6 + s6*z0,  s6*z1);
            v6 = std::clamp(_v6, 0.0F, 1.0F);
        }
    }
    v2 *= scale;
    v3 *= scale;
    v4 *= scale;
    v5 *= scale;
    v6 *= scale;
    return;
}





   