#include "AnaTreeWriter.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>
#include <fun4all/Fun4AllServer.h>

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

#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfo.h>

#include <mbd/MbdOut.h>
#include <mbd/MbdOutV1.h>
#include <mbd/MbdOutV2.h>

#include <ffarawobjects/Gl1Packet.h>
#include <ffarawobjects/Gl1Packetv2.h>

#include <jetbase/Jet.h>
#include <jetbase/Jetv2.h>
#include <jetbase/JetContainer.h>
#include <jetbase/JetContainerv1.h>

#include <jetbackground/TowerRho.h>
#include <jetbackground/TowerRhov1.h>

#include <TTree.h>

#include <cstdlib>
#include <iostream>
#include <cstdint>
#include <cmath>


AnaTreeWriter::AnaTreeWriter( const std::string  &outputfile ) 
  : SubsysReco("AnaTreeWriter") 
  , m_output_filename( outputfile ) 
{

}

int AnaTreeWriter::Init( PHCompositeNode * /*topNode*/ )
{
  
  // create output file
  PHTFileServer::get().open( m_output_filename, "RECREATE" );

  if ( Verbosity () > 0 ) {
    std::cout << "AnaTreeWriter::Init - opening file " << m_output_filename << std::endl;
  } 

  // Reset counters
  m_run_tree = new TTree( "RunTree", "RunTree" );
  m_run_tree -> Branch( "num_events", &m_nevents, "num_events/I" );
  m_run_tree -> Branch( "do_min_tower_energy", &m_do_min_tower_energy, "do_min_tower_energy/O" );
  m_run_tree -> Branch( "min_tower_energy", &m_min_tower_energy, "min_tower_energy/F" );
  m_run_tree -> Branch( "do_online_jet_eta_cut", &m_do_online_eta_cut, "do_online_jet_eta_cut/O" );
  m_run_tree -> Branch( "do_min_jet_energy", &m_do_min_jet_energy, "do_min_jet_energy/O" );
  m_run_tree -> Branch( "min_jet_energy", &m_min_jet_energy, "min_jet_energy/F" );
  if ( !m_gl1_node.empty() ) {
    m_run_tree -> Branch( "active_triggers", &m_run_trigger_status, "active_triggers[64]/O" );
    m_run_tree -> Branch( "run_scalars_live", &m_run_live_scalar, "run_scalars_live[64]/L" );
    m_run_tree -> Branch( "run_scalars_scaled", &m_run_scaled_scalar, "run_scalars_scaled[64]/L" );
    m_run_tree -> Branch( "run_scalars_raw", &m_run_raw_scalar, "run_scalars_raw[64]/L" );
    m_run_trigger_status.fill(false);
    m_run_live_scalar.fill(0);
    m_run_scaled_scalar.fill(0);
    m_run_raw_scalar.fill(0);
  }
  
  // create tree
  m_tree = new TTree( "EventTree", "EventTree" );
  m_tree -> Branch( "event_id", &m_event_id, "event_id/I" );
  m_event_id = -1;

  // zvtx
  if ( !m_zvrtx_node.empty() ) {
    m_tree -> Branch( "zvrtx", &m_zvtx, "zvrtx/F" );
    if ( Verbosity() > 0 ) {
      std::cout << "AnaTreeWriter::Init - Registered zvtx" << std::endl;
    }
  }

  // centrality
  if ( !m_cent_node.empty() ) {
    m_tree->Branch("centrality", &m_centrality, "centrality/I");
    if ( Verbosity() > 0 ) {
      std::cout << "AnaTreeWriter::Init - Registered centrality" << std::endl;
    }
  }
  
  // gl1 info
  if ( !m_gl1_node.empty() ) {
    m_gl1_tree = new TTree( m_gl1_node.c_str(), m_gl1_node.c_str() );
    m_gl1_tree -> Branch( "event_id", &m_event_id, "event_id/I" );
    m_gl1_tree -> Branch( "s_triggervec", &s_triggervec, "s_triggervec/l" );
    m_gl1_tree -> Branch( "l_triggervec", &l_triggervec, "l_triggervec/l" );
    m_gl1_tree -> Branch( "gl1_trigger_status", &m_gl1_trigger_status, "gl1_trigger_status[64]/O" );
    m_gl1_tree -> Branch( "gl1_live_scalar", &m_gl1_live_scalar, "gl1_live_scalar[64]/L" );
    m_gl1_tree -> Branch( "gl1_scaled_scalar", &m_gl1_scaled_scalar, "gl1_scaled_scalar[64]/L" );
    m_gl1_tree -> Branch( "gl1_raw_scalar", &m_gl1_raw_scalar, "gl1_raw_scalar[64]/L" );
    if ( Verbosity() > 0 ) {
      std::cout << "AnaTreeWriter::Init - GL1 node: " << m_gl1_node << std::endl;
    }

    m_gl1_tree_id = -1;
    m_tree -> Branch( Form("%s_entry", m_gl1_node.c_str()), &m_gl1_tree_id, Form("%s_entry/I", m_gl1_node.c_str()) );

  }

  // MBD
  if ( !m_mbd_node.empty() ) {
    m_mbd_tree = new TTree( m_mbd_node.c_str(), m_mbd_node.c_str() );
    m_mbd_tree -> Branch( "event_id", &m_event_id, "event_id/I" );
    m_mbd_tree -> Branch( "mbd_q_N", &m_mbd_q_N, "mbd_q_N/F" );
    m_mbd_tree -> Branch( "mbd_q_S", &m_mbd_q_S, "mbd_q_S/F" );
    m_mbd_tree -> Branch( "mbd_time_N", &m_mbd_time_N, "mbd_time_N/F" );
    m_mbd_tree -> Branch( "mbd_time_S", &m_mbd_time_S, "mbd_time_S/F" );
    m_mbd_tree -> Branch( "mbd_npmt_N", &m_mbd_npmt_N, "mbd_npmt_N/F" );
    m_mbd_tree -> Branch( "mbd_npmt_S", &m_mbd_npmt_S, "mbd_npmt_S/F" );
    if ( Verbosity() > 0 ) {
      std::cout << "AnaTreeWriter::Init - MBD node: " << m_mbd_node << std::endl;
    }

    m_mbd_tree_id = -1;
    m_tree -> Branch( Form("%s_entry", m_mbd_node.c_str()), &m_mbd_tree_id, Form("%s_entry/I", m_mbd_node.c_str()) );
  }

  // rho
  if ( m_rho_nodes.size() > 0 ) {
    m_rho_tree = new TTree( "RhoTree", "RhoTree" );
    m_rho_tree -> Branch( "event_id", &m_event_id, "event_id/I" );
    for ( unsigned int i = 0; i < m_rho_nodes.size(); ++i ) {
      m_rho_tree -> Branch( m_rho_nodes[i].c_str(), &m_rho_val[i], (m_rho_nodes[i] + "/F").c_str() );
      m_rho_tree -> Branch( (m_rho_nodes[i] + "_std").c_str(), &m_std_rho_val[i], (m_rho_nodes[i] + "_std/F").c_str() );
      if ( Verbosity() > 0 ) {
        std::cout << "AnaTreeWriter::Init - Registered rho nodes: " << m_rho_nodes[i] << std::endl;
      }
    }
    m_rho_tree_id = -1;
    m_tree -> Branch( "RhoTree_entry", &m_rho_tree_id, "RhoTree_entry/I" );
  }

  // calo nodes
  for ( unsigned int i = 0; i < m_calo_nodes.size(); ++i ) {
    m_calo_trees[i] = new TTree( m_calo_nodes[i].c_str(), m_calo_nodes[i].c_str() );
    m_calo_trees[i] -> Branch( "event_id", &m_event_id, "event_id/I" );
    m_calo_trees[i] -> Branch( "num_towers", &m_num_towers, "num_towers/i" );
    m_calo_trees[i] -> Branch( "num_towers_fired", &m_num_towers_fired, "num_towers_fired/i" );
    m_calo_trees[i] -> Branch( "num_towers_dead", &m_num_towers_dead, "num_towers_dead/i" );
    m_calo_trees[i] -> Branch( "tower_energy", &m_tower_energy );
    m_calo_trees[i] -> Branch( "tower_energy_time", &m_tower_energy_time );
    m_calo_trees[i] -> Branch( "tower_status_ieta_iphi", &m_tower_status_ieta_iphi );

    if ( Verbosity() > 0 ) {
      std::cout << "AnaTreeWriter::Init - Registered Calo nodes: " << m_calo_nodes[i] << std::endl;
    }

    m_calo_tree_id[i] = -1;
    m_tree -> Branch( Form("%s_entry", m_calo_nodes[i].c_str()), &m_calo_tree_id[i], Form("%s_entry/I", m_calo_nodes[i].c_str()) );
   

  }

  // jets
  for ( unsigned int i = 0; i < m_jet_nodes.size(); ++i ) {

    m_jet_trees[i] = new TTree( m_jet_nodes[i].c_str(), m_jet_nodes[i].c_str() );
    m_jet_trees[i] -> Branch( "event_id", &m_event_id, "event_id/I" );
    m_jet_trees[i] -> Branch( "num_jets", &m_num_jets, "num_jets/I" );
    m_jet_trees[i] -> Branch( "jet_energy", &m_jet_energy );
    m_jet_trees[i] -> Branch( "jet_eta", &m_jet_eta );
    m_jet_trees[i] -> Branch( "jet_phi", &m_jet_phi );
    m_jet_trees[i] -> Branch( "jet_pT", &m_jet_pT );
    m_jet_trees[i] -> Branch( "jet_eT", &m_jet_eT );
    m_jet_trees[i] -> Branch( "jet_ntowers", &m_jet_ntowers );
    m_jet_trees[i] -> Branch( "jet_ntowers_cemc", &m_jet_ntowers_cemc );
    m_jet_trees[i] -> Branch( "jet_ntowers_hcalin", &m_jet_ntowers_hcalin );
    m_jet_trees[i] -> Branch( "jet_ntowers_hcalout", &m_jet_ntowers_hcalout );
    m_jet_trees[i] -> Branch( "jet_energy_hcalin", &m_jet_energy_hcalin );
    m_jet_trees[i] -> Branch( "jet_energy_hcalout", &m_jet_energy_hcalout );
    m_jet_trees[i] -> Branch( "jet_energy_cemc", &m_jet_energy_cemc );

    if ( Verbosity() > 0 ) {
      std::cout << "AnaTreeWriter::Init - Registered Jet nodes: " << m_jet_nodes[i] << std::endl;
    }

    m_jet_tree_id[i] = -1;
    m_tree -> Branch( Form("%s_entry", m_jet_nodes[i].c_str()), &m_jet_tree_id[i], Form("%s_entry/I", m_jet_nodes[i].c_str()) );

  }


  if ( Verbosity () > 0 ){
    std::cout << "AnaTreeWriter::Init - done" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaTreeWriter::process_event( PHCompositeNode *topNode )
{

  m_event_id++; 

  if(Verbosity() > 1) {
    std::cout << "AnaTreeWriter::process_event - Process event " << m_event_id << std::endl;
  }

  if ( !m_cent_node.empty() ) { // get centrality
    auto res = GetCentInfo(topNode);
    if ( res != Fun4AllReturnCodes::EVENT_OK ) { return res; }
  }

  if ( !m_zvrtx_node.empty() ) { // get zvtx
    auto res = GetZvtx(topNode);
    if ( res != Fun4AllReturnCodes::EVENT_OK ) { return res; }
  }

  if ( !m_gl1_node.empty() ){  // get GL1
    auto res = GetGL1(topNode);
    if ( res != Fun4AllReturnCodes::EVENT_OK ) { return res; }
  }
 
  if ( !m_mbd_node.empty() ) { // get MBD
    auto res = GetMbdInfo(topNode);
    if ( res != Fun4AllReturnCodes::EVENT_OK ) { return res; }
  }

  for ( unsigned int i = 0; i < m_calo_nodes.size(); ++i ) { // get calo info
    auto res = GetCaloInfo(topNode, i);
    if ( res != Fun4AllReturnCodes::EVENT_OK ) { return res; }
  }

  for ( unsigned int i = 0; i < m_jet_nodes.size(); ++i ) {
    auto res = GetJetInfo(topNode, i);
    if ( res != Fun4AllReturnCodes::EVENT_OK ) { return res; }
  }

  if ( m_rho_nodes.size() > 0 ) {
    auto res = GetRhoInfo(topNode);
    if ( res != Fun4AllReturnCodes::EVENT_OK ) { return res; }
  }
    
  m_tree->Fill();
  
  return Fun4AllReturnCodes::EVENT_OK;

}

int AnaTreeWriter::End( PHCompositeNode * /*topNode*/ )
{
  
  if( Verbosity() > 0 ) {
    std::cout << "AnaTreeWriter::EndRun - End run " << std::endl;
    std::cout << "AnaTreeWriter::EndRun - Writing to " << m_output_filename << std::endl;
  }

  PHTFileServer::get().cd(m_output_filename); 
  m_tree->Write();
  for ( unsigned int i = 0; i < m_jet_nodes.size(); ++i ) {
    m_jet_trees[i]->Write();
  }
  for ( unsigned int i = 0; i < m_calo_nodes.size(); ++i ) {
    m_calo_trees[i]->Write();
  }
  if ( !m_gl1_node.empty() ) {
    m_gl1_tree->Write();
  }
  if ( !m_mbd_node.empty() ) {
    m_mbd_tree->Write();
  }
  if ( m_rho_nodes.size() > 0 ) {
    m_rho_tree->Write();
  }
  for (unsigned int i = 0; i < 64; ++i) {
    uint64_t live_i = m_run_live_scalar[i];
    uint64_t scaled_i = m_run_scaled_scalar[i];
    uint64_t raw_i = m_run_raw_scalar[i];
    uint64_t live_f = m_gl1_live_scalar[i];
    uint64_t scaled_f = m_gl1_scaled_scalar[i];
    uint64_t raw_f = m_gl1_raw_scalar[i];
    m_run_live_scalar[i] = live_f-live_i;
    m_run_scaled_scalar[i] = scaled_f-scaled_i;
    m_run_raw_scalar[i] = raw_f-raw_i;
  }
  m_nevents = m_event_id+1;
  m_run_tree->Fill();
  m_run_tree->Write();

  if ( Verbosity () > 0 ) {
    std::cout << "AnaTreeWriter::EndRun - done" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaTreeWriter::ResetEvent( PHCompositeNode * /*topNode*/ )
{ 
 
  // zvtx
  m_zvtx = 0;

  // gl1
  m_gl1_trigger_status.fill(false);
  m_gl1_raw_scalar.fill(0);
  m_gl1_scaled_scalar.fill(0);
  m_gl1_live_scalar.fill(0);

  // jets
  m_num_jets = 0;
  m_jet_eta.clear();
  m_jet_phi.clear();
  m_jet_eT.clear();
  m_jet_pT.clear();
  m_jet_energy.clear();
  m_jet_ntowers.clear();
  m_jet_energy_cemc.clear();
  m_jet_energy_hcalin.clear();
  m_jet_energy_hcalout.clear();
  m_jet_ntowers_cemc.clear();
  m_jet_ntowers_hcalin.clear();
  m_jet_ntowers_hcalout.clear();

  // calo
  m_num_towers = 0;
  m_num_towers_fired = 0;
  m_num_towers_dead = 0;
  m_tower_energy.clear();
  m_tower_energy_time.clear();
  m_tower_status_ieta_iphi.clear();

  // mbd
  m_rho_val.fill(0.0);
  m_std_rho_val.fill(0.0);

  return Fun4AllReturnCodes::EVENT_OK;

}

int AnaTreeWriter::GetGL1( PHCompositeNode *topNode )
{

  auto gl1 = findNode::getClass< Gl1Packetv2 >( topNode, m_gl1_node );
  if( !gl1 ) {
    std::cout << PHWHERE << " No GL1 packet found!" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // reset event trigger status
  m_gl1_trigger_status.fill(false);
  m_gl1_raw_scalar.fill(0);
  m_gl1_scaled_scalar.fill(0);
  m_gl1_live_scalar.fill(0);

  // get GL1 trigger vectors
  s_triggervec = gl1->getScaledVector();
  l_triggervec = gl1->getLiveVector();

  for (unsigned int i = 0; i < 64; ++i) {
    m_gl1_trigger_status[i] = (s_triggervec >> i) & 0x1;
    m_gl1_raw_scalar[i] = gl1->lValue(i, 0);
    m_gl1_scaled_scalar[i] = gl1->lValue(i, 1);
    m_gl1_live_scalar[i] = gl1->lValue(i, 2);
  }
  
  // update run trigger status
  for (unsigned int i = 0; i < 64; ++i) {
    if (m_gl1_trigger_status[i]) {
      m_run_trigger_status[i] = true;
    }
  }

  if ( m_event_id == 0 ) { // get first event scalars
    for (unsigned int i = 0; i < 64; ++i) {
      m_run_live_scalar[i] = m_gl1_live_scalar[i];
      m_run_scaled_scalar[i] = m_gl1_scaled_scalar[i];
      m_run_raw_scalar[i] = m_gl1_raw_scalar[i];
    }
  }

  m_gl1_tree->Fill();
  m_gl1_tree_id++;

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaTreeWriter::GetZvtx( PHCompositeNode *topNode )
{
    // get zvtx
    auto vertexmap = findNode::getClass<GlobalVertexMap>( topNode, m_zvrtx_node );
    if ( !vertexmap  || vertexmap->empty() ) {
      std::cout << PHWHERE << "" << m_zvrtx_node << " node missing, doing nothing." << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    auto vtx = vertexmap->begin()->second;
    m_zvtx = NAN;
    if ( vtx ) {
      m_zvtx = vtx->get_z();
    } 

    if ( std::isnan(m_zvtx) || m_zvtx > 1e3) {
      static bool once = true;
      if (once) {
        once = false;
        std::cout << PHWHERE << "vertex is " << m_zvtx << ". Drop all tower inputs (further vertex warning will be suppressed)." << std::endl;
      }
      return Fun4AllReturnCodes::ABORTEVENT;
    }


    if ( Verbosity() > 1 ) {
      std::cout << "AnaTreeWriter::GetZvtx - zvtx = " << m_zvtx << std::endl;
    }

    return Fun4AllReturnCodes::EVENT_OK;
}

int AnaTreeWriter::GetCentInfo( PHCompositeNode *topNode )
{
  // get centrality
  auto cent_node = findNode::getClass< CentralityInfo >( topNode, m_cent_node );
  if ( !cent_node ) {
    std::cout << PHWHERE << m_cent_node << " node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  
  if ( cent_node -> has_centrality_bin(CentralityInfo::PROP::mbd_NS) ) {
    m_centrality = (int)(cent_node->get_centrality_bin(CentralityInfo::PROP::mbd_NS));
  } else if ( cent_node -> has_centile(CentralityInfo::PROP::mbd_NS) ){
    m_centrality = (int)(cent_node->get_centile(CentralityInfo::PROP::mbd_NS));
  } 
  
  if (  m_centrality < 0 ) {
    std::cerr << PHWHERE << "CentralityInfo Node missing centrality information" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if ( Verbosity() > 1 ) {
    std::cout << "AnaTreeWriter::GetCentInfo - Centrality = " << m_centrality << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaTreeWriter::GetCaloInfo( PHCompositeNode *topNode ,  const int idx)
{
  const std::string & node_name = m_calo_nodes[idx];

  m_num_towers = 0;
  m_num_towers_fired = 0;
  m_num_towers_dead = 0;
  m_tower_energy.clear();
  m_tower_energy_time.clear();
  m_tower_status_ieta_iphi.clear();

  auto towers = findNode::getClass< TowerInfoContainer >( topNode, node_name );
  if ( !towers ) {
    std::cout << PHWHERE << " Input node "<< node_name << " Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN; // fatal error
  }

  m_num_towers = towers->size();
  for ( unsigned int channel = 0; channel < m_num_towers; channel++ ) {
      
    auto tower = towers -> get_tower_at_channel( channel );
    assert(tower);
    unsigned int key = towers -> encode_key(channel);
    int ieta = towers -> getTowerEtaBin(key);
    int iphi = towers -> getTowerPhiBin(key);
    uint8_t status = tower -> get_status();
    uint16_t ieta_iphi = ((ieta & 0xFF) << 8) | (iphi & 0xFF);
    uint32_t status_ieta_iphi = (status << 16) | ieta_iphi;
    float energy = tower -> get_energy();
    
    bool is_masked = tower->get_isHot() || tower->get_isNoCalib() || tower->get_isNotInstr() || tower->get_isBadChi2() || std::isnan(energy);
    if ( m_do_min_tower_energy && energy < m_min_tower_energy ) {
      is_masked = true;
    }
    if ( is_masked ) {
      m_num_towers_dead++;
    } else {      
      m_num_towers_fired++; 
      m_tower_status_ieta_iphi.push_back(status_ieta_iphi);
      m_tower_energy.push_back(energy);
      m_tower_energy_time.push_back(tower->get_time_float());
    }
    
  } 


  if ( Verbosity() > 1 ) {
    std::cout << "AnaTreeWriter::GetCaloInfo - Calo " << node_name << " - Num Towers: " << m_num_towers << ", Num Towers Fired: " << m_num_towers_fired << ", Num Towers Dead: " << m_num_towers_dead << std::endl;
  }
  
  m_calo_trees[idx]->Fill();
  m_calo_tree_id[idx]++;


  return Fun4AllReturnCodes::EVENT_OK;

}

int AnaTreeWriter::GetMbdInfo( PHCompositeNode *topNode )
{
  // get MBD info
  auto mbd = findNode::getClass<MbdOutV2>( topNode, m_mbd_node );
  if ( !mbd ) {
    static bool once = true;
    if ( once ) {
      std::cout << PHWHERE << m_mbd_node << " node missing, supressing further warnings." << std::endl;
      once = false;
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }

  m_mbd_q_N = mbd->get_q(k_mbd_north_code);
  m_mbd_q_S = mbd->get_q(k_mbd_south_code);
  m_mbd_time_N = mbd->get_time(k_mbd_north_code);
  m_mbd_time_S = mbd->get_time(k_mbd_south_code);
  m_mbd_npmt_N = mbd->get_npmt(k_mbd_north_code);
  m_mbd_npmt_S = mbd->get_npmt(k_mbd_south_code);

  if ( Verbosity() > 1 ) {
    std::cout << "AnaTreeWriter::GetMbdInfo - MBD info N:(q,t,n), S:(q,t,n) = (" << m_mbd_q_N << ", " << m_mbd_time_N << ", " << m_mbd_npmt_N << "), (" << m_mbd_q_S << ", " << m_mbd_time_S << ", " << m_mbd_npmt_S << ")" << std::endl;
  }

  m_mbd_tree->Fill();
  m_mbd_tree_id++;
  return Fun4AllReturnCodes::EVENT_OK;
}

int AnaTreeWriter::GetJetInfo( PHCompositeNode *topNode , const int idx)
{

  const std::string & node_name = m_jet_nodes[idx];

  m_num_jets = 0;
  m_jet_eta.clear();
  m_jet_phi.clear();
  m_jet_eT.clear();
  m_jet_pT.clear();
  m_jet_energy.clear();
  m_jet_ntowers.clear();
  m_jet_energy_cemc.clear();
  m_jet_energy_hcalin.clear();
  m_jet_energy_hcalout.clear();
  m_jet_ntowers_cemc.clear();
  m_jet_ntowers_hcalin.clear();
  m_jet_ntowers_hcalout.clear();


  // find "_Sub1" in the node name
  bool is_sub1 = std::string::npos != node_name.find("_Sub1");
  std::string jet_calo_cemc_node = "TOWERINFO_CALIB_CEMC_RETOWER";
  std::string jet_calo_hcalin_node = "TOWERINFO_CALIB_HCALIN";
  std::string jet_calo_hcalout_node = "TOWERINFO_CALIB_HCALOUT";
  if ( is_sub1 ) {
    jet_calo_cemc_node = "TOWERINFO_CALIB_CEMC_RETOWER_SUB1";
    jet_calo_hcalin_node = "TOWERINFO_CALIB_HCALIN_SUB1";
    jet_calo_hcalout_node = "TOWERINFO_CALIB_HCALOUT_SUB1";
  }


  auto towersEM3 = findNode::getClass< TowerInfoContainer >( topNode, jet_calo_cemc_node );
  auto towersIH3 = findNode::getClass< TowerInfoContainer >( topNode, jet_calo_hcalin_node);
  auto towersOH3 = findNode::getClass< TowerInfoContainer >( topNode, jet_calo_hcalout_node);
  bool do_frac = true;
  if( !towersEM3 || !towersIH3 || !towersOH3 ){
    std::cout
      <<"AnaTreeWriter::process_event - Error can not find raw tower node "
      << std::endl;
    do_frac = false;
  }

  
  auto jets = findNode::getClass<JetContainer>( topNode, node_name );
  if ( !jets ) {
    std::cout << PHWHERE << " Input node "<< node_name << " Node missing, doing nothing." << std::endl;
    exit(-1); // fatal error
  }

  m_jet_R = jets -> get_jetpar_R();

  // get jets
  m_num_jets = 0;
  for ( auto jet : *jets ) {

    if ( !accept_jet( jet->get_e(), jet->get_eta() ) ){ continue; }

    m_num_jets++;
    float eta = jet->get_eta();
    float phi = jet->get_phi();
    float pT = jet->get_pt();
    float eT = jet->get_et();
    float E = jet->get_e();

    int nconst = 0;
    int n_towers_cemc = 0, n_towers_hcalin = 0, n_towers_hcalout = 0;
    float cemc_energy = 0, hcalin_energy = 0, hcalout_energy = 0;
    if ( do_frac ){
      for (auto comp: jet->get_comp_vec()) {
        
        nconst++;
        if (comp.first == Jet::SRC::CEMC_TOWERINFO || comp.first == Jet::SRC::CEMC_TOWERINFO_RETOWER || comp.first == Jet::SRC::CEMC_TOWERINFO_SUB1) {
          auto tower = towersEM3 -> get_tower_at_channel(comp.second);
          if ( !tower  ) { continue; }
          cemc_energy += tower->get_energy();
          n_towers_cemc++;
        } else if (comp.first == Jet::SRC::HCALIN_TOWERINFO || comp.first == Jet::SRC::HCALIN_TOWERINFO_SUB1) {
          auto tower = towersIH3->get_tower_at_channel(comp.second);
          if ( !tower ) { continue; }
          hcalin_energy += tower->get_energy();
          n_towers_hcalin++;
        } else if (comp.first == Jet::SRC::HCALOUT_TOWERINFO || comp.first == Jet::SRC::HCALOUT_TOWERINFO_SUB1) {
          auto tower = towersOH3->get_tower_at_channel(comp.second);
          if ( !tower ) { continue; }
          hcalout_energy += tower->get_energy();
          n_towers_hcalout++;
        }
        
      }
    } else {
      for (auto comp: jet->get_comp_vec()) {
        nconst++;
        if (comp.first == Jet::SRC::CEMC_TOWERINFO || comp.first == Jet::SRC::CEMC_TOWERINFO_RETOWER || comp.first == Jet::SRC::CEMC_TOWERINFO_SUB1) {
          n_towers_cemc++;
        } else if (comp.first == Jet::SRC::HCALIN_TOWERINFO || comp.first == Jet::SRC::HCALIN_TOWERINFO_SUB1) {
          n_towers_hcalin++;
        } else if (comp.first == Jet::SRC::HCALOUT_TOWERINFO || comp.first == Jet::SRC::HCALOUT_TOWERINFO_SUB1) {
          n_towers_hcalout++;
        }
      }
    }

    m_jet_eta.push_back(eta);
    m_jet_phi.push_back(phi);
    m_jet_eT.push_back(eT);
    m_jet_pT.push_back(pT);
    m_jet_energy.push_back(E);
    m_jet_energy_cemc.push_back(cemc_energy);
    m_jet_energy_hcalin.push_back(hcalin_energy);
    m_jet_energy_hcalout.push_back(hcalout_energy);
    m_jet_ntowers.push_back(nconst);
    m_jet_ntowers_cemc.push_back(n_towers_cemc);
    m_jet_ntowers_hcalin.push_back(n_towers_hcalin);
    m_jet_ntowers_hcalout.push_back(n_towers_hcalout);

    if ( Verbosity() > 3 ) {
      std::cout << "AnaTreeWriter::GetJetInfo - Jet " << m_num_jets << " eta, phi, pT, eT, E = " << eta << ", " << phi << ", " << pT << ", " << eT << ", " << E << std::endl;
      std::cout << "AnaTreeWriter::GetJetInfo - Jet " << m_num_jets << " nconst, n_towers_cemc, n_towers_hcalin, n_towers_hcalout = " << nconst << ", " << n_towers_cemc << ", " << n_towers_hcalin << ", " << n_towers_hcalout << std::endl;
    }
   

  } // end of jet loop
     
  if ( Verbosity() > 0 ) {
    std::cout << "AnaTreeWriter::GetJetInfo - Number of jets = " << m_num_jets << std::endl;
  }


  m_jet_trees[idx]->Fill();
  m_jet_tree_id[idx]++;

  return Fun4AllReturnCodes::EVENT_OK;

}

int AnaTreeWriter::GetRhoInfo( PHCompositeNode *topNode )
{

  // get rho nodes
  for ( unsigned int i = 0; i < m_rho_nodes.size(); ++i ) {

    auto rho = findNode::getClass<TowerRhov1>(topNode, m_rho_nodes[i]);
    if ( !rho ) {
      std::cout << PHWHERE << " Input node " << m_rho_nodes[i] << " Node missing, doing nothing." << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    m_rho_val[i] = rho->get_rho();
    m_std_rho_val[i] = rho->get_sigma();
  }

  m_rho_tree->Fill();
  m_rho_tree_id++;
  return Fun4AllReturnCodes::EVENT_OK;
}



