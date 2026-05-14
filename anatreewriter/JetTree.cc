#include "JetTree.h"

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


int JetTree::Init( PHCompositeNode * /*topNode*/ )
{
  
  PHTFileServer::get().open( m_output_filename, "RECREATE" );

  if ( Verbosity () > 0 ) 
  {
    std::cout << "JetTree::Init - opening file " << m_output_filename << std::endl;
  } 

  m_event_id = -1;

  m_tree = new TTree( "T", "T" );
  m_tree -> Branch( "event_id", &m_event_id, "event_id/I" );
  if ( !m_zvrtx_node.empty() ) 
  {
    m_tree -> Branch( "zvrtx", &m_zvrtx, "zvrtx/F" );
  }
  if ( !m_cent_node.empty() ) 
  {
    m_tree -> Branch( "cent", &m_cent, "cent/I" );
  }
  if ( !m_header_node.empty() )
  {
    m_tree -> Branch( "b", &m_b, "b/F" );
    m_tree -> Branch( "ep_angle", &m_ep_angle, "ep_angle/F" );
    m_tree -> Branch( "ecc", &m_ecc, "ecc/F" );
    m_tree -> Branch( "ncoll", &m_ncoll, "ncoll/F" );
    m_tree -> Branch( "npart", &m_npart, "npart/F" );
    m_tree -> Branch( "psiN", m_psin, Form("psiN[%d]/F", (int)k_max_psin) );    
  }
  for ( unsigned int i = 0; i <m_ncalo_nodes; ++i ) 
  {
    const std::string & calo_node = m_calo_nodes[i];
    const std::string & calo_nick = m_calo_nicknames_map[calo_node];
    m_tree -> Branch( Form("%s_sumE", calo_nick.c_str()), &m_calo_sumE_map[calo_node], Form("%s_sumE/F", calo_nick.c_str()) );
    m_tree -> Branch( Form("%s_sumEt", calo_nick.c_str()), &m_calo_sumEt_map[calo_node], Form("%s_sumEt/F", calo_nick.c_str()) );
    m_tree -> Branch( Form("%s_ieta_avgE", calo_nick.c_str()), &m_calo_ieta_avgE_map[calo_node], Form("%s_ieta_avgE[%d]/F", calo_nick.c_str(), k_ieta) );
    m_tree -> Branch( Form("%s_iphi_avgE", calo_nick.c_str()), &m_calo_iphi_avgE_map[calo_node], Form("%s_iphi_avgE[%d]/F", calo_nick.c_str(), k_iphi) );
  }
  for ( unsigned int i = 0; i < m_nrho_nodes; ++i ) 
  {
    const std::string & rho_node = m_rho_nodes[i];
    const std::string & rho_nick = m_rho_nicknames_map[rho_node];
    m_tree -> Branch( Form("%s_mu", rho_nick.c_str()), &m_rho_map[rho_node], Form("%s_mu/F", rho_nick.c_str()) );
    m_tree -> Branch( Form("%s_sigma", rho_nick.c_str()), &m_rho_sigma_map[rho_node], Form("%s_sigma/F", rho_nick.c_str()) );
  }
  for ( unsigned int i = 0; i < m_ntowerbkgd_nodes; ++i ) 
  {
    const std::string & towerbkgd_node = m_towerbkgd_nodes[i];
    const std::string & towerbkgd_nick = m_towerbkgd_nicknames_map[towerbkgd_node];
    m_tree -> Branch( Form("%s_cemc", towerbkgd_nick.c_str()), &m_towerbkgd_cemc_map[towerbkgd_node].data(), Form("%s_cemc[%d]/F", towerbkgd_nick.c_str(), k_ieta) );
    m_tree -> Branch( Form("%s_hcalin", towerbkgd_nick.c_str()), &m_towerbkgd_hcalin_map[towerbkgd_node].data(), Form("%s_hcalin[%d]/F", towerbkgd_nick.c_str(), k_ieta) );
    m_tree -> Branch( Form("%s_hcalout", towerbkgd_nick.c_str()), &m_towerbkgd_hcalout_map[towerbkgd_node].data(), Form("%s_hcalout[%d]/F", towerbkgd_nick.c_str(), k_ieta) );
    m_tree -> Branch( Form("%s_v2", towerbkgd_nick.c_str()), &m_towerbkgd_v2_map[towerbkgd_node], Form("%s_v2/F", towerbkgd_nick.c_str()) );
    m_tree -> Branch( Form("%s_flowfail", towerbkgd_nick.c_str()), &m_towerbkgd_flowfail_map[towerbkgd_node], Form("%s_flowfail/I", towerbkgd_nick.c_str()) );
  }
  for ( unsigned int i = 0; i < m_njet_nodes ; ++i ) 
  {
    const std::string & jet_node = m_jet_nodes[i];
    const std::string & jet_nick = m_jet_nicknames_map[jet_node];
    const float jet_R = m_jetR_map[jet_node];
    const JET_TYPE jet_type = m_jet_type_map[jet_node];

    m_tree -> Branch( Form("%s_R", jet_nick.c_str()), &jet_R, Form("%s_R/F", jet_nick.c_str()) );
    m_tree -> Branch( Form("%s_pT", jet_nick.c_str()), &m_jet_pT[jet_node] );
    m_tree -> Branch( Form("%s_E", jet_nick.c_str()), &m_jet_E[jet_node] );
    m_tree -> Branch( Form("%s_eta", jet_nick.c_str()), &m_jet_eta[jet_node] );
    m_tree -> Branch( Form("%s_phi", jet_nick.c_str()), &m_jet_phi[jet_node] );
    if ( jet_type == JET_TYPE::TRUTH )
    {
      continue; // for truth jets, only fill basic kinematics for now
    }
    if ( jet_type == JET_TYPE::SEED )
    { 
      m_tree -> Branch( Form("%s_maxD", jet_nick.c_str()), &m_jet_maxD[jet_node] );
      m_tree -> Branch( Form("%s_avgD", jet_nick.c_str()), &m_jet_avgD[jet_node] );
      m_tree -> Branch( Form("%s_supercomp_eT", jet_nick.c_str()), &m_jet_supercomp_eT[jet_node] );
    }
    if ( jet_type == JET_TYPE::SUB1 )
    {
      m_tree -> Branch( Form("%s_unsub_pT", jet_nick.c_str()), &m_jet_unsub_pT[jet_node] );
      m_tree -> Branch( Form("%s_unsub_E", jet_nick.c_str()), &m_jet_unsub_E[jet_node] );
    }
    if ( jet_type == JET_TYPE::RAW )
    {
      m_tree -> Branch( Form("%s_area_layer", jet_nick.c_str()), &m_jet_area_layer[jet_node] );
      m_tree -> Branch( Form("%s_E_layer", jet_nick.c_str()), &m_jet_E_layer[jet_node] );
      m_tree -> Branch( Form("%s_N_layer", jet_nick.c_str()), &m_jet_N_layer[jet_node] );
    }
    m_tree -> Branch( Form("%s_comp_ieta", jet_nick.c_str()), &m_jet_comp_ieta[jet_node] );
    m_tree -> Branch( Form("%s_comp_iphi", jet_nick.c_str()), &m_jet_comp_iphi[jet_node] );
    m_tree -> Branch( Form("%s_comp_caloid", jet_nick.c_str()), &m_jet_comp_caloid[jet_node] );
    m_tree -> Branch( Form("%s_comp_E", jet_nick.c_str()), &m_jet_comp_E[jet_node] );
    m_tree -> Branch( Form("%s_comp_eta", jet_nick.c_str()), &m_jet_comp_eta[jet_node] );
    m_tree -> Branch( Form("%s_comp_phi", jet_nick.c_str()), &m_jet_comp_phi[jet_node] );
    if ( Verbosity() > 0 ) 
    {
      std::cout << "JetTree::Init - Registered jet node: " << jet_node << " with R = " << jet_R << " and type = " << jet_type << std::endl;
    }
  }

  
  if ( Verbosity () > 0 )
  {
    std::cout << "JetTree::Init - done" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int JetTree::process_event( PHCompositeNode *topNode )
{

  m_event_id++; 

  
  if(Verbosity() > 1) 
  {
    std::cout << "JetTree::process_event - Process event " << m_event_id << std::endl;
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

  GetCaloInfo(topNode);

  // fill tree
  m_tree->Fill();
  
  return Fun4AllReturnCodes::EVENT_OK;

}

int JetTree::End( PHCompositeNode * /*topNode*/ )
{
  
  if( Verbosity() > 0 ) 
  {
    std::cout << "JetTree::EndRun - End run " << std::endl;
    std::cout << "JetTree::EndRun - Writing to " << m_output_filename << std::endl;
  }

  PHTFileServer::get().cd(m_output_filename); 

  m_tree->Write();
  

  if ( Verbosity() > 0 ) 
  {
    std::cout << "JetTree::EndRun - Writing run tree" << std::endl;
  }
  PHTFileServer::get().close();
 
  if ( Verbosity () > 0 ) 
  {
    std::cout << "JetTree::EndRun - done" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int JetTree::GetZvtx( PHCompositeNode *topNode )
{
    // get zvtx
    ResetZvtx();
    auto vertexmap = findNode::getClass<GlobalVertexMap>( topNode, m_zvrtx_node );
    if ( !vertexmap  || vertexmap->empty() ) 
    {
      std::cout << PHWHERE << "" << m_zvrtx_node << " node missing, skipping event." << std::endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }
  
    
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

int JetTree::GetCentInfo( PHCompositeNode *topNode )
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

int JetTree::GetEventHeaderInfo( PHCompositeNode *topNode )
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

int JetTree::GetSub1JetInfo( PHCompositeNode *topNode , const int idx)
{

  const std::string & node_name = m_sub1jet_nodes[idx];
 
  ResetSub1Jet();
  ResetTowerBkgd();

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


  
  return Fun4AllReturnCodes::EVENT_OK;

}

int JetTree::GetTruthJetInfo( PHCompositeNode *topNode , const int idx)
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


  
  return Fun4AllReturnCodes::EVENT_OK;

}

void JetTree::SumCaloE( PHCompositeNode *topNode, const std::string &towerinfo_node, const std::string &geo_node, const float zvrtx, float &sum_e )
{
  sum_e = 0.0;
  auto towerinfos = findNode::getClass< TowerInfoContainer >( topNode, towerinfo_node );
  if ( !towerinfos ) 
  {
    return;
  }
  auto geocont = findNode::getClass<RawTowerGeomContainer>(topNode, geo_node);
  if ( !geocont )
  {
    std::cout << PHWHERE << " RawTowerGeomContainer node " << geo_node << " missing, skipping." << std::endl;
    return;
  }

  RawTowerDefs::CalorimeterId geocaloid{RawTowerDefs::CalorimeterId::NONE};

  const double CALO_RADII[3] = {93.5, 127.503, 225.87};
  int calo_id = 0;
  if ( towerinfo_node.find("CEMC") != std::string::npos ) 
  {
    calo_id = 0;
  }
  else if ( towerinfo_node.find("HCALIN") != std::string::npos ) 
  {
    calo_id = 1;
  }
  else if ( towerinfo_node.find("HCALOUT") != std::string::npos ) 
  {
    calo_id = 2;
  }
  
  const double this_R = CALO_RADII[calo_id];
  
  if ( geo_node.find("HCALIN") != std::string::npos ) 
  {
    geocaloid = RawTowerDefs::CalorimeterId::HCALIN;
  }
  else if ( geo_node.find("HCALOUT") != std::string::npos ) 
  {
    geocaloid = RawTowerDefs::CalorimeterId::HCALOUT;
  }
  else if ( geo_node.find("CEMC") != std::string::npos ) 
  {
    geocaloid = RawTowerDefs::CalorimeterId::CEMC;
  }
  else 
  {
    std::cout << PHWHERE << " Warning: cannot determine calo id from geo node name " << geo_node << ", skipping." << std::endl;
    return;
  }
  
  // std::cout << PHWHERE << " - Calo id = " << calo_id << ", R = " << this_R << std::endl;

  auto ntowers = towerinfos->size();
  for ( unsigned int ich = 0; ich < ntowers; ich++ ) 
  {
    auto tower = towerinfos->get_tower_at_channel(ich);
    assert(tower);

    unsigned int key = towerinfos -> encode_key(ich);
    int ieta = towerinfos -> getTowerEtaBin(key);
    int iphi = towerinfos -> getTowerPhiBin(key);
    const RawTowerDefs::keytype geokey = RawTowerDefs::encode_towerid(geocaloid, ieta, iphi);

    if ( ! tower->get_isGood() ) 
    {
      continue; // skip bad towers
    }
    if (std::isnan(tower->get_energy()) ) 
    {
      continue; // skip towers with bad energy
    }

    auto tower_geom = geocont -> get_tower_geometry(geokey);
    assert(tower_geom);

    double eta = tower_geom -> get_eta();
    double E = tower -> get_energy();
    double z0 = sinh(eta) * this_R;
    double z = z0 - zvrtx;
    double eta_shifted = asinh(z / this_R);  // eta after shift from
    double eT = E / cosh(eta_shifted);
    

    sum_e += eT;
  } 

  if ( Verbosity() > 1 ) 
  {
    std::cout << PHWHERE << " - " << towerinfo_node << ": sum_e = " << sum_e << std::endl;
  }

  return;
}

int JetTree::GetCaloInfo( PHCompositeNode *topNode )
{
  ResetCaloInfo();

  SumCaloE( topNode, "TOWERINFO_CALIB_CEMC_RETOWER", "TOWERGEOM_HCALIN", m_zvtx, m_sumeT_cemc );
  SumCaloE( topNode, "TOWERINFO_CALIB_HCALIN", "TOWERGEOM_HCALIN", m_zvtx, m_sumeT_hcalin );
  SumCaloE( topNode, "TOWERINFO_CALIB_HCALOUT", "TOWERGEOM_HCALOUT", m_zvtx, m_sumeT_hcalout );
  SumCaloE( topNode, "TOWERINFO_CALIB_CEMC_RETOWER_SUB1", "TOWERGEOM_HCALIN", m_zvtx, m_sumeT_cemc_sub1 );
  SumCaloE( topNode, "TOWERINFO_CALIB_HCALIN_SUB1", "TOWERGEOM_HCALIN", m_zvtx, m_sumeT_hcalin_sub1 );
  SumCaloE( topNode, "TOWERINFO_CALIB_HCALOUT_SUB1","TOWERGEOM_HCALOUT", m_zvtx,  m_sumeT_hcalout_sub1 );
  SumCaloE( topNode, "TOWERINFO_CALIB_ORIGINAL_CEMC", "TOWERGEOM_CEMC", m_zvtx, m_sumeT_cemc_org );
  SumCaloE( topNode, "TOWERINFO_CALIB_ORIGINAL_HCALIN", "TOWERGEOM_HCALIN", m_zvtx, m_sumeT_hcalin_org );
  SumCaloE( topNode, "TOWERINFO_CALIB_ORIGINAL_HCALOUT", "TOWERGEOM_HCALOUT", m_zvtx, m_sumeT_hcalout_org );

  if ( Verbosity() > 1 ) 
  {
    std::cout << PHWHERE << " - sumeT_cemc = " << m_sumeT_cemc << ", sumeT_hcalin = " << m_sumeT_hcalin << ", sumeT_hcalout = " << m_sumeT_hcalout << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}






   