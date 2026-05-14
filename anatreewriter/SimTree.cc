#include "SimTree.h"


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


#include <eventplaneinfo/Eventplaneinfo.h>
#include <eventplaneinfo/EventplaneinfoMap.h>


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


int SimTree::Init( PHCompositeNode * /*topNode*/ )
{
  
  PHTFileServer::get().open( m_output_filename, "RECREATE" );

  if ( Verbosity () > 0 ) 
  {
    std::cout << "SimTree::Init - opening file " << m_output_filename << std::endl;
  } 

  m_event_id = -1;

  m_tree = new TTree( "T", "T" );

  m_tree -> Branch( "event_id", &m_event_id, "event_id/I" );
  if ( !m_zvrtx_node.empty() )
  {
    m_tree -> Branch( "zvrtx", &m_zvtx, "zvrtx/F" );
  }
  if ( !m_cent_node.empty() )
  {
    m_tree -> Branch( "cent", &m_cent, "cent/I" );
  }

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
      std::cout << "SimTree::Init - Event header node: " << m_eventhead_node << std::endl;
    }
  }
  
  m_tree -> Branch( "sumeT_cemc", &m_sumeT_cemc, "sumeT_cemc/F" );
  m_tree -> Branch( "sumeT_hcalin", &m_sumeT_hcalin, "sumeT_hcalin/F" );
  m_tree -> Branch( "sumeT_hcalout", &m_sumeT_hcalout, "sumeT_hcalout/F" );
  
  m_tree -> Branch( "sumeT_cemc_sub1", &m_sumeT_cemc_sub1, "sumeT_cemc_sub1/F" );
  m_tree -> Branch( "sumeT_hcalin_sub1", &m_sumeT_hcalin_sub1, "sumeT_hcalin_sub1/F" );
  m_tree -> Branch( "sumeT_hcalout_sub1", &m_sumeT_hcalout_sub1, "sumeT_hcalout_sub1/F" );

  m_tree -> Branch( "sumeT_cemc_org", &m_sumeT_cemc_org, "sumeT_cemc_org/F" );
  m_tree -> Branch( "sumeT_hcalin_org", &m_sumeT_hcalin_org, "sumeT_hcalin_org/F" );
  m_tree -> Branch( "sumeT_hcalout_org", &m_sumeT_hcalout_org, "sumeT_hcalout_org/F" );

  if ( m_do_rho ) 
  {
    m_tree -> Branch( "rhoM_cemc", &m_rhoM_cemc, "rhoM_cemc/F" );
    m_tree -> Branch( "rhoM_hcalin", &m_rhoM_hcalin, "rhoM_hcalin/F" );
    m_tree -> Branch( "rhoM_hcalout", &m_rhoM_hcalout, "rhoM_hcalout/F" );
    m_tree -> Branch( "rhoM", &m_rhoM, "rhoM/F" );

    m_tree -> Branch( "rhoA_cemc", &m_rhoA_cemc, "rhoA_cemc/F" );
    m_tree -> Branch( "rhoA_hcalin", &m_rhoA_hcalin, "rhoA_hcalin/F" );
    m_tree -> Branch( "rhoA_hcalout", &m_rhoA_hcalout, "rhoA_hcalout/F" );
    m_tree -> Branch( "rhoA", &m_rhoA, "rhoA/F" );
  }
  if ( m_do_towerbkgd ) 
  {
    m_tree -> Branch( "sub1_towerbkgd_cemc", m_sub1_towerbkgd_cemc, Form("sub1_towerbkgd_cemc[%d]/F", k_ieta) );
    m_tree -> Branch( "sub1_towerbkgd_hcalin", m_sub1_towerbkgd_hcalin, Form("sub1_towerbkgd_hcalin[%d]/F", k_ieta) );
    m_tree -> Branch( "sub1_towerbkgd_hcalout", m_sub1_towerbkgd_hcalout, Form("sub1_towerbkgd_hcalout[%d]/F", k_ieta) );
    m_tree -> Branch( "sub2_towerbkgd_cemc", m_sub2_towerbkgd_cemc, Form("sub2_towerbkgd_cemc[%d]/F", k_ieta) );
    m_tree -> Branch( "sub2_towerbkgd_hcalin", m_sub2_towerbkgd_hcalin, Form("sub2_towerbkgd_hcalin[%d]/F", k_ieta) );
    m_tree -> Branch( "sub2_towerbkgd_hcalout", m_sub2_towerbkgd_hcalout, Form("sub2_towerbkgd_hcalout[%d]/F", k_ieta) );
    m_tree -> Branch( "sub1_v2", &m_sub1_v2, "sub1_v2/F" );
    m_tree -> Branch( "sub2_v2", &m_sub2_v2, "sub2_v2/F" );
    m_tree -> Branch( "sub1_flowfaliure", &m_sub1_flowfaliure, "sub1_flowfaliure/I" );
    m_tree -> Branch( "sub2_flowfaliure", &m_sub2_flowfaliure, "sub2_flowfaliure/I" );
    m_tree -> Branch( "sub1_psi2", &m_sub1_psi2, "sub1_psi2/F" );
    m_tree -> Branch( "sub2_psi2", &m_sub2_psi2, "sub2_psi2/F" );
  }
  if ( !m_sub1jet_node.empty() ) 
  {
    m_tree  -> Branch( "jet_E", &m_sub1_jet_E );
    m_tree  -> Branch( "jet_phi", &m_sub1_jet_phi );
    m_tree  -> Branch( "jet_eta", &m_sub1_jet_eta );
    m_tree  -> Branch( "jet_pT", &m_sub1_jet_pT );
    m_tree  -> Branch( "jet_unsub_pT", &m_sub1_jet_unsub_pT );
    m_tree  -> Branch( "jet_unsub_E", &m_sub1_jet_unsub_E );
  }
  if ( !m_truthjet_node.empty() ) 
  {
    m_tree  -> Branch("truth_jet_E", &m_truth_jet_E );
    m_tree  -> Branch("truth_jet_phi", &m_truth_jet_phi );
    m_tree  -> Branch("truth_jet_eta", &m_truth_jet_eta );
    m_tree  -> Branch("truth_jet_pT", &m_truth_jet_pT );
  }
  if ( !m_rawjet_node.empty() ) 
  {
    m_tree  -> Branch("raw_jet_E", &m_raw_jet_E );
    m_tree  -> Branch("raw_jet_phi", &m_raw_jet_phi );
    m_tree  -> Branch("raw_jet_eta", &m_raw_jet_eta );
    m_tree  -> Branch("raw_jet_pT", &m_raw_jet_pT );
  }
  if ( !m_multjet_node.empty() ) 
  {
    m_tree  -> Branch("mult_jet_E", &m_mult_jet_E );
    m_tree  -> Branch("mult_jet_phi", &m_mult_jet_phi );
    m_tree  -> Branch("mult_jet_eta", &m_mult_jet_eta );
    m_tree  -> Branch("mult_jet_pT", &m_mult_jet_pT ); 
  }
  if ( !m_areajet_node.empty() ) 
  {
    m_tree  -> Branch("area_jet_E", &m_area_jet_E );
    m_tree  -> Branch("area_jet_phi", &m_area_jet_phi );
    m_tree  -> Branch("area_jet_eta", &m_area_jet_eta );
    m_tree  -> Branch("area_jet_pT", &m_area_jet_pT ); 
  }

  if ( !m_eventplane_node.empty() ) 
  {
    m_tree -> Branch( "psi_shifted_NS", &m_psi_shifted_NS);
    m_tree -> Branch( "psi_shifted_S", &m_psi_shifted_S);
    m_tree -> Branch( "psi_shifted_N", &m_psi_shifted_N);
    m_tree -> Branch( "psi_NS", &m_psi_NS);
    m_tree -> Branch( "psi_S", &m_psi_S);
    m_tree -> Branch( "psi_N", &m_psi_N);
  }
  if ( !m_g4truth_node.empty() )
  {
    m_tree -> Branch( "g4truth_zvtx", &m_g4truth_zvtx, "g4truth_zvtx/F" );
    m_tree -> Branch( "m_g4truth_v2reco", &m_g4truth_v2reco );
    m_tree -> Branch( "m_g4truth_v3reco", &m_g4truth_v3reco );
    m_tree -> Branch( "m_g4truth_v4reco", &m_g4truth_v4reco );
    m_tree -> Branch( "m_g4truth_v5reco", &m_g4truth_v5reco );
    m_tree -> Branch( "m_g4truth_v6reco", &m_g4truth_v6reco );
  }
  


  if ( Verbosity () > 0 )
  {
    std::cout << "SimTree::Init - done" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int SimTree::process_event( PHCompositeNode *topNode )
{

  m_event_id++; 

  
  if(Verbosity() > 1) 
  {
    std::cout << "SimTree::process_event - Process event " << m_event_id << std::endl;
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

  if ( !m_sub1jet_node.empty() ) 
  { // get sub1 jet info
    auto res = GetSub1JetInfo(topNode);
    if ( res != Fun4AllReturnCodes::EVENT_OK ) { return res; }
  }

  if ( !m_truthjet_node.empty() ) 
  { // get truth jet info
    auto res = GetTruthJetInfo(topNode);
    if ( res != Fun4AllReturnCodes::EVENT_OK ) { return res; }
  }
  if ( !m_eventplane_node.empty() ) 
  { // get eventplane info
    auto res = GetEventPlaneInfo(topNode);
    if ( res != Fun4AllReturnCodes::EVENT_OK ) { return res; }
  }
  if ( !m_g4truth_node.empty() ) 
  { // get g4 truth info
    auto res = GetG4TruthInfo(topNode);
    if ( res != Fun4AllReturnCodes::EVENT_OK ) { return res; }
  }

  if ( !m_rawjet_node.empty() || !m_multjet_node.empty() || !m_areajet_node.empty() )
  { // get raw jet info
    auto res = GetRawJetInfo(topNode);
    if ( res != Fun4AllReturnCodes::EVENT_OK ) { return res; }
  }


  GetCaloInfo(topNode);

  // fill tree
  m_tree->Fill();
  
  return Fun4AllReturnCodes::EVENT_OK;

}

int SimTree::End( PHCompositeNode * /*topNode*/ )
{
  
  if( Verbosity() > 0 ) 
  {
    std::cout << "SimTree::EndRun - End run " << std::endl;
    std::cout << "SimTree::EndRun - Writing to " << m_output_filename << std::endl;
  }

  PHTFileServer::get().cd(m_output_filename); 

  m_tree->Write();
  

  if ( Verbosity() > 0 ) 
  {
    std::cout << "SimTree::EndRun - Writing run tree" << std::endl;
  }
  PHTFileServer::get().close();
 
  if ( Verbosity () > 0 ) 
  {
    std::cout << "SimTree::EndRun - done" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int SimTree::GetZvtx( PHCompositeNode *topNode )
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

int SimTree::GetCentInfo( PHCompositeNode *topNode )
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

int SimTree::GetEventHeaderInfo( PHCompositeNode *topNode )
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

int SimTree::GetSub1JetInfo( PHCompositeNode *topNode)
{

 
  ResetSub1Jet();
  ResetTowerBkgd();

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
  auto jets = findNode::getClass<JetContainer>( topNode, m_sub1jet_node );
  if ( !jets )
  {
    std::cout << PHWHERE << " Input node " << m_sub1jet_node << " Node missing, doing nothing." << std::endl;
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

  m_sub1_v2 = tower_background_sub1->get_v2();
  m_sub2_v2 = tower_background_sub2->get_v2();
  m_sub1_flowfaliure = tower_background_sub1->get_flow_failure_flag();
  m_sub2_flowfaliure = tower_background_sub2->get_flow_failure_flag();
  m_sub1_psi2 = tower_background_sub1->get_Psi2();
  m_sub2_psi2 = tower_background_sub2->get_Psi2();
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


  } 

  if ( Verbosity() > 0 ) 
  {
    std::cout << PHWHERE << " - Found " << jets->size() << " sub1 jets in node " << m_sub1jet_node << std::endl;
  }


  
  return Fun4AllReturnCodes::EVENT_OK;

}

int SimTree::GetTruthJetInfo( PHCompositeNode *topNode)
{

 
  ResetTruthJet();

  // get truth jets
  auto jets = findNode::getClass<JetContainer>( topNode, m_truthjet_node );
  if ( !jets )
  {
    std::cout << PHWHERE << " Input node " << m_truthjet_node << " Node missing, doing nothing." << std::endl;
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
    std::cout << PHWHERE << " - Found " << jets->size() << " truth jets in node " << m_truthjet_node << std::endl;
  }
  
  return Fun4AllReturnCodes::EVENT_OK;

}

int SimTree::GetRawJetInfo( PHCompositeNode *topNode)
{

 
  ResetRawJet();
  ResetMultJet();
  ResetAreaJet();
  ResetRho();
  
  auto rhoM_cemc = findNode::getClass<TowerRhov1>(topNode, "TowerRho_MULT_CEMC" );
  auto rhoM_hcalin = findNode::getClass<TowerRhov1>(topNode, "TowerRho_MULT_HCALIN" );
  auto rhoM_hcalout = findNode::getClass<TowerRhov1>(topNode, "TowerRho_MULT_HCALOUT" );
  auto rhoA_cemc = findNode::getClass<TowerRhov1>(topNode, "TowerRho_AREA_CEMC" );
  auto rhoA_hcalin = findNode::getClass<TowerRhov1>(topNode, "TowerRho_AREA_HCALIN" );
  auto rhoA_hcalout = findNode::getClass<TowerRhov1>(topNode, "TowerRho_AREA_HCALOUT" );
  auto rhoM = findNode::getClass<TowerRhov1>(topNode, "TowerRho_MULT" );
  auto rhoA = findNode::getClass<TowerRhov1>(topNode, "TowerRho_AREA" );
  if ( !rhoM_cemc || !rhoM_hcalin || !rhoM_hcalout || !rhoM ) {
    std::cout << PHWHERE << " One of the following nodes is missing: TowerRho_MULT_CEMC, TowerRho_MULT_HCALIN, TowerRho_MULT_HCALOUT, TowerRho_MULT. Skipping rho filling." << std::endl;
  }
  else {
    m_rhoM_cemc = rhoM_cemc->get_rho();
    m_rhoM_hcalin = rhoM_hcalin->get_rho();
    m_rhoM_hcalout = rhoM_hcalout->get_rho();
    m_rhoM = rhoM->get_rho();
  }
  if ( !rhoA_cemc || !rhoA_hcalin || !rhoA_hcalout || !rhoA ) {
    std::cout << PHWHERE << " One of the following nodes is missing: TowerRho_AREA_CEMC, TowerRho_AREA_HCALIN, TowerRho_AREA_HCALOUT, TowerRho_AREA. Skipping rho filling." << std::endl;
  }
  else {
    m_rhoA_cemc = rhoA_cemc->get_rho();
    m_rhoA_hcalin = rhoA_hcalin->get_rho();
    m_rhoA_hcalout = rhoA_hcalout->get_rho();
    m_rhoA = rhoA->get_rho();
  }


  // get raw jets
  auto jets = findNode::getClass<JetContainer>( topNode, m_rawjet_node );
  if ( !jets )
  {
    std::cout << PHWHERE << " Input node " << m_rawjet_node << " Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT; 
  }

  for ( auto jet : *jets )
  {

    float this_pt = jet->get_pt();
    float this_eta = jet->get_eta();
    float this_phi = jet->get_phi();
    float this_e = jet->get_e();

    m_raw_jet_E.push_back(this_e);
    m_raw_jet_eta.push_back(this_eta);
    m_raw_jet_phi.push_back(this_phi);
    m_raw_jet_pT.push_back(this_pt);

  } // end loop over jets

  if ( !m_multjet_node.empty() ) 
  {
    auto multjets = findNode::getClass<JetContainer>( topNode, m_multjet_node );
    if (!multjets )
    {
      std::cout << PHWHERE << " Input node " << m_multjet_node << " Node missing, skipping mult jets." << std::endl;
    }
    else 
    {
      for ( auto jet : *multjets )
      {

        float this_pt = jet->get_pt();
        float this_eta = jet->get_eta();
        float this_phi = jet->get_phi();
        float this_e = jet->get_e();

        m_mult_jet_E.push_back(this_e);
        m_mult_jet_eta.push_back(this_eta);
        m_mult_jet_phi.push_back(this_phi);
        m_mult_jet_pT.push_back(this_pt);

      } // end loop over mult jets
    }
  }
  if ( !m_areajet_node.empty() ) 
  {
    auto areajets = findNode::getClass<JetContainer>( topNode, m_areajet_node );
    if (!areajets )
    {
      std::cout << PHWHERE << " Input node " << m_areajet_node << " Node missing, skipping area jets." << std::endl;
    }
    else 
    {
      for ( auto jet : *areajets )
      {

        float this_pt = jet->get_pt();
        float this_eta = jet->get_eta();
        float this_phi = jet->get_phi();
        float this_e = jet->get_e();

        m_area_jet_E.push_back(this_e);
        m_area_jet_eta.push_back(this_eta);
        m_area_jet_phi.push_back(this_phi);
        m_area_jet_pT.push_back(this_pt);

      } // end loop over area jets
    }
  }

  if ( Verbosity() > 0 ) 
  {
    std::cout << PHWHERE << " - Found " << jets->size() << " raw jets in node " << m_rawjet_node << std::endl;
  }
  
  return Fun4AllReturnCodes::EVENT_OK;

}

void SimTree::SumCaloE( PHCompositeNode *topNode, const std::string &towerinfo_node, const std::string &geo_node, const float zvrtx, float &sum_e )
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

int SimTree::GetCaloInfo( PHCompositeNode *topNode )
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

int SimTree::GetG4TruthInfo( PHCompositeNode *topNode )
{
  ResetG4Truth();

  PHG4TruthInfoContainer * truthinfo = findNode::getClass<PHG4TruthInfoContainer>( topNode, m_g4truth_node );
  if ( !truthinfo )  {
    std::cout << PHWHERE << " Input node "<< m_g4truth_node << " Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  std::vector<float> mpsi{m_psi2, m_psi3, m_psi4, m_psi5, m_psi6};
  PHG4TruthInfoContainer::ConstRange range = truthinfo->GetSPHENIXPrimaryParticleRange();
  PHG4VtxPoint  * gvertex = truthinfo->GetPrimaryVtx(truthinfo->GetPrimaryVertexIndex());
  m_g4truth_zvtx = gvertex->get_z();
  std::vector<float> _sum_w_vn(5, 0.0);
  float _sum_w = 0.0;
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)  {
    auto g4particle = iter->second;
    if ( truthinfo-> isEmbeded ( g4particle->get_track_id() ) != 0) continue; // skip embedded particles
    auto t = new TLorentzVector( g4particle->get_px(), g4particle->get_py(), g4particle->get_pz(), g4particle->get_e() );
    if ( t -> Pt() < 0.4 ) continue; // apply pT cut to truth particles
    if ( std::abs(t -> Eta()) > 1.1 ) continue;
    for ( int i = 0; i < 5; i++ ) _sum_w_vn[i] += ( t->Pt() * TMath::Cos( static_cast<double>(i+2) * (t->Phi() - mpsi[i]) ) );
    _sum_w += t->Pt();
  }

  m_g4truth_v2reco = _sum_w_vn[0] / _sum_w;
  m_g4truth_v3reco = _sum_w_vn[1] / _sum_w;
  m_g4truth_v4reco = _sum_w_vn[2] / _sum_w;
  m_g4truth_v5reco = _sum_w_vn[3] / _sum_w;
  m_g4truth_v6reco = _sum_w_vn[4] / _sum_w; 

  if ( Verbosity() > 1 ) 
  {
    std::cout << PHWHERE << " - G4TruthInfo: zvtx = " << m_g4truth_zvtx << ", v2reco = " << m_g4truth_v2reco << ", v3reco = " << m_g4truth_v3reco << ", v4reco = " << m_g4truth_v4reco << ", v5reco = " << m_g4truth_v5reco << ", v6reco = " << m_g4truth_v6reco << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int SimTree::GetEventPlaneInfo( PHCompositeNode *topNode )
{
  ResetEventPlane();
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
              << "S: (" <<  m_psi_shifted_S[0] << ", " << m_psi_S[0] << "), "
              << "N: (" <<  m_psi_shifted_N[0] << ", " << m_psi_N[0] << ")" 
              << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}






   