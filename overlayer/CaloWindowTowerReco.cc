#include "CaloWindowTowerReco.h"
#include "CaloWindowMapv1.h"

#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>  
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv2.h>
#include <calobase/TowerInfoDefs.h>
#include <calobase/RawTowerDefs.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/recoConsts.h>

// standard includes
#include <cstdlib> 
#include <vector>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <iostream>


CaloWindowTowerReco::CaloWindowTowerReco(const std::string &name)
  : SubsysReco(name)
{
}

CaloWindowTowerReco::~CaloWindowTowerReco()
{
  m_inputs.clear();
  m_geom_names.clear();
  m_srcs.clear();
  m_window_names.clear();
}

int CaloWindowTowerReco::Init(PHCompositeNode *topNode)
{
  
  // check all user settings
  if ( m_inputs.empty() ) {
    std::cout << PHWHERE << "No input nodes are set, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if ( m_geom_names.empty() ) {
    std::cout << PHWHERE << "No geometry nodes are set, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if ( m_srcs.empty() ) {
    std::cout << PHWHERE << "No source nodes are set, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if ( m_window_names.empty() ) {
    for ( auto &src : m_srcs ) {
      std::string name = UEDefs::GetCaloTowerNode(src, "");
      m_window_names.push_back(m_window_prefix+name);
    }
  }

  if (Verbosity() > 0 ) {
    std::cout << "CaloWindowTowerReco::Init - settings" << std::endl;
    for ( unsigned int in = 0; in < m_inputs.size(); in++ ) {
      std::cout << "  input node = " << m_inputs[in] << " geometry node = " << m_geom_names[in] << " source = " << m_srcs[in] << " window node = " << m_window_names[in] << std::endl;
    }
  }
  
  return CreateNodes(topNode);

}

int CaloWindowTowerReco::process_event(PHCompositeNode *topNode)
{

  auto vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if ( !vertexmap  || vertexmap->empty() ) {
    std::cout << PHWHERE << "GlobalVertexMap node is missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
 
  auto vtx = vertexmap->begin()->second;
  float z_vrtx {NAN};
  if ( vtx ) {
    z_vrtx = vtx->get_z();
  } 
  if ( std::isnan(z_vrtx) || z_vrtx > 1e3) {
    static bool once = true;
    if (once) {
      once = false;
      std::cout << PHWHERE << "vertex is " << z_vrtx << ". Drop all tower inputs (further vertex warning will be suppressed)." << std::endl;
    }
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  
  for ( unsigned int in = 0; in < m_inputs.size(); in++ ) {
      
    TowerInfoContainer * towerinfos {nullptr};
    RawTowerGeomContainer * geom {nullptr};
    RawTowerDefs::CalorimeterId geocaloid {RawTowerDefs::CalorimeterId::NONE};

    towerinfos = findNode::getClass<TowerInfoContainer>(topNode, m_inputs[in]);
    geom = findNode::getClass<RawTowerGeomContainer>(topNode, m_geom_names[in]);
    if ( !towerinfos || !geom ) {
      std::cout << PHWHERE << "Can't find TowerInfoContainer or RawTowerGeomContainer node " << m_inputs[in] << " " << m_geom_names[in] << std::endl;
      exit(-1); // fatal error
    }

    if ( m_geom_names[in] == "TOWERGEOM_CEMC" ) {
      geocaloid = RawTowerDefs::CalorimeterId::CEMC;
    }
    else if ( m_geom_names[in] == "TOWERGEOM_HCALIN" ) {
      geocaloid = RawTowerDefs::CalorimeterId::HCALIN;
    }
    else if ( m_geom_names[in] == "TOWERGEOM_HCALOUT" ) {
      geocaloid = RawTowerDefs::CalorimeterId::HCALOUT;
    } else {
      std::cout << PHWHERE << "Unknown calorimeter name " << m_geom_names[in] << std::endl;
      exit(-1); // fatal error
    } // end of calorimeter id

    CaloWindowMap *window = findNode::getClass<CaloWindowMap>(topNode, m_window_names[in]);
    if ( !window ) {
      std::cout << PHWHERE << "CaloWindowMap node " << m_window_names[in] << " is missing, doing nothing." << std::endl;
      exit(-1); // fatal error
    }

    window -> Reset();
    window -> set_src(m_srcs[in]);
    if ( m_geom_names[in].find("CEMC") != std::string::npos ) {
      window -> set_nphi_neta(CaloWindowTowerReco::nphi_emcal , CaloWindowTowerReco::neta_emcal);
    } else {
      window -> set_nphi_neta(CaloWindowTowerReco::nphi_ihcal , CaloWindowTowerReco::neta_ihcal);
    }
    window -> clear_towers(); // sets length of m_towers to m_nphi*m_neta and sets all elements to 0

    unsigned int nchannels = towerinfos->size();
    for (unsigned int channel = 0; channel < nchannels; channel++) {
      auto tower = towerinfos->get_tower_at_channel(channel);
      assert(tower);
      unsigned int calokey = towerinfos->encode_key(channel);
      unsigned int ieta = towerinfos->getTowerEtaBin(calokey);
      unsigned int iphi = towerinfos->getTowerPhiBin(calokey);

      bool is_masked = tower->get_isHot() || tower->get_isNoCalib() || tower->get_isNotInstr() || tower->get_isBadChi2() || std::isnan(tower->get_energy());
      
      const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(geocaloid, ieta, iphi);
      auto tower_geom = geom->get_tower_geometry(key);
      assert(tower_geom);

      // find if input has RETOWER and CEMC in the name
      double r = tower_geom->get_center_radius();
      if (m_inputs[in].find("RETOWER") != std::string::npos && m_geom_names[in].find("CEMC") != std::string::npos) {
        auto EMCal_geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
        assert(EMCal_geom);
        const RawTowerDefs::keytype EMCal_key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, 0, 0);
        auto EMCal_tower_geom = EMCal_geom->get_tower_geometry(EMCal_key);
        assert(EMCal_tower_geom);
        r = EMCal_tower_geom->get_center_radius();
      }

      double towereta = tower_geom->get_eta();
      double z0 = sinh(towereta) * r;
      double z = z0 - z_vrtx;
      double eta = asinh(z / r);  
      double pt = tower->get_energy() / cosh(eta);

      window -> add_tower(iphi, ieta, pt, is_masked);
    }

  } // end of loop over inputs

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloWindowTowerReco::CreateNodes( PHCompositeNode *topNode )
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode) {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHCompositeNode *windowNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "CALOWINDOWS"));
  if ( !windowNode ){
    windowNode = new PHCompositeNode("CALOWINDOWS");
    dstNode->addNode(windowNode);
  }

  for ( auto &window_name : m_window_names ) {
    CaloWindowMap *window = findNode::getClass<CaloWindowMap>(topNode, window_name);
    if ( !window ) {
      window = new CaloWindowMapv1();
      PHIODataNode<PHObject> *node = new PHIODataNode<PHObject>(window, window_name, "PHObject");
      windowNode->addNode(node);
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}




