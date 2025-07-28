#include "RandomConeTowerReco.h"
#include "RandomConev1.h"

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

#include <TRandom3.h>
#include <TTimeStamp.h>

#include <jetbase/Jet.h>
#include <jetbase/JetContainer.h>
#include <jetbase/JetContainerv1.h>
#include <jetbase/JetMap.h>
#include <jetbase/JetMapv1.h>

// standard includes
#include <cstdlib> 
#include <vector>
#include <cmath>
#include <algorithm>
#include <cassert>



RandomConeTowerReco::RandomConeTowerReco(const std::string &name)
  : SubsysReco(name)
{
}

RandomConeTowerReco::~RandomConeTowerReco()
{
  m_inputs.clear();
  m_geom_names.clear();
  m_srcs.clear();
  if( m_random ) { delete m_random;}
}

int RandomConeTowerReco::Init(PHCompositeNode *topNode)
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

  if ( m_seed == 0 ) { // make sure random seed is set
    auto rc = recoConsts::instance(); // try to get random seed from flags
    if ( !rc->FlagExist("PPG04RANDOMSEED") ) { // if not set, use time
      std::cout << PHWHERE << "PPG04RANDOMSEED flag not set, using time." << std::endl;  
      auto t = new TTimeStamp();
      m_seed = static_cast<unsigned int>(t->GetNanoSec());
      delete t;
    } else { // if set, use the flag
      m_seed = static_cast<unsigned int>(recoConsts::instance()->get_IntFlag("PPG04RANDOMSEED"));
    }
  }

  m_random = new TRandom3(m_seed);
  if ( !m_random ) { // make sure random number generator is created
    std::cout << PHWHERE << "m_random is not set, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (m_output_node.empty()) {
    if (Verbosity()) { std::cout << PHWHERE << "output node is not set, setting default." << std::endl; }
    m_output_node = "RandomCone_r" + std::to_string(static_cast<int>(m_R*10));
    if ( m_avoid_lead_jet ) {
      m_output_node += "_avoidleadjet";
    }
  }

  if ( m_max_abs_eta < 0) {
    if( Verbosity()) { std::cout << PHWHERE << "max abs eta is not set, setting defualt 1.1 - R" << std::endl; }
    m_max_abs_eta = 1.1 - m_R;
  }

  if (Verbosity() > 0 ) {
    std::cout << "RandomConeTowerReco::Init - settings" << std::endl;
    std::cout << "  R = " << m_R << std::endl;
    std::cout << "  output node = " << m_output_node << std::endl;
    std::cout << "  max abs eta = " << m_max_abs_eta << std::endl;
    std::cout << "  masked threshold = " << m_masked_threshold << std::endl;
    std::cout << "  seed = " << m_seed << std::endl;
    if ( m_avoid_lead_jet ) {
      std::cout << "  avoid leading jet = " << m_lead_jet_node << " dR = " << m_lead_jet_dR << std::endl;
    }
  }
  
  return CreateNode(topNode);

}

int RandomConeTowerReco::process_event(PHCompositeNode *topNode)
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

  auto cone = findNode::getClass<RandomConev1>(topNode, m_output_node);
  if ( !cone ) {
    std::cout << PHWHERE << "RandomConev1 node " << m_output_node << " is missing, doing nothing." << std::endl;
    exit(-1); // fatal error
  }

  int ntries = 0;
  while (true) {

    float cone_eta = NAN, cone_phi = NAN;
    GetConeAxis(topNode, cone_eta, cone_phi);
    if ( std::isnan(cone_eta) || std::isnan(cone_phi) ) {
      std::cout << PHWHERE << "RandomConeTowerReco::process_event - cone axis is NAN, aborting event" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

    cone->Reset();
    cone->set_id(ntries);
    cone->set_R(m_R);
    cone->set_eta(cone_eta);
    cone->set_phi(cone_phi);
 

    for ( unsigned int in = 0; in < m_inputs.size(); in++){
      
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

      unsigned int nchannels = towerinfos->size();
      for (unsigned int channel = 0; channel < nchannels; channel++) {
        auto tower = towerinfos->get_tower_at_channel(channel);
        assert(tower);
        unsigned int calokey = towerinfos->encode_key(channel);
        int ieta = towerinfos->getTowerEtaBin(calokey);
        int iphi = towerinfos->getTowerPhiBin(calokey);
        const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(geocaloid, ieta, iphi);
        bool is_masked = tower->get_isHot() || tower->get_isNoCalib() || tower->get_isNotInstr() || tower->get_isBadChi2() || std::isnan(tower->get_energy());
        RawTowerGeom *tower_geom = geom->get_tower_geometry(key);
        assert(tower_geom);
        // find if input has RETOWER and CEMC in the name
        double r = tower_geom->get_center_radius();
        if (m_inputs[in].find("RETOWER") != std::string::npos && m_geom_names[in].find("CEMC") != std::string::npos) {
          auto EMCal_geom = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
          assert(EMCal_geom);
          const RawTowerDefs::keytype EMCal_key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, 0, 0);
          RawTowerGeom *EMCal_tower_geom = EMCal_geom->get_tower_geometry(EMCal_key);
          assert(EMCal_tower_geom);
          r = EMCal_tower_geom->get_center_radius();
        }

        double phi = atan2(tower_geom->get_center_y(), tower_geom->get_center_x());
        double towereta = tower_geom->get_eta();
        double z0 = sinh(towereta) * r;
        double z = z0 - z_vrtx;
        double eta = asinh(z / r);  // eta after shift from vertex
        double pt = tower->get_energy() / cosh(eta);

        if (std::isnan(pt)) {
          pt = 0;
        }

        double delta_eta = eta - cone_eta;
        double delta_phi = phi - cone_phi;
        // wrap around
        if (delta_phi > M_PI) {delta_phi -= 2 * M_PI;}
        if (delta_phi < -M_PI){ delta_phi += 2 * M_PI;}
        double dR = std::sqrt( (delta_eta * delta_eta) + (delta_phi * delta_phi) );
        
        if (dR < m_R) {
          cone->add_tower(m_srcs[in], pt, channel, is_masked);
        }
      }
    } // end of loop over inputs

    // std::cout << "RandomConeTowerReco::process_event : cone pt = " << cone->get_pt() << " eta = " << cone->get_eta() << " phi = " << cone->get_phi() << std::endl;

    // check masked threshold
    bool pass = true;
    for ( unsigned int in = 0; in < m_inputs.size(); in++) {
      if ( cone->masked_fraction(m_srcs[in]) > m_masked_threshold ) {
        pass = false;
      }
    } 

    if ( pass ) { break; }
    ntries++;
    if ( ntries > 10 ) {
      std::cout << "RandomConeTowerReco::process_event - cone failed masked threshold, trying again (" << ntries << ")" << std::endl;
      std::cout << "RandomConeTowerReco::process_event - aborting event" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    if (Verbosity() > 1) {
      std::cout << "RandomConeTowerReco::process_event - cone failed masked threshold, trying again (" << ntries << ")" << std::endl;
    }

  } // end of while loop

  // reset lead jet axis
  if ( m_avoid_lead_jet ) {
    m_lead_jet_eta = NAN;
    m_lead_jet_phi = NAN;
  }

  return Fun4AllReturnCodes::EVENT_OK;

}

int RandomConeTowerReco::CreateNode( PHCompositeNode *topNode )
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode) {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHCompositeNode *coneNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RANDOMCONE"));
  if ( !coneNode ){
    coneNode = new PHCompositeNode("RANDOMCONE");
    dstNode->addNode(coneNode);
  }

  RandomCone *cone = new RandomConev1();
  PHIODataNode<PHObject> *node = new PHIODataNode<PHObject>(cone, m_output_node, "PHObject");
  coneNode->addNode(node);
  return Fun4AllReturnCodes::EVENT_OK;
}

void RandomConeTowerReco::GetConeAxis(PHCompositeNode *topNode, float &cone_eta, float &cone_phi)
{

  if( m_avoid_lead_jet ) {
    if ( std::isnan(m_lead_jet_eta) || std::isnan(m_lead_jet_phi) ) {
      float lead_eta = NAN, lead_phi = NAN, lead_pt = -1; 
      auto jetcont = findNode::getClass<JetContainer>(topNode, m_lead_jet_node);
      if( !jetcont ) {
        // could be a jetmap
        auto jetmap = findNode::getClass<JetMap>(topNode, m_lead_jet_node);
        if( !jetmap ) {
          std::cout << PHWHERE << "Could not find leading jet node " << m_lead_jet_node << std::endl;
          std::cout << PHWHERE << "Disabling avoid leading jet mode for this event" << std::endl;
          m_avoid_lead_jet = false;          
          cone_eta = m_random->Uniform(-m_max_abs_eta, m_max_abs_eta);
          cone_phi = m_random->Uniform(-M_PI, M_PI);
          return;
        } // end of !jetmap
        
        for(JetMap::Iter iter = jetmap->begin(); iter != jetmap->end(); ++iter) {
          auto jet = iter->second;
          if(jet->get_pt() > lead_pt) {
            lead_pt = jet->get_pt();
            lead_eta = jet->get_eta();
            lead_phi = jet->get_phi();
          }
        } // end of loop over jets
        
      } else { // jet container
        for ( unsigned int ijet = 0; ijet < jetcont->size(); ijet++) {
          auto jet = jetcont->get_jet(ijet);
          if ( jet->get_pt() > lead_pt ) {
            lead_pt = jet->get_pt();
            lead_eta = jet->get_eta();
            lead_phi = jet->get_phi();
          }
        } // end of loop over jets
      } // end of jet container

      if ( lead_pt < 0 || std::isnan(lead_eta) || std::isnan(lead_phi) ) {
        std::cout << PHWHERE << "Could not find leading jet" << std::endl;
        cone_eta = m_random->Uniform(-m_max_abs_eta, m_max_abs_eta);
        cone_phi = m_random->Uniform(-M_PI, M_PI);
        return;
      }

      m_lead_jet_eta = lead_eta;
      m_lead_jet_phi = lead_phi;
    } // end of lead jet not set

    while(true) {
      // set jet_phi to 0 and select phi from dR - (2pi-dR) range
      float dphi = m_random->Uniform(m_lead_jet_dR, 2*M_PI - m_lead_jet_dR);
      cone_phi = m_lead_jet_phi + dphi;
      if(cone_phi > M_PI) cone_phi -= 2*M_PI;
      if(cone_phi < -M_PI) cone_phi += 2*M_PI;
      dphi = cone_phi - m_lead_jet_phi;
      if(dphi > M_PI) dphi -= 2*M_PI;
      if(dphi < -M_PI) dphi += 2*M_PI;

      cone_eta = m_random->Uniform(-m_max_abs_eta, m_max_abs_eta);

      float deta = cone_eta - m_lead_jet_eta;
      float dR = std::sqrt((deta*deta) + (dphi*dphi));
      if(dR > m_lead_jet_dR) { return; }
    } // end of while loop
  } else { // standard random cone
    cone_eta = m_random->Uniform(-m_max_abs_eta, m_max_abs_eta);
    cone_phi = m_random->Uniform(-M_PI, M_PI);
  }

  return;
}




