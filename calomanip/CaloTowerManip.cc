#include "CaloTowerManip.h"

#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>  // for TowerInfo
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv2.h>
#include <calobase/TowerInfoDefs.h>

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

#include <algorithm>
#include <vector>
#include <cstdlib>

#include <TRandom3.h>
#include <TTimeStamp.h>

CaloTowerManip::CaloTowerManip( const std::string &name )
  : SubsysReco(name)
{}

int CaloTowerManip::InitRun( PHCompositeNode * topNode )
{
  if ( m_input_node.empty() ) { // make sure input node is set
    std::cout << PHWHERE << "Input node not set, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  } // end of m_input_node

  if ( m_do_min_energy_cut && std::isnan(m_min_energy) ) { // make sure min energy is set if needed
    std::cout << PHWHERE << "Min energy not set, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  } // end of m_do_min_energy_cut

  if ( m_do_randomize_towers ) {// create random number generator if needed

    if ( m_random_seed == 0 ) { // make sure random seed is set
      auto rc = recoConsts::instance(); // try to get random seed from flags
      if ( !rc->FlagExist("PPG04RANDOMSEED") ) { // if not set, use time
        std::cout << PHWHERE << "PPG04RANDOMSEED flag not set, using time." << std::endl;  
        auto t = new TTimeStamp();
        m_random_seed = static_cast<unsigned int>(t->GetNanoSec());
        delete t;
      } else { // if set, use the flag
        m_random_seed = static_cast<unsigned int>(recoConsts::instance()->get_IntFlag("PPG04RANDOMSEED"));
      }
    }

    m_random = new TRandom3(m_random_seed);
    if ( !m_random ) { // make sure random number generator is created
      std::cout << PHWHERE << "m_random is not set, doing nothing." << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  } // end of m_random init

  if ( m_save_original_towers ) { // create output node if needed
    
    if ( m_output_node.empty() ) { // if output node is not set, create one
      // get calorimeter name
      m_output_node = m_input_node + "_ORIGINAL";
      std::cout << PHWHERE << "output node not set, using " << m_output_node << std::endl;
    }

    if ( CreateTowerInfoContainerCopy(topNode) != Fun4AllReturnCodes::EVENT_OK ) {
      std::cout << PHWHERE << "Failed to create tower info container copy, doing nothing." << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  } // end of m_save_original_towers

  if ( Verbosity() ) { // print out settings
    std::cout << PHWHERE << "input node = " << m_input_node << std::endl;
    if ( m_do_min_energy_cut ) {
      std::cout << PHWHERE << "Min energy = " << m_do_min_energy_cut << std::endl;
    }
    if ( m_do_randomize_towers ) {
      std::cout << PHWHERE << "Random seed = " << m_random_seed << std::endl;
      std::cout << PHWHERE << "Using new algo" << std::endl;
    }
    if ( m_save_original_towers ) {
      std::cout << PHWHERE << "output node = " << m_output_node << std::endl;
    }
  } // end of Verbosity()

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloTowerManip::process_event( PHCompositeNode * topNode )
{
  // get original tower info container
  auto manip_towerinfo = findNode::getClass<TowerInfoContainer>( topNode, m_input_node );
  if ( ! manip_towerinfo ) {
    std::cout << PHWHERE << " Input node " << m_input_node << " Node missing, doing nothing." << std::endl;
  }

  TowerInfoContainer * unmanip_towerinfo {nullptr};
  if ( m_save_original_towers ) {
    unmanip_towerinfo = findNode::getClass<TowerInfoContainer>( topNode, m_output_node );
    if ( ! unmanip_towerinfo ) {
      std::cout << PHWHERE << m_output_node << " Node missing, doing nothing." << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  } else {
    unmanip_towerinfo = (TowerInfoContainer *)manip_towerinfo->CloneMe();
  }
  unmanip_towerinfo->Reset();

  // copy original tower info to new tower info
  unsigned int ntowers = manip_towerinfo->size();
  // std::vector<unsigned int> tower_keys {};
  std::vector<unsigned int> unmasked_tower_keys {};
  std::vector<unsigned int> masked_tower_keys {};
  for (unsigned int channel = 0; channel < ntowers; channel++) {

    auto original_tower = manip_towerinfo->get_tower_at_channel(channel);
    auto copied_tower = unmanip_towerinfo->get_tower_at_channel(channel);
    copied_tower->copy_tower(original_tower);
    
    if (original_tower->get_isHot() 
        || original_tower->get_isNoCalib() 
        || original_tower->get_isNotInstr() 
        || original_tower->get_isBadChi2() 
        || std::isnan(original_tower->get_energy()))
    { // tower is masked
      masked_tower_keys.push_back(channel);
    } else {
      unmasked_tower_keys.push_back(channel);
    }
  }

  if ( Verbosity() > 1 ) {
    std::cout << PHWHERE << "input node = " << m_input_node << " output node = " << m_output_node << std::endl;
    std::cout << "Original tower info container: " << std::endl;
    manip_towerinfo->identify();
    std::cout << "Copied tower info container: " << std::endl;
    unmanip_towerinfo->identify();
  }


  // apply manipulations to tower ids
  if ( m_do_randomize_towers ) { // only randomize unmasked towers
    std::vector<unsigned int> unmasked_tower_keys_copy = unmasked_tower_keys;
    unsigned int ntowers_unmasked = unmasked_tower_keys.size();
    unmasked_tower_keys.clear();
    for (unsigned int i = 0; i < ntowers_unmasked; i++) {
      unsigned int j = static_cast<unsigned int>( m_random->Uniform( unmasked_tower_keys_copy.size() ) );
      unmasked_tower_keys.push_back(unmasked_tower_keys_copy.at(j));
      // remove tower from list
      unmasked_tower_keys_copy.erase(unmasked_tower_keys_copy.begin() + j);      
    }
    unmasked_tower_keys_copy.clear();
    unsigned int ntowers_now = unmasked_tower_keys.size();
    if ( ntowers_now != ntowers_unmasked ) {
      std::cout << PHWHERE << "Number of unmasked towers changed from " << ntowers_unmasked << " to " << ntowers_now << std::endl;
      exit(1);
    }
    // for (unsigned int iidx = unmasked_tower_keys.size() - 1; iidx > 0; --iidx) {
    //   unsigned int jidx = static_cast<unsigned int>(m_random->Uniform(iidx + 1));
    //   std::swap(unmasked_tower_keys[iidx], unmasked_tower_keys[jidx]);
    // }
  }

  unsigned int n_unasked_channels = 0;
  for (unsigned int channel = 0; channel < ntowers; channel++) {
    
    unsigned int target_channel;
    bool is_masked = std::find(masked_tower_keys.begin(), masked_tower_keys.end(), channel) != masked_tower_keys.end();
    if ( is_masked ) { // if tower is masked
      target_channel = channel; // do not modify masked towers
    } else {
      target_channel = unmasked_tower_keys.at(n_unasked_channels);
      n_unasked_channels++;
    }

    auto manip_tower = manip_towerinfo->get_tower_at_channel(channel);  // original tower
    auto unmanip_tower = unmanip_towerinfo->get_tower_at_channel(target_channel); // should be random if randomize is on
    manip_tower->copy_tower(unmanip_tower);    
    if ( is_masked ) {
      continue;
    }
    
    unsigned int key = manip_towerinfo->encode_key(channel);
    int ieta = manip_towerinfo->getTowerEtaBin(key);
    int iphi = manip_towerinfo->getTowerPhiBin(key);
    if ( m_do_min_energy_cut ) {
      float E = manip_tower->get_energy();
      if ( manip_tower->get_energy() < m_min_energy ) {
        if ( Verbosity() > 1 && E != 0 ) {
          std::cout << "removing tower with E = " << E << " at eta = " << ieta << " phi = " << iphi << std::endl;
        }
        manip_tower->set_energy(0);
      }
    } // end of m_do_min_energy_cut
  }  // end loop over channels

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloTowerManip::CreateTowerInfoContainerCopy( PHCompositeNode *topNode )
{
  
  // get original tower info container
  auto original_towerinfo = findNode::getClass<TowerInfoContainer>(topNode, m_input_node);
  if ( ! original_towerinfo ) {
    std::cout << PHWHERE << " Input node " << m_input_node << " Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // get copied tower info container
  auto copied_towerinfo = findNode::getClass<TowerInfoContainer>(topNode, m_output_node);
  if ( ! copied_towerinfo ) {
    // find parent node
    PHNodeIterator iter(topNode);
    auto dstNode = dynamic_cast< PHCompositeNode* >( iter.findFirst("PHCompositeNode", "DST") );
    if ( !dstNode ) {
      std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    PHNodeIterator dstiter(dstNode);
    PHNode * original_towerinfo_node = dstiter.findFirst("PHIODataNode", m_input_node);
    if ( !original_towerinfo_node ) { // probably redundant
      std::cout << PHWHERE << m_input_node << " Node missing, doing nothing." << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    auto DNode = dynamic_cast< PHDataNode<TowerInfoContainer>* >(original_towerinfo_node);
    if ( DNode ) {
      original_towerinfo = dynamic_cast<TowerInfoContainer*>(DNode->getData());
      if ( !original_towerinfo ) {
        std::cout << PHWHERE << "dynamic_cast<TowerInfoContainer*> failed, doing nothing." << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
      } // end of !original_towerinfo
    } // end of DNode

    auto parentNode = dynamic_cast<PHCompositeNode *>(original_towerinfo_node->getParent());
    if (!parentNode) {
      std::cout << PHWHERE << "parentNode is missing, doing nothing." << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    copied_towerinfo = dynamic_cast<TowerInfoContainer*>(original_towerinfo->CloneMe());
    PHIODataNode<PHObject> *newNode  = new PHIODataNode<PHObject>(copied_towerinfo, m_output_node, "PHObject");
    parentNode->addNode(newNode);
    
  } // end of !copied_towerinfo

  return Fun4AllReturnCodes::EVENT_OK;
} 

