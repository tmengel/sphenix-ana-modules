#include "OverlayToTTree.h"

#include "EmbedInfo.h"
#include "EmbedInfov1.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/PHTFileServer.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <phool/PHObject.h>

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeom.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/GlobalVertexMapv1.h>

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#include <HepMC/GenEvent.h>

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>

#include <jetbase/JetContainer.h>
#include <jetbase/Jet.h>

#include <TTree.h>
#include <TF1.h>
#include <TMath.h>

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cassert>


int OverlayToTTree::Init( PHCompositeNode * /*topNode*/ )
{

  PHTFileServer::get().open( _foutname, "RECREATE" );
  _tree = new TTree("T", "T");

  _tree -> Branch( "event", &_count, "event/I" );
  if ( _target_cent >= 0 )
  {
    _tree -> Branch( "target_cent", &_target_cent, "target_cent/I" );
  }
  _tree -> Branch( "jetpt_thres", &_jetpt_thres, "jetpt_thres/F" );
  _tree -> Branch( "cemc_scale", &_cemc_scale, "cemc_scale/F" );
  _tree -> Branch( "ihcal_scale", &_ihcal_scale, "ihcal_scale/F" );
  _tree -> Branch( "ohcal_scale", &_ohcal_scale, "ohcal_scale/F" );
  _tree -> Branch( "jetv2", &_jetv2, "jetv2/F" );
  _tree -> Branch( "v2_weight", &_v2_weight, "v2_weight/F" );
  _tree -> Branch( "dpsi2", &_dpsi2, "dpsi2/F" );
  _tree -> Branch( "emb_cent", &_emb_cent, "emb_cent/I" );
  _tree -> Branch( "emb_zvrtx", &_emb_zvrtx, "emb_zvrtx/F" );
  if ( !_is_data ) 
  {
    _tree -> Branch( "emb_b", &_emb_b, "emb_b/F" );
    _tree -> Branch( "emb_psi2", &_emb_psi2, "emb_psi2/F" );
  }
  _tree -> Branch( "emb_ohcal_sumet", &_emb_ohcal_sumet, "emb_ohcal_sumet/F" );
  _tree -> Branch( "emb_ihcal_sumet", &_emb_ihcal_sumet, "emb_ihcal_sumet/F" );
  _tree -> Branch( "emb_cemc_sumet", &_emb_cemc_sumet, "emb_cemc_sumet/F" );
  _tree -> Branch( "emb_ohcal_sumet_scaled", &_emb_ohcal_sumet_scaled, "emb_ohcal_sumet_scaled/F" );
  _tree -> Branch( "emb_ihcal_sumet_scaled", &_emb_ihcal_sumet_scaled, "emb_ihcal_sumet_scaled/F" );
  _tree -> Branch( "emb_cemc_sumet_scaled", &_emb_cemc_sumet_scaled, "emb_cemc_sumet_scaled/F" );
  _tree -> Branch( "pythia_z",&_pythia_z,"pythia_z/F");
  _tree -> Branch( "pythia_z_reco",&_pythia_z_reco,"pythia_z_reco/F");

  _tree -> Branch( "truth_jet_E",&_truth_jet_E);
  _tree -> Branch( "truth_jet_pt",&_truth_jet_pt);
  _tree -> Branch( "truth_jet_eta",&_truth_jet_eta);
  _tree -> Branch( "truth_jet_phi",&_truth_jet_phi);
  _tree -> Branch( "truth_jet_recoUE",&_truth_jet_recoUE);
  _tree -> Branch( "truth_jet_recoUE_cemc",&_truth_jet_recoUE_cemc);
  _tree -> Branch( "truth_jet_recoUE_ihcal",&_truth_jet_recoUE_ihcal);
  _tree -> Branch( "truth_jet_recoUE_ohcal",&_truth_jet_recoUE_ohcal);
  _tree -> Branch( "recoPythia_jet_E",&_recoPythia_jet_E);
  _tree -> Branch( "recoPythia_jet_pt",&_recoPythia_jet_pt);
  _tree -> Branch( "recoPythia_jet_eta",&_recoPythia_jet_eta);
  _tree -> Branch( "recoPythia_jet_phi",&_recoPythia_jet_phi);
  _tree -> Branch( "recoPythia_jet_EMfrac",&_recoPythia_jet_EMfrac);
  _tree -> Branch( "recoEmbed_jet_E",&_recoEmbed_jet_E);
  _tree -> Branch( "recoEmbed_jet_pt",&_recoEmbed_jet_pt);
  _tree -> Branch( "recoEmbed_jet_eta",&_recoEmbed_jet_eta);
  _tree -> Branch( "recoEmbed_jet_phi",&_recoEmbed_jet_phi);
  _tree -> Branch( "recoEmbed_jet_EMfrac",&_recoEmbed_jet_EMfrac);

  std::cout << " OverlayToTTree: initialized output tree in file " << _foutname << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int OverlayToTTree::process_event(PHCompositeNode *topNode)
{
  auto embedinfo = findNode::getClass<EmbedInfov1>(topNode, "EmbedInfo");
  if ( !embedinfo ) 
  {
    std::cout << " OverlayToTTree: ERROR - cannot find EmbedInfo node " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  _count = embedinfo->get_counter();
  _emb_cent = embedinfo->get_embed_cent();
  _emb_zvrtx = embedinfo->get_embed_zvrtx();
  if ( !_is_data ) 
  {
    _emb_b = embedinfo->get_embed_b();
    _emb_psi2 = embedinfo->get_embed_psi2();
  }
  _emb_ihcal_sumet = embedinfo->get_embed_ihcal_sumet();
  _emb_cemc_sumet = embedinfo->get_embed_cemc_sumet();
  _emb_ohcal_sumet = embedinfo->get_embed_ohcal_sumet();
  _emb_ihcal_sumet_scaled = embedinfo->get_embed_ihcal_sumet_scaled();
  _emb_cemc_sumet_scaled = embedinfo->get_embed_cemc_sumet_scaled();
  _emb_ohcal_sumet_scaled = embedinfo->get_embed_ohcal_sumet_scaled();
  _v2_weight = embedinfo->get_embed_v2_w();
  _dpsi2 = embedinfo->get_embed_dpsi2();
  _leadpt = embedinfo->get_embed_leadpt();
  _leadphi = embedinfo->get_embed_leadphi();
  
  if ( Verbosity() > 0 )
  {
    std::cout << " OverlayToTTree: processing event " << _count << std::endl;
    std::cout << " OverlayToTTree: retrieved from EmbedInfo: count / cent / zvrtx = " << _count << " / " << _emb_cent << " / " << _emb_zvrtx << std::endl;
    std::cout << " b / psi2 = " << _emb_b << " / " << _emb_psi2 << std::endl;
    std::cout << " v2_weight / dpsi2 = " << _v2_weight << " / " << _dpsi2 << std::endl;
    std::cout << " ihcal / cemc / ohcal SumET = " << _emb_ihcal_sumet << " / " << _emb_cemc_sumet << " / " << _emb_ohcal_sumet << std::endl;
    std::cout << " ihcal / cemc / ohcal SumET (scaled) = " << _emb_ihcal_sumet_scaled << " / " << _emb_cemc_sumet_scaled << " / " << _emb_ohcal_sumet_scaled << std::endl;
    std::cout << " lead jet pt / phi = " << _leadpt << " / " << _leadphi << std::endl;
  }


  auto genEventMap = findNode::getClass <PHHepMCGenEventMap> (topNode, "PHHepMCGenEventMap");
  if ( !genEventMap ) 
  {
    std::cout << " OverlayToTTree: ERROR - cannot find PHHepMCGenEventMap node " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  auto genEvt = genEventMap->get (1); // for primary simulated (?)
  _pythia_z = genEvt->get_collision_vertex().z();

  int ret = getZvertex(topNode);
  if ( ret != Fun4AllReturnCodes::EVENT_OK ) 
  {
    std::cout << " OverlayToTTree: ERROR - getZvertex() failed " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if ( Verbosity() > 0 )
  {
    std::cout << " OverlayToTTree: z (truth, reco, emb) = " << _pythia_z << " / " << _pythia_z_reco << " / " << _emb_zvrtx << std::endl;
  }

  auto *truth_jets = findNode::getClass<JetContainer>(topNode, "AntiKt_Truth_r03");
  if (!truth_jets)
  {
    std::cout << "OverlayFromTTree: ERROR - cannot find truth jet container node" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  } 
  _truth_jet_E.clear();
  _truth_jet_pt.clear();
  _truth_jet_eta.clear();
  _truth_jet_phi.clear();
  _truth_jet_recoUE.clear();
  _truth_jet_recoUE_cemc.clear();
  _truth_jet_recoUE_ihcal.clear();
  _truth_jet_recoUE_ohcal.clear();
  for (auto *this_jet : *truth_jets)
  {
    if (!this_jet) continue;

    const float pt  = this_jet->get_pt();
    const float eta = this_jet->get_eta();
    const float phi = this_jet->get_phi();
    const float e   = this_jet->get_e();
    const float UE = this_jet->get_property(truth_jets->property_index(Jet::PROPERTY::prop_JetCharge));
    const float UE_cemc = this_jet->get_property(truth_jets->property_index(Jet::PROPERTY::prop_BFrac));
    const float UE_ihcal = this_jet->get_property(truth_jets->property_index(Jet::PROPERTY::prop_area));
    const float UE_ohcal = this_jet->get_property(truth_jets->property_index(Jet::PROPERTY::prop_zg));
  
    _truth_jet_E.push_back(e);
    _truth_jet_pt.push_back(pt);
    _truth_jet_eta.push_back(eta);
    _truth_jet_phi.push_back(phi);
    _truth_jet_recoUE.push_back(UE);
    _truth_jet_recoUE_cemc.push_back(UE_cemc);
    _truth_jet_recoUE_ihcal.push_back(UE_ihcal);
    _truth_jet_recoUE_ohcal.push_back(UE_ohcal);
    if ( Verbosity() > 1 )
    {
      std::cout << " OverlayToTTree: truth jet E / pt / eta / phi = "
                << e << " / " << pt << " / "
                << eta << " / " << phi
                << ", reco UE = " << UE << std::endl;
    }

  } // end for truth_jets
  

  auto *reco_jets = findNode::getClass<JetContainer>(topNode, "AntiKt_Tower_r03_DVP");
  if (!reco_jets)
  {
    std::cout << "OverlayFromTTree: ERROR - cannot find reco jet container node" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  _recoPythia_jet_E.clear();
  _recoPythia_jet_pt.clear();
  _recoPythia_jet_eta.clear();
  _recoPythia_jet_phi.clear();
  _recoPythia_jet_EMfrac.clear();
  for (auto *this_jet : *reco_jets)
  {
    if (!this_jet) continue;

    const float pt = this_jet->get_pt();
    const float eta = this_jet->get_eta();
    const float phi = this_jet->get_phi();
    const float e = this_jet->get_e();
    if ( pt < _jetpt_thres ) continue;
    const float EMfrac = this_jet->get_property(reco_jets->property_index(Jet::PROPERTY::prop_JetCharge));
  
    _recoPythia_jet_E.push_back(e);
    _recoPythia_jet_pt.push_back(pt);
    _recoPythia_jet_eta.push_back(eta);
    _recoPythia_jet_phi.push_back(phi);
    _recoPythia_jet_EMfrac.push_back(EMfrac);

  } // end for reco_jets


  auto *towers_EM = getTowerInfos(topNode, "TOWERINFO_CALIB_CEMC_RETOWER_SUB1");
  auto *sub_jets = findNode::getClass<JetContainer>(topNode,"AntiKt_Tower_r03_Sub1");
  if (!sub_jets)
  {
    std::cout << "OverlayFromTTree: ERROR - cannot find Sub1 jet container node" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  _recoEmbed_jet_E.clear();  
  _recoEmbed_jet_pt.clear();
  _recoEmbed_jet_eta.clear();
  _recoEmbed_jet_phi.clear();
  _recoEmbed_jet_EMfrac.clear();
  for (auto *this_jet : *sub_jets)
  {
    if (!this_jet) continue;

    const float pt = this_jet->get_pt();
    const float eta = this_jet->get_eta();
    const float phi = this_jet->get_phi();
    const float e = this_jet->get_e();
    if ( pt < _jetpt_thres ) continue;
    float EMfrac = 0.0;
    for (const auto &comp_pair : this_jet->get_comp_vec())
    {
      if (comp_pair.first == Jet::SRC::CEMC_TOWERINFO_SUB1)
      {
        auto *tower = towers_EM->get_tower_at_channel(static_cast<int>(comp_pair.second));
        if (tower) EMfrac += tower->get_energy();
      }
    }
    if ( e > 0 ) EMfrac /= e;
  
    _recoEmbed_jet_E.push_back(e);
    _recoEmbed_jet_pt.push_back(pt);
    _recoEmbed_jet_eta.push_back(eta);
    _recoEmbed_jet_phi.push_back(phi);
    _recoEmbed_jet_EMfrac.push_back(EMfrac);

  } // end for sub_jets
   
  if ( Verbosity() > 1 )
  {
    std::cout << " OverlayToTTree: truth jets stored: " << _truth_jet_E.size() << std::endl;
    std::cout << " OverlayToTTree: recoPythia jets stored: " << _recoPythia_jet_E.size() << std::endl;
    std::cout << " OverlayToTTree: recoEmbed jets stored: " << _recoEmbed_jet_E.size() << std::endl;
  }

  _tree->Fill();
  if ( Verbosity() > 0 )
  {
    std::cout << " OverlayToTTree: filled tree for event " << _count << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
  
}

int OverlayToTTree::End(PHCompositeNode * /*topNode*/ )
{
  if( Verbosity() > 0 ) 
  {
    std::cout << " OverlayToTTree: writing output tree to file " << _foutname << std::endl;
    std::cout << " OverlayToTTree: total events processed: " << _count << std::endl;
  }

  PHTFileServer::get().cd( _foutname );
  _tree->Write();
  PHTFileServer::get().close();
 
  if ( Verbosity () > 0 ) 
  {
    std::cout << " OverlayToTTree: finished writing output file " << _foutname << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

inline TowerInfoContainer * OverlayToTTree::getTowerInfos(PHCompositeNode *topNode, const std::string &tower_node_name)
{
  auto *towerinfo = findNode::getClass<TowerInfoContainer>(topNode, tower_node_name);
  if (!towerinfo)
  {
    std::cout << PHWHERE << "Error: can't find TowerInfoContainer node " << tower_node_name << std::endl;
    std::exit(1);
  }
  return towerinfo;
}

int OverlayToTTree::getZvertex( PHCompositeNode *topNode )
{

  if ( Verbosity() > 0 )
  {
    std::cout << PHWHERE \
     << "Getting Z vertex from GlobalVertexMap node" << std::endl;
  }

  _pythia_z_reco = 0.0;

  auto vrtxmap = findNode::getClass< GlobalVertexMap >( topNode, "GlobalVertexMap");
  if ( !vrtxmap )
  {
    std::cerr << PHWHERE \
      << "GlobalVertexMap node not found, fatal error, exiting." << std::endl;
    exit(1);
  }

  if ( vrtxmap->empty() && Verbosity() > 0 )
  {
    std::cout << PHWHERE \
      << "empty vertex map, continuing as if zvtx = 0" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  auto vrtx = vrtxmap -> get_gvtxs_with_type( { GlobalVertex::MBD } );
  if ( vrtx.empty() && Verbosity() > 0 )
  {
    std::cout << PHWHERE \
      << "GlobalVertex is null, continuing as if zvtx = 0" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  if ( !vrtx.empty() )
  {
    if ( vrtx.at(0) )
    {
      _pythia_z_reco = vrtx.at(0) -> get_z();
    }
    else if ( Verbosity() > 0 )
    {
      std::cout << PHWHERE \
        << "GlobalVertex pointer within vector is null, continuing as if zvtx = 0" << std::endl;
    }
    
    
    if ( vrtx.size() > 1 && Verbosity() > 0 )
    {
      std::cout << PHWHERE \
        << "WARNING: more than one GlobalVertex with type " << GlobalVertex::MBD \
        << " found, using first one only." << std::endl;
    }

  }
  
  if ( std::isnan( _pythia_z_reco) )
  {
    static bool once = true;
    if ( once )
    {
      std::cout << PHWHERE \
        << "WARNING: retrieved NaN z vertex, resetting to 0" << std::endl;
      once = false;
    }
    _pythia_z_reco = 0.0;
  }

  if ( std::abs( _pythia_z_reco ) > 1e3 )
  {
    static bool once = true;
    if ( once )
    {
      std::cout << PHWHERE \
        << "WARNING: retrieved large z vertex = " << _pythia_z_reco \
        << ", resetting to 0" << std::endl;
      once = false;
    }
    _pythia_z_reco = 0.0;
  }

  return Fun4AllReturnCodes::EVENT_OK;

}