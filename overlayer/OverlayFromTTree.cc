#include "OverlayFromTTree.h"

#include "EmbedInfo.h"
#include "EmbedInfov1.h"

#include <fun4all/Fun4AllReturnCodes.h>

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

#include <jetbase/JetContainer.h>
#include <jetbase/Jet.h>

#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TMath.h>

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cassert>

double OverlayFromTTree::delta_psi2_phi(const double psi2, const double phi)
{
  double diffA = std::abs(psi2 - phi);
  if ( diffA < -TMath::Pi() ) { diffA += 2.0 * TMath::Pi(); }
  double diffB = TMath::Pi() - diffA;
  return std::min(std::abs(diffA), std::abs(diffB));
}

int OverlayFromTTree::Init(PHCompositeNode *topNode)
{
  _count = 0;

  _emb_file = TFile::Open(_emb_input_file.c_str(), "READ");
  if (!_emb_file || _emb_file->IsZombie())
  {
    std::cout << "OverlayFromTTree: ERROR - cannot open emb input file "
              << _emb_input_file << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _emb_tree = dynamic_cast<TTree*>(_emb_file->Get("T"));
  if (!_emb_tree)
  {
    std::cout << "OverlayFromTTree: ERROR - cannot find emb input tree in file "
              << _emb_input_file << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (Verbosity() > 0) _emb_tree->Show(0);

  _emb_tree->SetBranchStatus("*", false);

  if (!_is_data_input)
  {
    _emb_tree->SetBranchStatus("b", true);
    _emb_tree->SetBranchStatus("psi2", true);
    _emb_tree->SetBranchAddress("b", &_emb_b);
    _emb_tree->SetBranchAddress("psi2", &_emb_psi2);
  }

  _emb_tree->SetBranchStatus("sum_eT_all", true);
  _emb_tree->SetBranchAddress("sum_eT_all", &_emb_sumet);
  _emb_tree->SetBranchStatus("cent", true);
  _emb_tree->SetBranchStatus("zvrtx", true);
  _emb_tree->SetBranchAddress("cent", &_emb_cent);
  _emb_tree->SetBranchAddress("zvrtx", &_emb_zvrtx);

  _emb_tree->SetBranchStatus("sum_eT_cemc", true);
  _emb_tree->SetBranchStatus("cemc_E", true);
  // _emb_tree->SetBranchStatus("cemc_eta", true);
  // _emb_tree->SetBranchStatus("cemc_phi", true);
  _emb_tree->SetBranchStatus("cemc_isgood", true);
  _emb_tree->SetBranchAddress("sum_eT_cemc", &_emb_cemc_sumet);
  _emb_tree->SetBranchAddress("cemc_E", _emb_cemc_E);
  // _emb_tree->SetBranchAddress("cemc_eta", _emb_cemc_eta);
  // _emb_tree->SetBranchAddress("cemc_phi", _emb_cemc_phi);
  _emb_tree->SetBranchAddress("cemc_isgood", _emb_cemc_isgood);

  _emb_tree->SetBranchStatus("sum_eT_hcalin", true);
  _emb_tree->SetBranchStatus("hcalin_E", true);
  // _emb_tree->SetBranchStatus("hcalin_eta", true);
  // _emb_tree->SetBranchStatus("hcalin_phi", true);
  _emb_tree->SetBranchStatus("hcalin_isgood", true);
  _emb_tree->SetBranchAddress("sum_eT_hcalin", &_emb_ihcal_sumet);
  _emb_tree->SetBranchAddress("hcalin_E", _emb_ihcal_E);
  // _emb_tree->SetBranchAddress("hcalin_eta", _emb_ihcal_eta);
  // _emb_tree->SetBranchAddress("hcalin_phi", _emb_ihcal_phi);
  _emb_tree->SetBranchAddress("hcalin_isgood", _emb_ihcal_isgood);

  _emb_tree->SetBranchStatus("sum_eT_hcalout", true);
  _emb_tree->SetBranchStatus("hcalout_E", true);
  // _emb_tree->SetBranchStatus("hcalout_eta", true);
  // _emb_tree->SetBranchStatus("hcalout_phi", true);
  _emb_tree->SetBranchStatus("hcalout_isgood", true);
  _emb_tree->SetBranchAddress("sum_eT_hcalout", &_emb_ohcal_sumet);
  _emb_tree->SetBranchAddress("hcalout_E", _emb_ohcal_E);
  // _emb_tree->SetBranchAddress("hcalout_eta", _emb_ohcal_eta);
  // _emb_tree->SetBranchAddress("hcalout_phi", _emb_ohcal_phi);
  _emb_tree->SetBranchAddress("hcalout_isgood", _emb_ohcal_isgood);

  std::cout << "OverlayFromTTree: opened emb input file " << _emb_input_file << std::endl;

  // init TF1
  if (_jetv2_func) { delete _jetv2_func; _jetv2_func = nullptr; }
  _jetv2_func = new TF1("jetv2_func", OverlayFromTTree::vn_function, 0.0, TMath::Pi()/2.0, 1);
  _jetv2_func->SetParameter(0, _jetv2);
  _jetv2_func->SetNpx(1000);

  return CreateNode(topNode);
}

int OverlayFromTTree::process_event(PHCompositeNode *topNode)
{
  static bool first = true;
  if (first)
  {
    _cemc_R  = getCaloRadius(topNode, RawTowerDefs::CalorimeterId::CEMC);
    _ihcal_R = getCaloRadius(topNode, RawTowerDefs::CalorimeterId::HCALIN);
    _ohcal_R = getCaloRadius(topNode, RawTowerDefs::CalorimeterId::HCALOUT);
    first = false;

    if (Verbosity() > 0)
    {
      std::cout << "OverlayFromTTree: calo radii - CEMC / IHCal / OHCal = "
                << _cemc_R << " / " << _ihcal_R << " / " << _ohcal_R << std::endl;
    }
  }

  _emb_tree->GetEntry(_count);

  if (_target_cent >= 0)
  {
    while (_emb_cent > _target_cent)
    {
      ++_count;
      _emb_tree->GetEntry(_count);
    }
  }
  ++_count;

  auto *laudered_info = findNode::getClass<EmbedInfo>(topNode, "EmbedInfo");
  if (!laudered_info)
  {
    std::cout << "OverlayFromTTree: ERROR - cannot find EmbedInfo node" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  laudered_info -> set_counter(_count);
  laudered_info -> set_embed_cent(_emb_cent);
  laudered_info -> set_embed_zvrtx(_emb_zvrtx);

  if (!_is_data_input)
  {
    laudered_info->set_embed_b(_emb_b);
    laudered_info->set_embed_psi2(_emb_psi2);
  }

  std::memset(_emb_ihcal_E_scaled, 0, sizeof(_emb_ihcal_E_scaled));
  std::memset(_emb_cemc_E_scaled, 0, sizeof(_emb_cemc_E_scaled));
  std::memset(_emb_ohcal_E_scaled, 0, sizeof(_emb_ohcal_E_scaled));

  _emb_ihcal_sumet_scaled = _emb_cemc_sumet_scaled = _emb_ohcal_sumet_scaled = 0.0f;

  // for (unsigned int ieta = 0; ieta < 24; ++ieta)
  // {

  //   for (unsigned int iphi = 0; iphi < 64; ++iphi)
  //   { 
  //     if (_emb_ihcal_isgood[ieta][iphi] > 0)
  //     {
  //       // _emb_ihcal_E_scaled[ieta][iphi] = _emb_ihcal_E[ieta][iphi] * _ihcal_scale;
  //       // const double eta = correct_calo_eta(_emb_ihcal_eta[ieta][iphi], _ihcal_R, (double)_emb_zvrtx);
  //       // _emb_ihcal_sumet        += _emb_ihcal_E[ieta][iphi]        / std::cosh(eta);
  //       // _emb_ihcal_sumet_scaled += _emb_ihcal_E_scaled[ieta][iphi] / std::cosh(eta);
  //       _emb_ihcal_E_scaled[ieta][iphi] = _emb_ihcal_E[ieta][iphi];
  //     }

  //     if (_emb_cemc_isgood[ieta][iphi] > 0)
  //     {
  //       // _emb_cemc_E_scaled[ieta][iphi] = _emb_cemc_E[ieta][iphi] * _cemc_scale;
  //       // const double eta = correct_calo_eta(_emb_cemc_eta[ieta][iphi], _cemc_R, (double)_emb_zvrtx);
  //       // _emb_cemc_sumet        += _emb_cemc_E[ieta][iphi]        / std::cosh(eta);
  //       // _emb_cemc_sumet_scaled += _emb_cemc_E_scaled[ieta][iphi] / std::cosh(eta);

  //     }

  //     if (_emb_ohcal_isgood[ieta][iphi] > 0)
  //     {
  //       _emb_ohcal_E_scaled[ieta][iphi] = _emb_ohcal_E[ieta][iphi] * _ohcal_scale;
  //       const double eta = correct_calo_eta(_emb_ohcal_eta[ieta][iphi], _ohcal_R, (double)_emb_zvrtx);
  //       _emb_ohcal_sumet        += _emb_ohcal_E[ieta][iphi]        / std::cosh(eta);
  //       _emb_ohcal_sumet_scaled += _emb_ohcal_E_scaled[ieta][iphi] / std::cosh(eta);
  //     }
  //   }

  // }

  auto *truth_jets = findNode::getClass<JetContainer>(topNode, "AntiKt_Truth_r03");
  if (!truth_jets)
  {
    std::cout << "OverlayFromTTree: ERROR - cannot find truth jet container node" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  float lead_jet_phi = 0.0f;
  float lead_jet_eta = 0.0f;
  float lead_jet_pT  = -1.0f;

  for (auto *this_jet : *truth_jets)
  {
    if (!this_jet) continue;

    const float pt  = this_jet->get_pt();
    const float eta = this_jet->get_eta();
    const float phi = this_jet->get_phi();

    if (pt > lead_jet_pT)
    {
      lead_jet_pT  = pt;
      lead_jet_eta = eta;
      lead_jet_phi = phi;
    }


    double sum_eTreco_cemc = 0.0, sum_eTreco_ihcal = 0.0, sum_eTreco_ohcal = 0.0;

    for (unsigned int ieta2 = 0; ieta2 < 24; ++ieta2)
    for (unsigned int iphi2 = 0; iphi2 < 64; ++iphi2)
    {
      if (_emb_cemc_isgood[ieta2][iphi2] > 0 &&
          deltaR(eta, phi, _emb_cemc_eta[ieta2][iphi2], _emb_cemc_phi[ieta2][iphi2]) < 0.3)
      {
        const double et = _emb_cemc_E_scaled[ieta2][iphi2] /
                          std::cosh(correct_calo_eta(_emb_cemc_eta[ieta2][iphi2], _cemc_R, (double)_emb_zvrtx));
        sum_eTreco_cemc += et;
      }

      if (_emb_ihcal_isgood[ieta2][iphi2] > 0 &&
          deltaR(eta, phi, _emb_ihcal_eta[ieta2][iphi2], _emb_ihcal_phi[ieta2][iphi2]) < 0.3)
      {
        const double et = _emb_ihcal_E_scaled[ieta2][iphi2] /
                          std::cosh(correct_calo_eta(_emb_ihcal_eta[ieta2][iphi2], _ihcal_R, (double)_emb_zvrtx));
        sum_eTreco_ihcal += et;
      }

      if (_emb_ohcal_isgood[ieta2][iphi2] > 0 &&
          deltaR(eta, phi, _emb_ohcal_eta[ieta2][iphi2], _emb_ohcal_phi[ieta2][iphi2]) < 0.3)
      {
        const double et = _emb_ohcal_E_scaled[ieta2][iphi2] /
                          std::cosh(correct_calo_eta(_emb_ohcal_eta[ieta2][iphi2], _ohcal_R, (double)_emb_zvrtx));
        sum_eTreco_ohcal += et;
      }
    }

    const double SumET_reco = sum_eTreco_cemc + sum_eTreco_ihcal + sum_eTreco_ohcal;

    if (Verbosity() > 1)
    {
      std::cout << "OverlayFromTTree: truth jet E / pt / eta / phi = "
                << this_jet->get_e() << " / " << pt << " / "
                << this_jet->get_eta() << " / " << this_jet->get_phi()
                << ", SumET reco = " << SumET_reco << std::endl;
    }

    this_jet->set_property(truth_jets->property_index(Jet::PROPERTY::prop_JetCharge), SumET_reco);
    this_jet->set_property(truth_jets->property_index(Jet::PROPERTY::prop_BFrac),     sum_eTreco_cemc);
    this_jet->set_property(truth_jets->property_index(Jet::PROPERTY::prop_area),     sum_eTreco_ihcal);
    this_jet->set_property(truth_jets->property_index(Jet::PROPERTY::prop_zg),       sum_eTreco_ohcal);
  }

  laudered_info->set_embed_ihcal_sumet(_emb_ihcal_sumet);
  laudered_info->set_embed_cemc_sumet(_emb_cemc_sumet);
  laudered_info->set_embed_ohcal_sumet(_emb_ohcal_sumet);
  laudered_info->set_embed_ihcal_sumet_scaled(_emb_ihcal_sumet_scaled);
  laudered_info->set_embed_cemc_sumet_scaled(_emb_cemc_sumet_scaled);
  laudered_info->set_embed_ohcal_sumet_scaled(_emb_ohcal_sumet_scaled);

  const double dpsi2 = (!_is_data_input) ? delta_psi2_phi(_emb_psi2, lead_jet_phi) : 0.0;
  const double v2_w  = (_jetv2_func && !_is_data_input) ? _jetv2_func->Eval(dpsi2) : 1.0;

  if (Verbosity() > 1)
  {
    std::cout << "OverlayFromTTree: lead truth jet pT / eta / phi = "
              << lead_jet_pT << " / " << lead_jet_eta << " / " << lead_jet_phi
              << ", dpsi2 = " << dpsi2 << " , v2 weight = " << v2_w << std::endl;
  }

  laudered_info->set_embed_v2_w(v2_w);
  laudered_info->set_embed_dpsi2(dpsi2);
  laudered_info->set_embed_leadpt(lead_jet_pT);
  laudered_info->set_embed_leadphi(lead_jet_phi);

  auto *towers_EM = getTowerInfos(topNode, _cemc_tower_node);
  auto *towers_IH = getTowerInfos(topNode, _ihcal_tower_node);
  auto *towers_OH = getTowerInfos(topNode, _ohcal_tower_node);

  auto *reco_jets = findNode::getClass<JetContainer>(topNode, "AntiKt_Tower_r03_DVP");
  if (!reco_jets)
  {
    std::cout << "OverlayFromTTree: ERROR - cannot find reco jet container node" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  for (auto *this_jet : *reco_jets)
  {
    if (!this_jet) continue;

    const float jet_pt = this_jet->get_pt();
    if (jet_pt < _jetpt_thres) continue;

    double em_frac = 0.0;
    for (const auto &comp_pair : this_jet->get_comp_vec())
    {
      if (comp_pair.first == Jet::SRC::CEMC_TOWERINFO_RETOWER)
      {
        auto *tower = towers_EM->get_tower_at_channel(static_cast<int>(comp_pair.second));
        if (tower) em_frac += tower->get_energy();
      }
    }

    const double jet_e = this_jet->get_e();
    if (jet_e > 0) em_frac /= jet_e;

    if (Verbosity() > 1)
    {
      std::cout << "OverlayFromTTree: reco jet E / pt / eta / phi = "
                << this_jet->get_e() << " / " << jet_pt << " / "
                << this_jet->get_eta() << " / " << this_jet->get_phi()
                << ", EM frac = " << em_frac << std::endl;
    }

    this_jet->set_property(reco_jets->property_index(Jet::PROPERTY::prop_JetCharge), em_frac);
  }

  for (unsigned int ich = 0; ich < towers_EM->size(); ++ich)
  {
    auto *tower = towers_EM->get_tower_at_channel(ich);
    assert(tower);

    const int ieta = towers_EM->getTowerEtaBin(ich);
    const int iphi = towers_EM->getTowerPhiBin(ich);

    tower->set_energy(tower->get_energy() + _emb_cemc_E[ieta][iphi]);
  }

  for (unsigned int ich = 0; ich < towers_OH->size(); ++ich)
  {
    auto *tower = towers_OH->get_tower_at_channel(ich);
    assert(tower);

    const int ieta = towers_OH->getTowerEtaBin(ich);
    const int iphi = towers_OH->getTowerPhiBin(ich);

    tower->set_energy(tower->get_energy() + _emb_ohcal_E[ieta][iphi]);
  }

  for (unsigned int ich = 0; ich < towers_IH->size(); ++ich)
  {
    auto *tower = towers_IH->get_tower_at_channel(ich);
    assert(tower);

    const int ieta = towers_IH->getTowerEtaBin(ich);
    const int iphi = towers_IH->getTowerPhiBin(ich);

    tower->set_energy(tower->get_energy() + _emb_ihcal_E[ieta][iphi]);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int OverlayFromTTree::End(PHCompositeNode * /*topNode*/)
{
  if (_jetv2_func) { delete _jetv2_func; _jetv2_func = nullptr; }

  if (_emb_file)
  {
    _emb_file->Close();
    delete _emb_file;
    _emb_file = nullptr;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

inline double OverlayFromTTree::getCaloRadius(PHCompositeNode *topNode, const RawTowerDefs::CalorimeterId caloid)
{
  if (Verbosity() > 2)
  {
    std::cout << PHWHERE << "Getting calo radius for calo id "
              << RawTowerDefs::convert_caloid_to_name(caloid) << std::endl;
  }

  auto *towgeo = findNode::getClass<RawTowerGeomContainer>(
      topNode, "TOWERGEOM_" + RawTowerDefs::convert_caloid_to_name(caloid));

  if (!towgeo)
  {
    std::cout << PHWHERE << "Error: can't find RawTowerGeomContainer node for calo id "
              << RawTowerDefs::convert_caloid_to_name(caloid) << std::endl;
    std::exit(1);
  }

  auto *geo = towgeo->get_tower_geometry(RawTowerDefs::encode_towerid(caloid, 0, 0));
  if (!geo)
  {
    std::cout << PHWHERE << "Error: can't find tower geometry for calo id "
              << RawTowerDefs::convert_caloid_to_name(caloid) << std::endl;
    std::exit(1);
  }

  return geo->get_center_radius();
}

inline TowerInfoContainer* OverlayFromTTree::getTowerInfos(PHCompositeNode *topNode, const std::string &tower_node_name)
{
  auto *towerinfo = findNode::getClass<TowerInfoContainer>(topNode, tower_node_name);
  if (!towerinfo)
  {
    std::cout << PHWHERE << "Error: can't find TowerInfoContainer node " << tower_node_name << std::endl;
    std::exit(1);
  }
  return towerinfo;
}

int OverlayFromTTree::CreateNode(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  auto *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  auto *bkgNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "LAURDERED"));
  if (!bkgNode)
  {
    bkgNode = new PHCompositeNode("LAURDERED");
    dstNode->addNode(bkgNode);
  }

  auto *embedinfo = findNode::getClass<EmbedInfo>(topNode, "EmbedInfo");
  if (!embedinfo)
  {
    embedinfo = new EmbedInfov1();
    auto *node = new PHIODataNode<PHObject>(embedinfo, "EmbedInfo", "PHObject");
    bkgNode->addNode(node);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

inline RawTowerGeomContainer* OverlayFromTTree::getTowerGeoms(PHCompositeNode *topNode, const std::string &geom_node_name)
{
  auto *towgeom = findNode::getClass<RawTowerGeomContainer>(topNode, geom_node_name);
  if (!towgeom)
  {
    std::cout << PHWHERE << "Error: can't find RawTowerGeomContainer node " << geom_node_name << std::endl;
    std::exit(1);
  }
  return towgeom;
}

inline double OverlayFromTTree::deltaR(const double eta1, const double phi1, const double eta2, const double phi2)
{
  double dphi = phi1 - phi2;
  if (dphi > M_PI)      dphi -= 2 * M_PI;
  else if (dphi < -M_PI) dphi += 2 * M_PI;

  const double deta = eta1 - eta2;
  return std::sqrt(deta * deta + dphi * dphi);
}

inline double OverlayFromTTree::correct_calo_eta(const double eta, const double R, const double zvrtx)
{
  const double z0 = std::sinh(eta) * R;
  const double z1 = z0 - zvrtx;
  return std::asinh(z1 / R);
}

inline float OverlayFromTTree::get_centdeb_cemc(int cent)
{
  const float scale_factors[10] = {1.05, 1.04, 1.03, 1.02, 1.01, 1.00, 0.99, 0.98, 0.97, 0.96};
  return scale_factors[cent];
}
inline float OverlayFromTTree::get_centdeb_ihcal(int cent)
{
  const float scale_factors[10] = {1.10, 1.08, 1.06, 1.04, 1.02, 1.00, 0.98, 0.96, 0.94, 0.92};
  return scale_factors[cent];
}
inline float OverlayFromTTree::get_centdeb_ohcal(int cent)
{
  const float scale_factors[10] = {1.15, 1.12, 1.09, 1.06, 1.03, 1.00, 0.97, 0.94, 0.91, 0.88};
  return scale_factors[cent];
}