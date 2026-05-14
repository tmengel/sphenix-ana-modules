#ifndef _OVERLAYTOTTREE_H_
#define _OVERLAYTOTTREE_H_

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

class TFile;
class TTree;
class TF1;

class TowerInfoContainer;

class PHCompositeNode;

class OverlayToTTree : public SubsysReco
{
 public:

  OverlayToTTree(const std::string & fname = "OverlayFromTTree_output.root")
    : SubsysReco("OverlayFromTTree")
    , _foutname( fname )
  {}
  ~OverlayToTTree() override = default;

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void is_data(bool val) { _is_data = val; }

  void set_pT_threshold(const float pT) { _jetpt_thres = pT; }

  void set_jetv2(float v2) { _jetv2 = v2; }

  void set_emcal_scale(float scale) { _cemc_scale = scale; }
  void set_ihcal_scale(float scale) { _ihcal_scale = scale; }
  void set_ohcal_scale(float scale) { _ohcal_scale = scale; }

  float get_jetv2() const { return _jetv2; }
  float get_emcal_scale() const { return _cemc_scale; }
  float get_ihcal_scale() const { return _ihcal_scale; }
  float get_ohcal_scale() const { return _ohcal_scale; }

 private:

  std::string _foutname {""};
  bool _is_data {false};

  TTree *_tree {nullptr};

  float _jetpt_thres {5.0};
  int _target_cent {-1};

  TF1* _jetv2_func {nullptr};
  float _jetv2 {0.0};
  float _v2_weight {0.0};
  float _dpsi2 {0.0};
  float _leadpt {0.0};
  float _leadphi {0.0};

  float _cemc_scale {1.0};
  float _ihcal_scale {1.0};
  float _ohcal_scale {1.0};

  int   _emb_cent { -1 };
  int   _emb_zvrtx { 0 };
  float _emb_b {0.0f};
  float _emb_psi2 {0.0f};

  float _emb_ohcal_sumet {0.0f};
  float _emb_ohcal_sumet_scaled {0.0f};

  float _emb_ihcal_sumet {0.0f};
  float _emb_ihcal_sumet_scaled {0.0f};

  float _emb_cemc_sumet {0.0f};
  float _emb_cemc_sumet_scaled {0.0f};

  int _count {0};

  float _pythia_z {0.0f};
  float _pythia_z_reco {0.0f};
  std::vector<float> _truth_jet_E{};
  std::vector<float> _truth_jet_pt{};
  std::vector<float> _truth_jet_eta{};
  std::vector<float> _truth_jet_phi{};
  std::vector<float> _truth_jet_recoUE{};
  std::vector<float> _truth_jet_recoUE_cemc{};
  std::vector<float> _truth_jet_recoUE_ihcal{};
  std::vector<float> _truth_jet_recoUE_ohcal{};

  std::vector<float> _recoPythia_jet_E{};
  std::vector<float> _recoPythia_jet_pt{};
  std::vector<float> _recoPythia_jet_eta{};
  std::vector<float> _recoPythia_jet_phi{};
  std::vector<float> _recoPythia_jet_EMfrac{};

  std::vector<float> _recoEmbed_jet_E{};
  std::vector<float> _recoEmbed_jet_pt{};
  std::vector<float> _recoEmbed_jet_eta{};
  std::vector<float> _recoEmbed_jet_phi{};
  std::vector<float> _recoEmbed_jet_EMfrac{};

  inline TowerInfoContainer* getTowerInfos(PHCompositeNode *topNode, const std::string &tower_node_name);
  int getZvertex(PHCompositeNode *topNode);


};

#endif // _OVERLAYTOTTREE_H_
