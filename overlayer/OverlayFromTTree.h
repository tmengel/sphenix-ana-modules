#ifndef _OVERLAYFORMTTREE_H_
#define _OVERLAYFORMTTREE_H_

#include <fun4all/SubsysReco.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

// sPHENIX / calo
#include <calobase/RawTowerDefs.h>

// ROOT math (used in static helpers in this header)
#include <TMath.h>

class TFile;
class TTree;
class TF1;

class TowerInfoContainer;
class RawTowerGeomContainer;

class PHCompositeNode;

class OverlayFromTTree : public SubsysReco
{
 public:
  OverlayFromTTree(const std::string & fname = "OverlayFromTTree_output.root")
    : SubsysReco("OverlayFromTTree")
    , _emb_input_file(fname)
  {}

  ~OverlayFromTTree() override = default;

  int Init(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void set_emb_input(const std::string &emb_input, bool is_data_input = false)
  {
    _is_data_input = is_data_input;
    _emb_input_file = emb_input;
  }

  void set_target_cent(const int cent) { _target_cent = cent; }
  void set_pT_threshold(const float pT) { _jetpt_thres = pT; }

  void set_jetv2(float v2) { _jetv2 = v2; }
  void set_emcal_scale(float scale) { _cemc_scale = scale; }
  void set_ihcal_scale(float scale) { _ihcal_scale = scale; }
  void set_ohcal_scale(float scale) { _ohcal_scale = scale; }
  // void set_centdep_scalling(bool val) { _centdeb_scalling = val; }

  float get_jetv2() const { return _jetv2; }
  float get_emcal_scale() const { return _cemc_scale; }
  float get_ihcal_scale() const { return _ihcal_scale; }
  float get_ohcal_scale() const { return _ohcal_scale; }
  // void get_centdep_scalling() const { return _centdeb_scalling; }


  // Optional: let caller override node names if your DST naming differs
  void set_cemc_tower_node(const std::string& s) { _cemc_tower_node = s; }
  void set_ihcal_tower_node(const std::string& s) { _ihcal_tower_node = s; }
  void set_ohcal_tower_node(const std::string& s) { _ohcal_tower_node = s; }

 private:
  int CreateNode(PHCompositeNode *topNode);

  std::string _emb_input_file {""};
  bool _is_data_input {false};

  TFile *_emb_file {nullptr};
  TTree *_emb_tree {nullptr};

  float _jetpt_thres {5.0};
  int _target_cent {-1};

  TF1* _jetv2_func {nullptr};
  float _jetv2 {0.0};

  float _cemc_scale {1.0};
  float _ihcal_scale {1.0};
  float _ohcal_scale {1.0};
  bool _centdeb_scalling {false};
  inline float get_centdeb_cemc(int cent);
  inline float get_centdeb_ihcal(int cent);
  inline float get_centdeb_ohcal(int cent);

  // tower node names (defaults; adjust if your DST uses different ones)
  std::string _cemc_tower_node {"TOWERINFO_CALIB_CEMC_RETOWER"};
  std::string _ihcal_tower_node {"TOWERINFO_CALIB_HCALIN"};
  std::string _ohcal_tower_node {"TOWERINFO_CALIB_HCALOUT"};

  int   _emb_cent { -1 };
  float   _emb_zvrtx { 0 };
  float _emb_b {0.0f};
  float _emb_psi2 {0.0f};
  float _emb_sumet {0.0f};
  float _emb_sumet_scaled {0.0f};

  float _emb_ohcal_sumet {0.0f};
  float _emb_ohcal_sumet_scaled {0.0f};
  float _emb_ohcal_E[24][64] {};
  float _emb_ohcal_E_scaled[24][64] {};
  float _emb_ohcal_eta[24][64] {};
  float _emb_ohcal_phi[24][64] {};
  int   _emb_ohcal_isgood[24][64] {};

  float _emb_ihcal_sumet {0.0f};
  float _emb_ihcal_sumet_scaled {0.0f};
  float _emb_ihcal_E[24][64] {};
  float _emb_ihcal_E_scaled[24][64] {};
  float _emb_ihcal_eta[24][64] {};
  float _emb_ihcal_phi[24][64] {};
  int   _emb_ihcal_isgood[24][64] {};

  float _emb_cemc_sumet {0.0f};
  float _emb_cemc_sumet_scaled {0.0f};
  float _emb_cemc_E[24][64] {};
  float _emb_cemc_E_scaled[24][64] {};
  float _emb_cemc_eta[24][64] {};
  float _emb_cemc_phi[24][64] {};
  int   _emb_cemc_isgood[24][64] {};

  double _cemc_R {0.0};
  double _ihcal_R {0.0};
  double _ohcal_R {0.0};

  int _count {0};

  static double vn_function(double *x, double *par)
  {
    const double v2  = par[0];
    const double phi = x[0];
    return 1.0 + 2.0 * v2 * TMath::Cos(2.0 * phi);
  }

  static double delta_psi2_phi(const double psi2, const double phi);

  inline double deltaR(const double eta1, const double phi1, const double eta2, const double phi2);
  inline double correct_calo_eta(const double eta, const double R, const double zvrtx);

  inline double getCaloRadius(PHCompositeNode *topNode, const RawTowerDefs::CalorimeterId caloid);
  inline TowerInfoContainer* getTowerInfos(PHCompositeNode *topNode, const std::string &tower_node_name);
  inline RawTowerGeomContainer* getTowerGeoms(PHCompositeNode *topNode, const std::string &geom_node_name);
};

#endif // _OVERLAYFORMTTREE_H_
