#ifndef SimTree_H
#define SimTree_H

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <string>
#include <vector>
#include <array>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <cstdio>
#include <cstring>

class PHCompositeNode;
class TTree;
class TRandom3;

class SimTree : public SubsysReco
{
 public:

  // super simple tree maker
  SimTree( const std::string & outputfile = "output.root" ) 
    : SubsysReco("SimTree")
    , m_output_filename( outputfile )
  {}

  ~SimTree() override
  {
    if ( Verbosity() > 0 ) 
    {
      std::cout << "SimTree::~SimTree - done" << std::endl;
    }
    return;
  }


  int Init( PHCompositeNode * /*topNode*/) override;
  
  int process_event( PHCompositeNode * topNode ) override;
  
  int End( PHCompositeNode * /*topNode*/ ) override;
  
  int ResetEvent( PHCompositeNode * /*topNode*/ ) override 
  {
    ResetZvtx();
    ResetCent();
    ResetEventHeader();
    ResetSub1Jet();
    ResetTruthJet();
    ResetCaloInfo();
    ResetRawJet();
    ResetG4Truth();
    return Fun4AllReturnCodes::EVENT_OK;
  }

  void add_zvrtx_node ( const std::string & name = "GlobalVertexMap" ){ m_zvrtx_node = name; }
  void add_cent_node ( const std::string & name = "CentralityInfo" ) { m_cent_node = name;  }
  void add_sub1jet_node ( const std::string & name ) { m_sub1jet_node = name; }
  void add_event_header( const std::string & name =  "EventHeader"){ m_eventhead_node = name; }
  void add_truthjet_node ( const std::string & name ) { m_truthjet_node = name; }

  void add_rawjet_node ( const std::string & name ) { m_rawjet_node = name; }
  void add_multjet_node ( const std::string & name ) { m_multjet_node = name; }
  void add_areajet_node ( const std::string & name ) { m_areajet_node = name; }

  void do_towerbkgd( const bool do_bkgd = true ) { m_do_towerbkgd = do_bkgd; }
  void do_rho( const bool do_rho = true ) { m_do_rho = do_rho; }

  void add_truthnode( const std::string & name ) { m_g4truth_node = name; }
  void add_ep_info ( const std::string &name = "EventPlaneInfo" ) { m_eventplane_node = name; }


 private:
    
  // output file name
  std::string m_output_filename { "" };


  TTree * m_tree {nullptr};
  int m_event_id {-1};

  float m_sumeT_cemc { 0.0 };
  float m_sumeT_hcalin { 0.0 };
  float m_sumeT_hcalout { 0.0 };
  float m_sumeT_cemc_sub1 { 0.0 };
  float m_sumeT_hcalin_sub1 { 0.0 };
  float m_sumeT_hcalout_sub1 { 0.0 };
  float m_sumeT_cemc_org { 0.0 };
  float m_sumeT_hcalin_org { 0.0 };
  float m_sumeT_hcalout_org { 0.0 };

  void ResetCaloInfo()
  {
    m_sumeT_cemc = 0.0;
    m_sumeT_hcalin = 0.0;
    m_sumeT_hcalout = 0.0;
    m_sumeT_cemc_sub1 = 0.0;
    m_sumeT_hcalin_sub1 = 0.0;
    m_sumeT_hcalout_sub1 = 0.0;
    m_sumeT_cemc_org = 0.0;
    m_sumeT_hcalin_org = 0.0;
    m_sumeT_hcalout_org = 0.0;    
  }

  std::string m_zvrtx_node { "" };
  float m_zvtx { 0.0 };
  void ResetZvtx()
  {
    m_zvtx = -999;
  }

  // centrality info
  std::string m_cent_node { "" };
  int m_cent {-1};
  void ResetCent()
  {
    m_cent = -1;
  }

  // event header info
  std::string m_eventhead_node { "" };
  float m_b { 0.0 };
  float m_ep_angle { 0.0 };
  float m_ecc { 0.0 };
  float m_psi1 { 0.0 };
  float m_psi2 { 0.0 };
  float m_psi3 { 0.0 };
  float m_psi4 { 0.0 };
  float m_psi5 { 0.0 };
  float m_psi6 { 0.0 };
  float m_ncoll { 0.0 };
  float m_npart { 0.0 };
  void ResetEventHeader()
  {
    m_b = -999;
    m_ep_angle = -999;
    m_ecc = -999;
    m_psi1 = -999;
    m_psi2 = -999;
    m_psi3 = -999;
    m_psi4 = -999;
    m_psi5 = -999;
    m_psi6 = -999;
    m_ncoll = -999;
    m_npart = -999;
  }

  // sub1 jet info
  std::string m_sub1jet_node { "" };
  std::vector < float > m_sub1_jet_E {};
  std::vector < float > m_sub1_jet_phi {};
  std::vector < float > m_sub1_jet_eta {};
  std::vector < float > m_sub1_jet_pT {};
  std::vector < float > m_sub1_jet_unsub_pT {};
  std::vector < float > m_sub1_jet_unsub_E {};
  std::vector < std::vector < int > > m_sub1_jet_comp_ieta {};
  std::vector < std::vector < int > > m_sub1_jet_comp_iphi {};
  std::vector < std::vector < int > > m_sub1_jet_comp_caloid {};
  std::vector < std::vector < int > > m_sub1_jet_comp_status {};
  std::vector < std::vector < float > > m_sub1_jet_comp_E {};
  std::vector < std::vector < float > > m_sub1_jet_comp_eta {};
  std::vector < std::vector < float > > m_sub1_jet_comp_phi {};
  void ResetSub1Jet()
  {
    m_sub1_jet_E.clear();
    m_sub1_jet_phi.clear();
    m_sub1_jet_eta.clear();
    m_sub1_jet_pT.clear();
    m_sub1_jet_unsub_pT.clear();
    m_sub1_jet_unsub_E.clear();
    m_sub1_jet_comp_ieta.clear();
    m_sub1_jet_comp_iphi.clear();
    m_sub1_jet_comp_caloid.clear();
    m_sub1_jet_comp_status.clear();
    m_sub1_jet_comp_E.clear();
    m_sub1_jet_comp_eta.clear();
    m_sub1_jet_comp_phi.clear();
  }

  static const int k_ieta = 24;
  static const int k_iphi = 64;
 
  bool m_do_towerbkgd = false;
  bool m_do_rho = false;

  float m_sub1_towerbkgd_cemc[k_ieta] {};
  float m_sub1_towerbkgd_hcalin[k_ieta] {};
  float m_sub1_towerbkgd_hcalout[k_ieta] {};
  float m_sub2_towerbkgd_cemc[k_ieta] {};
  float m_sub2_towerbkgd_hcalin[k_ieta] {};
  float m_sub2_towerbkgd_hcalout[k_ieta] {}; 
  float m_sub1_v2 { 0.0 };
  float m_sub2_v2 { 0.0 };
  int m_sub1_flowfaliure { 0 };
  int m_sub2_flowfaliure { 0 };
  float m_sub1_psi2 { 0.0 };
  float m_sub2_psi2 { 0.0 };

  float m_rhoM_cemc { 0.0 };
  float m_rhoM_hcalin { 0.0 };
  float m_rhoM_hcalout { 0.0 };
  float m_rhoM { 0.0 };
  
  float m_rhoA_cemc { 0.0 };
  float m_rhoA_hcalin { 0.0 };
  float m_rhoA_hcalout { 0.0 };
  float m_rhoA { 0.0 };
  
  float m_rhoM_cemc_std { 0.0 };
  float m_rhoM_hcalin_std { 0.0 };
  float m_rhoM_hcalout_std { 0.0 };
  float m_rhoM_std { 0.0 };

  float m_rhoA_cemc_std { 0.0 };
  float m_rhoA_hcalin_std { 0.0 };
  float m_rhoA_hcalout_std { 0.0 };
  float m_rhoA_std { 0.0 };
  void ResetTowerBkgd()
  {
    memset(m_sub1_towerbkgd_cemc, 0, sizeof(m_sub1_towerbkgd_cemc));
    memset(m_sub1_towerbkgd_hcalin, 0, sizeof(m_sub1_towerbkgd_hcalin));
    memset(m_sub1_towerbkgd_hcalout, 0, sizeof(m_sub1_towerbkgd_hcalout));
    memset(m_sub2_towerbkgd_cemc, 0, sizeof(m_sub2_towerbkgd_cemc));
    memset(m_sub2_towerbkgd_hcalin, 0, sizeof(m_sub2_towerbkgd_hcalin));
    memset(m_sub2_towerbkgd_hcalout, 0, sizeof(m_sub2_towerbkgd_hcalout));
    m_sub1_v2 = 0.0;
    m_sub2_v2 = 0.0;
    m_sub1_flowfaliure = 0;
    m_sub2_flowfaliure = 0;
  }
  void ResetRho()
  {
    m_rhoM_cemc = 0.0;
    m_rhoM_hcalin = 0.0;
    m_rhoM_hcalout = 0.0;
    m_rhoM = 0.0;
    
    m_rhoA_cemc = 0.0;
    m_rhoA_hcalin = 0.0;
    m_rhoA_hcalout = 0.0;
    m_rhoA = 0.0;

    m_rhoM_cemc_std = 0.0;
    m_rhoM_hcalin_std = 0.0;
    m_rhoM_hcalout_std = 0.0;
    m_rhoM_std = 0.0;

    m_rhoA_cemc_std = 0.0;
    m_rhoA_hcalin_std = 0.0;
    m_rhoA_hcalout_std = 0.0;
    m_rhoA_std = 0.0;    
  }

  // truth jet info
  std::string m_truthjet_node { "" };
  std::vector < float > m_truth_jet_E {};
  std::vector < float > m_truth_jet_phi {};
  std::vector < float > m_truth_jet_eta {};
  std::vector < float > m_truth_jet_pT {};
  void ResetTruthJet()
  {
    m_truth_jet_E.clear();
    m_truth_jet_phi.clear();
    m_truth_jet_eta.clear();
    m_truth_jet_pT.clear();
  }

  std::string m_rawjet_node { "" };
  std::vector < float > m_raw_jet_E {};
  std::vector < float > m_raw_jet_phi {};
  std::vector < float > m_raw_jet_eta {};
  std::vector < float > m_raw_jet_pT {};
  
  std::string m_multjet_node { "" };
  std::vector < float > m_mult_jet_E {};
  std::vector < float > m_mult_jet_phi {};
  std::vector < float > m_mult_jet_eta {};
  std::vector < float > m_mult_jet_pT {};

  std::string m_areajet_node { "" };
  std::vector < float > m_area_jet_E {};
  std::vector < float > m_area_jet_phi {};
  std::vector < float > m_area_jet_eta {};
  std::vector < float > m_area_jet_pT {};
  void ResetRawJet()
  {
    m_raw_jet_E.clear();
    m_raw_jet_phi.clear();
    m_raw_jet_eta.clear();
    m_raw_jet_pT.clear();
    
  }
  void ResetMultJet()
  {
    m_mult_jet_E.clear();
    m_mult_jet_phi.clear();
    m_mult_jet_eta.clear();
    m_mult_jet_pT.clear();
  }
  void ResetAreaJet()
  {
    m_area_jet_E.clear();
    m_area_jet_phi.clear();
    m_area_jet_eta.clear();
    m_area_jet_pT.clear();
  }
  
  std::string m_g4truth_node { "" };
  float m_g4truth_zvtx { 0.0 };
  float m_g4truth_v2reco { 0.0 };
  float m_g4truth_v3reco { 0.0 };
  float m_g4truth_v4reco { 0.0 };
  float m_g4truth_v5reco { 0.0 };
  float m_g4truth_v6reco { 0.0 };
  void ResetG4Truth()
  {
    m_g4truth_v2reco = -999;
    m_g4truth_v3reco = -999;
    m_g4truth_v4reco = -999;
    m_g4truth_v5reco = -999;
    m_g4truth_v6reco = -999;
    m_g4truth_zvtx = -999;
  }

  std::string m_eventplane_node { "" };
  std::vector< float > m_psi_shifted_NS {};
  std::vector< float > m_psi_NS {};
  std::vector< float > m_psi_shifted_S {};
  std::vector< float > m_psi_S {};
  std::vector< float > m_psi_shifted_N {};
  std::vector< float > m_psi_N {};
  void ResetEventPlane()
  {
    m_psi_shifted_NS.clear();
    m_psi_NS.clear();
    m_psi_shifted_S.clear();
    m_psi_S.clear();
    m_psi_shifted_N.clear();
    m_psi_N.clear();
  }

 
  int GetZvtx( PHCompositeNode *topNode );
  int GetCentInfo( PHCompositeNode *topNode );
  int GetEventHeaderInfo( PHCompositeNode *topNode );
  int GetSub1JetInfo( PHCompositeNode *topNode );
  int GetTruthJetInfo( PHCompositeNode *topNode );
  int GetRawJetInfo( PHCompositeNode *topNode );
  int GetRhoInfo( PHCompositeNode *topNode );
  int GetEventPlaneInfo( PHCompositeNode *topNode );
  int GetCaloInfo( PHCompositeNode *topNode );
  void SumCaloE( PHCompositeNode *topNode, const std::string &towerinfo_node, const std::string &geo_node, const float zvrtx, float &sum_e );
  
  int GetG4TruthInfo( PHCompositeNode *topNode );  

};


#endif
