#ifndef TreeWriter_H
#define TreeWriter_H

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

class TreeWriter : public SubsysReco
{
 public:

  // super simple tree maker

  TreeWriter( const std::string & outputfile = "output.root" ) 
    : SubsysReco("TreeWriter")
    , m_output_filename( outputfile )
  {}

  ~TreeWriter() override
  {
    if ( Verbosity() > 0 ) 
    {
      std::cout << "TreeWriter::~TreeWriter - done" << std::endl;
    }
    return;
  }


  int Init( PHCompositeNode * /*topNode*/) override;
  
  int process_event( PHCompositeNode * topNode ) override;
  
  int End( PHCompositeNode * /*topNode*/ ) override;
  
  int ResetEvent( PHCompositeNode * /*topNode*/ ) override 
  {
    ResetGL1();
    ResetZvtx();
    ResetCent();
    ResetEventHeader();
    ResetMbd();
    ResetCaloArrays();
    ResetCaloSub1Arrays();
    ResetRhoArrays();
    ResetTowerBkgd();
    ResetRawSeed();
    ResetRawJet();
    ResetSub1Jet();
    ResetTruthJet();
    ResetG4Truth();
    return Fun4AllReturnCodes::EVENT_OK;
  }

  void add_gl1_node ( const std::string & name = "GL1Packet"){  m_gl1_node = name; }
  void add_zvrtx_node ( const std::string & name = "GlobalVertexMap" ){ m_zvrtx_node = name; }
  void add_cent_node ( const std::string & name = "CentralityInfo" ) { m_cent_node = name;  }
  void add_mbd_node ( const std::string & name = "MbdOut" ) { m_mbd_node = name; }

  
  void add_sepd_node( const std::string &name = "TOWERINFO_CALIB_EPD" ) { m_sepd_node = name; }
  void add_ep_info ( const std::string &name = "EventPlaneInfo" ) { m_eventplane_node = name; }


  void add_cemc_node ( const std::string & name = "TOWERINFO_CALIB_CEMC_RETOWER" ) { m_cemc_node = name; }
  void add_hcalin_node ( const std::string & name = "TOWERINFO_CALIB_HCALIN" ) { m_hcalin_node = name; }
  void add_hcalout_node ( const std::string & name = "TOWERINFO_CALIB_HCALOUT" ) { m_hcalout_node = name; }
  void add_cemc_sub1_node ( const std::string & name = "TOWERINFO_CALIB_CEMC_RETOWER_SUB1" ) { m_cemc_sub1_node = name; }
  void add_hcalin_sub1_node ( const std::string & name = "TOWERINFO_CALIB_HCALIN_SUB1" ) { m_hcalin_sub1_node = name; }
  void add_hcalout_sub1_node ( const std::string & name = "TOWERINFO_CALIB_HCALOUT_SUB1" ) { m_hcalout_sub1_node = name; }
   
  void add_rho_node ( const std::string & name ) { m_rho_nodes.push_back(name); }
  void do_towerbkgd_nodes( const bool b ) { m_do_towerbkgd = b; }

  void add_rawseed_node ( const std::string & name = "AntiKt_TowerInfo_HIRecoSeedsRaw_r02" ) { m_rawseed_node = name; }
  void add_rawjet_node ( const std::string & name ) { m_rawjet_nodes.push_back(name); }
  void add_sub1jet_node ( const std::string & name ) { m_sub1jet_nodes.push_back(name); }


  void add_event_header( const std::string & name =  "EventHeader"){ m_eventhead_node = name; }
  void add_weight ( const float w = 1.0 ) { m_weight = w; }
  void add_g4truth_node( const std::string &name ) { m_g4truth_node = name; }
  void do_flow ( const bool doflow , const bool do_fluc, const float scale  ) { m_do_flow = doflow;  m_do_fluc = do_fluc;  m_flow_scale = scale; }
  void add_truthjet_node ( const std::string & name ) { m_truthjet_nodes.push_back(name); }





 private:
    
  // output file name
  std::string m_output_filename { "" };

  TRandom3 * m_rand { nullptr };

  TTree * m_run_tree {nullptr};
  int m_nevents {0};
  float m_weight {1.0};

  // event tree
  TTree * m_tree {nullptr};
  int m_event_id {-1};

  // gl1
  std::string m_gl1_node {""};
  uint64_t s_triggervec { 0 };
  uint64_t l_triggervec { 0 };
  void ResetGL1()
  {
    s_triggervec = 0;
    l_triggervec = 0;
  }

  // z vertex info
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

  // sepd
  std::string m_sepd_node { "" };
  std::vector< float > m_sepd_energy {};
  std::vector< int  > m_sepd_isgood {};
  std::vector< int  > m_sepd_arm {};
  std::vector< float > m_sepd_radius {};
  std::vector< float > m_sepd_phi {};
  void ResetSepd()
  {
    m_sepd_energy.clear();
    m_sepd_isgood.clear();
    m_sepd_arm.clear();
    m_sepd_radius.clear();
    m_sepd_phi.clear();
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

  // mbd
  std::string m_mbd_node { "" };
  const int k_mbd_south_code = 0;
  const int k_mbd_north_code = 1;
  float m_mbd_q_N { 0 };
  float m_mbd_q_S { 0 };
  float m_mbd_t_N { 0 };
  float m_mbd_t_S { 0 };
  float m_mbd_n_N { 0 };
  float m_mbd_n_S { 0 };
  void ResetMbd()
  {
    m_mbd_q_N = -999;
    m_mbd_q_S = -999;
    m_mbd_t_N = -999;
    m_mbd_t_S = -999;
    m_mbd_n_N = -999;
    m_mbd_n_S = -999;
  }

  // calo info
  static const int k_ieta = 24;
  static const int k_iphi = 64;
  std::string m_cemc_node { "" };
  std::string m_hcalin_node { "" };
  std::string m_hcalout_node { "" };

  float m_cemc_E[k_ieta][k_iphi] {};
  float m_hcalin_E[k_ieta][k_iphi] {};
  float m_hcalout_E[k_ieta][k_iphi] {};
  float m_cemc_isgood[k_ieta][k_iphi] {};
  float m_hcalin_isgood[k_ieta][k_iphi] {};
  float m_hcalout_isgood[k_ieta][k_iphi] {};
  float m_cemc_time[k_ieta][k_iphi] {};
  float m_hcalin_time[k_ieta][k_iphi] {};
  float m_hcalout_time[k_ieta][k_iphi] {};
  float m_cemc_eta[k_ieta][k_iphi] {};
  float m_hcalin_eta[k_ieta][k_iphi] {};
  float m_hcalout_eta[k_ieta][k_iphi] {};
  float m_cemc_phi[k_ieta][k_iphi] {};
  float m_hcalin_phi[k_ieta][k_iphi] {};
  float m_hcalout_phi[k_ieta][k_iphi] {};
  void ResetCaloArrays()
  {
    memset(m_cemc_E, 0, sizeof(m_cemc_E));
    memset(m_cemc_isgood, 0, sizeof(m_cemc_isgood));
    memset(m_cemc_time, 0, sizeof(m_cemc_time));
    memset(m_cemc_eta, 0, sizeof(m_cemc_eta));
    memset(m_cemc_phi, 0, sizeof(m_cemc_phi));
    memset(m_hcalin_E, 0, sizeof(m_hcalin_E));
    memset(m_hcalin_isgood, 0, sizeof(m_hcalin_isgood));
    memset(m_hcalin_time, 0, sizeof(m_hcalin_time));
    memset(m_hcalin_eta, 0, sizeof(m_hcalin_eta));
    memset(m_hcalin_phi, 0, sizeof(m_hcalin_phi));
    memset(m_hcalout_E, 0, sizeof(m_hcalout_E));
    memset(m_hcalout_isgood, 0, sizeof(m_hcalout_isgood));
    memset(m_hcalout_time, 0, sizeof(m_hcalout_time));
    memset(m_hcalout_eta, 0, sizeof(m_hcalout_eta));
    memset(m_hcalout_phi, 0, sizeof(m_hcalout_phi));
  }

  std::string m_cemc_sub1_node { "" };
  std::string m_hcalin_sub1_node { "" };
  std::string m_hcalout_sub1_node { "" };
  float m_cemc_sub1_E[k_ieta][k_iphi] {};
  float m_hcalin_sub1_E[k_ieta][k_iphi] {};
  float m_hcalout_sub1_E[k_ieta][k_iphi] {};
  float m_cemc_sub1_isgood[k_ieta][k_iphi] {};
  float m_hcalin_sub1_isgood[k_ieta][k_iphi] {};
  float m_hcalout_sub1_isgood[k_ieta][k_iphi] {};
  float m_cemc_sub1_time[k_ieta][k_iphi] {};
  float m_hcalin_sub1_time[k_ieta][k_iphi] {};
  float m_hcalout_sub1_time[k_ieta][k_iphi] {};
  float m_cemc_sub1_eta[k_ieta][k_iphi] {};
  float m_hcalin_sub1_eta[k_ieta][k_iphi] {};
  float m_hcalout_sub1_eta[k_ieta][k_iphi] {};
  float m_cemc_sub1_phi[k_ieta][k_iphi] {};
  float m_hcalin_sub1_phi[k_ieta][k_iphi] {};
  float m_hcalout_sub1_phi[k_ieta][k_iphi] {};
  void ResetCaloSub1Arrays()
  {
    memset(m_cemc_sub1_E, 0, sizeof(m_cemc_sub1_E));
    memset(m_cemc_sub1_isgood, 0, sizeof(m_cemc_sub1_isgood));
    memset(m_cemc_sub1_time, 0, sizeof(m_cemc_sub1_time));
    memset(m_cemc_sub1_eta, 0, sizeof(m_cemc_sub1_eta));
    memset(m_cemc_sub1_phi, 0, sizeof(m_cemc_sub1_phi));
    memset(m_hcalin_sub1_E, 0, sizeof(m_hcalin_sub1_E));
    memset(m_hcalin_sub1_isgood, 0, sizeof(m_hcalin_sub1_isgood));
    memset(m_hcalin_sub1_time, 0, sizeof(m_hcalin_sub1_time));
    memset(m_hcalin_sub1_eta, 0, sizeof(m_hcalin_sub1_eta));
    memset(m_hcalin_sub1_phi, 0, sizeof(m_hcalin_sub1_phi));
    memset(m_hcalout_sub1_E, 0, sizeof(m_hcalout_sub1_E));
    memset(m_hcalout_sub1_isgood, 0, sizeof(m_hcalout_sub1_isgood));
    memset(m_hcalout_sub1_time, 0, sizeof(m_hcalout_sub1_time));
    memset(m_hcalout_sub1_eta, 0, sizeof(m_hcalout_sub1_eta));
    memset(m_hcalout_sub1_phi, 0, sizeof(m_hcalout_sub1_phi));
  }

  // rho info
  static const int k_maxrho = 24;
  std::vector< std::string > m_rho_nodes {};
  float m_rho_val[k_maxrho] {};
  float m_std_rho_val[k_maxrho] {};
  void ResetRhoArrays()
  {
    memset(m_rho_val, -1.0, sizeof(m_rho_val));
    memset(m_std_rho_val, -1.0, sizeof(m_std_rho_val));
  }


  // tower background nodes
  std::string m_towerbkgd_node_sub1 { "TowerInfoBackground_Sub1" };
  std::string m_towerbkgd_node_sub2 { "TowerInfoBackground_Sub2" };
  bool m_do_towerbkgd = false;
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

  
  // seed info
  std::string m_rawseed_node { "" };
  std::vector < float > m_raw_seed_E {};
  std::vector < float > m_raw_seed_phi {};
  std::vector < float > m_raw_seed_eta {};
  std::vector < float > m_raw_seed_pT {};
  std::vector< float > m_mean_super_tower_E {};
  std::vector< float > m_max_super_tower_E {};
  std::vector < std::vector < int > >   m_raw_seed_comp_ieta {};
  std::vector < std::vector < int > >   m_raw_seed_comp_iphi {};
  std::vector < std::vector < int > >   m_raw_seed_comp_caloid {};
  std::vector < std::vector < int > >   m_raw_seed_comp_status {};
  std::vector < std::vector < float > > m_raw_seed_comp_E {};
  std::vector < std::vector < float > > m_raw_seed_comp_eta {};
  std::vector < std::vector < float > > m_raw_seed_comp_phi {};
  std::vector < std::vector < float > > m_raw_seed_super_tower_E {};
  void ResetRawSeed()
  {
    m_raw_seed_E.clear();
    m_raw_seed_phi.clear();
    m_raw_seed_eta.clear();
    m_raw_seed_pT.clear();
    m_mean_super_tower_E.clear();
    m_max_super_tower_E.clear();
    m_raw_seed_comp_ieta.clear();
    m_raw_seed_comp_iphi.clear();
    m_raw_seed_comp_caloid.clear();
    m_raw_seed_comp_status.clear();
    m_raw_seed_comp_E.clear();
    m_raw_seed_comp_eta.clear();
    m_raw_seed_comp_phi.clear();
    m_raw_seed_super_tower_E.clear();
  }

  // raw jet info
  std::vector< std::string > m_rawjet_nodes {};
  TTree * m_raw_jet_trees[16] = { nullptr };
  std::vector < float > m_raw_jet_E {};
  std::vector < float > m_raw_jet_phi {};
  std::vector < float > m_raw_jet_eta {};
  std::vector < float > m_raw_jet_pT {};
  std::vector < std::vector < int > > m_raw_jet_comp_ieta {};
  std::vector < std::vector < int > > m_raw_jet_comp_iphi {};
  std::vector < std::vector < int > > m_raw_jet_comp_caloid {};
  std::vector < std::vector < int > > m_raw_jet_comp_status {};
  std::vector < std::vector < float > > m_raw_jet_comp_E {};
  std::vector < std::vector < float > > m_raw_jet_comp_eta {};
  std::vector < std::vector < float > > m_raw_jet_comp_phi {};
  void ResetRawJet()
  {
    m_raw_jet_E.clear();
    m_raw_jet_phi.clear();
    m_raw_jet_eta.clear();
    m_raw_jet_pT.clear();
    m_raw_jet_comp_ieta.clear();
    m_raw_jet_comp_iphi.clear();
    m_raw_jet_comp_caloid.clear();
    m_raw_jet_comp_status.clear();
    m_raw_jet_comp_E.clear();
    m_raw_jet_comp_eta.clear();
    m_raw_jet_comp_phi.clear();
  }
 
  // sub1 jet info
  std::vector< std::string > m_sub1jet_nodes {};
  TTree * m_sub1_jet_trees[16] = { nullptr };
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

  // truth jet info
  std::vector< std::string > m_truthjet_nodes {};
  TTree * m_truth_jet_trees[16] = { nullptr };
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

  // truth info
  std::string m_g4truth_node { "" };
  bool m_do_flow { false };
  bool m_do_fluc { false };
  float m_flow_scale { 1.0 };

  float m_g4truth_zvtx { 0.0 };
  float m_g4truth_v2reco { 0.0 };
  float m_g4truth_v3reco { 0.0 };
  float m_g4truth_v4reco { 0.0 };
  float m_g4truth_v5reco { 0.0 };
  float m_g4truth_v6reco { 0.0 };
  std::vector < int > m_g4truth_id {};
  std::vector < int > m_g4truth_status {};
  std::vector < float > m_g4truth_px {};
  std::vector < float > m_g4truth_py {};
  std::vector < float > m_g4truth_pz {};
  std::vector < float > m_g4truth_E {};
  std::vector < float > m_g4truth_eta {};
  std::vector < float > m_g4truth_phi {};
  std::vector < float > m_g4truth_pT {};
  std::vector < float > m_g4truth_v2 {};
  std::vector < float > m_g4truth_v3 {};
  std::vector < float > m_g4truth_v4 {};
  std::vector < float > m_g4truth_v5 {};
  std::vector < float > m_g4truth_v6 {};
  void ResetG4Truth()
  {
    m_g4truth_id.clear();
    m_g4truth_status.clear();
    m_g4truth_px.clear();
    m_g4truth_py.clear();
    m_g4truth_pz.clear();
    m_g4truth_E.clear();
    m_g4truth_eta.clear();
    m_g4truth_phi.clear();
    m_g4truth_pT.clear();
    m_g4truth_v2.clear();
    m_g4truth_v3.clear();
    m_g4truth_v4.clear();
    m_g4truth_v5.clear();
    m_g4truth_v6.clear();
    m_g4truth_v2reco = -999;
    m_g4truth_v3reco = -999;
    m_g4truth_v4reco = -999;
    m_g4truth_v5reco = -999;
    m_g4truth_v6reco = -999;
    m_g4truth_zvtx = -999;
  }


  int GetGL1( PHCompositeNode *topNode );
  int GetZvtx( PHCompositeNode *topNode );
  int GetCentInfo( PHCompositeNode *topNode );
  int GetEventHeaderInfo( PHCompositeNode *topNode );
  int GetMbdInfo( PHCompositeNode *topNode );
  int GetRawCaloInfo( PHCompositeNode *topNode );
  int GetSubCaloInfo( PHCompositeNode *topNode );
  int GetRhoInfo( PHCompositeNode *topNode );
  int GetTowerBkgdInfo( PHCompositeNode *topNode );
  int GetSeedInfo( PHCompositeNode *topNode );
  int GetRawJetInfo( PHCompositeNode *topNode , const int idx );
  int GetSub1JetInfo( PHCompositeNode *topNode , const int idx );
  int GetTruthJetInfo( PHCompositeNode *topNode , const int idx );
  int GetG4TruthInfo( PHCompositeNode *topNode );
  int GetSepdInfo( PHCompositeNode *topNode );
  int GetEventPlaneInfo( PHCompositeNode *topNode );
  
  static void CalcFlow( float b , float eta, float pt, 
    float &v2 , float &v3, float &v4, float &v5, float &v6 ,
    TRandom3 * engine, bool do_fluc = false, float scale = 1.0);
 
};


#endif
