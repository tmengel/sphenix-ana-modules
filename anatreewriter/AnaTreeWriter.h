#ifndef ANATREEWRITER_H
#define ANATREEWRITER_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>
#include <array>
#include <cstdint>
#include <cmath>

class PHCompositeNode;
class TTree;

class AnaTreeWriter : public SubsysReco
{
 public:

  // super simple tree maker

  AnaTreeWriter( const std::string & outputfile = "output.root" );
  
  ~AnaTreeWriter() override {};

  int Init( PHCompositeNode * /*topNode*/) override;
  
  int process_event( PHCompositeNode * topNode ) override;
  
  int End( PHCompositeNode * /*topNode*/ ) override;
  
  int ResetEvent( PHCompositeNode * /*topNode*/ ) override;

  // add info to the tree ! default are not added
  void add_gl1_node ( const std::string & name = "GL1Packet"){  m_gl1_node = name; }
  void add_mbd_node ( const std::string & name = "MbdOut" ) { m_mbd_node = name; }
  void add_zvrtx_node ( const std::string & name = "GlobalVertexMap" ){ m_zvrtx_node = name; }
  void add_cent_node ( const std::string & name = "CentralityInfo" ) { m_cent_node = name; }
  
  void add_calo_node ( const std::string & name ) { m_calo_nodes.push_back(name); }  
  void set_min_tower_energy ( const float e ) { m_min_tower_energy = e; m_do_min_tower_energy = true; }

  void add_jet_node ( const std::string & name ) { m_jet_nodes.push_back(name); }
  void do_jet_calo_eta_cut( const bool b ) { m_do_online_eta_cut = b; }
  void set_min_jet_energy (const float e ) {  m_min_jet_energy = e, m_do_min_jet_energy = true; }

  void add_rho_node ( const std::string & name ) { m_rho_nodes.push_back(name); }

 private:
    
  // output file name
  std::string m_output_filename { "" };

  // node names
  std::string m_gl1_node {""};
  std::string m_mbd_node { "" };
  std::string m_zvrtx_node { "" };
  std::string m_cent_node { "" };
  std::vector< std::string > m_calo_nodes {};
  std::vector< std::string > m_jet_nodes {};
  std::vector< std::string > m_rho_nodes {};

  TTree * m_run_tree {nullptr};
  // event counters
  int m_nevents {0};
  std::array< uint64_t, 64 > m_run_trigger_status {};
  std::array< uint64_t, 64 > m_run_live_scalar {};
  std::array< uint64_t, 64 > m_run_scaled_scalar {};
  std::array< uint64_t, 64 > m_run_raw_scalar {};

  // event tree
  TTree * m_tree {nullptr};
  int m_event_id {-1};
  int m_gl1_tree_id {-1};
  int m_mbd_tree_id {-1};
  int m_rho_tree_id {-1};
  std::array<int, 8> m_calo_tree_id {};
  std::array<int, 8> m_jet_tree_id {};

  // centrality info
  int m_centrality {-1};

  // z vertex info
  float m_zvtx { 0.0 };
    
  // gl1
  TTree * m_gl1_tree { nullptr };
  uint64_t s_triggervec { 0 };
  uint64_t l_triggervec { 0 };
  std::array< bool, 64 > m_gl1_trigger_status {};
  std::array< uint64_t, 64 > m_gl1_live_scalar {};
  std::array< uint64_t, 64 > m_gl1_scaled_scalar {};
  std::array< uint64_t, 64 > m_gl1_raw_scalar {}; 

  // mbd
  TTree * m_mbd_tree { nullptr };
  const int k_mbd_south_code = 0;
  const int k_mbd_north_code = 1;
  float m_mbd_q_N { 0 };
  float m_mbd_q_S { 0 };
  float m_mbd_time_N { 0 };
  float m_mbd_time_S { 0 };
  float m_mbd_npmt_N { 0 };
  float m_mbd_npmt_S { 0 };

  // rho info
  TTree * m_rho_tree { nullptr };
  std::array< float, 20 > m_rho_val {};
  std::array< float, 20 > m_std_rho_val {};

  // calo info
  TTree * m_calo_trees[16] = { nullptr };
  float m_min_tower_energy { -10.0 };
  bool m_do_min_tower_energy { false };
  unsigned int m_num_towers { 0 };
  unsigned int m_num_towers_fired { 0 };
  unsigned int m_num_towers_dead { 0 };
  std::vector< float > m_tower_energy {};
  std::vector< float > m_tower_energy_time {};
  std::vector< uint32_t > m_tower_status_ieta_iphi {};

  // jet info
  bool m_do_min_jet_energy { false };
  float m_min_jet_energy { -10.0 };
  bool m_do_online_eta_cut { false };
  float m_jet_R { 0.0 };

  TTree * m_jet_trees[16] = { nullptr };
  int m_num_jets { 0 };
  std::vector < float > m_jet_energy {};
  std::vector < float > m_jet_phi {};
  std::vector < float > m_jet_eta {};
  std::vector < float > m_jet_pT {};
  std::vector < float > m_jet_eT {};
  std::vector < int > m_jet_ntowers {};
  std::vector < int > m_jet_ntowers_cemc {};
  std::vector < int > m_jet_ntowers_hcalin {};
  std::vector < int > m_jet_ntowers_hcalout {};
  std::vector < float > m_jet_energy_hcalin {};
  std::vector < float > m_jet_energy_hcalout {};
  std::vector < float > m_jet_energy_cemc {};

  int GetZvtx( PHCompositeNode *topNode );
  int GetGL1( PHCompositeNode *topNode );
  int GetCaloInfo( PHCompositeNode *topNode , const int idx );
  int GetJetInfo( PHCompositeNode *topNode , const int idx);
  int GetRhoInfo( PHCompositeNode *topNode );
  int GetMbdInfo( PHCompositeNode *topNode );
  int GetCentInfo( PHCompositeNode *topNode );
  
  const float radius_EM = 93.5;
  const float minz_EM = -130.23;
  const float maxz_EM = 130.23;

  const float radius_IH = 127.503;
  const float minz_IH = -170.299;
  const float maxz_IH = 170.299;

  const float radius_OH = 225.87;
  const float minz_OH = -301.683;
  const float maxz_OH = 301.683;

  float get_emcal_mineta_zcorrected( float zvertex ) {
    float z = minz_EM - zvertex;
    float eta_zcorrected = asinh(z / (float)radius_EM);
    return eta_zcorrected;
  }

  float get_emcal_maxeta_zcorrected(float zvertex) {
    float z = maxz_EM - zvertex;
    float eta_zcorrected = asinh(z / (float)radius_EM);
    return eta_zcorrected;
  }

  float get_ihcal_mineta_zcorrected(float zvertex) {
    float z = minz_IH - zvertex;
    float eta_zcorrected = asinh(z / (float)radius_IH);
    return eta_zcorrected;
  }

  float get_ihcal_maxeta_zcorrected(float zvertex) {
    float z = maxz_IH - zvertex;
    float eta_zcorrected = asinh(z / (float)radius_IH);
    return eta_zcorrected;
  }

  float get_ohcal_mineta_zcorrected(float zvertex) {
    float z = minz_OH - zvertex;
    float eta_zcorrected = asinh(z / (float)radius_OH);
    return eta_zcorrected;
  }

  float get_ohcal_maxeta_zcorrected(float zvertex) {
    float z = maxz_OH - zvertex;
    float eta_zcorrected = asinh(z / (float)radius_OH);
    return eta_zcorrected;
  }

  bool is_bad_jet_eta(float jet_eta, float zertex, float jet_radius) {
    float emcal_mineta = get_emcal_mineta_zcorrected(zertex);
    float emcal_maxeta = get_emcal_maxeta_zcorrected(zertex);
    float ihcal_mineta = get_ihcal_mineta_zcorrected(zertex);
    float ihcal_maxeta = get_ihcal_maxeta_zcorrected(zertex);
    float ohcal_mineta = get_ohcal_mineta_zcorrected(zertex);
    float ohcal_maxeta = get_ohcal_maxeta_zcorrected(zertex);
    float minlimit = emcal_mineta;
    if (ihcal_mineta > minlimit) { minlimit = ihcal_mineta; }
    if (ohcal_mineta > minlimit) { minlimit = ohcal_mineta; }
    float maxlimit = emcal_maxeta;
    if (ihcal_maxeta < maxlimit) { maxlimit = ihcal_maxeta; }
    if (ohcal_maxeta < maxlimit) { maxlimit = ohcal_maxeta; }
    minlimit += jet_radius;
    maxlimit -= jet_radius;
    return (jet_eta < minlimit || jet_eta > maxlimit);
  }

  bool accept_jet( float jet_energy, float jet_eta ) {
    if ( m_do_min_jet_energy && jet_energy < m_min_jet_energy ) { return false; }
    if ( m_do_online_eta_cut && is_bad_jet_eta(jet_eta, m_zvtx, m_jet_R) ) { return false; }
    return true;
  }


};


#endif
