#ifndef _JET_TREE_H_
#define _JET_TREE_H_

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <string>
#include <vector>
#include <map>
#include <utility>

class PHCompositeNode;
class TTree;

class JetTree : public SubsysReco
{

 public:

  enum JET_TYPE
  {
    RAW = 0,
    SUB1 = 1,
    TRUTH = 2,
    SEED = 3
  };
  

  JetTree( const std::string & outputfile = "output.root" ) : SubsysReco("JetTree"), m_output_filename( outputfile ) {}
  ~JetTree() override {}

  int Init( PHCompositeNode * /*topNode*/) override;
  int process_event( PHCompositeNode * topNode ) override;
  int End( PHCompositeNode * /*topNode*/ ) override;
  int ResetEvent( PHCompositeNode * /*topNode*/ ) override;

  void set_zvrtx_node( const std::string & name = "GlobalVertexMap" ){ m_zvrtx_node = name; }
  
  void set_cent_node( const std::string & name = "CentralityInfo" ) { m_cent_node = name;  }
  
  void set_header_node( const std::string & name =  "EventHeader"){ m_header_node = name; }

  void add_calo_node( const std::string & name, const std::string & nickname ) 
  { 
    m_calo_nodes.push_back(name); 
    m_calo_nicknames_map[name] = nickname;
    m_calo_sumE_map[name] = 0.0;
    m_calo_sumEt_map[name] = 0.0;
    m_calo_ieta_avgE_map[name] = std::vector < float > (k_ieta, 0.0);
    m_calo_iphi_avgE_map[name] = std::vector < float > (k_iphi, 0.0);
    m_ncalo_node = m_calo_nodes.size();
  }

  void add_rho_node( const std::string & name , const std::string & nickname )
  {
    m_rho_nodes.push_back(name);
    m_rho_nicknames_map[name] = nickname;
    m_rho_map[name] = 0.0;
    m_rho_sigma_map[name] = 0.0;
    m_nrho_nodes = m_rho_nodes.size();
  }

  void add_towerbkgd( const std::string & name, const std::string & nickname )
  {
    m_towerbkgd_nodes.push_back(name);
    m_towerbkgd_nicknames_map[name] = nickname;
    m_towerbkgd_v2_map[name] = 0.0;
    m_towerbkgd_flowfail_map[name] = 0;
    m_towerbkgd_cemc_map[name] = std::vector < float > (k_ieta, 0.0);
    m_towerbkgd_hcalin_map[name] = std::vector < float > (k_ieta, 0.0);
    m_towerbkgd_hcalout_map[name] = std::vector < float > (k_ieta, 0.0);
    m_ntowerbkgd_nodes = m_towerbkgd_nodes.size();
  }    

  void add_jet_node( const std::string & name, const std::string & nickname, const JET_TYPE type, const float R )
  {
    m_jet_nodes.push_back(name);
    m_jet_nicknames_map[name] = nickname;
    m_jet_type_map[name] = type;
    m_jetR_map[name] = R;

    m_jet_minpT_map[name] = 0.0;
    m_jet_minE_map[name] = 0.0;
    m_jet_doetacut_map[name] = true;

    if ( type == SEED ) 
    {
      m_jet_minE_map[name] = -999.0; // for seeds, do not apply energy cut by default since they can be very low energy
      m_jet_doetacut_map[name] = false; // for seeds, do not apply eta cut by default since they can be outside of acceptance
    }
    
    // initialize jet info maps
    m_jet_E[name] = std::vector < float > {};
    m_jet_eta[name] = std::vector < float > {};
    m_jet_phi[name] = std::vector < float > {};
    m_jet_pT[name] = std::vector < float > {};

    m_jet_unsub_pT[name] = std::vector < float > {};
    m_jet_unsub_E[name] = std::vector < float > {};
    
    m_jet_area_layer[name] = std::vector < std::vector < float > > {};
    m_jet_E_layer[name] = std::vector < std::vector < float > > {};
    m_jet_N_layer[name] = std::vector < std::vector < int > > {};


    m_jet_maxD[name] = std::vector < float > {};
    m_jet_avgD[name] = std::vector < float > {};
    m_jet_supercomp_eT[name] = std::vector < std::vector < float > > {};
    
    m_jet_comp_ieta[name] = std::vector < std::vector < int > > {};
    m_jet_comp_iphi[name] = std::vector < std::vector < int > > {};
    m_jet_comp_caloid[name] = std::vector < std::vector < int > > {};

    m_jet_comp_E[name] = std::vector < std::vector < float > > {};
    m_jet_comp_eta[name] = std::vector < std::vector < float > > {};
    m_jet_comp_phi[name] = std::vector < std::vector < float > > {};
    m_njet_nodes = m_jet_nodes.size();
  }
  void set_minpt_jet_node( const std::string & name, const float minpT )
  {
    m_jet_minpT_map[name] = minpT;
  }
  void set_minE_jet_node( const std::string & name, const float minE )
  {
    m_jet_minE_map[name] = minE;
  }
  void set_etacut_jet_node( const std::string & name, const bool doetacut )
  {
    m_jet_doetacut_map[name] = doetacut;
  }

 private:
    
  std::string m_output_filename { "" };

  std::string m_zvrtx_node { "GlobalVertexMap" };
  float m_zvrtx { 0.0 };
  void _reset_zvrtx()
  {
    m_zrvtx = -999;
  }

  std::string m_cent_node { "" };
  int m_cent {-1};
  void _reset_cent()
  {
    m_cent = -1;
  }

  std::string m_header_node { "" };
  float m_b { 0.0 };
  float m_ep_angle { 0.0 };
  float m_ncoll { 0.0 };
  float m_npart { 0.0 };
  static const int k_max_psin = 6;
  float m_psin[k_max_psin];
  void _reset_header()
  {
    m_b = -999;
    m_ep_angle = -999;
    m_ncoll = -999;
    m_npart = -999;
    memset(m_psin, -999, sizeof(m_psin));
  }

  static const int k_ieta = 24;
  static const int k_iphi = 64;
  unsigned int m_ncalo_nodes { 0 };
  std::vector< std::string > m_calo_nodes {};
  std::map< std::string, std::string > m_calo_nicknames_map {};
  std::map< std::string, float > m_calo_sumE_map {};
  std::map< std::string, float > m_calo_sumEt_map {};
  std::map< std::string, std::vector < float > > m_calo_ieta_avgE_map {};
  std::map< std::string, std::vector < float > > m_calo_iphi_avgE_map {};
  void _reset_calo_maps()
  {
    for ( auto & calo : m_calo_nodes )
    {
      m_calo_sumE_map[calo] = 0.0;
      m_calo_sumEt_map[calo] = 0.0;
      m_calo_ieta_avgE_map[calo] = std::vector < float > (k_ieta, 0.0);
      m_calo_iphi_avgE_map[calo] = std::vector < float > (k_iphi, 0.0);
    }
  }

  unsigned int m_nrho_nodes { 0 };
  std::vector< std::string > m_rho_nodes {};
  std::map< std::string, std::string > m_rho_nicknames_map {};
  std::map< std::string, float > m_rho_map {};
  std::map< std::string, float > m_rho_sigma_map {};
  void _reset_rho_maps()
  {
    for ( auto & rho : m_rho_nodes )
    {
      m_rho_map[rho] = 0.0;
      m_rho_sigma_map[rho] = 0.0;
    }
  }

  unsigned int m_ntowerbkgd_nodes { 0 };
  std::vector< std::string > m_towerbkgd_nodes {};
  std::map< std::string, std::string > m_towerbkgd_nicknames_map {};
  std::map< std::string, std::vector < float > > m_towerbkgd_cemc_map {};
  std::map< std::string, std::vector < float > > m_towerbkgd_hcalin_map {};
  std::map< std::string, std::vector < float > > m_towerbkgd_hcalout_map {};
  std::map< std::string, float > m_towerbkgd_v2_map {};
  std::map< std::string, int > m_towerbkgd_flowfail_map {};
  void _reset_towerbkgd_maps()
  {
    for ( auto & towerbkgd : m_towerbkgd_nodes )
    {
      m_towerbkgd_cemc_map[towerbkgd] = std::vector < float > (k_ieta, 0.0);
      m_towerbkgd_hcalin_map[towerbkgd] = std::vector < float > (k_ieta, 0.0);
      m_towerbkgd_hcalout_map[towerbkgd] = std::vector < float > (k_ieta, 0.0);
      m_towerbkgd_v2_map[towerbkgd] = 0.0;
      m_towerbkgd_flowfail_map[towerbkgd] = 0;
    }
  }

  unsigned int m_njet_nodes { 0 };
  std::vector< std::string > m_jet_nodes {};
  std::map< std::string, std::string > m_jet_nicknames_map {};
  std::map< std::string, JetTree::JET_TYPE > m_jet_type_map {};
  std::map< std::string, float > m_jetR_map {};
  // kinematic cuts for jets
  std::map< std::string, float > m_jet_minpT_map {};
  std::map< std::string, float > m_jet_minE_map {};
  std::map< std::string, bool >  m_jet_doetacut_map {};
  // variables to be filled for jets
  std::map< std::string, std::vector < float > > m_jet_E {};
  std::map< std::string, std::vector < float > > m_jet_phi {};
  std::map< std::string, std::vector < float > > m_jet_eta {};
  std::map< std::string, std::vector < float > > m_jet_pT {};
  std::map< std::string, std::vector < float > > m_jet_unsub_pT {};
  std::map< std::string, std::vector < float > > m_jet_unsub_E {};
  std::map< std::string, std::vector < std::vector < float > > > m_jet_area_layer {};
  std::map< std::string, std::vector < std::vector < float > > > m_jet_E_layer {};
  std::map< std::string, std::vector < std::vector < int > > > m_jet_N_layer {};
  std::map< std::string, std::vector < float > > m_jet_maxD {};
  std::map< std::string, std::vector < float > > m_jet_avgD {};
  std::map< std::string, std::vector < std::vector < float > > > m_jet_supercomp_eT {};

  std::map< std::string, std::vector < std::vector < int > > > m_jet_comp_ieta {};
  std::map< std::string, std::vector < std::vector < int > > > m_jet_comp_iphi {};
  std::map< std::string, std::vector < std::vector < int > > > m_jet_comp_caloid {};

  std::map< std::string, std::vector < std::vector < float > > > m_jet_comp_E {};
  std::map< std::string, std::vector < std::vector < float > > > m_jet_comp_eta {};
  std::map< std::string, std::vector < std::vector < float > > > m_jet_comp_phi {};
  void _reset_jet_maps_for_node( std::string & jet_node )
  {
    m_jet_E[jet_node].clear();
    m_jet_phi[jet_node].clear();
    m_jet_eta[jet_node].clear();
    m_jet_pT[jet_node].clear();

    m_jet_unsub_pT[jet_node].clear();
    m_jet_unsub_E[jet_node].clear();
    
    m_jet_area_layer[jet_node].clear();
    m_jet_E_layer[jet_node].clear();
    m_jet_N_layer[jet_node].clear();

    m_jet_maxD[jet_node].clear();
    m_jet_avgD[jet_node].clear();
    m_jet_supercomp_eT[jet_node].clear();

    m_jet_comp_ieta[jet_node].clear();
    m_jet_comp_iphi[jet_node].clear();
    m_jet_comp_caloid[jet_node].clear();
    
    m_jet_comp_E[jet_node].clear();
    m_jet_comp_eta[jet_node].clear();
    m_jet_comp_phi[jet_node].clear();
  }
  void _reset_jet_maps()
  {
    for ( auto & jet_node : m_jet_nodes )
    {
      _reset_jet_maps_for_node( jet_node );
    }
  }

  TTree * m_tree { nullptr };
  int m_event_id {-1};
  void _reset_all()
  {
    m_event_id++;
    _reset_zvrtx();
    _reset_cent();
    _reset_header();
    _reset_calo_maps();
    _reset_rho_maps();
    _reset_towerbkgd_maps();
    _reset_jet_maps();
  }
 


  int GetZvtx( PHCompositeNode *topNode );
  int GetCentInfo( PHCompositeNode *topNode );
  int GetEventHeaderInfo( PHCompositeNode *topNode );
  int GetSub1JetInfo( PHCompositeNode *topNode , const int idx );
  int GetTruthJetInfo( PHCompositeNode *topNode , const int idx );
  int GetEventPlaneInfo( PHCompositeNode *topNode );
  int GetCaloInfo( PHCompositeNode *topNode );
  void SumCaloE( PHCompositeNode *topNode, const std::string &towerinfo_node, const std::string &geo_node, const float zvrtx, float &sum_e );
  
  

};


#endif
