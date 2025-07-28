#ifndef UNDERLYINGEVENT_RANDOMCONETOWERRECO_H
#define UNDERLYINGEVENT_RANDOMCONETOWERRECO_H

//===========================================================
/// \file RandomConeTowerReco.h
/// \brief UE background rho calculator
/// \author Tanner Mengel
//===========================================================

#include "UEDefs.h"

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>
#include <iostream> 
#include <cmath>
#include <cstdlib>

#include <jetbase/Jet.h>

class PHCompositeNode;
class TRandom3;

class RandomConeTowerReco : public SubsysReco
{
  public:
    
    RandomConeTowerReco(const std::string &name = "RandomConeTowerReco");
    ~RandomConeTowerReco() override;

    // standard Fun4All methods
    int Init(PHCompositeNode * topNode) override;
    int process_event(PHCompositeNode * topNode) override;

    void add_input( Jet::SRC src, const std::string & prefix = "TOWERINFO_CALIB" ) {
      std::string name = UEDefs::GetCaloTowerNode(src, prefix);
      std::string geom_name = UEDefs::GetCaloGeomNode(src);
      m_inputs.push_back(name);
      m_geom_names.push_back(geom_name);
      m_srcs.push_back(src);
    }
    void add_input( const std::string &name, const std::string &geom_name, Jet::SRC src ) {
      m_inputs.push_back(name);
      m_geom_names.push_back(geom_name);
      m_srcs.push_back(src);
    }
    void set_output_node( const std::string &name ) { m_output_node = name; }

    void set_R(const float R) { m_R = R; }
    void set_abs_eta(const float max_abs_eta){ m_max_abs_eta = max_abs_eta; }
    void set_masked_threshold(const float threshold){ m_masked_threshold = threshold; }
    void set_user_seed(const unsigned int seed){ m_seed = seed; }

    void set_avoid_lead_jet( const std::string &lead_jet_node, const float dR = -1.0 ) { 
      m_lead_jet_node = lead_jet_node; 
      m_avoid_lead_jet = true; 
      if (dR > 0.0) {
        m_lead_jet_dR = dR;
      } else {
        m_lead_jet_dR = 1.0 + m_R;
      }
    }


  private:

    std::vector<std::string> m_inputs {}; // input node names
    std::vector<std::string> m_geom_names {}; // geometry node names
    std::vector<Jet::SRC> m_srcs {}; // source of input
    std::string m_output_node;

    float m_R {0.4}; // cone radius
    
    float m_max_abs_eta{-1.0}; // max abs eta for cone
    float m_masked_threshold{0.0}; // threshold for masking towers

    bool m_avoid_lead_jet {false}; // avoid leading jet
    float m_lead_jet_dR {1.4}; // minimum dR between leading jet and random cone
    std::string m_lead_jet_node{""}; // node name of leading jet to avoid
    float m_lead_jet_eta {NAN}; // eta of leading jet
    float m_lead_jet_phi {NAN}; // phi of leading jet
    
    unsigned int m_seed{0}; 
    TRandom3 *m_random {nullptr};

    int CreateNode(PHCompositeNode *topNode);
    void GetConeAxis(PHCompositeNode *topNode, float &cone_eta, float &cone_phi); // get random cone axis
};

#endif // UNDERLYINGEVENT_RANDOMCONETOWERRECO_H
