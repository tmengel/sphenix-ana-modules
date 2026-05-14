#ifndef UNDERLYINGEVENT_CALOWINDOWTOWERRECO_H
#define UNDERLYINGEVENT_CALOWINDOWTOWERRECO_H

//===========================================================
/// \file CaloWindowTowerReco.h 
/// \brief SubsysReco module to create a group of calo windows
/// \author Tanner Mengel
//===========================================================

#include "UEDefs.h"

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>
#include <iostream> 

#include <jetbase/Jet.h>

class PHCompositeNode;

class CaloWindowTowerReco : public SubsysReco
{
  public:
    
    CaloWindowTowerReco(const std::string &name = "CaloWindowTowerReco");
    ~CaloWindowTowerReco() override;

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

    void set_window_prefix( const std::string &prefix ) { m_window_prefix = prefix; }

  private:

    std::vector<std::string> m_inputs {}; // input node names
    std::vector<std::string> m_geom_names {}; // geometry node names
    std::vector<Jet::SRC> m_srcs {}; // source of input
    std::string m_window_prefix {"CaloWindowMap"}; // prefix for window nodes
    std::vector<std::string> m_window_names {}; // window names

    static const int neta_ihcal = 24;
    static const int neta_emcal = 96;
    static const int nphi_ihcal = 64;
    static const int nphi_emcal = 256;

    int CreateNodes(PHCompositeNode *topNode);

};

#endif // UNDERLYINGEVENT_CALOWINDOWTOWERRECO_H
