#ifndef UNDERLYINGEVENT_CALOWINDOWMAP_H
#define UNDERLYINGEVENT_CALOWINDOWMAP_H

//===========================================================
/// \file CaloWindowMap.h
/// \brief PHObject for Calo Windows. Just a vector of calo energies dynamically sectioned into windows
/// \author Tanner Mengel
//===========================================================

#include <phool/PHObject.h>

#include <iostream>
#include <vector>
#include <utility>

#include <jetbase/Jet.h>

class CaloWindowMap : public PHObject
{
  public:

    CaloWindowMap() {};
    ~CaloWindowMap() override{};

    // const expression for masked tower energy
    constexpr static const float kMASK_ENERGY {-99999.0};

    void identify(std::ostream &os = std::cout) const override { os << "CaloWindowMap base class" << std::endl; };
    int isValid() const override { return 0; }

    virtual void set_src(Jet::SRC /*src*/) { return; }
    virtual void set_nphi_neta (unsigned int /*nphi*/, unsigned int /*neta*/) { return; }

    virtual Jet::SRC get_src() const { return Jet::SRC::VOID; }
    virtual unsigned int get_nphi() const { return 0; }
    virtual unsigned int get_neta() const { return 0; }
    virtual unsigned int get_ntowers() const { return 0; }

    virtual void clear_towers() { return; }
    virtual void add_tower(const unsigned int /*iphi*/, const unsigned int /*ieta*/, const float /*pt*/, bool /*is_masked*/) { return; }

    virtual std::vector<float> get_calo_windows(const unsigned int /*dphi*/, const unsigned int /*deta*/) const { return {}; }
    virtual std::vector< std::pair<unsigned int, unsigned int> > get_window_comps_phieta(const unsigned int /*dphi*/, const unsigned int /*deta*/, const unsigned int /*iwindow*/) const { return {}; }
    virtual std::vector< std::pair<float, unsigned int> > get_window_comps_energy_key(const unsigned int /*dphi*/, const unsigned int /*deta*/, const unsigned int /*iwindow*/) const { return {}; }  

    ClassDefOverride(CaloWindowMap, 1);
};

#endif // UNDERLYINGEVENT_CALOWINDOW_H
