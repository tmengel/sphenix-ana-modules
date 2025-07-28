#ifndef UNDERLYINGEVENT_CALOWINDOWMAPV1_H
#define UNDERLYINGEVENT_CALOWINDOWMAPV1_H

//===========================================================
/// \file CaloWindowMapv1.h
/// \brief PHObject for Calo Windows. Just a vector of calo energies dynamically sectioned into windows
/// \author Tanner Mengel
//===========================================================

#include "CaloWindowMap.h"

class CaloWindowMapv1 : public CaloWindowMap
{
  public:

    CaloWindowMapv1();
    ~CaloWindowMapv1() override {};

    void identify(std::ostream &os = std::cout) const override;
    void Reset() override { clear_towers(); }
    int isValid() const override;
    
    void set_src(Jet::SRC src) override { m_src = src; }
    void set_nphi_neta(unsigned int nphi, unsigned int neta) override { m_nphi = nphi; m_neta = neta; clear_towers(); }

    Jet::SRC get_src() const override { return m_src; }
    unsigned int get_nphi() const override { return m_nphi; }
    unsigned int get_neta() const override { return m_neta; }
    unsigned int get_ntowers() const override { return m_nphi*m_neta; }

    void clear_towers() override { m_towers.clear(); m_towers.resize(m_nphi*m_neta, 0); }
    void add_tower(const unsigned int iphi, const unsigned int ieta, const float pt, bool is_masked) override;

    std::vector<float> get_calo_windows(const unsigned int dphi, const unsigned int deta) const override;
    std::vector< std::pair<unsigned int, unsigned int> > get_window_comps_phieta(const unsigned int dphi, const unsigned int deta, const unsigned int iwindow) const override;
    std::vector< std::pair<float, unsigned int> > get_window_comps_energy_key(const unsigned int dphi, const unsigned int deta, const unsigned int iwindow) const override;

  private:

    Jet::SRC m_src {Jet::SRC::VOID};
    unsigned int m_nphi {0};
    unsigned int m_neta {0};

    std::vector<float> m_towers {};

    ClassDefOverride(CaloWindowMapv1, 1);
};


#endif // UNDERLYINGEVENT_CALOWINDOWMAPV1_H
