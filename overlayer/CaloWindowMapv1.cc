#include "CaloWindowMapv1.h"

#include <jetbase/Jet.h>

#include <iostream>
#include <utility> 
#include <cmath>
#include <algorithm>


CaloWindowMapv1::CaloWindowMapv1()
  : m_src(Jet::SRC::VOID)
  , m_nphi(0)
  , m_neta(0)
  , m_towers()
{
}

void CaloWindowMapv1::identify(std::ostream &os) const
{
    os << "CaloWindowMapv1: (nphi, neta) = (" << m_nphi << ", " << m_neta << ") with " << m_towers.size() << " towers" << std::endl;
}

int CaloWindowMapv1::isValid() const
{
    if (m_src == Jet::SRC::VOID || m_nphi == 0 || m_neta == 0) { return 0; }
    return 1;
}

void CaloWindowMapv1::add_tower(const unsigned int iphi, const unsigned int ieta, const float pt, bool is_masked)
{
    if (iphi >= m_nphi || ieta >= m_neta) {
        std::cerr << "CaloWindowMapv1::add_tower - invalid iphi or ieta: " 
                    << iphi << ", " << ieta << std::endl;
        return;
    }

    // Ensure the tower vector has the correct size
    if (m_towers.size() < m_nphi * m_neta) {
        m_towers.resize(m_nphi * m_neta, 0);
    }

    // Compute index and update tower value
    unsigned int index = iphi + ieta * m_nphi;
    m_towers[index] = is_masked ? kMASK_ENERGY : pt;
    return;
}

std::vector<float> CaloWindowMapv1::get_calo_windows(const unsigned int dphi, const unsigned int deta) const
{
    

    if (dphi == 0 || dphi > m_nphi || deta == 0 || deta > m_neta) {
        std::cout << "CaloWindowMapv1::get_calowindows - invalid dphi or deta: " << dphi << ", " << deta << std::endl;
        return std::vector<float>();
    }

    unsigned int num_windows = (m_neta - deta + 1) * m_nphi;
    std::vector<float> calo_windows(num_windows, 0);

    for (unsigned int window_idx = 0; window_idx < num_windows; ++window_idx) {
        unsigned int eta_start = window_idx / m_nphi;
        unsigned int phi_start = window_idx % m_nphi;

        bool is_masked = false;
        float window_sum = 0;

        for (unsigned int deta_idx = 0; deta_idx < deta; ++deta_idx) {
            unsigned int eta = eta_start + deta_idx;

            for (unsigned int dphi_idx = 0; dphi_idx < dphi; ++dphi_idx) {
                unsigned int phi = (phi_start + dphi_idx) % m_nphi;
                unsigned int tower_idx = phi + eta * m_nphi;

                if (m_towers[tower_idx] == kMASK_ENERGY) {
                    is_masked = true;
                }

                window_sum += m_towers[tower_idx];

            }
        }

        calo_windows[window_idx] = is_masked ? kMASK_ENERGY : window_sum;
    }

    return calo_windows;
}

std::vector< std::pair<unsigned int, unsigned int> > CaloWindowMapv1::get_window_comps_phieta(const unsigned int dphi, const unsigned int deta, const unsigned int iwindow) const 
{
    if (dphi == 0 || dphi > m_nphi || deta == 0 || deta > m_neta) {
        std::cerr << "CaloWindowMapv1::get_window_comps_phieta - invalid dphi or deta: " 
                  << dphi << ", " << deta << std::endl;
        return {};
    }

    unsigned int nwindows = (m_neta - deta + 1) * m_nphi;
    if (iwindow >= nwindows) {
        std::cerr << "CaloWindowMapv1::get_window_comps_phieta - invalid iwindow: " 
                  << iwindow << std::endl;
        return {};
    }

    // Calculate starting eta and phi indices
    unsigned int eta_start = iwindow / m_nphi;
    unsigned int phi_start = iwindow % m_nphi;

    // Populate window components
    std::vector<std::pair<unsigned int, unsigned int>> window_comps;
    window_comps.reserve(dphi * deta); // Reserve space to avoid reallocations

    for (unsigned int deta_idx = 0; deta_idx < deta; ++deta_idx) {
        unsigned int eta = eta_start + deta_idx;

        for (unsigned int dphi_idx = 0; dphi_idx < dphi; ++dphi_idx) {
            unsigned int phi = (phi_start + dphi_idx) % m_nphi;
            window_comps.emplace_back(phi, eta);
        }
    }

    return window_comps;
}

std::vector< std::pair<float, unsigned int> > CaloWindowMapv1::get_window_comps_energy_key(const unsigned int dphi, const unsigned int deta, const unsigned int iwindow) const 
{
        if (dphi == 0 || dphi > m_nphi || deta == 0 || deta > m_neta) {
        std::cerr << "CaloWindowMapv1::get_window_comps_phieta - invalid dphi or deta: " 
                  << dphi << ", " << deta << std::endl;
        return {};
    }

    unsigned int nwindows = (m_neta - deta + 1) * m_nphi;
    if (iwindow >= nwindows) {
        std::cerr << "CaloWindowMapv1::get_window_comps_phieta - invalid iwindow: " 
                  << iwindow << std::endl;
        return {};
    }

    // Calculate starting eta and phi indices
    unsigned int eta_start = iwindow / m_nphi;
    unsigned int phi_start = iwindow % m_nphi;

    // Populate window components
    std::vector<std::pair<float, unsigned int>> window_comps;
    window_comps.reserve(dphi * deta); // Reserve space to avoid reallocations

    for (unsigned int deta_idx = 0; deta_idx < deta; ++deta_idx) {
        unsigned int eta = eta_start + deta_idx;

        for (unsigned int dphi_idx = 0; dphi_idx < dphi; ++dphi_idx) {
            unsigned int phi = (phi_start + dphi_idx) % m_nphi;
            unsigned int key = phi + (eta << 16U);
            window_comps.emplace_back(m_towers[key], key);
        }
    }

    return window_comps;
}


