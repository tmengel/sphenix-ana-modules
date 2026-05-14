#ifndef UNDERLYINGEVENT_RANDOMCONEV1_H
#define UNDERLYINGEVENT_RANDOMCONEV1_H

#include "RandomCone.h"

#include <cmath>
#include <map>

class RandomConev1 : public RandomCone
{
  public:

    RandomConev1() {};
    ~RandomConev1() override {}

    void identify(std::ostream &os = std::cout) const override;
    void Reset() override; 
    int isValid() const override;

    void set_id(unsigned int id) override { m_id = id; }
    void set_R(float R) override { m_R = R; }
    void set_fixed_eta(float eta) override { m_feta = eta; }
    void set_fixed_phi(float phi) override { m_fphi = phi; }

    unsigned int get_id() const override { return m_id; }
    float get_R() const override { return m_R; }
    float get_fixed_eta() const override { return m_feta; }
    float get_fixed_phi() const override { return m_fphi; }

    void add_comp( Jet * particle , bool is_masked = false) override;
    void add_tower( Jet::SRC src, float E, float eta, float phi, unsigned int ch, bool is_masked = false ) override;

    float get_pt(Jet::SRC src = Jet::SRC::VOID) const override;
    float get_eta(Jet::SRC src = Jet::SRC::VOID) const override;
    float get_phi(Jet::SRC src = Jet::SRC::VOID) const override;
    float get_px(Jet::SRC src = Jet::SRC::VOID) const override;
    float get_py(Jet::SRC src = Jet::SRC::VOID) const override;
    float get_pz(Jet::SRC src = Jet::SRC::VOID) const override;
    float get_e(Jet::SRC src = Jet::SRC::VOID) const override;
    float get_et(Jet::SRC src = Jet::SRC::VOID) const override;

    float sum_px() const override;
    float sum_py() const override;
    float sum_pz() const override;
    float sum_e() const override;
    float sum_pt() const override;
    float sum_phi() const override;
    float sum_eta() const override;
    float sum_et() const override;
    float sum_p() const override;

    unsigned int n_clustered( Jet::SRC src = Jet::SRC::VOID) const override;
    unsigned int n_masked( Jet::SRC src = Jet::SRC::VOID) const override;
    float masked_fraction( Jet::SRC src = Jet::SRC::VOID) const override;
    std::vector<Jet::SRC> get_src_vec() override;
    std::vector < std::pair < Jet::SRC, unsigned int > > get_comp_src_vec( Jet::SRC src = Jet::SRC::VOID) const override;

    
  private:

    unsigned int m_id {0};
    std::map < Jet::SRC, float > m_px_map {};
    std::map < Jet::SRC, float > m_py_map {};
    std::map < Jet::SRC, float > m_pz_map {};
    std::map < Jet::SRC, float > m_e_map {};
    float m_feta {NAN};
    float m_fphi {NAN};
    float m_R {NAN};

    std::map< Jet::SRC, unsigned int> m_n_comps {};
    std::map< Jet::SRC, unsigned int> m_n_masked_comps {};
    std::vector < std::pair < Jet::SRC, unsigned int > > m_comp_src_vec {};

    ClassDefOverride(RandomConev1, 1);
};

#endif // UNDERLYINGEVENT_RANDOMCONEV1_H
