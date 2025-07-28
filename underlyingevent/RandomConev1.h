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
    void set_eta(float eta) override { m_eta = eta; }
    void set_phi(float phi) override { m_phi = phi; }
    void set_R(float R) override { m_R = R; }

    unsigned int get_id() const override { return m_id; }
    float get_pt(Jet::SRC src = Jet::SRC::VOID) const override;
    float get_eta() const override { return m_eta; }
    float get_phi() const override { return m_phi; }
    float get_R() const override { return m_R; }

    void add_comp( Jet *particle , bool is_masked = false) override;
    void add_tower( Jet::SRC src, float pt, unsigned int ch, bool is_masked = false ) override;

    unsigned int n_clustered( Jet::SRC src = Jet::SRC::VOID) const override;
    unsigned int n_masked( Jet::SRC src = Jet::SRC::VOID) const override;
    float masked_fraction( Jet::SRC src = Jet::SRC::VOID) const override;
    std::vector<Jet::SRC> get_src_vec() override;
    std::vector < std::pair < Jet::SRC, unsigned int > > get_comp_src_vec( Jet::SRC src = Jet::SRC::VOID) const override;

    
  private:

    unsigned int m_id {0};
    float m_pt {0};
    std::map < Jet::SRC, float > m_pt_map {};
    float m_eta {NAN};
    float m_phi {NAN};
    float m_R {NAN};

    std::map< Jet::SRC, unsigned int> n_comps {};
    std::map< Jet::SRC, unsigned int> n_masked_comps {};
    std::vector < std::pair < Jet::SRC, unsigned int > > m_comp_src_vec {};

    ClassDefOverride(RandomConev1, 1);
};

#endif // UNDERLYINGEVENT_RANDOMCONEV1_H
