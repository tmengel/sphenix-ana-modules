#ifndef UNDERLYINGEVENT_RANDOMCONE_H
#define UNDERLYINGEVENT_RANDOMCONE_H

//===========================================================
/// \file RandomCone.h
/// \brief PHObject for Random Cone
/// \author Tanner Mengel
//===========================================================

#include <phool/PHObject.h>

#include <iostream>
#include <vector>
#include <utility>

#include <jetbase/Jet.h>

class RandomCone : public PHObject
{
  public:

    ~RandomCone() override{};

    void identify(std::ostream &os = std::cout) const override { os << "RandomCone base class" << std::endl; };
    int isValid() const override { return 0; }

    virtual void set_id (unsigned int /*id*/) { return; }
    virtual void set_R(float /*R*/) { return; }
    virtual void set_eta(float /*eta*/) { return; }
    virtual void set_phi(float /*phi*/) { return; }
    
    virtual unsigned int get_id() const { return 0; }
    virtual float get_pt(Jet::SRC /*src*/ = Jet::SRC::VOID) const { return 0; }
    virtual float get_eta() const { return 0; }
    virtual float get_phi() const { return 0; }
    virtual float get_R() const { return 0; }

    virtual void add_comp (Jet * /*jet*/, bool /*is_masked*/ = false) { return; }
    virtual void add_tower(Jet::SRC /*src*/, float /*pt*/, unsigned int /*ch*/,  bool /*is_masked*/ = false) { return; }

    virtual unsigned int n_clustered(Jet::SRC /*src*/ = Jet::SRC::VOID) const { return 0; }
    virtual unsigned int n_masked(Jet::SRC /*src*/ = Jet::SRC::VOID) const { return 0; }
    virtual float masked_fraction(Jet::SRC /*src*/ = Jet::SRC::VOID) const { return 0; }
    virtual std::vector<Jet::SRC> get_src_vec() { return {}; }
    virtual std::vector<std::pair<Jet::SRC, unsigned int>> get_comp_src_vec(Jet::SRC /*src*/ = Jet::SRC::VOID) const { return {}; }
      
  protected:
    
    RandomCone() {} // ctor

  private:
    
    ClassDefOverride(RandomCone, 1);
};

#endif // UNDERLYINGEVENT_RANDOMCONE_H
