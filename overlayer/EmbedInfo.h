#ifndef _EMBEDINFO_H_
#define _EMBEDINFO_H_


#include <phool/PHObject.h>

#include <iostream>

class EmbedInfo : public PHObject
{
 public:


  ~EmbedInfo() override {};

  void identify(std::ostream &os = std::cout) const override { os << "EmbedInfo base class" << std::endl; };
  int isValid() const override { return 0; }

  // setters
  virtual void set_counter(int /*counter*/) {}
  virtual void set_embed_cent(int /*cent*/) {}

  virtual void set_embed_zvrtx(float /*zvrtx*/) {}
  virtual void set_embed_b(float /*b*/) {}

  virtual void set_embed_v2_w(float /*v2_w*/) {}
  virtual void set_embed_dpsi2(float /*dpsi2*/) {}
  virtual void set_embed_psi2(float /*phi*/) {}

  virtual void set_embed_leadpt(float /*pt*/) {}
  virtual void set_embed_leadphi(float /*phi*/) {}

  virtual void set_embed_ihcal_sumet(float /*sumet*/) {}
  virtual void set_embed_cemc_sumet(float /*sumet*/) {}
  virtual void set_embed_ohcal_sumet(float /*sumet*/) {}
  virtual void set_embed_ihcal_sumet_scaled(float /*sumet*/) {}
  virtual void set_embed_cemc_sumet_scaled(float /*sumet*/) {}
  virtual void set_embed_ohcal_sumet_scaled(float /*sumet*/) {}

  // getters
  virtual int get_counter() const { return -1; }
  virtual int get_embed_cent() const { return -1; }
  virtual float get_embed_zvrtx() const { return -9999; }
  virtual float get_embed_b() const { return -9999; }
  virtual float get_embed_psi2() const { return -9999; }  
  virtual float get_embed_v2_w() const { return -9999; }
  virtual float get_embed_dpsi2() const { return -9999; }
  virtual float get_embed_leadpt() const { return -9999; }
  virtual float get_embed_leadphi() const { return -9999; }
  virtual float get_embed_ihcal_sumet() const { return -9999; }
  virtual float get_embed_cemc_sumet() const { return -9999; }
  virtual float get_embed_ohcal_sumet() const { return -9999; }
  virtual float get_embed_ihcal_sumet_scaled() const { return -9999; }
  virtual float get_embed_cemc_sumet_scaled() const { return -9999; }
  virtual float get_embed_ohcal_sumet_scaled() const { return -9999; }

  

 protected:
  EmbedInfo() {}  // ctor

 private:
  ClassDefOverride(EmbedInfo, 1);
};

#endif  // _EMBEDINFO_H_
