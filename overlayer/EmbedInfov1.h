#ifndef _EMBEDINFOV1_H_
#define _EMBEDINFOV1_H_


#include "EmbedInfo.h"

class EmbedInfov1 : public EmbedInfo
{
 public:
  EmbedInfov1();
  ~EmbedInfov1() override = default;

  void identify(std::ostream &os = std::cout) const override;
  void Reset() override {}  // new in v1
  int isValid() const override { return 1; }

  // // setters
  void set_counter(int counter) override { _counter = counter; }
  void set_embed_cent(int cent) override { _embed_cent = cent; }
  void set_embed_zvrtx(float zvrtx) override { _embed_zvrtx = zvrtx; }
  void set_embed_b(float b) override { _embed_b = b; }
  void set_embed_psi2(float psi2) override { _embed_psi2 = psi2; }  
  void set_embed_v2_w(float v2_w) override { _embed_v2_w = v2_w; }
  void set_embed_dpsi2(float dpsi2) override { _embed_dpsi2 = dpsi2; }
  void set_embed_leadpt(float pt) override { _embed_leadpt = pt; }
  void set_embed_leadphi(float phi) override { _embed_leadphi = phi; }
  void set_embed_ihcal_sumet(float sumet) override { _embed_ihcal_sumet = sumet; }
  void set_embed_cemc_sumet(float sumet) override { _embed_cemc_sumet = sumet; }
  void set_embed_ohcal_sumet(float sumet) override { _embed_ohcal_sumet = sumet; }
  void set_embed_ihcal_sumet_scaled(float sumet) override { _embed_ihcal_sumet_scaled = sumet; }
  void set_embed_cemc_sumet_scaled(float sumet) override { _embed_cemc_sumet_scaled = sumet; }
  void set_embed_ohcal_sumet_scaled(float sumet) override { _embed_ohcal_sumet_scaled = sumet; }
  // // getters
  int get_counter() const override { return _counter; }
  int get_embed_cent() const override { return _embed_cent; }
  float get_embed_zvrtx() const override { return _embed_zvrtx; }
  float get_embed_b() const override { return _embed_b; }
  float get_embed_psi2() const override { return _embed_psi2; }  
  float get_embed_v2_w() const override { return _embed_v2_w; }
  float get_embed_dpsi2() const override { return _embed_dpsi2; }
  float get_embed_leadpt() const override { return _embed_leadpt; }
  float get_embed_leadphi() const override { return _embed_leadphi; }
  float get_embed_ihcal_sumet() const override { return _embed_ihcal_sumet; }
  float get_embed_cemc_sumet() const override { return _embed_cemc_sumet; }
  float get_embed_ohcal_sumet() const override { return _embed_ohcal_sumet; }
  float get_embed_ihcal_sumet_scaled() const override { return _embed_ihcal_sumet_scaled; }
  float get_embed_cemc_sumet_scaled() const override { return _embed_cemc_sumet_scaled; }
  float get_embed_ohcal_sumet_scaled() const override { return _embed_ohcal_sumet_scaled; } 
  private:
  int _counter;
  int _embed_cent;
  float _embed_zvrtx;
  float _embed_b;
  float _embed_psi2;  
  float _embed_v2_w;
  float _embed_dpsi2;
  float _embed_leadpt;
  float _embed_leadphi;
  float _embed_ihcal_sumet;
  float _embed_cemc_sumet;
  float _embed_ohcal_sumet;
  float _embed_ihcal_sumet_scaled;
  float _embed_cemc_sumet_scaled;
  float _embed_ohcal_sumet_scaled;

  

  ClassDefOverride(EmbedInfov1, 1);
};

#endif  // _EMBEDINFOV1_H_
