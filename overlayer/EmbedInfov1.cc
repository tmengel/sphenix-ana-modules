#include "EmbedInfov1.h"

#include <cstdlib>

EmbedInfov1::EmbedInfov1()
{
}

void EmbedInfov1::identify(std::ostream& os) const
{
  os << "EmbedInfov1: " << std::endl;
  os << " counter: " << _counter << std::endl;
  os << " embed_cent: " << _embed_cent << std::endl;
  os << " embed_zvrtx: " << _embed_zvrtx << std::endl;
  os << " embed_b: " << _embed_b << std::endl;
  os << " embed_psi2: " << _embed_psi2 << std::endl;
  os << " embed_v2_w: " << _embed_v2_w << std::endl;
  os << " embed_dpsi2: " << _embed_dpsi2 << std::endl;
  os << " embed_leadpt: " << _embed_leadpt << std::endl;
  os << " embed_leadphi: " << _embed_leadphi << std::endl;
  os << " embed_ihcal_sumet: " << _embed_ihcal_sumet << std::endl;
  os << " embed_cemc_sumet: " << _embed_cemc_sumet << std::endl;
  os << " embed_ohcal_sumet: " << _embed_ohcal_sumet << std::endl;
  os << " embed_ihcal_sumet_scaled: " << _embed_ihcal_sumet_scaled << std::endl;
  os << " embed_cemc_sumet_scaled: " << _embed_cemc_sumet_scaled << std::endl;
  os << " embed_ohcal_sumet_scaled: " << _embed_ohcal_sumet_scaled << std::endl;
}
