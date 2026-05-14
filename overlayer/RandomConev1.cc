#include "RandomConev1.h"

#include <jetbase/Jet.h>

#include <iostream>
#include <vector>
#include <utility> 

void RandomConev1::identify(std::ostream& os) const
{
  os << "RandomConesv1: (eta, phi, R) = (" << m_feta << ", " << m_fphi << ", " << m_R << ")" << " n_comp = " << m_n_comps.size() << std::endl; 
  return ;
}

void RandomConev1::Reset()
{
  m_id = 0;
  m_feta = NAN;
  m_fphi = NAN;
  m_R = NAN;
  m_px_map.clear();
  m_py_map.clear();
  m_pz_map.clear();
  m_e_map.clear();
  m_n_comps.clear();
  m_n_masked_comps.clear();
  m_comp_src_vec.clear();
  return ;
}

int RandomConev1::isValid() const
{
  if ( std::isnan(m_feta) || std::isnan(m_fphi) || std::isnan(m_R) || std::isnan(m_R) || m_R <= 0 )
  {
    return 0;
  }
  return 1;
}

void RandomConev1::add_comp( Jet* particle , bool is_masked)
{
  if ( !particle )
  {
    std::cout << "Random Cone: add_comp - empty particle" << std::endl;
    return ;
  }


  auto comps = particle->get_comp_vec();
  for ( auto comp : comps ) 
  {

    if ( m_n_comps.find(comp.first) == m_n_comps.end() )
    { 
      m_n_comps[comp.first] = 0; 
      m_n_masked_comps[comp.first] = 0;
      m_px_map[comp.first] = 0;
      m_py_map[comp.first] = 0;
      m_pz_map[comp.first] = 0;
      m_e_map[comp.first] = 0;
    }

    float _px = particle->get_px();
    float _py = particle->get_py();
    float _pz = particle->get_pz();
    float _e = particle->get_e();
    if ( std::isnan(_px) || std::isnan(_py) || std::isnan(_pz) || std::isnan(_e) ) 
    {
      std::cout << "Random Cone: add_comp - particle with NaN in momentum or energy" << std::endl;
      continue;
    }
    if ( is_masked ) 
    { 
      _px = 0; 
      _py = 0; 
      _pz = 0; 
      _e = 0; 
      m_n_masked_comps[comp.first]++;
    }
    else 
    {
      m_comp_src_vec.push_back(comp); 
    }

    m_px_map[comp.first] += _px;
    m_py_map[comp.first] += _py;
    m_pz_map[comp.first] += _pz;
    m_e_map[comp.first] += _e;
    m_n_comps[comp.first]++;
    
  }

  return ;

}

void RandomConev1::add_tower(Jet::SRC src, float E, float eta, float phi, unsigned int ch, bool is_masked )
{
  if ( std::isnan(E) || std::isnan(eta) || std::isnan(phi) ) 
  {
    std::cout << "Random Cone: add_tower - tower with NaN in E, eta, or phi" << std::endl;
    return ;
  }
  if ( m_n_comps.find(src) == m_n_comps.end() )
  { 
    m_n_comps[src] = 0; 
    m_n_masked_comps[src] = 0;
    m_px_map[src] = 0;
    m_py_map[src] = 0;
    m_pz_map[src] = 0;
    m_e_map[src] = 0;
  }

  float pt = E / cosh(eta);
  float px = pt * cos(phi);
  float py = pt * sin(phi);
  float pz = pt * sinh(eta);
  if ( std::isnan(px) || std::isnan(py) || std::isnan(pz) ) 
  {
    std::cout << "Random Cone: add_tower - tower with NaN in px, py, or pz" << std::endl;
    return ;
  }

  if ( is_masked ) 
  { 
    px = 0; 
    py = 0; 
    pz = 0; 
    E = 0; 
    m_n_masked_comps[src]++;
  }
  else 
  {
    m_comp_src_vec.push_back(std::make_pair(src, ch)); 
  }

  m_px_map[src] += px;
  m_py_map[src] += py;
  m_pz_map[src] += pz;
  m_e_map[src] += E;
  m_n_comps[src]++;


  return ;
}

float RandomConev1::sum_px() const
{
  float sum = 0;
  for ( auto comp : m_px_map ) 
  {
    sum += comp.second;
  }
  return sum;
}

float RandomConev1::sum_py() const
{
  float sum = 0;
  for ( auto comp : m_py_map ) 
  {
    sum += comp.second;
  }
  return sum;
}

float RandomConev1::sum_pz() const
{
  float sum = 0;
  for ( auto comp : m_pz_map ) 
  {
    sum += comp.second;
  }
  return sum;
}

float RandomConev1::sum_e() const
{
  float sum = 0;
  for ( auto comp : m_e_map ) 
  {
    sum += comp.second;
  }
  return sum;
}

float RandomConev1::sum_p() const
{
  return std::sqrt( (sum_px() * sum_px()) + (sum_py() * sum_py()) + (sum_pz() * sum_pz()) );
}

float RandomConev1::sum_pt() const
{
  return std::sqrt( (sum_px() * sum_px()) + (sum_py() * sum_py()) );
}

float RandomConev1::sum_et() const
{
  float p = sum_p();
  if ( p == 0 ) 
  {
    return 0;
  }
  return sum_pt() / p * sum_e();
}

float RandomConev1::sum_eta() const
{
  float pt = sum_pt();
  if ( pt == 0 ) 
  {
    return NAN;
  }
  return std::asinh( sum_pz() / pt );
}

float RandomConev1::sum_phi() const
{
  return std::atan2( sum_py(), sum_px() );
}

float RandomConev1::get_pt(Jet::SRC src) const
{
  if ( src == Jet::SRC::VOID ) 
  {
    return sum_pt();
  } 
  else 
  {
    if ( m_px_map.find(src) != m_px_map.end() && m_py_map.find(src) != m_py_map.end() ) 
    {
      return std::sqrt( (m_px_map.at(src) * m_px_map.at(src)) + (m_py_map.at(src) * m_py_map.at(src)) );
    }
  }
  return 0;
}

float RandomConev1::get_eta(Jet::SRC src) const
{
  if ( src == Jet::SRC::VOID ) 
  {
    return sum_eta();
  } 
  else 
  {
    if ( m_px_map.find(src) != m_px_map.end() && m_py_map.find(src) != m_py_map.end() && m_pz_map.find(src) != m_pz_map.end() ) 
    {
      float pt = std::sqrt( (m_px_map.at(src) * m_px_map.at(src)) + (m_py_map.at(src) * m_py_map.at(src)) );
      if ( pt == 0 ) 
      {
        return NAN;
      }
      return std::asinh( m_pz_map.at(src) / pt );
    }
  }
  return NAN;
}

float RandomConev1::get_phi(Jet::SRC src) const
{
  if ( src == Jet::SRC::VOID ) 
  {
    return sum_phi();
  } 
  else 
  {
    if ( m_px_map.find(src) != m_px_map.end() && m_py_map.find(src) != m_py_map.end() ) 
    {
      return std::atan2( m_py_map.at(src), m_px_map.at(src) );
    }
  }
  return NAN;
}

float RandomConev1::get_px(Jet::SRC src) const
{
  if ( src == Jet::SRC::VOID ) 
  {
    return sum_px();
  } 
  else 
  {
    if ( m_px_map.find(src) != m_px_map.end() ) 
    {
      return m_px_map.at(src);
    }
  }
  return 0;
}

float RandomConev1::get_py(Jet::SRC src) const
{
  if ( src == Jet::SRC::VOID ) 
  {
    return sum_py();
  } 
  else 
  {
    if ( m_py_map.find(src) != m_py_map.end() ) 
    {
      return m_py_map.at(src);
    }
  }
  return 0;
}

float RandomConev1::get_pz(Jet::SRC src) const
{
  if ( src == Jet::SRC::VOID ) 
  {
    return sum_pz();
  } 
  else 
  {
    if ( m_pz_map.find(src) != m_pz_map.end() ) 
    {
      return m_pz_map.at(src);
    }
  }
  return 0;
}

float RandomConev1::get_e(Jet::SRC src) const
{
  if ( src == Jet::SRC::VOID ) 
  {
    return sum_e();
  } 
  else 
  {
    if ( m_e_map.find(src) != m_e_map.end() ) 
    {
      return m_e_map.at(src);
    }
  }
  return 0;
}

float RandomConev1::get_et(Jet::SRC src) const
{
  if ( src == Jet::SRC::VOID ) 
  {
    return sum_et();
  } 
  else 
  {
    float p = get_pz(src);
    if ( p == 0 ) 
    {
      return 0;
    }
    float pt = get_pt(src);
    float e = get_e(src);
    return pt / p * e;
  }
  return 0;
}

unsigned int RandomConev1::n_clustered(Jet::SRC src) const
{
  unsigned int n = 0;
  if ( src == Jet::SRC::VOID ) 
  {
    for (auto comp : m_n_comps)
    {
      n += comp.second;
    }
  } 
  else 
  {
    if ( m_n_comps.find(src) != m_n_comps.end() ) 
    {
      n = m_n_comps.at(src);
    }
  }
  return n;
}

unsigned int RandomConev1::n_masked(Jet::SRC src) const
{
  unsigned int n = 0;
  if ( src == Jet::SRC::VOID ) 
  {
    for (auto comp : m_n_masked_comps)
    {
      n += comp.second;
    }
  } 
  else 
  {
    if ( m_n_masked_comps.find(src) != m_n_masked_comps.end() ) 
    {
      n = m_n_masked_comps.at(src);
    }
  }
  return n;
}

float RandomConev1::masked_fraction(Jet::SRC src) const
{
  unsigned int nclus = n_clustered(src);
  if ( nclus == 0 ) 
  {
    return 0;
  }
  return static_cast<float>(n_masked(src)) / static_cast<float>(nclus);

}

std::vector<Jet::SRC> RandomConev1::get_src_vec()
{
  std::vector<Jet::SRC> srcs {};
  for (auto comp : m_n_comps)
  {
    srcs.push_back(comp.first);
  }
  return srcs;
}


std::vector<std::pair<Jet::SRC, unsigned int>> RandomConev1::get_comp_src_vec(Jet::SRC src) const
{
  std::vector<std::pair<Jet::SRC, unsigned int>> comp_srcs {};
  if ( src == Jet::SRC::VOID ) 
  {
    comp_srcs = m_comp_src_vec;
  } 
  else {
    for ( auto comp : m_comp_src_vec ) 
    {
      if ( comp.first == src ) 
      {
        comp_srcs.push_back(comp);
      }
    }
  }
  return comp_srcs;
}








