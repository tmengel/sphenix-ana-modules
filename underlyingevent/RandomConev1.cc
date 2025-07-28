#include "RandomConev1.h"

#include <jetbase/Jet.h>

#include <iostream>
#include <vector>
#include <utility> 

void RandomConev1::identify(std::ostream& os) const
{
  os << "RandomConesv1: (pt, eta, phi, R) = (" << m_pt << ", " << m_eta << ", " << m_phi << ", " << m_R << ")" << std::endl;
  return ;
}

void RandomConev1::Reset()
{
  m_id = 0;
  m_pt = 0;
  m_eta = NAN;
  m_phi = NAN;
  m_R = NAN;
  n_comps.clear();
  n_masked_comps.clear();
  m_pt_map.clear();
  m_comp_src_vec.clear();
  return ;
}

int RandomConev1::isValid() const
{
  if ( std::isnan(m_pt) || std::isnan(m_eta) || std::isnan(m_phi) || std::isnan(m_R) ) {
    return 0;
  }
  return 1;
}

void RandomConev1::add_comp( Jet* particle , bool is_masked)
{
  if ( !particle ) {
    std::cout << "Random Cone: add_comp - empty particle" << std::endl;
    return ;
  }

  float pt = particle->get_e()/cosh(particle->get_eta());
  m_pt += is_masked ? 0 : pt;

  auto comps = particle->get_comp_vec();
  for ( auto comp : comps ) {
    if ( n_comps.find(comp.first) == n_comps.end() ) { n_comps[comp.first] = 0; }
    if ( m_pt_map.find(comp.first) == m_pt_map.end() ) { m_pt_map[comp.first] = 0; }
    if ( n_masked_comps.find(comp.first) == n_masked_comps.end() ) { n_masked_comps[comp.first] = 0; }
    n_comps[comp.first]++;
    m_pt_map[comp.first] += is_masked ? 0 : pt;
    if ( is_masked ) { n_masked_comps[comp.first]++; }
    else { m_comp_src_vec.push_back(comp); }
  }
  return ;

}

void RandomConev1::add_tower(Jet::SRC src, float pt, unsigned int ch, bool is_masked)
{
  if ( n_comps.find(src) == n_comps.end() ) { n_comps[src] = 0; }
  if ( m_pt_map.find(src) == m_pt_map.end() ) { m_pt_map[src] = 0; }
  if ( n_masked_comps.find(src) == n_masked_comps.end() ) { n_masked_comps[src] = 0; }

  n_comps[src]++;
  m_pt += is_masked ? 0 : pt;
  m_pt_map[src] += is_masked ? 0 : pt;
  if ( is_masked ) { n_masked_comps[src]++; }
  else { m_comp_src_vec.push_back(std::make_pair(src, ch)); }
  return ;
}


float RandomConev1::get_pt(Jet::SRC src) const
{
  if ( src == Jet::SRC::VOID ) {
    return m_pt;
  } else {
    if ( m_pt_map.find(src) != m_pt_map.end() ) {
      return m_pt_map.at(src);
    }
  }
  return 0;
}

unsigned int RandomConev1::n_clustered(Jet::SRC src) const
{
  unsigned int n = 0;
  if ( src == Jet::SRC::VOID ) {
    for (auto comp : n_comps) {
      n += comp.second;
    }
  } else {
    if ( n_comps.find(src) != n_comps.end() ) {
      n = n_comps.at(src);
    }
  }
  return n;
}

unsigned int RandomConev1::n_masked(Jet::SRC src) const
{
  unsigned int n = 0;
  if ( src == Jet::SRC::VOID ) {
    for (auto comp : n_masked_comps) {
      n += comp.second;
    }
  } else {
    if ( n_masked_comps.find(src) != n_masked_comps.end() ) {
      n = n_masked_comps.at(src);
    }
  }
  return n;
}

float RandomConev1::masked_fraction(Jet::SRC src) const
{
  unsigned int n_masked = 0;
  unsigned int n_total = 0;
  if ( src == Jet::SRC::VOID ) {
    for (auto comp : n_masked_comps) {
      n_masked += comp.second;
    }
    for (auto comp : n_comps) {
      n_total += comp.second;
    }
  } else {
    if ( n_masked_comps.find(src) != n_masked_comps.end() ) {
      n_masked = n_masked_comps.at(src);
    }
    if ( n_comps.find(src) != n_comps.end() ) {
      n_total = n_comps.at(src);
    }
  }
  if ( n_total == 0 ) {
    return 0;
  }
  return static_cast<float>(n_masked)/static_cast<float>(n_total);
}

std::vector<Jet::SRC> RandomConev1::get_src_vec()
{
  std::vector<Jet::SRC> srcs {};
  for (auto comp : n_comps) {
    srcs.push_back(comp.first);
  }
  return srcs;
}


std::vector<std::pair<Jet::SRC, unsigned int>> RandomConev1::get_comp_src_vec(Jet::SRC src) const
{
  std::vector<std::pair<Jet::SRC, unsigned int>> comp_srcs {};
  if ( src == Jet::SRC::VOID ) {
    comp_srcs = m_comp_src_vec;
  } else {
    for ( auto comp : m_comp_src_vec ) {
      if ( comp.first == src ) {
        comp_srcs.push_back(comp);
      }
    }
  }
  return comp_srcs;
}








