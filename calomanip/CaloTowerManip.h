#ifndef CALOMANIP_CALOTOWERMANIP_H
#define CALOMANIP_CALOTOWERMANIP_H

#include <fun4all/SubsysReco.h>

#include <iostream>
#include <string>
#include <cmath>

class PHCompositeNode;
class TRandom3;

class CaloTowerManip : public SubsysReco
{
 public:

  CaloTowerManip(const std::string &name = "CaloTowerManip");
  ~CaloTowerManip() override {}

  int InitRun(PHCompositeNode * topNode) override;
  int process_event(PHCompositeNode * topNode) override;

  void SetInputNode( const std::string &name ) { m_input_node = name; }

  void SaveCopyOutputNode( const std::string &name = "" ) {
    m_output_node = name;
    m_save_original_towers = true;
    return;
  }

  void SetMinEnergy( float val ) {
    m_min_energy = val;
    m_do_min_energy_cut = true;
    return;
  }

  void RandomizeTowers( bool val ) { m_do_randomize_towers = val; }
  void SetRandomSeed( unsigned int seed ) {  m_random_seed = seed; }

 private:

  std::string m_input_node;

  bool m_save_original_towers{ false };
  std::string m_output_node;

  bool m_do_min_energy_cut{ false };
  float m_min_energy{ NAN };

  bool m_do_randomize_towers{ false };
  TRandom3 * m_random{ nullptr };
  unsigned int m_random_seed{ 0 };

  int CreateTowerInfoContainerCopy( PHCompositeNode *topNode );

};

#endif  // PPG04BASE_CALOTOWERMANIP_H
