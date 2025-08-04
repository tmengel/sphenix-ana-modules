#include "MissingSebFilter.h"

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

// tower info
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>

#include <cstdlib>

MissingSebFilter::MissingSebFilter() : EventCut("MissingSebFilter") 
{ 
    SetNodeNames({ "TOWERINFO_CALIB_CEMC", "TOWERINFO_CALIB_HCALIN", "TOWERINFO_CALIB_HCALOUT" }); // Set the default node names
}


void MissingSebFilter::identify(std::ostream &os) const
{
    os << Name() + "::identify: " << std::endl;
    for (const auto &node_name : GetNodeNames())
    {
        os << "  Node name: " << node_name << std::endl;
    }
    return;
}

bool MissingSebFilter::operator()(PHCompositeNode *topNode)
{

    for (const auto &node_name : GetNodeNames())
    {
        auto * towers = findNode::getClass<TowerInfoContainer>(topNode, node_name);
        if (!towers) {
            std::cerr << Name() + "::operator(PHCompositeNode *topNode) Could not find " << node_name << " node" << std::endl;
            exit(-1); // this is a fatal error
        }

        // std::vector<int> dead_channels = ; 
        unsigned int nChannels = towers->size();
        unsigned int nDead = 0;
        for(unsigned int channel = 0; channel < nChannels; channel++){
            auto * tower = towers->get_tower_at_channel(channel);
            float energy = tower->get_energy();
            if ( tower->get_isHot() || tower->get_isNoCalib() || tower->get_isNotInstr() || tower->get_isBadChi2() || std::isnan(energy) ) {
                nDead++;
            }
        }
        Passed(nDead < (nChannels/16)); // if more than 1/16th of the channels are dead, we fail the cut
        if (!Passed()){
            if(Verbosity()){
                std::cout << Name() + "::operator(PHCompositeNode *topNode) Node " << node_name << " failed, " 
                          << nDead << " dead channels out of " << nChannels << std::endl;
            }
            return Passed();
        } 
    }

    return Passed();
    
}
