#ifndef UNDERLYINGEVENT_UEDEFS_H
#define UNDERLYINGEVENT_UEDEFS_H

#include <jetbase/Jet.h>

#include <map>
#include <string>

namespace UEDefs {

    const std::map<Jet::SRC, std::string> src_to_towernode_suffix = {
        {Jet::SRC::CEMC_TOWERINFO, "CEMC"},
        {Jet::SRC::HCALIN_TOWERINFO, "HCALIN"},
        {Jet::SRC::HCALOUT_TOWERINFO, "HCALOUT"},
        {Jet::SRC::CEMC_TOWERINFO_RETOWER, "CEMC_RETOWER"},
        {Jet::SRC::CEMC_TOWERINFO_SUB1, "CEMC_RETOWER_SUB1"},
        {Jet::SRC::HCALIN_TOWERINFO_SUB1, "HCALIN_SUB1"},
        {Jet::SRC::HCALOUT_TOWERINFO_SUB1, "HCALOUT_SUB1"}, 
    };

    const std::map<Jet::SRC, std::string> src_to_geomnode_suffix = {
        {Jet::SRC::CEMC_TOWERINFO, "CEMC"},
        {Jet::SRC::HCALIN_TOWERINFO, "HCALIN"},
        {Jet::SRC::HCALOUT_TOWERINFO, "HCALOUT"},
        {Jet::SRC::CEMC_TOWERINFO_RETOWER, "HCALIN"},
        {Jet::SRC::CEMC_TOWERINFO_SUB1, "HCALIN"},
        {Jet::SRC::HCALIN_TOWERINFO_SUB1, "HCALIN"},
        {Jet::SRC::HCALOUT_TOWERINFO_SUB1, "HCALOUT"}, 
    };

    const std::string ERROR_MSG = "ERROR: SRC not found in UEDefs map. Supply a known SRC. or include the geom_node and Jet::Src";
    bool ExistsInUEMap(const Jet::SRC src);
    void CheckUEMap(const Jet::SRC src);
    std::string GetCaloGeomNode(const Jet::SRC src);
    std::string GetCaloTowerNode(const Jet::SRC src, const std::string prefix = "TOWERINFO_CALIB");

}

#endif // UNDERLYINGEVENT_UEDEFS_H