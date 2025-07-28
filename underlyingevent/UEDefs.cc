#include "UEDefs.h"

#include <iostream>
#include <cstdlib>

bool UEDefs::ExistsInUEMap(const Jet::SRC src) {
    return UEDefs::src_to_towernode_suffix.find(src) != UEDefs::src_to_towernode_suffix.end();
}

void UEDefs::CheckUEMap(const Jet::SRC src) {
    if (!UEDefs::ExistsInUEMap(src)) {
        std::cerr << UEDefs::ERROR_MSG << std::endl;
        exit(-1);
    }
}

std::string UEDefs::GetCaloGeomNode(const Jet::SRC src) {
    UEDefs::CheckUEMap(src);
    return "TOWERGEOM_" + UEDefs::src_to_geomnode_suffix.at(src);
}

std::string UEDefs::GetCaloTowerNode(const Jet::SRC src, const std::string prefix) {
    UEDefs::CheckUEMap(src);
    return prefix + "_" + UEDefs::src_to_towernode_suffix.at(src);
}



