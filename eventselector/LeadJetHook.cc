#include "LeadJetHook.h"

// phool includes
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

// jet reco
#include <jetbase/Jetv1.h>
#include <jetbase/JetContainerv1.h>

#include <cstdlib>

void LeadJetHook::identify(std::ostream &os) const
{
    os << Name() + "::identify: " << std::endl;
    os << "  Node name: " << GetNodeName() << std::endl;
    os << "  Min pT: " << GetMin() << std::endl;
    return;
}

bool LeadJetHook::operator()(PHCompositeNode *topNode)
{

    // get truth jet nodes
    auto * jets = findNode::getClass<JetContainer>(topNode, GetNodeName());
    if( !jets ) {
        std::cerr << Name() + "::operator(PHCompositeNode *topNode) Could not find " + GetNodeName() << std::endl;
        exit(-1); // this is a fatal error
    }

    // get leading truth jet pT
    bool found = false;
    for ( auto jet : *jets ) {
        if ( jet->get_et() > GetMin() ) {
            found = true;
            break;
        }
    }

    Passed(found);

    if(!Passed() && Verbosity()){
       std::cout << Name() + "::operator(PHCompositeNode *topNode) Event failed" + Name() + ". Lead jet pT in " + GetNodeName() + " > " + std::to_string(GetMin()) << std::endl;
    }

    return Passed();
}
