/*!
 * \file MissingSebFilter.h
 * \brief Filters events based on the presence of a all SEBs 
 * \author Tanner Mengel <tmengel@bnl.gov>
 * \version $Verison: 1.0.1 $
 * \date $Date: 09/10/2024 $
 */

#ifndef EVENTSELECTION_MISSINGSEBFILETER_H
#define EVENTSELECTION_MISSINGSEBFILETER_H

#include "EventCut.h"

class PHCompositeNode;

class MissingSebFilter  : public EventCut
{
  public:

    MissingSebFilter();
    ~MissingSebFilter() override {}

    void identify(std::ostream &os = std::cout) const override;
    bool operator()(PHCompositeNode* topNode) override;

};

#endif // EVENTSELECTION_MISSINGSEBFILETER_H