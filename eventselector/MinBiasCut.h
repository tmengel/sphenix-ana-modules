/*!
 * \file MinBiasCut.h
 * \brief Cuts events based on minimum bias information
 * \author Tanner Mengel <tmengel@bnl.gov>
 * \version $Verison: 1.0.1 $
 * \date $Date: 09/10/2024 $
 */

#ifndef EVENTSELECTION_MINBIASCUT_H
#define EVENTSELECTION_MINBIASCUT_H

#include "EventCut.h"

class PHCompositeNode;

class MinBiasCut  : public EventCut
{
  public:

    MinBiasCut();
    ~MinBiasCut() override {} 

    void identify(std::ostream &os = std::cout) const override;
    bool operator()(PHCompositeNode* topNode) override;

};

#endif // EVENTSELECTION_MINBIASCUT_H