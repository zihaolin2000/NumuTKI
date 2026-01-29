// File: ElasticFSIReweighter.h
// Brief: Boosted-Decision-Tree reweighter that reweights MINERvA medium energy
// CCQE-like numu-Carbon events
// from source MC
//    GENIE v2.12.6
// to target MC
//    GENIE 3.0.4 AR23.
// Author: Zihao Lin zlin22@ur.rochester.edu

#ifndef PLOTUTILS_ELASTICFSIREWEIGHTER_H
#define PLOTUTILS_ELASTICFSIREWEIGHTER_H

// //PlotUtils includes
// #include "utilities/NSFDefaults.h"
// #include "universes/GenieSystematics.cxx" //IsNonResPi()
// #include "universes/MnvTuneSystematics.cxx" //IsCCRes()

// Vector class to pass event features
#include <vector>

// Reweighter includes
#include "weighters/Reweighter.h"


namespace PlotUtils
{
    template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
    class ElasticFSIReweighter: public Reweighter<UNIVERSE, EVENT>
    {
        public:
            ElasticFSIReweighter(): Reweighter<UNIVERSE, EVENT> ()
            {
            }

            virtual ~ElasticFSIReweighter() = default;

            double GetWeight(const UNIVERSE& univ, const EVENT& /*event*/) const override
            {
                // only check for CCQE-like 1p0n final states
                if((univ.GetCCQELikeCategory(50,10) == 2) && univ.GetIsElasticFSIBugFate())
                    return 0.0;
                else
                    return 1.0;
            }

            std::string GetName() const override { return "Elastic FSI bug reweight"; }

            bool DependsReco() const override { return false; }

    };
}

#endif //PLOTUTILS_ELASTICFSIREWEIGHTER_H