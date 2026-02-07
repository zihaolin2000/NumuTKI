// File: CCQEnuBDTReweighter.h
// Brief: Boosted-Decision-Tree reweighter that reweights MINERvA medium energy
// CCQE-like numu-Carbon events
// from source MC
//    GENIE v2.12.6
// to target MC
//    GENIE 3.0.4 AR23.
// Author: Zihao Lin zlin22@ur.rochester.edu

#ifndef PLOTUTILS_CCQENUBDTREWEIGHTER_H
#define PLOTUTILS_CCQENUBDTREWEIGHTER_H

// //PlotUtils includes
// #include "utilities/NSFDefaults.h"
// #include "universes/GenieSystematics.cxx" //IsNonResPi()
// #include "universes/MnvTuneSystematics.cxx" //IsCCRes()

// Vector class to pass event features
#include <vector>

// Reweighter includes
#include "weighters/Reweighter.h"

// Python includes
#include <pybind11/embed.h>
#include <pybind11/stl.h>
// *** Also Need to set up embed in environment's python ***

// Python utilities namespace
namespace py = pybind11;

namespace PlotUtils
{
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class CCQELikeBDTReweighter: public Reweighter<UNIVERSE, EVENT>
  {
    public:
      CCQELikeBDTReweighter(): Reweighter<UNIVERSE, EVENT>(),
      fSumWeights{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
      fNEvents{0, 0, 0, 0, 0, 0, 0}
      {
        py::initialize_interpreter();
        // sys insert to find the api
        py::module sys = py::module::import("sys");
        sys.attr("path").attr("insert")(0, "/exp/minerva/app/users/zihaolin/MINERvA101_2025/MINERvA-101-Cross-Section/reweight");
        py::module api = py::module::import("BDTReweight_api");
        predict_0p0n = api.attr("predict_weight_0p0n");
        predict_0pNn = api.attr("predict_weight_0pNn");
        predict_1p0n = api.attr("predict_weight_1p0n");
        predict_1pNn = api.attr("predict_weight_1pNn");
        predict_2p0n = api.attr("predict_weight_2p0n");
        predict_2pNn = api.attr("predict_weight_2pNn");
        predict_others = api.attr("predict_weight_others");
      }

      ~CCQELikeBDTReweighter()
      {
        py::finalize_interpreter();
      }

      double GetWeight(const UNIVERSE& univ, const EVENT& /*event*/) const override
      {
        py::gil_scoped_acquire gil;
        int category = univ.GetCCQELikeCategory(50, 10); // proton, neutron detecting KE threshold: 50 MeV, 10 MeV
        std::vector<double> features = univ.GetReactionFrameReweightFeatures(category);
        double weight(0.0);
        if (category == 0)
        {
          weight = predict_0p0n(features).cast<double>();
          fSumWeights[0] += weight;
          fNEvents[0] += 1;
        }
        else if (category == 1)
        {
          weight = predict_0pNn(features).cast<double>();
          fSumWeights[1] += weight;
          fNEvents[1] += 1;
        }
        else if (category == 2)
        {
          weight = predict_1p0n(features).cast<double>();
          fSumWeights[2] += weight;
          fNEvents[2] += 1;
        }
        else if (category == 3)
        {
          weight = predict_1pNn(features).cast<double>();
          fSumWeights[3] += weight;
          fNEvents[3] += 1;
        }
        else if (category == 4)
        {
          weight = predict_2p0n(features).cast<double>();
          fSumWeights[4] += weight;
          fNEvents[4] += 1;
        }
        else if (category == 5)
        {
          weight = predict_2pNn(features).cast<double>();
          fSumWeights[5] += weight;
          fNEvents[5] += 1;
        }
        else if (category == 6)
        {
          weight = predict_others(features).cast<double>();
          fSumWeights[6] += weight;
          fNEvents[6] += 1;
        }
        else
        {
          weight = 1.0;
        }

        if(weight > 100) return 0.0; // Drop large weight events
        else return weight;
      }

      std::string GetName() const override { return "CCQELikeBDTGENIE"; }

      bool DependsReco() const override { return false; }

      std::vector<double> GetTotalWeights() const
      {
        return fSumWeights;
      }

      std::vector<int> GetnCCQELikeEvents() const
      {
        return fNEvents;
      }
      
    private:
      py::object predict_0p0n;
      py::object predict_0pNn;
      py::object predict_1p0n;
      py::object predict_1pNn;
      py::object predict_2p0n;
      py::object predict_2pNn;
      py::object predict_others;
      mutable std::vector<int> fNEvents;
      mutable std::vector<double> fSumWeights;
  };
}

#endif //PLOTUTILS_CCQENUBDTREWEIGHTER_H