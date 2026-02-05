
// File: NumuTKISignal.h
// This file is forked from Carlos Pernas's NuETKI/cuts/NuETKISignal.h, forked
// from Andrew's CCInclusiveSignal.h.
// Some Nue quantities are switched to Numu.
// -- Zihao Lin zlin22@ur.rochester.edu


//File: CCInclusiveSignal.h
//Brief: Signal definition constraints for a CCInclusive analysis.
//       Everything you need to reproduce Dan's 2D inclusive signal
//       definition.
//Author: Andrew Olivier aolivier@ur.rochester.edu

namespace truth
{
  template <class UNIVERSE>
  class IsNeutrino: public PlotUtils::SignalConstraint<UNIVERSE>
  {
    public:
      IsNeutrino(): PlotUtils::SignalConstraint<UNIVERSE>("IsNeutrino")
      {
      }

    private:
      bool checkConstraint(const UNIVERSE& univ) const override
      {
        return univ.GetTruthNuPDG() == 14;
      }
  };


  template <class UNIVERSE>
  class IsAntiNeutrino: public PlotUtils::SignalConstraint<UNIVERSE>
  {
    public:
      IsAntiNeutrino(): PlotUtils::SignalConstraint<UNIVERSE>("IsAntiNeutrino")
      {
      }

    private:
      bool checkConstraint(const UNIVERSE& univ) const override
      {
        return univ.GetTruthNuPDG() == -14;
      }
  };

  // Checks muon flavor, will accept both nu and anti nu
  template <class UNIVERSE>
  class IsNumu: public PlotUtils::SignalConstraint<UNIVERSE>
  {
    public:
      IsNumu(): PlotUtils::SignalConstraint<UNIVERSE>("IsNumu or AntiNumu")
      {
      }

    private:
      bool checkConstraint(const UNIVERSE& univ) const override
      {
        return abs(univ.GetTruthNuPDG()) == 14;
      }
  };



  template <class UNIVERSE>
  class IsCC: public PlotUtils::SignalConstraint<UNIVERSE>
  {
    public:
      IsCC(): PlotUtils::SignalConstraint<UNIVERSE>("IsCC")
      {
      }

    private:
      bool checkConstraint(const UNIVERSE& univ) const override
      {
        return univ.GetCurrent() == 1;
      }
  };

// TODO: Eavail and TrueElectronEnergy switch to numu, muon quantities
/*
  //electron neutrino truth definitions here (Carlos P, Nov 2023)
  template <class UNIVERSE>
    class Eavail: public PlotUtils::SignalConstraint<UNIVERSE>
    {
    public:
    Eavail(const double max): PlotUtils::SignalConstraint<UNIVERSE>(std::string("Eavail < ") + std::to_string(max)), fMax(max)
      {
      }

    private:
      bool checkConstraint(const UNIVERSE& univ) const override
      {
        return univ.GetEavailTrue() <= fMax;
      }

      const double fMax;
    };

  template <class UNIVERSE>
    class TrueElectronEnergy: public PlotUtils::SignalConstraint<UNIVERSE>
    {
    public:
    TrueElectronEnergy(const double min): PlotUtils::SignalConstraint<UNIVERSE>(std::string("True electron E > ") + std::to_string(min)), fMin(min)
      {
      }

    private:
      bool checkConstraint(const UNIVERSE& univ) const override
      {
        return univ.GetElectronEnergyTrue() >= fMin;
      }

      const double fMin;
    };
*/
  template <class UNIVERSE>
  class Is2p2hExtendQ0Cut: public PlotUtils::SignalConstraint<UNIVERSE>
  {
    public:
      Is2p2hExtendQ0Cut(): PlotUtils::SignalConstraint<UNIVERSE>("Is 2p2h q3 < 1200 MeV")
      {
      }

    private:
      bool checkConstraint(const UNIVERSE& univ) const override
      {
        return univ.GetIsq3Below1200MeV();
      }
  };

  template <class UNIVERSE>
  class HasAbove50MeVProton: public PlotUtils::SignalConstraint<UNIVERSE>
  {
    public:
      HasAbove50MeVProton(): PlotUtils::SignalConstraint<UNIVERSE>("Has FS Proton with Tp >= 50 MeV")
      {
      }

    private:
      bool checkConstraint(const UNIVERSE& univ) const override
      {
        return univ.GetHasAbove50MeVProton();
      }
  };

  template <class UNIVERSE>
  class Is1p0nTopology: public PlotUtils::SignalConstraint<UNIVERSE>
  {
    public:
      Is1p0nTopology(): PlotUtils::SignalConstraint<UNIVERSE>("Is 1p0n Topology")
      {
      }

    private:
      bool checkConstraint(const UNIVERSE& univ) const override
      {
        return (univ.GetCCQELikeCategory(50,10) == 2);
      }
  };

  template <class UNIVERSE>
  class Is1pNnTopology: public PlotUtils::SignalConstraint<UNIVERSE>
  {
    public:
      Is1pNnTopology(): PlotUtils::SignalConstraint<UNIVERSE>("Is 1pNn Topology")
      {
      }

    private:
      bool checkConstraint(const UNIVERSE& univ) const override
      {
        return (univ.GetCCQELikeCategory(50,10) == 3);
      }
  };

  template <class UNIVERSE>
  class Is2p0nTopology: public PlotUtils::SignalConstraint<UNIVERSE>
  {
    public:
      Is2p0nTopology(): PlotUtils::SignalConstraint<UNIVERSE>("Is 2p0n Topology")
      {
      }

    private:
      bool checkConstraint(const UNIVERSE& univ) const override
      {
        return (univ.GetCCQELikeCategory(50,10) == 4);
      }
  };

  template <class UNIVERSE>
  class Is2pNnTopology: public PlotUtils::SignalConstraint<UNIVERSE>
  {
    public:
      Is2pNnTopology(): PlotUtils::SignalConstraint<UNIVERSE>("Is 2pNn Topology")
      {
      }

    private:
      bool checkConstraint(const UNIVERSE& univ) const override
      {
        return (univ.GetCCQELikeCategory(50,10) == 5);
      }
  };

    template <class UNIVERSE>
  class Is3pOthersTopology: public PlotUtils::SignalConstraint<UNIVERSE>
  {
    public:
      Is3pOthersTopology(): PlotUtils::SignalConstraint<UNIVERSE>("Is >=3p Topology")
      {
      }

    private:
      bool checkConstraint(const UNIVERSE& univ) const override
      {
        return (univ.GetCCQELikeCategory(50,10) == 6);
      }
  };

  template <class UNIVERSE>
  class IsTargetCarbon: public PlotUtils::SignalConstraint<UNIVERSE>
  {
    public:
      IsTargetCarbon(): PlotUtils::SignalConstraint<UNIVERSE>("Is Target Carbon")
      {
      }

    private:
      bool checkConstraint(const UNIVERSE& univ) const override
      {
        return univ.GetIsTargetCarbon();
      }
  };


  // Copying classes from Carlos's NuETKI/cuts/NuETKISignal.h.
  // -- Ziggy
  template <class UNIVERSE>
  class HasSignalProton: public PlotUtils::SignalConstraint<UNIVERSE>
  {
    public:
      HasSignalProton(): PlotUtils::SignalConstraint<UNIVERSE>("HasSignalProton (450 < P_p < 1200 MeV/c, (Theta_p < 70 deg OR > 110 deg))")
        {
        }

    private:
      bool checkConstraint(const UNIVERSE& univ) const override
      {
        return (univ.GetHasSignalFSProton()==1);
      }
  };

  template <class UNIVERSE>
  class HasNoMeson: public PlotUtils::SignalConstraint<UNIVERSE>
  {
    public:
      HasNoMeson(): PlotUtils::SignalConstraint<UNIVERSE>("HasNoMesons")
        {
        }

    private:
      bool checkConstraint(const UNIVERSE& univ) const override
      {
        // return (univ.GetHasFSMeson()==0);
        return univ.GetHasFSMeson() == false;
      }
  };

  template <class UNIVERSE>
  class HasNoPhoton: public PlotUtils::SignalConstraint<UNIVERSE>
  {
    public:
      HasNoPhoton(): PlotUtils::SignalConstraint<UNIVERSE>("HasNoPhoton")
        {
        }

    private:
      bool checkConstraint(const UNIVERSE& univ) const override
      {
        // return (univ.GetHasFSPhoton()==0);
        return univ.GetHasFSPhoton() == false;
      }
  };

  template <class UNIVERSE>
  class MuonAngle: public PlotUtils::SignalConstraint<UNIVERSE>
    {
      public:
        MuonAngle(const double angleMax): PlotUtils::SignalConstraint<UNIVERSE>(std::string("Muon Angle ") + std::to_string(angleMax)), fMax(angleMax*M_PI/180.)
        {
        }

      private:
        bool checkConstraint(const UNIVERSE& univ) const override
        {
          return univ.GetThetalepTrue() <= fMax;
        }

        const double fMax;
    };

  template <class UNIVERSE>
  class PZMuMin: public PlotUtils::SignalConstraint<UNIVERSE>
  {
    public:
      PZMuMin(const double min): PlotUtils::SignalConstraint<UNIVERSE>(std::string("PzMu > ") + std::to_string(min)), fMin(min)
      {
      }

    private:
      bool checkConstraint(const UNIVERSE& univ) const override
      {
        return univ.GetPlepTrue() * cos(univ.GetThetalepTrue()) >= fMin;
      }

      const double fMin;
  };

  template <class UNIVERSE>
  class Apothem: public PlotUtils::SignalConstraint<UNIVERSE>
  {
    public:
      Apothem(const double apothem): PlotUtils::SignalConstraint<UNIVERSE>(std::string("Apothem ") + std::to_string(apothem)), fApothem(apothem)
      {
      }

    private:
      bool checkConstraint(const UNIVERSE& univ) const override
      {
        const auto vertex = univ.GetTrueVertex();
        return (fabs(vertex.y()) < fSlope*fabs(vertex.x()) + 2.*fApothem/sqrt(3.))
               && (fabs(vertex.x()) < fApothem);
      }

      const double fApothem;
      const double fSlope = -1./sqrt(3.); //A regular hexagon has angles of 2*M_PI/3, so I can find this is 1/tan(M_PI/3.)
  };

  template <class UNIVERSE>
  class ZRange: public PlotUtils::SignalConstraint<UNIVERSE>
  {
    public:
      ZRange(const std::string& name, const double zMin, const double zMax): PlotUtils::SignalConstraint<UNIVERSE>(name), fMin(zMin), fMax(zMax)
      {
      }

    private:
      bool checkConstraint(const UNIVERSE& univ) const override
      {
        return univ.GetTrueVertex().z() >= fMin && univ.GetTrueVertex().z() <= fMax;
      }

      const double fMin;
      const double fMax;
  };

  template <class UNIVERSE>
  PlotUtils::constraints_t<UNIVERSE> GetCCInclusive2DPhaseSpace()
  {
    PlotUtils::constraints_t<UNIVERSE> signalDef;

    signalDef.emplace_back(new ZRange<UNIVERSE>("Tracker", 5980, 8422));
    signalDef.emplace_back(new Apothem<UNIVERSE>(850.));
    signalDef.emplace_back(new MuonAngle<UNIVERSE>(20.));
    signalDef.emplace_back(new PZMuMin<UNIVERSE>(1500.));

    return signalDef;
  }

  template <class UNIVERSE>
  PlotUtils::constraints_t<UNIVERSE> GetCCInclusive2DSignal()
  {
    PlotUtils::constraints_t<UNIVERSE> signalDef;

    signalDef.emplace_back(new IsNeutrino<UNIVERSE>());
    signalDef.emplace_back(new IsCC<UNIVERSE>());

    return signalDef;
  }
}