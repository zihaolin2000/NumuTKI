// =============================================================================
// Base class for an un-systematically shifted (i.e. CV) universe. Implement
// "Get" functions for all the quantities that you need for your analysis.
//
// This class inherits from PU::MinervaUniverse, which in turn inherits from
// PU::BaseUniverse. PU::BU defines the interface with anatuples.
// 
// Within the class, "WeightFunctions" and "MuonFunctions" are included to gain
// access to standardized weight and muon variable getters. See:
// https://cdcvs.fnal.gov/redmine/projects/minerva-sw/wiki/MinervaUniverse_Structure_
// for a full list of standardized functions you can use. In general, if a
// standard version of a function is available, you should be using it.
// =============================================================================
#ifndef CVUNIVERSE_H
#define CVUNIVERSE_H

#include <iostream>
#include "PlotUtils/MinervaUniverse.h"

// ROOT includes
#include "Math/RotationX.h"
#include "Math/Vector3D.h"
#include "Math/Vector2D.h"
#include "TMath.h"

class CVUniverse : public PlotUtils::MinervaUniverse
{

  public:
  #include "PlotUtils/MuonFunctions.h" // GetMinosEfficiencyWeight
  #include "PlotUtils/TruthFunctions.h" //Getq3True
  // ========================================================================
  // Constructor/Destructor
  // ========================================================================
  CVUniverse(PlotUtils::ChainWrapper* chw, double nsigma = 0)
      : PlotUtils::MinervaUniverse(chw, nsigma) {}

  virtual ~CVUniverse() {}

  // ========================================================================
  // Quantities defined here as constants for the sake of below. Definition
  // matched to Dan's CCQENuInclusiveME variables from:
  // `/minerva/app/users/drut1186/cmtuser/Minerva_v22r1p1_OrigCCQENuInc/Ana/CCQENu/ana_common/include/CCQENuUtils.h`
  // ========================================================================
  static constexpr double M_n = 939.56536;
  static constexpr double M_p = 938.272013;
  static constexpr double M_nucleon = (1.5*M_n+M_p)/2.5;

  static constexpr int PDG_n = 2112;
  static constexpr int PDG_p = 2212;
  static constexpr double pi = 3.141592653589793;

  // ========================================================================
  // Write a "Get" function for all quantities access by your analysis.
  // For composite quantities (e.g. Enu) use a calculator function.
  //
  // In order to properly calculate muon variables and systematics use the
  // various functions defined in MinervaUniverse.
  // E.g. GetPmu, GetEmu, etc.
  // ========================================================================

  // Quantities only needed for cuts
  // Although unlikely, in principle these quanties could be shifted by a
  // systematic. And when they are, they'll only be shifted correctly if we
  // write these accessor functions.


  //============================================================================
  // CCQEnu selection custom CV universe Get functions
  // -- Ziggy
  // to get phys_n_dead_discr_pair_upstream_prim_track_proj: use GetTDead()
  // to get muon_theta: use GetThetamu()
  //============================================================================

  int GetNEventExtraTrackPID() const
  {
    return GetInt("event_extra_track_PID_sz");
  }
 
  std::vector<double> GetEventExtraTrackPID() const //FIXME: return type
  {
    return GetVecDouble("event_extra_track_PID");
  }

  int GetHasInteractionVertex() const
  {
    return GetInt("has_interaction_vertex");
  }

  double GetProtonEndX() const
  {
    return GetDouble("proton_track_endx");
  }
  
  double GetProtonEndY() const
  {
    return GetDouble("proton_track_endy");
  }

  virtual int GetNImprovedMichel() const
  {
    return GetInt("improved_michel_vertex_type_sz");
  }

  std::vector<double> GetBlobsStartZ() const
  {
    return GetVecDouble("nonvtx_iso_blobs_start_position_z_in_prong");
  }
  
  int GetNBlobsStartZ() const
  {
    return GetInt("nonvtx_iso_blobs_start_position_z_in_prong_sz");
  }

  std::vector<double> GetProtonNodesNormE() const
  {
    return GetVecDouble((GetAnaToolName()+"_proton_nodes_nodesNormE").c_str());
  }

  int GetNProtonNodesNormE() const
  {
    return GetInt((GetAnaToolName()+"_proton_nodes_nodesNormE_sz").c_str());
  }

  int GetNuHelicity() const
  {
    return GetInt((GetAnaToolName()+"_nuHelicity").c_str());
  }

  double GetCCQEnuQ2() const
  {
    return GetDouble((GetAnaToolName()+"_Q2_CCQE").c_str());
  }

  double GetMuonTheta() const // reco
  {
    return GetThetamu();
  }

  double GetMultiplicity() const
  {
    return GetInt("multiplicity");
  }

  double GetCCQELikeTotalTpReco() const // reco, GeV
  {
    return GetDouble("recoil_energy_nonmuon_nonvtx0mm")/1000;
  }

  // Get CCQE-like category based on nubmer of above-threshold final state
  // protons and neutrons.  
  // KE thresholds in MeV.
  int GetCCQELikeCategory(double proton_threshold, double neutron_threshold) const
  {
    int n_protons(0);
    int n_neutrons(0);
    std::vector<int> FSParticles = GetVecInt("mc_FSPartPDG");
    std::vector<double> energies = GetVecDouble("mc_FSPartE");
    for (size_t i = 0; i < FSParticles.size(); i++)
    {
      if (FSParticles[i] == 2212 && (energies[i] - M_p) >= proton_threshold)
        n_protons++;
      else if (FSParticles[i] == 2112 && (energies[i] - M_n) >= neutron_threshold)
        n_neutrons++;
    }
    if (n_protons == 0 && n_neutrons == 0)
      return 0; //0p0n
    else if (n_protons == 0 && n_neutrons > 0)
      return 1; //0pNn
    else if (n_protons == 1 && n_neutrons == 0)
      return 2; //1p0n
    else if (n_protons == 1 && n_neutrons > 0)
      return 3; //1pNn
    else if (n_protons == 2 && n_neutrons == 0)
      return 4; //2p0n
    else if (n_protons == 2 && n_neutrons > 0)
      return 5; //2pNn
    else
      return 6; //others (>=3p)
  }

  // Get truth muon index
  int GetMuonIndex() const
  {
    int index(-999);
    std::vector<int> FSParticles = GetVecInt("mc_FSPartPDG");
    for (size_t i = 0; i < FSParticles.size(); i++)
    {
      if (FSParticles[i] == 13)
      {
          index = i;
          break;
      }
    }
    return index;
  }

  int GetLeadingProtonIndex() const
  {
    double leadingE(-999.0);
    int leading_i(-999);
    std::vector<int> FSParticles = GetVecInt("mc_FSPartPDG");
    std::vector<double> energies = GetVecDouble("mc_FSPartE");
    for (size_t i = 0; i < FSParticles.size(); i++)
    {
      if (FSParticles[i] == 2212 && energies[i] > leadingE)
      {
        leadingE = energies[i];
        leading_i = i;
      }
    }
    return leading_i;
  }

  int GetLeadingNeutronIndex() const
  {
    double leadingE(-999.0);
    int leading_i(-999);
    std::vector<int> FSParticles = GetVecInt("mc_FSPartPDG");
    std::vector<double> energies = GetVecDouble("mc_FSPartE");
    for (size_t i = 0; i < FSParticles.size(); i++)
    {
      if (FSParticles[i] == 2112 && energies[i] > leadingE)
      {
        leadingE = energies[i];
        leading_i = i;
      }
    }
    return leading_i;
  }

  std::vector<int> Get2HighestKEProtonIndices() const
  {
    double leadingE(-999.0);
    double subleadingE(-9999.0);
    int leading_i(-999);
    int subleading_i(-999);
    std::vector<int> FSParticles = GetVecInt("mc_FSPartPDG");
    std::vector<double> energies = GetVecDouble("mc_FSPartE");
    for (size_t i = 0; i < FSParticles.size(); i++)
    {
      if (FSParticles[i] == 2212 && energies[i] > leadingE)
      {
        subleading_i = leading_i;
        subleadingE = leadingE;
        leading_i = i;
        leadingE = energies[i];
      }
      else if (FSParticles[i] == 2212 && energies[i] < leadingE && energies[i] > subleadingE)
      {
        subleading_i = i;
        subleadingE = energies[i];
      }
    }
    return std::vector<int>{leading_i, subleading_i};
  }

  // Helper function to convert from lab frame to beam frame.
  // Beam points into earth by 3.3 degree, so momentum P in beam frame is given by:
  // rotating lab frame P by -3.3 degree around positive x-axis.
  // Neutrino direction is [0, 0, 1] in beam frame.
  ROOT::Math::XYZVector ConvertToBeamFrame(ROOT::Math::XYZVector labVec) const
  {
    ROOT::Math::RotationX r(-3.3 * (pi / 180.));
    return r(labVec);
  }

  // Get **BEAM FRAME** truth momentum 3-vector in GeV.
  ROOT::Math::XYZVector GetParticlePVec(int index) const // truth
  {
    if (index == -999)
      return ROOT::Math::XYZVector(-999.0, -999.0, -999.0);
    else
    {  
      ROOT::Math::XYZVector labP(GetVecElem("mc_FSPartPx", index)/1000, GetVecElem("mc_FSPartPy", index)/1000, GetVecElem("mc_FSPartPz", index)/1000);
      return ConvertToBeamFrame(labP);
    }
  }

  // Convert beam frame particle p vec to reaction frame given beam frame muon p vec
  ROOT::Math::XYZVector ConvertToReactionFrame(const ROOT::Math::XYZVector beam_Pp, const ROOT::Math::XYZVector beam_Pmu) const
  {
    // treat muon transverse direction as transverse plane y vector
    ROOT::Math::XYVector plane_Pp(beam_Pp.X(), beam_Pp.Y());
    ROOT::Math::XYVector plane_y(- beam_Pmu.X(), - beam_Pmu.Y());
    // Normalize y vector
    plane_y = plane_y.Unit();
    // x vector is simply y vector rotated clockwise by 90
    ROOT::Math::XYVector plane_x(- plane_y.Y(), plane_y.X());
    double new_px = plane_Pp.Dot(plane_x);
    double new_py = plane_Pp.Dot(plane_y);
    return ROOT::Math::XYZVector(new_px, new_py, beam_Pp.Z());
  }

  // Get *BEAM FRAME* momentum 3-vector in GeV.
  ROOT::Math::XYZVector GetTotalProtonPvec() const // truth
  {
    ROOT::Math::XYZVector totalP(0.0, 0.0, 0.0);
    std::vector<int> FSParticles = GetVecInt("mc_FSPartPDG");
    for (size_t i = 0; i < FSParticles.size(); i++)
    {
      if (FSParticles[i] == 2212)
      {
        ROOT::Math::XYZVector protonP = GetParticlePVec(i);
        totalP += protonP;
      }
    }
    return totalP;
  }

  // Get *BEAM FRAME* momentum 3-vector in GeV.
  ROOT::Math::XYZVector GetTotalNeutronPvec() const // truth
  {
    ROOT::Math::XYZVector totalP(0.0, 0.0, 0.0);
    std::vector<int> FSParticles = GetVecInt("mc_FSPartPDG");
    for (size_t i = 0; i < FSParticles.size(); i++)
    {
      if (FSParticles[i] == 2112)
      {
        ROOT::Math::XYZVector neutronP = GetParticlePVec(i);
        totalP += neutronP;
      }
    }
    return totalP;
  
  }

  double GetTotalProtonTp() const // truth, GeV
  {
    double totalTp(0.0);
    std::vector<double> energies = GetVecDouble("mc_FSPartE");
    std::vector<int> FSParticles = GetVecInt("mc_FSPartPDG");
    for (size_t i = 0; i < FSParticles.size(); i++)
    {
      if (FSParticles[i] == 2212)
        totalTp += (energies[i] - M_p)/1000;
    }
    return totalTp;
  }

  double GetTotalNeutronTn() const // truth, GeV
  {
    double totalTn(0.0);
    std::vector<double> energies = GetVecDouble("mc_FSPartE");
    std::vector<int> FSParticles = GetVecInt("mc_FSPartPDG");
    for (size_t i = 0; i < FSParticles.size(); i++)
    {
      if (FSParticles[i] == 2112)
        totalTn += (energies[i] - M_n)/1000;
    }
    return totalTn;
  }

  // Get features for CCQELikeBDTReweighter to predict weight of an event.
  // Return a vector of values of reweight variables.
  std::vector<double> GetReactionFrameReweightFeatures(const int category) const
  {
    std::vector<double> features = {-999, -999, -999};
    ROOT::Math::XYZVector beamPmu = GetParticlePVec(GetMuonIndex());
    double totalTp = GetTotalProtonTp(), muonPy = - GetMuonPTTrue()/1000, muonPz = GetMuonPzTrue()/1000; // Use truth muon PT Pz
    if (category == 0) // 0p0n
    {
      ROOT::Math::XYZVector totalPp = ConvertToReactionFrame(GetTotalProtonPvec(), beamPmu);
      features = {totalPp.X(), totalPp.Y(), totalPp.Z(), totalTp, muonPy, muonPz};
    }
    else if (category == 1) // 0pNn
    {
      ROOT::Math::XYZVector leadingPn = ConvertToReactionFrame(GetParticlePVec(GetLeadingNeutronIndex()), beamPmu);
      ROOT::Math::XYZVector totalPp = ConvertToReactionFrame(GetTotalProtonPvec(), beamPmu);
      features =
      {
        leadingPn.X(), leadingPn.Y(), leadingPn.Z(),
        totalPp.X(), totalPp.Y(), totalPp.Z(),
        totalTp, muonPy, muonPz
      };
    }
    else if (category == 2  || category == 6) // 1p0n or others
    {
      ROOT::Math::XYZVector leadingPp = ConvertToReactionFrame(GetParticlePVec(GetLeadingProtonIndex()), beamPmu);
      features = {leadingPp.X(), leadingPp.Y(), leadingPp.Z(), totalTp, muonPy, muonPz};
    }
    else if (category == 3) // 1pNn
    {
      ROOT::Math::XYZVector leadingPp = ConvertToReactionFrame(GetParticlePVec(GetLeadingProtonIndex()), beamPmu);
      ROOT::Math::XYZVector leadingPn = ConvertToReactionFrame(GetParticlePVec(GetLeadingNeutronIndex()), beamPmu);
      features =
      {
        leadingPp.X(), leadingPp.Y(), leadingPp.Z(),
        totalTp, muonPy, muonPz,
        leadingPn.X(), leadingPn.Y(), leadingPn.Z()
      };
    }
    else if (category == 4) // 2p0n
    {
      std::vector<int> ids = Get2HighestKEProtonIndices();
      int leading_i(ids[0]);
      int subleading_i(ids[1]);
      ROOT::Math::XYZVector leadingPp = ConvertToReactionFrame(GetParticlePVec(leading_i), beamPmu);
      ROOT::Math::XYZVector subLeadingPp = ConvertToReactionFrame(GetParticlePVec(subleading_i), beamPmu);
      features =
      {
        leadingPp.X(), leadingPp.Y(), leadingPp.Z(),
        totalTp, muonPy, muonPz,
        subLeadingPp.X(), subLeadingPp.Y(), subLeadingPp.Z()
      };
    }
    else if (category == 5) // 2pNn
    {
      std::vector<int> ids = Get2HighestKEProtonIndices();
      int leading_i(ids[0]);
      int subleading_i(ids[1]);
      ROOT::Math::XYZVector leadingPp = ConvertToReactionFrame(GetParticlePVec(leading_i), beamPmu);
      ROOT::Math::XYZVector subLeadingPp = ConvertToReactionFrame(GetParticlePVec(subleading_i), beamPmu);
      ROOT::Math::XYZVector leadingPn = ConvertToReactionFrame(GetParticlePVec(GetLeadingNeutronIndex()), beamPmu);
      features =
      {
        leadingPp.X(), leadingPp.Y(), leadingPp.Z(),
        totalTp, muonPy, muonPz,
        leadingPn.X(), leadingPn.Y(), leadingPn.Z(),
        subLeadingPp.X(), subLeadingPp.Y(), subLeadingPp.Z()
      };
    }
    return features;
  }

  ROOT::Math::XYZVector GetLeadingProtonReactionFramePvecTrue() const
  {
    ROOT::Math::XYZVector muonP = GetParticlePVec(GetMuonIndex());
    ROOT::Math::XYZVector protonP = GetParticlePVec(GetLeadingProtonIndex());
    return ConvertToReactionFrame(protonP, muonP);
  }

  double GetLeadingProtonReactionFramePxTrue() const
  {
    return GetLeadingProtonReactionFramePvecTrue().X();
  }

  double GetLeadingProtonReactionFramePyTrue() const
  {
    return GetLeadingProtonReactionFramePvecTrue().Y();
  }

  double GetLeadingProtonReactionFramePzTrue() const
  {
    return GetLeadingProtonReactionFramePvecTrue().Z();
  }

  ROOT::Math::XYZVector GetLeadingProtonReactionFramePvecReco() const
  {
    ROOT::Math::XYZVector protonP(GetDouble("MasterAnaDev_proton_Px_fromdEdx")/1000., GetDouble("MasterAnaDev_proton_Py_fromdEdx")/1000., GetDouble("MasterAnaDev_proton_Pz_fromdEdx")/1000.);
    protonP = ConvertToBeamFrame(protonP);
    ROOT::Math::PxPyPzEVector muon4V = GetMuon4V(); // MeV, in beam Frame
    ROOT::Math::XYZVector muonP(muon4V.Px()/1000, muon4V.Py()/1000, muon4V.Pz()/1000); // GeV/c
    return ConvertToReactionFrame(protonP, muonP);
  }

  double GetLeadingProtonReactionFramePxReco() const
  {
    return GetLeadingProtonReactionFramePvecReco().X();
  }

  double GetLeadingProtonReactionFramePyReco() const
  {
    return GetLeadingProtonReactionFramePvecReco().Y();
  }

  double GetLeadingProtonReactionFramePzReco() const
  {
    return GetLeadingProtonReactionFramePvecReco().Z();
  }

  // A sinple check of signal >=50 MeV proton -- Ziggy's CCQE-like selection  
  bool GetHasAbove50MeVProton() const
  {
    std::vector<int> FSParticles = GetVecInt("mc_FSPartPDG");
    std::vector<double> energies = GetVecDouble("mc_FSPartE");
    for (size_t i = 0; i < FSParticles.size(); i++)
    {
      if (FSParticles[i] == 2212 && energies[i] >= (M_p + 50.0))
        return true;
    }
    return false;
  }

  // check if the event has elastic FSI bug fate particle, so to reweight them to 0 -- Ziggy
  bool GetIsElasticFSIBugFate() const
  {
    int mc_incoming = GetInt("mc_incoming");
    int mc_primaryLepton = GetInt("mc_primaryLepton");
    int mc_charm = GetInt("mc_charm");
    int mc_intType = GetInt("mc_intType");
    int mc_targetA = GetInt("mc_targetA");
    int mc_targetZ = GetInt("mc_targetZ");
    int mc_er_nPart = GetInt("mc_er_nPart");
    std::vector<int> mc_er_ID     = GetVecInt("mc_er_ID");
    std::vector<int> mc_er_status = GetVecInt("mc_er_status");
    std::vector<int> mc_er_FD     = GetVecInt("mc_er_FD");
    std::vector<int> mc_er_LD     = GetVecInt("mc_er_LD");
    std::vector<int> mc_er_mother = GetVecInt("mc_er_mother");
    std::vector<double> mc_er_Px = GetVecDouble("mc_er_Px");
    std::vector<double> mc_er_Py = GetVecDouble("mc_er_Py");
    std::vector<double> mc_er_Pz = GetVecDouble("mc_er_Pz");
    std::vector<double> mc_er_E  = GetVecDouble("mc_er_E");

    std::vector<int> fates = calcFates(
      mc_incoming, mc_primaryLepton, mc_charm, mc_intType, mc_targetA,
      mc_targetZ, mc_er_nPart, mc_er_ID, mc_er_status, mc_er_FD,
      mc_er_LD, mc_er_mother, mc_er_Px, mc_er_Py, mc_er_Pz, mc_er_E
    );

    for (size_t i = 0; i < fates.size(); i++)
    if (fates[i] == 3) return true;
    return false;
  }

  // calculate fates of each particle in event record
  // Copied from Andrew's weighters/weight_fsi.cxx
  std::vector<Int_t> calcFates(
    int mc_incoming,
    int mc_primaryLepton,
    int mc_charm,
    int mc_intType,
    int mc_targetA,
    int mc_targetZ,
    int mc_er_nPart,
    const std::vector<int>& mc_er_ID,
    const std::vector<int>& mc_er_status,
    const std::vector<int>& mc_er_FD,
    const std::vector<int>& mc_er_LD,
    const std::vector<int>& mc_er_mother,
    const std::vector<double>& mc_er_Px,
    const std::vector<double>& mc_er_Py,
    const std::vector<double>& mc_er_Pz,
    const std::vector<double>& mc_er_E) const
  {
    // int vector initialized with fate = -1 for all particles. -- Ziggy 2026/1/26
    std::vector<Int_t> tempFates(mc_er_nPart, -1); 
    // the fate codes are what GENIE hA in GENIEv2 used.
    // the fate codes for hN and GENIEv3 have different numerology.
    // beware if you are forward porting this code, or it will hurt.
    for(Int_t i=0; i < mc_er_nPart; ++i)
    {
      if(mc_er_status[i] == 14)
      {
        Int_t tempfate = -1;
        // Only track fates for proton, neutron, pizero, piplus, piminus.
        if (mc_er_ID[i] == 2212 || mc_er_ID[i] == 2112 || mc_er_ID[i] == 111 || TMath::Abs(mc_er_ID[i])==211)
        {
          // The distinguishing feature is first daughter last daughter
          // get an easy handle to these for use later
          Int_t fd = mc_er_FD[i];
          Int_t ld = mc_er_LD[i];

          // is one of these a pion
          Int_t isaPion = 0;
          if(mc_er_ID[i] == 111 || TMath::Abs(mc_er_ID[i] == 211)) isaPion = 1;
          
          // how many daughters are pions
          Int_t FSpions = 0;
          for(int jj = mc_er_FD[i]; jj <= mc_er_LD[i]; ++jj)
          if(TMath::Abs(mc_er_ID[jj]) == 211 || mc_er_ID[jj] == 111) FSpions++;
          
          if(mc_er_FD[i] == mc_er_LD[i])
          {
          // only one daughter.
          //fate is either 1 = no interaction or 3 = elastic
            // or 5 = absorption with a 3xxxxxxxx hadblob 
          // This test of 25.00000 MeV is to separate elastic
          // but the actual number changes with nucleus
          // but most nuclei I know to six sig figs.
          // unknown nuclei get an offset like 8 MeV
          Double_t offset =  getGenieBEinMeV(mc_targetA);  //25.000000;
          Double_t tolerance = 0.0001; //within machine precision
          // unknown nuclei need a larger tolerance
          if(TMath::Abs(offset - 8.0) < 0.1)tolerance = 1.2;
          if(offset < 6.0) tolerance = 1.2;
          if(mc_intType != 1) offset = 0.0;

          if(mc_er_ID[fd] > 2000000000)
            tempfate = 5;   // absorption on 3 nucleons
          else if(TMath::Abs(mc_er_E[i] - mc_er_E[fd] - offset) < tolerance)
            tempfate = 1;   // no scattering
          else if(mc_er_ID[fd] != mc_er_ID[i])
            tempfate = 2;   // charge exchange ?
          else
            tempfate = 3;   // elastic fate
          }
          else if(FSpions > 0 && !isaPion)
            tempfate = 8;   // nucleon goes to pions
          else if(isaPion && FSpions == 0)
          {
            tempfate = 5;
            // was this two body absorption?
            if(ld-fd == 1)
            {
              int p = 0;
              int n = 0;
              if(mc_er_ID[fd] == 2212)p++;
              if(mc_er_ID[fd] == 2112)n++;
              if(mc_er_ID[ld] == 2212)p++;
              if(mc_er_ID[ld] == 2112)n++;
              // use special codes for producing nn, pn, pp
              if(p==1 && n==1)tempfate = 51;
              if(p==2)tempfate = 52;
              if(n==2)tempfate = 50;
            }
          }
          else if(isaPion && FSpions >= 2)
            tempfate = 8;
          else if(isaPion && FSpions == 1)
          {
            tempfate = 4;
          for(int jj = mc_er_FD[i]; jj <= mc_er_LD[i]; ++jj)
            if((mc_er_ID[jj] == 111 || TMath::Abs(mc_er_ID[jj]) == 211) && mc_er_ID[jj] != mc_er_ID[i]) tempfate = 2;
          }
          else 
          {
            // is not a pion, so is a nucleon.
            // if (verbose) std::cout << " Got unknown fate " << tempfate << std::endl;
            // can't tell 2 CEX from 4 for nucleons.  Just assign to 4.
            tempfate = 4;  // or four, can't tell.
          }
        } // if nucleons and pions
        tempFates[i] = tempfate;
      }   // end status14
    } // loop over mc_er_nPart
  
    return tempFates;
  }

  // Get GENIE binding energy in MeV
  double getGenieBEinMeV(const int A) const
  {
    // From UserPhysicsOptions.xml file
    // these are hard coded, so we can get an exact match
    // but beware if you try to port this code to any new GENIE.

    if (A == 1) return 0; // hydrogen
    if (A == 6) return 17.0; // lithium
    if (A == 12) return 25.0; // carbon
    if (A == 16) return 27.0; // oxygen
    if (A == 24) return 32.0; // magnesium
    if (A == 40) return 29.5; // argon
    if (A == 48) return 30.0; // Ti48
    if (A == 56) return 36.0; // 56 iron
    if (A == 58) return 36.0; // 58 nickel
    if (A >= 206 && A <= 208) return 44.0; // 208 lead

    // Specialty in MINERvA,picked off numbers by hand.
    if (A == 28) return 8.219751;  // silicon
    if (A == 27) return 8.115287;  // aluminum
    if (A == 14) return 7.185166;  // nitrogen
    if (A == 55) return 8.653063;  // iron55 or manganese55
    if (A == 35) return 8.347164;  // chlorine  

    if (A == 4) return 5.0;  // this is rough, 

    // else
    // this is a problem, all other numbers come from a semi-empirical binding energy formula
    // need to back off the QE precision for those.
    return 8.0;   // GENIE defaults to something like this.  Needs to be exact.  Check it.
  }

  // Following functions are copied from Carlos's NuETKI/event/CVUniverse.h.
  // Change Carlos's electron naming to lepton.
  // -- Ziggy
  bool GetHasSignalFSProton() const
  {
    bool hasProton = false;
    int i = GetHighestEnergySignalProtonIndex();
    if (i > -1)
    {
      hasProton = true;
    }
    return hasProton;
  }

  int GetHighestEnergySignalProtonIndex() const
  {
    double highestEnergy = -999;
    int index = -999;
    std::vector<int> FSParticles = GetVecInt("mc_FSPartPDG");
    //fs particle energies in MeV
    std::vector<double> energies = GetVecDouble("mc_FSPartE");
    for (size_t i = 0; i < FSParticles.size(); i++){
      //So I can choose to count any protons, or only protons above our reco threshold. 
      //if (FSParticles[i] == 2212){
      if (FSParticles[i] == 2212 && energies[i] > highestEnergy){
        //require momentum between 450 and 1200 MeV/C , and angle under 70 degrees (same as other TKI analyses)
        
        double protonP = sqrt(pow(GetVecElem("mc_FSPartE", i),2) - pow(M_p, 2)); //in MeV/C
        ROOT::Math::XYZVector p(GetVecElem("mc_FSPartPx", i), GetVecElem("mc_FSPartPy", i), GetVecElem("mc_FSPartPz", i)); 
        ROOT::Math::RotationX r(-3.3 * (pi / 180.));
        double protonTheta = (r(p)).Theta()*(180/pi); //in degrees
        //if (protonP>450 && protonP<1200 && protonTheta<70){
        if (protonP>450 && protonP<1200 && (protonTheta<70 || protonTheta>110))
        { //Testing allowing backwards protons??
          highestEnergy = energies[i];
          index = i;
        }
      }
    }
    return index;
  }

  // Checks for pions & kaons
  bool GetHasFSMeson() const
  {
    std::vector<int> FSParticles = GetVecInt("mc_FSPartPDG");
    bool hasMeson = false;
    for (size_t i = 0; i < FSParticles.size(); i++)
    {
      //std::cout << "final state particle: " << i << ": " << FSParticles[i] << std::endl;
      if (abs(FSParticles[i]) == 211 || abs(FSParticles[i]) == 321 || abs(FSParticles[i]) == 311 || abs(FSParticles[i]) == 130 || abs(FSParticles[i]) == 111)
      {
      	hasMeson = true;
      }
    }
    //std::cout << "hasMeson: " << hasMeson<< std::endl;
    return hasMeson;
  }

  bool GetHasFSPhoton() const
  {
    std::vector<int> FSParticles = GetVecInt("mc_FSPartPDG");
    std::vector<double> energies = GetVecDouble("mc_FSPartE");
    bool hasPhoton = false;
    for (size_t i = 0; i < FSParticles.size(); i++)
    {
      //std::cout << "final state particle: " << i << ": " << FSParticles[i] << std::endl;
      if (abs(FSParticles[i]) == 22 && energies[i] > 10)
      {
      	hasPhoton = true;
      }
    }
    //std::cout << "hasPhoton: " << hasPhoton<< std::endl;
    return hasPhoton;
  }

  // Reco proton Pt vector
  ROOT::Math::XYZVector GetProtonPtVec() const
  {
    ROOT::Math::XYZVector protonP_vec(GetDouble("MasterAnaDev_proton_Px_fromdEdx")/1000., GetDouble("MasterAnaDev_proton_Py_fromdEdx")/1000., GetDouble("MasterAnaDev_proton_Pz_fromdEdx")/1000.);
    ROOT::Math::RotationX r(-3.3 * (pi / 180.)); 
    ROOT::Math::RotationX r2(3.3 * (pi / 180.));

    ROOT::Math::XYZVector protonPt_vec = r(protonP_vec);
    protonPt_vec.SetZ(0);
    protonPt_vec = r2(protonPt_vec);

    return protonPt_vec;
  }

  //Returns a root XYZVector object containing the sum of the proton transverse 3 momentum
  //and the lepton (muon for my CCQEnu study -- Ziggy) transverse 3 momentum
  //which is then used to calculate TKI variables
  //Remember: z direction != beam direction so transverse doesn't exactly mean z components are zero, although they should be small
  ROOT::Math::XYZVector GetDeltaPtVec() const
  {
    ROOT::Math::XYZVector leptonPt_vec = GetLeptonPtVec();
    ROOT::Math::XYZVector protonPt_vec = GetProtonPtVec();
        
    //ROOT::Math::XYZVector deltaP_vec = leptonP_vec + protonP_vec; //sum of the full 3 momenta. Not sure if I need this so commenting it out for now
    ROOT::Math::XYZVector deltaPt_vec = leptonPt_vec + protonPt_vec;
    
    //std::cout << "delta P total (kinda useless?): " << sqrt(deltaP_vec.Mag2()) << std::endl;
    //std::cout << "delta Ptx: " << deltaPt_vec.X() << std::endl;
    //std::cout << "delta Pty: " << deltaPt_vec.Y() << std::endl;
    //std::cout << "delta Pt: " << sqrt(deltaPt_vec.Mag2()) << std::endl;
    //std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n" << std::endl;
    return deltaPt_vec;
  }

  //delta pT (magnitude of the vector), in GeV
  double GetDeltaPt() const
  {
    ROOT::Math::XYZVector deltaPt_vec = GetDeltaPtVec();
    //std::cout << "delta Pt RECO: " << sqrt(deltaPt_vec.Mag2()) << std::endl;
    return sqrt(deltaPt_vec.Mag2()); 
  }

  //Returns a root XYZVector object containing the lepton (muon for my CCQEnu study -- Ziggy) transverse 3 momentum
  ROOT::Math::XYZVector GetLeptonPtVec() const
  {
    std::vector<std::vector<double> > leptonP = GetVecOfVecDouble("prong_part_E");
    ROOT::Math::XYZVector leptonP_vec(leptonP[0][0]/1000., leptonP[0][1]/1000., leptonP[0][2]/1000.);

    ROOT::Math::RotationX r(-3.3 * (pi / 180.)); //This object represents a slight rotation about the x axis so that the theta we get is wrt to the beam direction and not the z axis, which points down to the ground 3.3 degrees. we can then apply r to the two vectors above
    ROOT::Math::RotationX r2(3.3 * (pi / 180.)); //the reverse rotation to get back to lab coords (3.3 is positive now)

    //perform the rotation, basically transforming to new beam based coord system
    ROOT::Math::XYZVector leptonPt_vec = r(leptonP_vec);

    //Now set z to zero in the beam frame to get only the transverse components, then rotate back
    leptonPt_vec.SetZ(0); 
    leptonPt_vec = r2(leptonPt_vec);

    return leptonPt_vec;
  }

  //Phi_t 
  double GetPhiT() const
  {
    ROOT::Math::XYZVector leptonPt_vec = GetLeptonPtVec();
    ROOT::Math::XYZVector protonPt_vec = GetProtonPtVec();

    //angle between the two vectors 
    double numerator = ((-1*leptonPt_vec).Dot(protonPt_vec));
    double denominator = ( sqrt(leptonPt_vec.Mag2()) * sqrt(protonPt_vec.Mag2()) );
    double phi = std::acos(numerator/denominator);

    //std::cout << "TKI phi: " << phi * 180/pi << std::endl;
    return phi * 180/pi;  //return in degrees cause that's how I've set up my bins for now
  }
  
  //Alpha_t, the TKI boosting angle. it's the angle between inverted lepton pT and delta pT
  double GetAlphaT() const
  {
    ROOT::Math::XYZVector leptonPt_vec = GetLeptonPtVec();
    ROOT::Math::XYZVector deltaPt_vec = GetDeltaPtVec();

    //angle between the two vectors
    double numerator = ((-1*leptonPt_vec).Dot(deltaPt_vec));
    double denominator = ( sqrt(leptonPt_vec.Mag2()) * sqrt(deltaPt_vec.Mag2()) );
    double alpha = std::acos(numerator/denominator);

    //std::cout << "boosting angle alpha: " << alpha * 180/pi << std::endl;
    return alpha * 180/pi;  //return alpha in degrees cause that's how I've set up my bins for now
  }

  double GetDeltaPtTrue() const
  {
    int i = GetHighestEnergySignalProtonIndex();
    if (i > -1){
      ROOT::Math::XYZVector leptonP_vec(GetVecElem("mc_primFSLepton", 0), GetVecElem("mc_primFSLepton", 1), GetVecElem("mc_primFSLepton", 2));
      ROOT::Math::XYZVector protonP_vec(GetVecElem("mc_FSPartPx",i), GetVecElem("mc_FSPartPy",i), GetVecElem("mc_FSPartPz",i));

      ROOT::Math::RotationX r(-3.3 * (pi / 180.)); //rotation into beam frame
      ROOT::Math::RotationX r2(3.3 * (pi / 180.)); //rotation back into lab frame

      ROOT::Math::XYZVector leptonPt_vec = r(leptonP_vec);
      ROOT::Math::XYZVector protonPt_vec = r(protonP_vec);
      leptonPt_vec.SetZ(0);
      protonPt_vec.SetZ(0);

      //do I need to rotate back for this? I don't think so but double check
      ROOT::Math::XYZVector deltaPt_vec = leptonPt_vec + protonPt_vec;      
      return sqrt(deltaPt_vec.Mag2())/1000.;
    }
    else {
      return -999; //what do I return for delta pt true if no true protons?
    }
  }

  double GetAlphaTTrue() const
  {
    int i = GetHighestEnergySignalProtonIndex();
    if (i > -1)
    {
      ROOT::Math::XYZVector leptonP_vec(GetVecElem("mc_primFSLepton", 0), GetVecElem("mc_primFSLepton", 1), GetVecElem("mc_primFSLepton", 2));
      ROOT::Math::XYZVector protonP_vec(GetVecElem("mc_FSPartPx",i), GetVecElem("mc_FSPartPy",i), GetVecElem("mc_FSPartPz",i));

      ROOT::Math::RotationX r(-3.3 * (pi / 180.)); //rotation into beam frame
      ROOT::Math::RotationX r2(3.3 * (pi / 180.)); //rotation back into lab frame
      
      ROOT::Math::XYZVector leptonPt_vec = r(leptonP_vec);
      ROOT::Math::XYZVector protonPt_vec = r(protonP_vec);
      leptonPt_vec.SetZ(0);
      protonPt_vec.SetZ(0);

      //Still in beam frame
      ROOT::Math::XYZVector deltaPt_vec = leptonPt_vec + protonPt_vec;

      double numerator = ((-1*leptonPt_vec).Dot(deltaPt_vec));
      double denominator = ( sqrt(leptonPt_vec.Mag2()) * sqrt(deltaPt_vec.Mag2()) );
      double alpha = std::acos(numerator/denominator);
      
      return alpha * 180/pi;
    }
    else
    {
      return -999;
    }
  }

  double GetPhiTTrue() const
  {
    int i = GetHighestEnergySignalProtonIndex();
    if (i > -1)
    {
      ROOT::Math::XYZVector leptonP_vec(GetVecElem("mc_primFSLepton", 0), GetVecElem("mc_primFSLepton", 1), GetVecElem("mc_primFSLepton", 2));
      ROOT::Math::XYZVector protonP_vec(GetVecElem("mc_FSPartPx",i), GetVecElem("mc_FSPartPy",i), GetVecElem("mc_FSPartPz",i));

      ROOT::Math::RotationX r(-3.3 * (pi / 180.)); //rotation into beam frame
      ROOT::Math::RotationX r2(3.3 * (pi / 180.)); //rotation back into lab frame
      
      ROOT::Math::XYZVector leptonPt_vec = r(leptonP_vec);
      ROOT::Math::XYZVector protonPt_vec = r(protonP_vec);
      leptonPt_vec.SetZ(0);
      protonPt_vec.SetZ(0);

      double numerator = ((-1*leptonPt_vec).Dot(protonPt_vec));
      double denominator = ( sqrt(leptonPt_vec.Mag2()) * sqrt(protonPt_vec.Mag2()) );
      double phi = std::acos(numerator/denominator);
      
      return phi * 180/pi;
    }
    else
    {
      return -999;
    }
  }

  //============================================================================
  //============================================================================


  //Muon kinematics
  double GetMuonPT() const //GeV/c
  {
    return GetPmu()/1000. * sin(GetThetamu());
  }

  double GetMuonPz() const //GeV/c
  {
    return GetPmu()/1000. * cos(GetThetamu());
  }

  double GetMuonPTTrue() const //GeV/c
  {
    return GetPlepTrue()/1000. * sin(GetThetalepTrue());
  }

  double GetMuonPzTrue() const //GeV/c
  {
    return GetPlepTrue()/1000. * cos(GetThetalepTrue());
  }

  double GetEmuGeV() const //GeV
  {
    return GetEmu()/1000.;
  }

  double GetElepTrueGeV() const //GeV
  {
    return GetElepTrue()/1000.;
  }

  int GetInteractionType() const {
    return GetInt("mc_intType");
  }

  int GetTargetNucleon() const {
    return GetInt("mc_targetNucleon");
  }
  
  double GetBjorkenXTrue() const {
    return GetDouble("mc_Bjorkenx");
  }

  double GetBjorkenYTrue() const {
    return GetDouble("mc_Bjorkeny");
  }

  virtual bool IsMinosMatchMuon() const {
    return GetInt("has_interaction_vertex") == 1;
  }
  
  ROOT::Math::XYZTVector GetVertex() const
  {
    ROOT::Math::XYZTVector result;
    result.SetCoordinates(GetVec<double>("vtx").data());
    return result;
  }

  ROOT::Math::XYZTVector GetTrueVertex() const
  {
    ROOT::Math::XYZTVector result;
    result.SetCoordinates(GetVec<double>("mc_vtx").data());
    return result;
  }

  virtual int GetTDead() const
  {
    return GetInt("phys_n_dead_discr_pair_upstream_prim_track_proj");
  }
  
  //TODO: If there was a spline correcting Eavail, it might not really be Eavail.
  //      Our energy correction spline, one of at least 2 I know of, corrects q0
  //      so that we get the right neutrino energy in an inclusive sample.  So,
  //      this function could be correcting for neutron energy which Eavail should
  //      not do.
  virtual double GetEavail() const
  {
    return GetDouble("recoilE_SplineCorrected");
  }
  
  virtual double GetQ2Reco() const
  {
    return GetDouble("qsquared_recoil");
  }

  //GetRecoilE is designed to match the NSF validation suite
  virtual double GetRecoilE() const
  {
    return GetVecElem("recoil_summed_energy", 0);
  }
  
  virtual double Getq3() const
  {
    double eavail = GetEavail()/pow(10,3);
    double q2 = GetQ2Reco() / pow(10,6);
    double q3mec = sqrt(eavail*eavail + q2);
    return q3mec;
  }
   
  virtual int GetCurrent() const { return GetInt("mc_current"); }

  virtual int GetTruthNuPDG() const { return GetInt("mc_incoming"); }

  virtual double GetMuonQP() const
  {
    return GetDouble((GetAnaToolName() + "_minos_trk_qp").c_str());
  }

  //Some functions to match CCQENuInclusive treatment of DIS weighting. Name matches same Dan area as before.
  virtual double GetTrueExperimentersQ2() const
  {
    double Enu = GetEnuTrue(); //MeV
    double Emu = GetElepTrue(); //MeV
    double thetaMu = GetThetalepTrue();
    return 4.0*Enu*Emu*pow(sin(thetaMu/2.0),2.0);//MeV^2
  }

  virtual double CalcTrueExperimentersQ2(double Enu, double Emu, double thetaMu) const
  {
    return 4.0*Enu*Emu*pow(sin(thetaMu/2.0),2.0);//MeV^2
  }

  virtual double GetTrueExperimentersW() const
  {
    double nuclMass = M_nucleon;
    int struckNucl = GetTargetNucleon();
    if (struckNucl == PDG_n)
    {
      nuclMass=M_n;
    }
    else if (struckNucl == PDG_p)
    {
      nuclMass=M_p;
    }
    double Enu = GetEnuTrue();
    double Emu = GetElepTrue();
    double thetaMu = GetThetalepTrue();
    double Q2 = CalcTrueExperimentersQ2(Enu, Emu, thetaMu);
    return TMath::Sqrt(pow(nuclMass,2) + 2.0*(Enu-Emu)*nuclMass - Q2);
  }

  //Still needed for some systematics to compile, but shouldn't be used for reweighting anymore.
  protected:
  #include "PlotUtils/WeightFunctions.h" // Get*Weight
};

#endif
