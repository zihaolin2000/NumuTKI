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

class CVUniverse : public PlotUtils::MinervaUniverse {

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

  double GetMuonTheta() const
  {
    return GetThetamu();
  }

  double GetMultiplicity() const
  {
    return GetInt("multiplicity");
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
    for (int i = 0; i < FSParticles.size(); i++){
      //So I can choose to count any protons, or only protons above our reco threshold. 
      //if (FSParticles[i] == 2212){
      if (FSParticles[i] == 2212 && energies[i] > highestEnergy){
        //require momentum between 450 and 1200 MeV/C , and angle under 70 degrees (same as other TKI analyses)
        
        double protonP = sqrt(pow(GetVecElem("mc_FSPartE", i),2) - pow(M_p, 2)); //in MeV/C
        ROOT::Math::XYZVector p(GetVecElem("mc_FSPartPx", i), GetVecElem("mc_FSPartPy", i), GetVecElem("mc_FSPartPz", i)); 
        ROOT::Math::RotationX r(-3.3 * (pi / 180.));
        double protonTheta = (r(p)).Theta()*(180/pi); //in degrees
        //if (protonP>450 && protonP<1200 && protonTheta<70){
        if (protonP>450 && protonP<1200 && (protonTheta<70 || protonTheta>110)){ //Testing allowing backwards protons??
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
    for (int i = 0; i < FSParticles.size(); i++)
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
    for (int i = 0; i < FSParticles.size(); i++)
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

  //Returns a root XYZVector object containing the lepton (muon for my CCQEnu study -- Ziggy) transverse 3 momentum
  ROOT::Math::XYZVector GetProtonPtVec() const{
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
  ROOT::Math::XYZVector GetDeltaPtVec() const{
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
  ROOT::Math::XYZVector GetLeptonPtVec() const{
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

  virtual int GetTDead() const {
    return GetInt("phys_n_dead_discr_pair_upstream_prim_track_proj");
  }
  
  //TODO: If there was a spline correcting Eavail, it might not really be Eavail.
  //      Our energy correction spline, one of at least 2 I know of, corrects q0
  //      so that we get the right neutrino energy in an inclusive sample.  So,
  //      this function could be correcting for neutron energy which Eavail should
  //      not do.
  virtual double GetEavail() const {
    return GetDouble("recoilE_SplineCorrected");
  }
  
  virtual double GetQ2Reco() const{
    return GetDouble("qsquared_recoil");
  }

  //GetRecoilE is designed to match the NSF validation suite
  virtual double GetRecoilE() const {
    return GetVecElem("recoil_summed_energy", 0);
  }
  
  virtual double Getq3() const{
    double eavail = GetEavail()/pow(10,3);
    double q2 = GetQ2Reco() / pow(10,6);
    double q3mec = sqrt(eavail*eavail + q2);
    return q3mec;
  }
   
  virtual int GetCurrent() const { return GetInt("mc_current"); }

  virtual int GetTruthNuPDG() const { return GetInt("mc_incoming"); }

  virtual double GetMuonQP() const {
    return GetDouble((GetAnaToolName() + "_minos_trk_qp").c_str());
  }

  //Some functions to match CCQENuInclusive treatment of DIS weighting. Name matches same Dan area as before.
  virtual double GetTrueExperimentersQ2() const {
    double Enu = GetEnuTrue(); //MeV
    double Emu = GetElepTrue(); //MeV
    double thetaMu = GetThetalepTrue();
    return 4.0*Enu*Emu*pow(sin(thetaMu/2.0),2.0);//MeV^2
  }

  virtual double CalcTrueExperimentersQ2(double Enu, double Emu, double thetaMu) const{
    return 4.0*Enu*Emu*pow(sin(thetaMu/2.0),2.0);//MeV^2
  }

  virtual double GetTrueExperimentersW() const {
    double nuclMass = M_nucleon;
    int struckNucl = GetTargetNucleon();
    if (struckNucl == PDG_n){
      nuclMass=M_n;
    }
    else if (struckNucl == PDG_p){
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
