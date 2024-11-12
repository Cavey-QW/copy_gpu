#pragma once
#include <map>
#include <array>
#include <vector>

class Domain {
private:
#pragma region 变量
  //! rank of the local process
  int _localRank;

  //! Potential of the local process
  double _localUpot;
  //! by Stefan Becker: component specific potential of the local process (fluid-fluid and fluid-solid, but not solid-solid)
  double _localUpotCspecif{};
  //! by Stefan Becker: number of fluid components specified in the config file
  //! --> used to determine the average component specific potential energy
  unsigned _numFluidComponent{};

  //! Virial of the local process
  double _localVirial;
  //! global Potential
  double _globalUpot;
  //! global component specific potential (fluid-fluid and fluid-solid but NOT solid-solid)
  double _globalUpotCspecif{};
  //! global virial
  double _globalVirial;
  //! global density
  double _globalRho;
  //! global Number of Molecules
  //! @todo redundancy?
  unsigned long _globalNumMolecules{};
  std::array<double,3> offset{};
  typedef unsigned long int _uint64_t;
  typedef _uint64_t uint64_t;
  template<typename T>       //TODO 先这样
  class CommVar
  {
  public:
    T local;
    T global;
  };
  CommVar<uint64_t> _maxMoleculeID{};

  //! side length of the cubic simulation box
  double _globalLength[3]{};

  //! does a componentwise thermostat apply?
  bool _componentwiseThermostat;
  //! thermostat IDs. negative: no thermostat, 0: global, positive: componentwise
  //! in the case of a componentwise thermostat, all components are assigned
  //! a thermostat ID different from zero.
  std::map<int, int> _componentToThermostatIdMap;

  //! _localThermostatN[0] and _universalThermostatN[0] are always the total number
  //! of particles in the subdomain and, respectively, the entire domain
  std::map<int, unsigned long> _localThermostatN;
  std::map<int, unsigned long> _universalThermostatN;
  std::map<int, unsigned long> _localRotationalDOF;
  std::map<int, unsigned long> _universalRotationalDOF;
  //! _globalTemperatureMap[0] is always the temperature of the whole system,
  //! including components to which no thermostat is applied.
  //! The _globalTemperatureMap stores actual CURRENT temperatures, whereas
  //! the temperature objective of the thermostat is stored in _universalTargetTemperature
  std::map<int, double> _globalTemperatureMap;
  std::map<int, double> _universalTargetTemperature;
  std::map<int, double> _universalBTrans;
  std::map<int, double> _universalBRot;
  //! should the directed movement be subtracted when calculating the temperature?
  std::map<int, bool> _universalUndirectedThermostat;
  //! stores the velocity of the directed movement
  std::map<int, std::array<double, 3>> _universalThermostatDirectedVelocity;
  std::map<int, std::array<double, 3>> _localThermostatDirectedVelocity;

  /* FIXME: This info should go into an ensemble class */
  bool _universalNVE;

  //! computation of the isochoric heat capacity
  unsigned _globalUSteps;
  double _globalSigmaU;
  double _globalSigmaUU;
  //! which components should be considered?
  std::map<unsigned, bool> _universalProfiledComponents;
  double _universalProfiledComponentMass{};  // set from outside
  double _universalLambda{};  // set from outside
  float _globalDecisiveDensity{};  // set from outside

  int _universalSelectiveThermostatCounter;
  int _universalSelectiveThermostatWarning;
  int _universalSelectiveThermostatError;

  //! local sum (over all molecules) of the mass multiplied with the squared velocity
  std::map<int, double> _local2KETrans;
  //! local sum (over all molecules) of the moment of inertia
  //! multiplied with the squared  rotational velocity
  std::map<int, double> _local2KERot;

  //! reaction field
  //!
  //! This is neither "local" nor "global" but a parameter of the reaction field method.
  //! (So one might regard it as "global" formally.)
  //! For a description of the reaction field method cf. the dissertation of J. Stoll.
  //! It was introduced by J. A. Barker and R. O. Watts, Mol. Phys. 26: 789-792 (1973).
  double _epsilonRF{};

  //! Global potential correction for the error made by the cutoff
  double _UpotCorr{};
  //! Global virial correction for the error made by the cutoff
  double _VirialCorr{};

  //! parameter streams for each possible pair of molecule-types
  // Comp2Param _comp2params;   TODO 暂时不要
  //! modified Lorentz-Berthelot mixing rule parameters
  //! @todo more explanation
  std::vector<double> _mixcoeff;

  // explosion heuristics, NOTE: turn off when using slab thermostat
  bool _bDoExplosionHeuristics;

  int _global_rbe_p_number = -1;

  double  _global_rbe_alpha_number = -0.1;


#pragma endregion


  Domain();
  Domain(Domain &domain);

  Domain& operator=(Domain &domain);

public:
  //! The constructor sets _localRank to rank and initializes all member variables
  Domain(int rank);

  double GetGlobalLength(int d) const;

  void SetGlobalLength(int dim,double dim_global_length);
  void SetGlobaMinMax(std::array<double,3> min,std::array<double,3> max);

  void  SetGlobalRbePnumber(int rbe_p_number);

  int  GetGlobalRbePnumber();

  void  SetGlobalRbeAlphaPnumber(double rbe_alpha_number);

  double  GetGlobalRbeAlphaPnumber();

  double GetOffset(int dim);




  //TODO  setGlobalTemperature



};
