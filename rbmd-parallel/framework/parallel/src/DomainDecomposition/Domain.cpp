#include "Domain.h"

#include <DomainDecompBase.h>

Domain::Domain(int rank) {
    _localRank = rank;
    _localUpot = 0;
    _localVirial = 0;
    _globalUpot = 0;
    _globalVirial = 0;
    _globalRho = 0;

    this->_componentToThermostatIdMap = std::map<int, int>();
    this->_localThermostatN = std::map<int, unsigned long>();
    this->_localThermostatN[-1] = 0;
    this->_localThermostatN[0] = 0;
    this->_universalThermostatN = std::map<int, unsigned long>();
    this->_universalThermostatN[-1] = 0;
    this->_universalThermostatN[0] = 0;
    this->_localRotationalDOF = std::map<int, unsigned long>();
    this->_localRotationalDOF[-1] = 0;
    this->_localRotationalDOF[0] = 0;
    this->_universalRotationalDOF = std::map<int, unsigned long>();
    this->_universalRotationalDOF[-1] = 0;
    this->_universalRotationalDOF[0] = 0;
    this->_globalLength[0] = 0;
    this->_globalLength[1] = 0;
    this->_globalLength[2] = 0;
    this->_universalBTrans = std::map<int, double>();
    this->_universalBTrans[0] = 1.0;
    this->_universalBRot = std::map<int, double>();
    this->_universalBRot[0] = 1.0;
    this->_universalTargetTemperature = std::map<int, double>();
    this->_universalTargetTemperature[0] = 1.0;
    this->_globalTemperatureMap = std::map<int, double>();
    this->_globalTemperatureMap[0] = 1.0;
    this->_local2KETrans[0] = 0.0;
    this->_local2KERot[0] = 0.0;

    this->_universalNVE = false;
    this->_globalUSteps = 0;
    this->_globalSigmaU = 0.0;
    this->_globalSigmaUU = 0.0;
    this->_componentwiseThermostat = false;
#ifdef COMPLEX_POTENTIAL_SET    //TODO 暂时不知道
    this->_universalUndirectedThermostat = std::map<int, bool>();
    this->_universalThermostatDirectedVelocity = std::map<int, std::array<double,3> >();
    this->_localThermostatDirectedVelocity = std::map<int, std::array<double,3> >();
#endif
    this->_universalSelectiveThermostatCounter = 0;
    this->_universalSelectiveThermostatWarning = 0;
    this->_universalSelectiveThermostatError = 0;
    // explosion heuristics, NOTE: turn off when using slab thermostat
    _bDoExplosionHeuristics = true;

}

double Domain::GetGlobalLength(int d) const {
    return _globalLength[d];
}

void Domain::SetGlobalLength(int dim,double dim_global_length) {
    _globalLength[dim] = dim_global_length;
}

void Domain::SetGlobaMinMax(std::array<double, 3> min, std::array<double, 3> max) {
    for (int i = 0; i < DIM; ++i) {
        _globalLength[i] = max[i]-min[i];
        offset[i] = min[i];
    }
}

void Domain::SetGlobalRbePnumber(int rbe_p_number) {
    this->_global_rbe_p_number = rbe_p_number;
}

int Domain::GetGlobalRbePnumber() {
    return this->_global_rbe_p_number;
}

void Domain::SetGlobalRbeAlphaPnumber(double rbe_alpha_number) {
    this->_global_rbe_alpha_number = rbe_alpha_number;

}

double Domain::GetGlobalRbeAlphaPnumber() {
    return this->_global_rbe_alpha_number;
}

double Domain::GetOffset(int dim) {
    return offset[dim];
}






