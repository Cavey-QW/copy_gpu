#include "berendsen_controller.h"

#include <thrust/device_ptr.h>

#include <cmath>

#include "device_types.h"
#include "unit_factor.h"
#include "update_temperature_op.h"

rbmd::Real _test_temperature;
BerendsenController::BerendsenController() {
  CHECK_RUNTIME(MALLOC(&_d_temp_contrib, sizeof(rbmd::Real)));
}
BerendsenController::~BerendsenController() {
  CHECK_RUNTIME(FREE(_d_temp_contrib));
};

void BerendsenController::Init() {
  auto temperature_array=DataManager::getInstance().getConfigData()->
GetArray<rbmd::Real>("temperature", "execution"); //[1.0,1.0,0.1]
  _temperature_start = temperature_array[0];
  _temperature_stop = temperature_array[1];
  _temperature_damp = temperature_array[2];

  _dt =  DataManager::getInstance().getConfigData()->Get<rbmd::Real>(
          "timestep", "execution");//0.001

  auto unit = DataManager::getInstance().getConfigData()->Get
    <std::string>("unit", "init_configuration", "read_data");

  UNIT unit_factor = unit_factor_map[unit];
  switch (unit_factor) {
    case UNIT::LJ:
      _mvv2e = UnitFactor<UNIT::LJ>::_mvv2e;
      _kB = UnitFactor<UNIT::LJ>::_kb;
      break;
    case UNIT::REAL:
      _mvv2e = UnitFactor<UNIT::REAL>::_mvv2e;
      _kB = UnitFactor<UNIT::REAL>::_kb;
      break;
    default:
      break;
  }
}

void BerendsenController::Update() {
  //ComputeTemp();

  UpdataVelocity();
}

void BerendsenController::ComputeTemp() {
  extern int test_current_step;
  rbmd::Id num_atoms = *(_structure_info_data->_num_atoms);

  CHECK_RUNTIME(MEMSET(_d_temp_contrib, 0, sizeof(rbmd::Real)));

  op::ComputeTemperatureOp<device::DEVICE_GPU> compute_temperature_op;
  compute_temperature_op(num_atoms, _mvv2e,
      thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
      thrust::raw_pointer_cast(_device_data->_d_mass.data()),
      thrust::raw_pointer_cast(_device_data->_d_vx.data()),
      thrust::raw_pointer_cast(_device_data->_d_vy.data()),
      thrust::raw_pointer_cast(_device_data->_d_vz.data()), _d_temp_contrib);

  CHECK_RUNTIME(MEMCPY(&_temp_sum, _d_temp_contrib, sizeof(rbmd::Real), D2H));

  bool available_shake = true; //TODO: need to judge

  if (available_shake)  // H2O / NACl / EAM ...
  {
    bool shake = DataManager::getInstance().getConfigData()->GetJudge<bool>( "fix_shake", "hyper_parameters", "extend");
    if (shake) {
      _temp = 0.5 * _temp_sum / ((3 * num_atoms - num_atoms - 3) * _kB / 2.0);
    } else {
      _temp = 0.5 * _temp_sum / ((3 * num_atoms - 3) * _kB / 2.0);
    }
  } else  // PEO
  {
    _temp = 0.5 * _temp_sum / ((3 * num_atoms - 3) * _kB / 2.0);
  }
  _test_temperature = _temp;

  std::cout << "_temp=" << _temp << std::endl;
  // out
  std::ofstream outfile("temp.txt", std::ios::app);
  outfile << test_current_step << " " << _temp << std::endl;
  outfile.close();

  // CHECK_RUNTIME(FREE(temp_contrib));
}

void BerendsenController::UpdataVelocity() {
  rbmd::Real coeff_Berendsen =
      SQRT(1.0 + (_dt / _temperature_damp) * (_temperature_start/ _temp - 1.0));

  op::UpdataVelocityRescaleOp<device::DEVICE_GPU> updata_velocity_op;
  updata_velocity_op(*(_structure_info_data->_num_atoms), coeff_Berendsen,
                     thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                     thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                     thrust::raw_pointer_cast(_device_data->_d_vz.data()));
}
