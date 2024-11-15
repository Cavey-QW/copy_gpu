#include "langevin_controller.h"
#include <thrust/device_ptr.h>
#include "device_types.h"
#include "neighbor_list/include/linked_cell/linked_cell_locator.h"
#include "rbmd_define.h"
#include "unit_factor.h"
#include "update_temperature_op.h"
#include <random>
#include <numeric>
#include <cmath>
LangevinController::LangevinController() {
  CHECK_RUNTIME(MALLOC(&_d_temp_contrib, sizeof(rbmd::Real)));
}
LangevinController::~LangevinController() {
  CHECK_RUNTIME(FREE(_d_temp_contrib));
}

void LangevinController::Init() {
  _dt = DataManager::getInstance().getConfigData()->Get<rbmd::Real>(
      "timestep", "execution");    // TODO: read values from json file;
  auto unit = "LJ";      // TODO: the conditions of LangevinController (only FarForce?)                    
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

void LangevinController::Update() {
  //ComputeTemp();

  UpdataForce();
}

void LangevinController::ComputeTemp() {
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

    bool available_shake = false;

    if (available_shake)  // H2O / NACl / EAM ...
    {
        bool shake = true;
        if (shake) {
            _temp = 0.5 * _temp_sum / ((3 * num_atoms - num_atoms - 3) * _kB / 2.0);
        }
        else {
            _temp = 0.5 * _temp_sum / ((3 * num_atoms - 3) * _kB / 2.0);
        }
    }
    else  // PEO
    {
        _temp = 0.5 * _temp_sum / ((3 * num_atoms - 3) * _kB / 2.0);
    }

    std::cout << "_temp=" << _temp << std::endl;
    // out
    std::ofstream outfile("temp.txt", std::ios::app);
    outfile << test_current_step << " " << _temp << std::endl;
    outfile.close();

    // CHECK_RUNTIME(FREE(temp_contrib));
}

void LangevinController::UpdataForce() {

  bool random = true; // TODO: Related to the of velocity_type  
  auto temp_value = FetchSample_1D(random, rbmd::Real(0.0), rbmd::Real(1.0));
  std::vector<rbmd::Real> gaussian(3, temp_value);
  
  rbmd::Real kbT = 1;
  rbmd::Real gamma = 100.0;
  op::UpdataForceLangevinOp<device::DEVICE_GPU> updata_force_op;
  updata_force_op(*(_structure_info_data->_num_atoms), gaussian[0], gaussian[1], gaussian[2], kbT, gamma,_dt,
                  thrust::raw_pointer_cast(_device_data->_d_mass.data()),
                  thrust::raw_pointer_cast(_device_data->_d_vx.data()),
                  thrust::raw_pointer_cast(_device_data->_d_vy.data()),
                  thrust::raw_pointer_cast(_device_data->_d_vz.data()),
                  thrust::raw_pointer_cast(_device_data->_d_fx.data()),
                  thrust::raw_pointer_cast(_device_data->_d_fy.data()),
                  thrust::raw_pointer_cast(_device_data->_d_fz.data()));
}

rbmd::Real LangevinController:: FetchSample_1D(const bool& random, const rbmd::Real& mu,const rbmd::Real& sigma) // Fetch 1D sample from Gaussion contribution
{
    rbmd::Real U1, U2, epsilon;
    epsilon = 1e-6;
    if (random)
    {
        do
        {
          U1 = RandomValue(0.0, 1.0);
        } while (U1 < epsilon);
        U2 = RandomValue(0.0, 1.0);
        std::vector<rbmd::Real> ChooseSample{ 0.0, 0.0 };
        ChooseSample[0] = sigma * std::sqrt(-2.0 * std::log(U1)) * std::cos(2 * M_PI * U2) + mu;
        //ChooseSample[1] = sigma * vtkm::Sqrt(-2.0 * vtkm::Log(U1)) * sin(2 * vtkm::Pi() * U2) + mu;
        return ChooseSample[0];
    }
    else
    {
        U1 = 0.5;
        U2 = 0.5;
        std::vector<rbmd::Real> ChooseSample{ 0.0, 0.0 };
        ChooseSample[0] = sigma * std::sqrt(-2.0 * std::log(U1)) * std::cos(2 * M_PI * U2) + mu;
        //ChooseSample[1] = sigma * vtkm::Sqrt(-2.0 * vtkm::Log(U1)) * sin(2 * vtkm::Pi() * U2) + mu;
        return ChooseSample[0];
    }
}

rbmd::Real LangevinController::RandomValue(const rbmd::Real& Min, const rbmd::Real& Max)
{
    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<rbmd::Real> dis(Min, Max);
    return dis(gen);
};
