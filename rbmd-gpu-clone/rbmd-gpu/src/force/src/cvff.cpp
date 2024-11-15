#include "cvff.h"

#include <thrust/device_ptr.h>
#include "thrust/sort.h"

#include "../../common/device_types.h"
#include "../../common/rbmd_define.h"
#include "../../common/types.h"
#include "../common/unit_factor.h"
#include "ljforce_op/ljforce_op.h"
#include "../common/RBEPSample.h"
#include "../common/erf_table.h"
#include "neighbor_list/include/linked_cell/linked_cell_locator.h"
#include "neighbor_list/include/neighbor_list_builder/half_neighbor_list_builder.h"
#include "neighbor_list/include/neighbor_list_builder/full_neighbor_list_builder.h"
#include "neighbor_list/include/neighbor_list_builder/rbl_full_neighbor_list_builder.h"
// #include <hipcub/hipcub.hpp>
// #include <hipcub/backend/rocprim/block/block_reduce.hpp>

extern int test_current_step;
extern std::map<std::string, UNIT> unit_factor_map;

CVFF::CVFF()
{
  _rbl_neighbor_list_builder = std::make_shared<RblFullNeighborListBuilder>();
  _neighbor_list_builder = std::make_shared<FullNeighborListBuilder>();

  CHECK_RUNTIME(MALLOC(&_d_total_evdwl, sizeof(rbmd::Real)));
  CHECK_RUNTIME(MALLOC(&_d_total_ecoul, sizeof(rbmd::Real)));
  CHECK_RUNTIME(MALLOC(&_d_total_e_specialcoul, sizeof(rbmd::Real)));
  CHECK_RUNTIME(MALLOC(&_d_total_ebond, sizeof(rbmd::Real)));
  CHECK_RUNTIME(MALLOC(&_d_total_eangle, sizeof(rbmd::Real)));
  CHECK_RUNTIME(MALLOC(&_d_total_edihedral, sizeof(rbmd::Real)));
}

CVFF::~CVFF()
{
  CHECK_RUNTIME(FREE(_d_total_evdwl));
  CHECK_RUNTIME(FREE(_d_total_ecoul));
  CHECK_RUNTIME(FREE(_d_total_e_specialcoul));
  CHECK_RUNTIME(FREE(_d_total_ebond));
  CHECK_RUNTIME(FREE(_d_total_eangle));
  CHECK_RUNTIME(FREE(_d_total_edihedral));

}

void CVFF::Init()
{
  auto unit = DataManager::getInstance().getConfigData()->Get
    <std::string>("unit", "init_configuration", "read_data");
  UNIT unit_factor = unit_factor_map[unit];

  switch (unit_factor) {
    case UNIT::LJ:
      _qqr2e = UnitFactor<UNIT::LJ>::_qqr2e;
    break;

    case UNIT::REAL:
      _qqr2e = UnitFactor<UNIT::REAL>::_qqr2e;
    break;

    default:
      break;
  }

   _cut_off = DataManager::getInstance().getConfigData()->Get
    <rbmd::Real>("cut_off", "hyper_parameters", "neighbor");
   _neighbor_type = DataManager::getInstance().getConfigData()->Get
      <std::string>("type", "hyper_parameters", "neighbor");

   _coulomb_type =DataManager::getInstance().getConfigData()->Get<std::string>(
        "type", "hyper_parameters", "coulomb");
    _RBE_P = DataManager::getInstance().getConfigData()->Get<rbmd::Id>(
      "coulomb_sample_num", "hyper_parameters", "coulomb");
    _alpha = DataManager::getInstance().getConfigData()->Get<rbmd::Real>(
  "alpha", "hyper_parameters", "coulomb");
    _Kmax = DataManager::getInstance().getConfigData()->Get<rbmd::Id>(
"kmax", "hyper_parameters", "coulomb");
    _num_k =  POW(2 * _Kmax + 1,3.0) - 1;
    ERFInit();
    RBEInit(_device_data->_d_box,_alpha,_RBE_P);
}

void CVFF::Execute()
{
  ComputeLJCutCoulForce();
  ComputeSpecialCoulForce();
  ComputeKspaceForce();

  ComputeBondForce();
  ComputeAngleForce();
  //ComputeDihedralForce();
  SumForces();
}

void CVFF::ComputeLJCutCoulForce()
{
  if ("RBL" ==_neighbor_type)
  {
    ComputeLJRBL();
  }
  else
  {
    ComputeLJVerlet();
  }
}

void CVFF::ComputeLJRBL()
{
    // rbl_neighbor_list_build
    auto start = std::chrono::high_resolution_clock::now();
    _rbl_list = _rbl_neighbor_list_builder->Build();

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<rbmd::Real> duration = end - start;
    std::cout << "构建RBL邻居列表耗时" << duration.count() << "秒" << std::endl;

    // compute force
    const auto r_core =
      DataManager::getInstance().getConfigData()->Get<rbmd::Real>(
          "r_core", "hyper_parameters", "neighbor");

     const auto neighbor_sample_num =
     DataManager::getInstance().getConfigData()->Get<rbmd::Id>(
         "neighbor_sample_num", "hyper_parameters", "neighbor");

    auto num_atoms = *(_structure_info_data->_num_atoms);
    op::SpecialLJCutCoulRBLForceOp<device::DEVICE_GPU> lj_cut_coul_rbl_force_op;
    lj_cut_coul_rbl_force_op(_device_data->_d_box,_device_data->_d_erf_table,
      r_core, _cut_off, num_atoms,neighbor_sample_num,
      _rbl_list->_selection_frequency,_alpha,_qqr2e,
        thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
        thrust::raw_pointer_cast(_device_data->_d_atoms_id.data()),
        thrust::raw_pointer_cast(_device_data->_d_sigma.data()),
        thrust::raw_pointer_cast(_device_data->_d_eps.data()),
        thrust::raw_pointer_cast(_rbl_list->_start_idx.data()),
        thrust::raw_pointer_cast(_rbl_list->_end_idx.data()),
        thrust::raw_pointer_cast(_rbl_list->_d_neighbors.data()),
        thrust::raw_pointer_cast(_rbl_list->_d_random_neighbor.data()),
        thrust::raw_pointer_cast(_rbl_list->_d_random_neighbor_num.data()),
        thrust::raw_pointer_cast(_device_data->_d_special_ids.data()),
        thrust::raw_pointer_cast(_device_data->_d_special_weights.data()),
        thrust::raw_pointer_cast(_device_data->_d_special_offsets.data()),
        thrust::raw_pointer_cast(_device_data->_d_special_count.data()),
        thrust::raw_pointer_cast(_device_data->_d_charge.data()),
        thrust::raw_pointer_cast(_device_data->_d_px.data()),
        thrust::raw_pointer_cast(_device_data->_d_py.data()),
        thrust::raw_pointer_cast(_device_data->_d_pz.data()),
        thrust::raw_pointer_cast(_device_data->_d_force_ljcoul_x.data()),
        thrust::raw_pointer_cast(_device_data->_d_force_ljcoul_y.data()),
        thrust::raw_pointer_cast(_device_data->_d_force_ljcoul_z.data()));

    _corr_value_x =
        thrust::reduce(_device_data->_d_force_ljcoul_x.begin(), _device_data->_d_force_ljcoul_x.end(),
                       0.0f, thrust::plus<rbmd::Real>()) /num_atoms;
    _corr_value_y =
        thrust::reduce(_device_data->_d_force_ljcoul_y.begin(), _device_data->_d_force_ljcoul_y.end(),
                       0.0f, thrust::plus<rbmd::Real>()) /num_atoms;
    _corr_value_z =
        thrust::reduce(_device_data->_d_force_ljcoul_z.begin(), _device_data->_d_force_ljcoul_z.end(),
                       0.0f, thrust::plus<rbmd::Real>()) /num_atoms;

    // fix RBL:   rbl_force = f - corr_value
    op::FixRBLForceOp<device::DEVICE_GPU> fix_rbl_force_op;
    fix_rbl_force_op(num_atoms, _corr_value_x, _corr_value_y, _corr_value_z,
                        thrust::raw_pointer_cast(_device_data->_d_force_ljcoul_x.data()),
                        thrust::raw_pointer_cast(_device_data->_d_force_ljcoul_y.data()),
                        thrust::raw_pointer_cast(_device_data->_d_force_ljcoul_z.data()));

    //energy
    ComputeLJCoulEnergy();

    // //out
    // std::vector<rbmd::Real> h_force_ljcoul_x(num_atoms);
    // std::vector<rbmd::Real> h_force_ljcoul_y(num_atoms);
    // std::vector<rbmd::Real> h_force_ljcoul_z(num_atoms);
    // thrust::copy(_device_data->_d_force_ljcoul_x.begin(),
    //   _device_data->_d_force_ljcoul_x.end(), h_force_ljcoul_x.begin());
    // thrust::copy(_device_data->_d_force_ljcoul_y.begin(),
    // _device_data->_d_force_ljcoul_y.end(), h_force_ljcoul_y.begin());
    // thrust::copy(_device_data->_d_force_ljcoul_z.begin(),
    // _device_data->_d_force_ljcoul_z.end(), h_force_ljcoul_z.begin());
    //
    // thrust::host_vector<rbmd::Real> h_atoms_id = _device_data->_d_atoms_id;
    // std::ofstream output_file("output_force_special_lj_RBL.txt");
    // for (size_t i = 0; i < h_force_ljcoul_x.size(); ++i)
    // {
    //   output_file << h_atoms_id[i]<< " " << h_force_ljcoul_x[i] << " "<<h_force_ljcoul_y[i] <<" "
    //   << h_force_ljcoul_z[i] << std::endl;
    // }
    // output_file.close();
}

void CVFF::ComputeLJVerlet()
{
  //neighbor_list_build
  auto start = std::chrono::high_resolution_clock::now();
  _list = _neighbor_list_builder->Build();

  auto end = std::chrono::high_resolution_clock::now();

  std::chrono::duration<rbmd::Real> duration = end - start;
  std::cout << "构建verlet-list耗时" << duration.count() << "秒" << std::endl;

  //

  rbmd::Real h_total_evdwl = 0.0;
  rbmd::Real h_total_ecoul = 0.0;

  CHECK_RUNTIME(MEMSET(_d_total_evdwl, 0, sizeof(rbmd::Real)));
  CHECK_RUNTIME(MEMSET(_d_total_ecoul, 0, sizeof(rbmd::Real)));

  //
  auto num_atoms = *(_structure_info_data->_num_atoms);
  op::SpecialLJCutCoulForceOp<device::DEVICE_GPU> lj_cut_coul_force_op;
  lj_cut_coul_force_op(_device_data->_d_box,_device_data->_d_erf_table, _cut_off, num_atoms,_alpha,_qqr2e,
                  thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
                  thrust::raw_pointer_cast(_device_data->_d_atoms_id.data()),
                  thrust::raw_pointer_cast(_device_data->_d_sigma.data()),
                  thrust::raw_pointer_cast(_device_data->_d_eps.data()),
                  thrust::raw_pointer_cast(_list->_start_idx.data()),
                  thrust::raw_pointer_cast(_list->_end_idx.data()),
                  thrust::raw_pointer_cast(_list->_d_neighbors.data()),
                  thrust::raw_pointer_cast(_device_data->_d_special_ids.data()),
                  thrust::raw_pointer_cast(_device_data->_d_special_weights.data()),
                  thrust::raw_pointer_cast(_device_data->_d_special_offsets.data()),
                  thrust::raw_pointer_cast(_device_data->_d_special_count.data()),
                  thrust::raw_pointer_cast(_device_data->_d_charge.data()),
                  thrust::raw_pointer_cast(_device_data->_d_px.data()),
                  thrust::raw_pointer_cast(_device_data->_d_py.data()),
                  thrust::raw_pointer_cast(_device_data->_d_pz.data()),
                  thrust::raw_pointer_cast(_device_data->_d_force_ljcoul_x.data()),
                  thrust::raw_pointer_cast(_device_data->_d_force_ljcoul_y.data()),
                  thrust::raw_pointer_cast(_device_data->_d_force_ljcoul_z.data()),
                  _d_total_evdwl,_d_total_ecoul);

  CHECK_RUNTIME(MEMCPY(&h_total_evdwl,_d_total_evdwl , sizeof(rbmd::Real), D2H));
  CHECK_RUNTIME(MEMCPY(&h_total_ecoul,_d_total_ecoul , sizeof(rbmd::Real), D2H));

  // 打印累加后的总能量
  _ave_evdwl = h_total_evdwl/num_atoms;
  _ave_ecoul = h_total_ecoul/num_atoms;


  std::cout << "test_current_step:" << test_current_step <<  " ,"
  << "average_energy_vdwl:" << _ave_evdwl << " " << "average_coul_energy_old:" <<
    _ave_ecoul  << std::endl;

  //out
  std::ofstream outfile("ave_energy_lj.txt", std::ios::app);
  outfile << test_current_step << " " << _ave_evdwl << std::endl;
  outfile.close();


    // //out
    // std::vector<rbmd::Real> h_force_ljcoul_x(num_atoms);
    // std::vector<rbmd::Real> h_force_ljcoul_y(num_atoms);
    // std::vector<rbmd::Real> h_force_ljcoul_z(num_atoms);
    // thrust::copy(_device_data->_d_force_ljcoul_x.begin(),
    //   _device_data->_d_force_ljcoul_x.end(), h_force_ljcoul_x.begin());
    // thrust::copy(_device_data->_d_force_ljcoul_y.begin(),
    // _device_data->_d_force_ljcoul_y.end(), h_force_ljcoul_y.begin());
    //
    // thrust::copy(_device_data->_d_force_ljcoul_z.begin(),
    // _device_data->_d_force_ljcoul_z.end(), h_force_ljcoul_z.begin());
    //
    // thrust::host_vector<rbmd::Real> h_atoms_id = _device_data->_d_atoms_id;
    // std::ofstream output_file("output_force_special_lj_verlet.txt");
    // for (size_t i = 0; i < h_force_ljcoul_x.size(); ++i)
    // {
    //   output_file << h_atoms_id[i] << " " << h_force_ljcoul_x[i] << " "<<h_force_ljcoul_y[i] <<" "
    //   << h_force_ljcoul_z[i] << std::endl;
    // }
    // output_file.close();
}

void CVFF::ComputeSpecialCoulForce()
{
  rbmd::Real h_total_e_specialcoul = 0.0;
  CHECK_RUNTIME(MEMSET(_d_total_e_specialcoul, 0, sizeof(rbmd::Real)));

  auto _atom_id_to_idx =
    LinkedCellLocator::GetInstance().GetLinkedCell()->_atom_id_to_idx;
  //
  auto num_atoms = *(_structure_info_data->_num_atoms);
  op::ComputeSpecialCoulForceOp<device::DEVICE_GPU> special_coul_force_op;
  special_coul_force_op(_device_data->_d_box,num_atoms,_qqr2e,
  thrust::raw_pointer_cast(_device_data->_d_atoms_id.data()),
  thrust::raw_pointer_cast(_atom_id_to_idx.data()),
    thrust::raw_pointer_cast(_device_data->_d_atoms_vec.data()),
    thrust::raw_pointer_cast(_device_data->_d_atoms_offset.data()),
    thrust::raw_pointer_cast(_device_data->_d_atoms_count.data()),
    thrust::raw_pointer_cast(_device_data->_d_special_ids.data()),
    thrust::raw_pointer_cast(_device_data->_d_special_weights.data()),
    thrust::raw_pointer_cast(_device_data->_d_special_offsets.data()),
    thrust::raw_pointer_cast(_device_data->_d_special_count.data()),
    thrust::raw_pointer_cast(_device_data->_d_charge.data()),
    thrust::raw_pointer_cast(_device_data->_d_px.data()),
    thrust::raw_pointer_cast(_device_data->_d_py.data()),
    thrust::raw_pointer_cast(_device_data->_d_pz.data()),
    thrust::raw_pointer_cast(_device_data->_d_force_specialcoul_x.data()),
    thrust::raw_pointer_cast(_device_data->_d_force_specialcoul_y.data()),
    thrust::raw_pointer_cast(_device_data->_d_force_specialcoul_z.data()),
    _d_total_e_specialcoul);

  CHECK_RUNTIME(MEMCPY(&h_total_e_specialcoul,_d_total_e_specialcoul , sizeof(rbmd::Real), D2H));

  _ave_e_specialcoul = h_total_e_specialcoul/num_atoms;
  _ave_ecoul = _ave_ecoul -  _ave_e_specialcoul;

  // 打印累加后的总能量
  std::cout << "test_current_step:" << test_current_step <<  " ,"<<
    "average_energy_specialcoul:" << _ave_e_specialcoul << " ,"
  <<"average_energy_ecoul:" << _ave_ecoul << std::endl;


  //out
  std::ofstream outfile("ave_energy_special_coul.txt", std::ios::app);
  outfile << test_current_step << " "
  << _ave_e_specialcoul << " " << _ave_ecoul<<  std::endl;
  outfile.close();

  // //out
  // std::vector<rbmd::Real> h_specialcoulx(num_atoms);
  // std::vector<rbmd::Real> h_specialcouly(num_atoms);
  // std::vector<rbmd::Real> h_specialcoulz(num_atoms);
  //
  // thrust::copy(_device_data->_d_force_specialcoul_x.begin(),
  //   _device_data->_d_force_specialcoul_x.end(), h_specialcoulx.begin());
  // thrust::copy(_device_data->_d_force_specialcoul_y.begin(),
  // _device_data->_d_force_specialcoul_y.end(), h_specialcouly.begin());
  // thrust::copy(_device_data->_d_force_specialcoul_z.begin(),
  // _device_data->_d_force_specialcoul_z.end(), h_specialcoulz.begin());
  //
  // thrust::host_vector<rbmd::Real> h_atoms_id = _device_data->_d_atoms_id;
  // std::ofstream output_file1("output_force_special_coul.txt");
  // for (size_t i = 0; i < h_specialcoulx.size(); ++i)
  // {
  //   output_file1  << h_atoms_id[i] << " "
  //   << h_specialcoulx[i] << " " << h_specialcouly[i]  << " " << h_specialcoulz[i]
  //   << std::endl;
  // }
  // output_file1.close();

}

void CVFF::ComputeKspaceForce()
{
  if("RBE" == _coulomb_type)
  {
      ComputeRBE();
  }
  else
  {
     ComputeEwlad();
  }
}

void CVFF::SumForces()
{
  TransformForces(_device_data->_d_fx,_device_data->_d_force_ljcoul_x,
                  _device_data->_d_force_specialcoul_x,_device_data->_d_force_ewald_x,
                  _device_data->_d_force_bond_x,_device_data->_d_force_angle_x,
                  _device_data->_d_force_dihedral_x);

  TransformForces(_device_data->_d_fy,_device_data->_d_force_ljcoul_y,
                _device_data->_d_force_specialcoul_y,_device_data->_d_force_ewald_y,
                _device_data->_d_force_bond_y,_device_data->_d_force_angle_y,
                _device_data->_d_force_dihedral_y);

  TransformForces(_device_data->_d_fz,_device_data->_d_force_ljcoul_z,
                _device_data->_d_force_specialcoul_z,_device_data->_d_force_ewald_z,
                _device_data->_d_force_bond_z,_device_data->_d_force_angle_z,
                _device_data->_d_force_dihedral_z);
}

void CVFF::ComputeChargeStructureFactorEwald(
    Box* box,
    rbmd::Id num_atoms,
    rbmd::Id Kmax,
    rbmd::Real alpha,
    rbmd::Real qqr2e,
    rbmd::Real* value_Re_array,
    rbmd::Real* value_Im_array)
{
    //thrust::fill(density_real.begin(), density_real.end(), 0.0f);
    //thrust::fill(density_imag.begin(), density_imag.end(), 0.0f);
    thrust::device_vector<rbmd::Real> density_real_atom;
    thrust::device_vector<rbmd::Real> density_imag_atom;
    density_real_atom.resize(num_atoms);
    density_imag_atom.resize(num_atoms);

    rbmd::Real total_energy_kspace= 0;
    rbmd::Id index = 0;
    for (rbmd::Id i = -Kmax; i <= Kmax; i++)
    {
        for (rbmd::Id j = -Kmax; j <= Kmax; j++)
        {
            for (rbmd::Id k = -Kmax; k <= Kmax; k++)
            {
                if (!(i == 0 && j == 0 && k == 0))
                {
                    Real3 K = make_Real3(rbmd::Real(2 * M_PI * i / box->_length[0]),
                                           rbmd::Real(2 * M_PI * j / box->_length[1]),
                                           rbmd::Real(2 * M_PI * k / box->_length[2]));
                    rbmd::Real Range_K = SQRT(K.x * K.x + K.y * K.y + K.z * K.z);
                    rbmd::Real Range_K2 = Range_K*Range_K;
                    rbmd::Real alpha_inv =  1 / alpha;

                    op::ComputeChargeStructureFactorOp<device::DEVICE_GPU> charge_structure_factor_op;
                    charge_structure_factor_op(num_atoms, K,
                        thrust::raw_pointer_cast(_device_data->_d_charge.data()),
                        thrust::raw_pointer_cast(_device_data->_d_px.data()),
                        thrust::raw_pointer_cast(_device_data->_d_py.data()),
                        thrust::raw_pointer_cast(_device_data->_d_pz.data()),
                        thrust::raw_pointer_cast(density_real_atom.data()),
                        thrust::raw_pointer_cast(density_imag_atom.data()));

                    rbmd::Real value_Re = thrust::reduce(density_real_atom.begin(), density_real_atom.end(), 0.0f, thrust::plus<rbmd::Real>());
                    rbmd::Real value_Im = thrust::reduce(density_imag_atom.begin(), density_imag_atom.end(), 0.0f, thrust::plus<rbmd::Real>());
                    rbmd::Real Range_density2 = POW(value_Re, 2.0) + POW(value_Im, 2.0);

                    total_energy_kspace +=
                      EXP(-0.25 * Range_K2 * alpha_inv) * Range_density2 / Range_K2;

                    value_Re_array[index] = value_Re;
                    value_Im_array[index] = value_Im;
                    index++;
                }
            }
        }
    }

  //energy

  //charge self energy//
  ComputeSelfEnergy(alpha,qqr2e,_ave_self_energy);

  //compute Kspace energy//
  rbmd::Real volume = box->_length[0] * box->_length[1]*box->_length[2];
  total_energy_kspace = qqr2e * (2 * M_PI / volume) * total_energy_kspace;
  _ave_ekspace = total_energy_kspace / num_atoms;

  _ave_ekspace = _ave_ekspace + _ave_self_energy;

  //out
   std::cout << "test_current_step:" << test_current_step <<  " ,"
   << "average_energy_ewald:" << _ave_ekspace << std::endl;

  std::ofstream outfile("ave_energy_ewald.txt", std::ios::app);
  outfile << test_current_step << " "<< _ave_ekspace << std::endl;
  outfile.close();
}

void CVFF::ComputeEwlad()
{
    auto num_atoms = *(_structure_info_data->_num_atoms);
    rbmd::Real* value_Re_array;
    rbmd::Real* value_Im_array;

    //EwaldForce//
    CHECK_RUNTIME(MALLOC(&value_Re_array, _num_k * sizeof(rbmd::Real)));
    CHECK_RUNTIME(MALLOC(&value_Im_array, _num_k * sizeof(rbmd::Real)));
    CHECK_RUNTIME(MEMSET(value_Re_array, 0, _num_k *sizeof(rbmd::Real)));
    CHECK_RUNTIME(MEMSET(value_Im_array, 0, _num_k *sizeof(rbmd::Real)));

    ComputeChargeStructureFactorEwald(_device_data->_d_box, num_atoms, _Kmax,
      _alpha,_qqr2e, value_Re_array, value_Im_array);


    op::ComputeEwaldForceOp<device::DEVICE_GPU> ewlad_force_op;
    ewlad_force_op(_device_data->_d_box,num_atoms, _Kmax, _alpha,_qqr2e,
        value_Re_array,value_Im_array,
        thrust::raw_pointer_cast(_device_data->_d_charge.data()),
        thrust::raw_pointer_cast(_device_data->_d_px.data()),
        thrust::raw_pointer_cast(_device_data->_d_py.data()),
        thrust::raw_pointer_cast(_device_data->_d_pz.data()),
        thrust::raw_pointer_cast(_device_data->_d_force_ewald_x.data()),
        thrust::raw_pointer_cast(_device_data->_d_force_ewald_y.data()),
        thrust::raw_pointer_cast(_device_data->_d_force_ewald_z.data()));

    CHECK_RUNTIME(FREE(value_Re_array));
    CHECK_RUNTIME(FREE(value_Im_array));

  // //out
  // std::vector<rbmd::Real> h_force_ewald_x(num_atoms);
  // std::vector<rbmd::Real> h_force_ewald_y(num_atoms);
  // std::vector<rbmd::Real> h_force_ewald_z(num_atoms);
  // thrust::copy(_device_data->_d_force_ewald_x.begin(),
  //   _device_data->_d_force_ewald_x.end(), h_force_ewald_x.begin());
  // thrust::copy(_device_data->_d_force_ewald_y.begin(),
  // _device_data->_d_force_ewald_y.end(), h_force_ewald_y.begin());
  //
  // thrust::copy(_device_data->_d_force_ewald_z.begin(),
  // _device_data->_d_force_ewald_z.end(), h_force_ewald_z.begin());
  //
  // thrust::host_vector<rbmd::Real> h_atoms_id = _device_data->_d_atoms_id;
  // std::ofstream output_file("output_force_ewald.txt");
  // for (size_t i = 0; i < h_force_ewald_x.size(); ++i)
  // {
  //   output_file << h_atoms_id[i]<< " " << h_force_ewald_x[i] << " "<<h_force_ewald_y[i] <<" "
  //   << h_force_ewald_z[i] << std::endl;
  // }
  // output_file.close();
}

void CVFF::ERFInit()
{
  _device_data->_d_erf_table->init();
}

void CVFF::RBEInit(Box* box,rbmd::Real alpha,rbmd::Id RBE_P)
{
  Real3 sigma = { rbmd::Real((SQRT(alpha / 2.0) * box->_length[0]/M_PI)),
                  rbmd::Real((SQRT(alpha / 2.0) * box->_length[1]/M_PI)),
                  rbmd::Real((SQRT(alpha / 2.0) * box->_length[2]/M_PI))};
  auto random = true;
  RBEPSAMPLE rbe_presolve_psample = { alpha, box, RBE_P, random};

  _P_Sample_x.resize(RBE_P);
  _P_Sample_y.resize(RBE_P);
  _P_Sample_z.resize(RBE_P);

  rbe_presolve_psample.Fetch_P_Sample(0.0, sigma,
    thrust::raw_pointer_cast(_P_Sample_x.data()),
    raw_pointer_cast(_P_Sample_y.data()),
    raw_pointer_cast(_P_Sample_z.data()));

  //index key
  auto num_atoms = *(_structure_info_data->_num_atoms);
  _psample_key.resize(num_atoms * RBE_P);
  op::GenerateIndexArrayOp<device::DEVICE_GPU> generate_index_array_op;
  generate_index_array_op(num_atoms,RBE_P,
    thrust::raw_pointer_cast(_psample_key.data()));
}

void CVFF::ComputeChargeStructureFactorRBE(
   Box* box,
   rbmd::Id num_atoms,
   rbmd::Id Kmax,
   rbmd::Real alpha,
   rbmd::Id RBE_P,
   rbmd::Real qqr2e,
   thrust::device_vector<rbmd::Real> rhok_real_redue,
   thrust::device_vector<rbmd::Real> rhok_image_redue)
{
  //
  thrust::device_vector<rbmd::Real>  rhok_real_atom;
  thrust::device_vector<rbmd::Real>  rhok_image_atom;

  rhok_real_atom.resize(num_atoms* RBE_P);
  rhok_image_atom.resize(num_atoms* RBE_P);
  auto p_number= RBE_P;

  //Charge Structure Factor
  op::ComputePnumberChargeStructureFactorOp<device::DEVICE_GPU> pnumber_charge_structure_factor_op;
  pnumber_charge_structure_factor_op(_device_data->_d_box, num_atoms, p_number,
      thrust::raw_pointer_cast(_device_data->_d_charge.data()),
      thrust::raw_pointer_cast(_P_Sample_x.data()),
      raw_pointer_cast(_P_Sample_y.data()),
      raw_pointer_cast(_P_Sample_z.data()),
      thrust::raw_pointer_cast(_device_data->_d_px.data()),
      thrust::raw_pointer_cast(_device_data->_d_py.data()),
      thrust::raw_pointer_cast(_device_data->_d_pz.data()),
      thrust::raw_pointer_cast(rhok_real_atom.data()),
      thrust::raw_pointer_cast(rhok_image_atom.data()));

  // reduce_by_key for rhok
  thrust::device_vector<rbmd::Id>  psamplekey_out;
  psamplekey_out.resize(RBE_P);

   reduce_by_key(_psample_key.begin(), _psample_key.end(),
    rhok_real_atom.begin(),psamplekey_out.begin(), rhok_real_redue.begin(),
    thrust::equal_to<rbmd::Id>(),thrust::plus<rbmd::Real>());

   reduce_by_key(_psample_key.begin(), _psample_key.end(),
    rhok_image_atom.begin(),psamplekey_out.begin(), rhok_image_redue.begin(),
    thrust::equal_to<rbmd::Id>(),thrust::plus<rbmd::Real>());

  //energy

  //charge self energy//
  ComputeSelfEnergy(alpha,qqr2e,_ave_self_energy);

  //kspace energy
  ComputeKspaceEnergy(_device_data->_d_box, num_atoms, Kmax,
      alpha, qqr2e ,_ave_ekspace);
  _ave_ekspace = _ave_ekspace +_ave_self_energy;

    //out
   std::cout << "test_current_step:" << test_current_step <<  " ,"
   << "average_energy_rbe:" << _ave_ekspace << std::endl;

  std::ofstream outfile("ave_energy_rbe.txt", std::ios::app);
  outfile << test_current_step << " "<< _ave_ekspace << std::endl;
  outfile.close();
}

void CVFF::ComputeRBE()
{
  //
  auto num_atoms = *(_structure_info_data->_num_atoms);
  _rhok_real_redue.resize(_RBE_P);
  _rhok_image_redue.resize(_RBE_P);

  ComputeChargeStructureFactorRBE(_device_data->_d_box, num_atoms, _Kmax,
      _alpha,_RBE_P,_qqr2e,_rhok_real_redue,_rhok_image_redue);

   //RBE Force
  op::ComputeRBEForceOp<device::DEVICE_GPU> rbe_force_op;
  rbe_force_op(_device_data->_d_box,num_atoms, _RBE_P,_alpha,_qqr2e,
        thrust::raw_pointer_cast(_rhok_real_redue.data()),
        thrust::raw_pointer_cast(_rhok_image_redue.data()),
        thrust::raw_pointer_cast(_device_data->_d_charge.data()),
        thrust::raw_pointer_cast(_P_Sample_x.data()),
        raw_pointer_cast(_P_Sample_y.data()),
        raw_pointer_cast(_P_Sample_z.data()),
        thrust::raw_pointer_cast(_device_data->_d_px.data()),
        thrust::raw_pointer_cast(_device_data->_d_py.data()),
        thrust::raw_pointer_cast(_device_data->_d_pz.data()),
        thrust::raw_pointer_cast(_device_data->_d_force_ewald_x.data()),
        thrust::raw_pointer_cast(_device_data->_d_force_ewald_y.data()),
        thrust::raw_pointer_cast(_device_data->_d_force_ewald_z.data()));

  // //out
  // std::vector<rbmd::Real> h_force_ewald_x(num_atoms);
  // std::vector<rbmd::Real> h_force_ewald_y(num_atoms);
  // std::vector<rbmd::Real> h_force_ewald_z(num_atoms);
  // thrust::copy(_device_data->_d_force_ewald_x.begin(),
  //   _device_data->_d_force_ewald_x.end(), h_force_ewald_x.begin());
  // thrust::copy(_device_data->_d_force_ewald_y.begin(),
  // _device_data->_d_force_ewald_y.end(), h_force_ewald_y.begin());
  //
  // thrust::copy(_device_data->_d_force_ewald_z.begin(),
  // _device_data->_d_force_ewald_z.end(), h_force_ewald_z.begin());
  //
  // thrust::host_vector<rbmd::Real> h_atoms_id = _device_data->_d_atoms_id;
  // std::ofstream output_file("output_force_rbe.txt");
  // for (size_t i = 0; i < h_force_ewald_x.size(); ++i)
  // {
  //   output_file  << h_atoms_id[i] << " " << h_force_ewald_x[i] << " "<<h_force_ewald_y[i] <<" "
  //   << h_force_ewald_z[i] << std::endl;
  // }
  // output_file.close();
}

void CVFF::ComputeLJCoulEnergy()
{
  // energy
  //neighbor_list_build
  auto start = std::chrono::high_resolution_clock::now();
  _list = _neighbor_list_builder->Build();

  auto end = std::chrono::high_resolution_clock::now();

  std::chrono::duration<rbmd::Real> duration = end - start;
  std::cout << "后处理---构建verlet-list耗时---" << duration.count() << "秒" << std::endl;


  rbmd::Real h_total_evdwl = 0.0;
  rbmd::Real h_total_ecoul = 0.0;

  CHECK_RUNTIME(MEMSET(_d_total_evdwl, 0, sizeof(rbmd::Real)));
  CHECK_RUNTIME(MEMSET(_d_total_ecoul, 0, sizeof(rbmd::Real)));

  auto num_atoms = *(_structure_info_data->_num_atoms);
  op::SpeciaLJCutCoulEnergyOp<device::DEVICE_GPU> lj_cut_coul_energy_op;
  lj_cut_coul_energy_op(_device_data->_d_box,_device_data->_d_erf_table,_cut_off, num_atoms,_alpha,_qqr2e,
                thrust::raw_pointer_cast(_device_data->_d_atoms_type.data()),
                thrust::raw_pointer_cast(_device_data->_d_atoms_id.data()),
                thrust::raw_pointer_cast(_device_data->_d_sigma.data()),
                thrust::raw_pointer_cast(_device_data->_d_eps.data()),
                thrust::raw_pointer_cast(_list->_start_idx.data()),
                thrust::raw_pointer_cast(_list->_end_idx.data()),
                thrust::raw_pointer_cast(_list->_d_neighbors.data()),
                thrust::raw_pointer_cast(_device_data->_d_special_ids.data()),
                thrust::raw_pointer_cast(_device_data->_d_special_weights.data()),
                thrust::raw_pointer_cast(_device_data->_d_special_offsets.data()),
                thrust::raw_pointer_cast(_device_data->_d_special_count.data()),
                thrust::raw_pointer_cast(_device_data->_d_charge.data()),
                thrust::raw_pointer_cast(_device_data->_d_px.data()),
                thrust::raw_pointer_cast(_device_data->_d_py.data()),
                thrust::raw_pointer_cast(_device_data->_d_pz.data()),
                _d_total_evdwl,_d_total_ecoul);
  CHECK_RUNTIME(MEMCPY(&h_total_evdwl,_d_total_evdwl , sizeof(rbmd::Real), D2H));
  CHECK_RUNTIME(MEMCPY(&h_total_ecoul,_d_total_ecoul , sizeof(rbmd::Real), D2H));

  // 打印累加后的总能量
  _ave_evdwl = h_total_evdwl/num_atoms;
  _ave_ecoul = h_total_ecoul/num_atoms;

  std::cout << "test_current_step:" << test_current_step <<  " ,"
  << "average_energy_vdwl:" << _ave_evdwl  << std::endl;

  //out
  std::ofstream outfile("ave_energy_specia_lj_rbl.txt", std::ios::app);
  outfile << test_current_step << " " << _ave_evdwl << std::endl;
  outfile.close();
}

void CVFF::ComputeSelfEnergy(
  rbmd::Real  alpha,
  rbmd::Real  qqr2e,
  rbmd::Real& ave_self_energy)
{
  //compute self_energy
  auto num_atoms = *(_structure_info_data->_num_atoms);
  thrust::device_vector<rbmd::Real> sq_charge;
  sq_charge.resize(num_atoms);

  op::SqchargeOp<device::DEVICE_GPU> sq_charge_op;
  sq_charge_op(num_atoms,
    thrust::raw_pointer_cast(_device_data->_d_charge.data()),
    thrust::raw_pointer_cast(sq_charge.data()));

  rbmd::Real sum_sq_charge = thrust::reduce(sq_charge.begin(), sq_charge.end(), 0.0f, thrust::plus<rbmd::Real>());
  rbmd::Real total_self_energy = qqr2e * (- sqrt(alpha / M_PI) *sum_sq_charge);

  ave_self_energy =  total_self_energy / num_atoms;
}

void CVFF::ComputeKspaceEnergy(
    Box* box,
    rbmd::Id num_atoms,
    rbmd::Id Kmax,
    rbmd::Real alpha,
    rbmd::Real qqr2e,
    rbmd::Real&  ave_ekspace)
{
    thrust::device_vector<rbmd::Real> density_real_atom;
    thrust::device_vector<rbmd::Real> density_imag_atom;
    density_real_atom.resize(num_atoms);
    density_imag_atom.resize(num_atoms);

    rbmd::Real total_energy_ewald = 0;
    for (rbmd::Id i = -Kmax; i <= Kmax; i++)
    {
        for (rbmd::Id j = -Kmax; j <= Kmax; j++)
        {
            for (rbmd::Id k = -Kmax; k <= Kmax; k++)
            {
                if (!(i == 0 && j == 0 && k == 0))
                {
                    Real3 K = make_Real3(rbmd::Real(2 * M_PI * i / box->_length[0]),
                                           rbmd::Real(2 * M_PI * j / box->_length[1]),
                                           rbmd::Real(2 * M_PI * k / box->_length[2]));
                    rbmd::Real Range_K = SQRT(K.x * K.x + K.y * K.y + K.z * K.z);
                    rbmd::Real Range_K2 = Range_K*Range_K;
                    rbmd::Real alpha_inv =  1 / alpha;

                    op::ComputeChargeStructureFactorOp<device::DEVICE_GPU> charge_structure_factor_op;
                    charge_structure_factor_op(num_atoms, K,
                        thrust::raw_pointer_cast(_device_data->_d_charge.data()),
                        thrust::raw_pointer_cast(_device_data->_d_px.data()),
                        thrust::raw_pointer_cast(_device_data->_d_py.data()),
                        thrust::raw_pointer_cast(_device_data->_d_pz.data()),
                        thrust::raw_pointer_cast(density_real_atom.data()),
                        thrust::raw_pointer_cast(density_imag_atom.data()));

                    rbmd::Real value_Re = thrust::reduce(density_real_atom.begin(), density_real_atom.end(), 0.0f, thrust::plus<rbmd::Real>());
                    rbmd::Real value_Im = thrust::reduce(density_imag_atom.begin(), density_imag_atom.end(), 0.0f, thrust::plus<rbmd::Real>());
                    rbmd::Real Range_density2 = POW(value_Re, 2.0) + POW(value_Im, 2.0);

                    total_energy_ewald +=
                      EXP(-0.25 * Range_K2 * alpha_inv) * Range_density2 / Range_K2;
                }
            }
        }
    }

  rbmd::Real volume = box->_length[0] * box->_length[1]*box->_length[2];
  total_energy_ewald = qqr2e * (2 * M_PI / volume) * total_energy_ewald;
  ave_ekspace = total_energy_ewald / num_atoms;
}

void CVFF::ComputeBondForce()
{
  auto _atom_id_to_idx =
    LinkedCellLocator::GetInstance().GetLinkedCell()->_atom_id_to_idx;

  rbmd::Real h_energy_bond = 0.0;
  CHECK_RUNTIME(MEMSET(_d_total_ebond, 0, sizeof(rbmd::Real)));

  thrust::fill(_device_data->_d_force_bond_x.begin(),
    _device_data->_d_force_bond_x.end(), 0.0f);
  thrust::fill(_device_data->_d_force_bond_y.begin(),
    _device_data->_d_force_bond_y.end(), 0.0f);
  thrust::fill(_device_data->_d_force_bond_z.begin(),
    _device_data->_d_force_bond_z.end(), 0.0f);

  auto num_bonds = *(_structure_info_data->_num_bonds);
  op::ComputeBondForceOp<device::DEVICE_GPU> bond_force_op;
  bond_force_op(_device_data->_d_box,num_bonds,thrust::raw_pointer_cast(_atom_id_to_idx.data()),
    thrust::raw_pointer_cast(_device_data->_d_bond_coeffs_k.data()),
    thrust::raw_pointer_cast(_device_data->_d_bond_coeffs_equilibrium.data()),
    thrust::raw_pointer_cast(_device_data->_d_bond_type.data()),
    thrust::raw_pointer_cast(_device_data->_d_bond_id0.data()),
    thrust::raw_pointer_cast(_device_data->_d_bond_id1.data()),
    thrust::raw_pointer_cast(_device_data->_d_px.data()),
    thrust::raw_pointer_cast(_device_data->_d_py.data()),
    thrust::raw_pointer_cast(_device_data->_d_pz.data()),
    thrust::raw_pointer_cast(_device_data->_d_force_bond_x.data()),
    thrust::raw_pointer_cast(_device_data->_d_force_bond_y.data()),
    thrust::raw_pointer_cast(_device_data->_d_force_bond_z.data()),
    thrust::raw_pointer_cast(_device_data->_d_temp_atom_ids.data()),
    _d_total_ebond);

  CHECK_RUNTIME(MEMCPY(&h_energy_bond,_d_total_ebond , sizeof(rbmd::Real), D2H));

  // 打印累加后的总能量
  _ave_ebond = h_energy_bond/num_bonds;

  std::cout << "test_current_step:" << test_current_step <<  " ,"
  << "average_energy_bond:" << _ave_ebond  << std::endl;

  //out
  std::ofstream outfile("ave_energy_bond.txt", std::ios::app);
  outfile << test_current_step << " " << _ave_ebond << std::endl;
  outfile.close();

  // //append force bond
  // std::vector<rbmd::Real> h_force_bondx(num_atoms);
  // std::vector<rbmd::Real> h_force_bondy(num_atoms);
  // std::vector<rbmd::Real> h_force_bondz(num_atoms);
  //
  // thrust::copy(_device_data->_d_force_bond_x.begin(),
  //   _device_data->_d_force_bond_x.end(), h_force_bondx.begin());
  // thrust::copy(_device_data->_d_force_bond_y.begin(),
  // _device_data->_d_force_bond_y.end(), h_force_bondy.begin());
  // thrust::copy(_device_data->_d_force_bond_z.begin(),
  // _device_data->_d_force_bond_z.end(), h_force_bondz.begin());
  //
  //  thrust::host_vector<rbmd::Real> h_atoms_id = _device_data->_d_atoms_id;
  // std::ofstream output_file("output_force_bond_atomadd.txt");
  // for (size_t i = 0; i < h_force_bondx.size(); ++i)
  // {
  //   output_file  << h_atoms_id[i] << " "<< h_force_bondx[i] << " "
  //   << h_force_bondy[i]  << " " << h_force_bondz[i] << std::endl;
  // }
  // output_file.close();

}

void CVFF::ComputeAngleForce()
{
  auto atom_id_to_idx =
    LinkedCellLocator::GetInstance().GetLinkedCell()->_atom_id_to_idx;

  rbmd::Real h_energy_bond = 0.0;
  CHECK_RUNTIME(MEMSET(_d_total_eangle, 0, sizeof(rbmd::Real)));

  // thrust::fill(_device_data->_d_force_angle_x.begin(),
  // _device_data->_d_force_angle_x.end(), 0.0f);
  // thrust::fill(_device_data->_d_force_angle_y.begin(),
  //   _device_data->_d_force_angle_y.end(), 0.0f);
  // thrust::fill(_device_data->_d_force_angle_z.begin(),
  //   _device_data->_d_force_angle_z.end(), 0.0f);

  auto num_angles = *(_structure_info_data->_num_angles);
  op::ComputeAngleForceOp<device::DEVICE_GPU> angle_force_op;
  angle_force_op(_device_data->_d_box,num_angles,
    thrust::raw_pointer_cast(atom_id_to_idx.data()),
    thrust::raw_pointer_cast(_device_data->_d_angle_coeffs_k.data()),
    thrust::raw_pointer_cast(_device_data->_d_angle_coeffs_equilibrium.data()),
    thrust::raw_pointer_cast(_device_data->_d_angle_type.data()),
    thrust::raw_pointer_cast(_device_data->_d_angle_id0.data()),
    thrust::raw_pointer_cast(_device_data->_d_angle_id1.data()),
    thrust::raw_pointer_cast(_device_data->_d_angle_id2.data()),
    thrust::raw_pointer_cast(_device_data->_d_px.data()),
    thrust::raw_pointer_cast(_device_data->_d_py.data()),
    thrust::raw_pointer_cast(_device_data->_d_pz.data()),
    thrust::raw_pointer_cast(_device_data->_d_force_angle_x.data()),
    thrust::raw_pointer_cast(_device_data->_d_force_angle_y.data()),
    thrust::raw_pointer_cast(_device_data->_d_force_angle_z.data()),
    _d_total_eangle);

  CHECK_RUNTIME(MEMCPY(&h_energy_bond,_d_total_eangle , sizeof(rbmd::Real), D2H));

  // 打印累加后的总能量
  _ave_eangle = h_energy_bond/num_angles;

  std::cout << "test_current_step:" << test_current_step <<  " ," <<
    "average_energy_angle:" << _ave_eangle << std::endl;

  //out
  std::ofstream outfile("ave_energy_angle.txt", std::ios::app);
  outfile << test_current_step << " " << _ave_eangle << std::endl;
  outfile.close();

  // //append force angle
  // std::vector<rbmd::Real> h_force_anglex(num_atoms);
  // std::vector<rbmd::Real> h_force_angley(num_atoms);
  // std::vector<rbmd::Real> h_force_anglez(num_atoms);
  //
  // thrust::copy(_device_data->_d_force_angle_x.begin(),
  //   _device_data->_d_force_angle_x.end(), h_force_anglex.begin());
  // thrust::copy(_device_data->_d_force_angle_y.begin(),
  // _device_data->_d_force_angle_y.end(), h_force_angley.begin());
  // thrust::copy(_device_data->_d_force_angle_z.begin(),
  // _device_data->_d_force_angle_z.end(), h_force_anglez.begin());
  //
  // thrust::host_vector<rbmd::Real> h_atoms_id = _device_data->_d_atoms_id;
  // std::ofstream output_file("output_force_angle.txt");
  // for (size_t i = 0; i < h_force_anglex.size(); ++i)
  // {
  //   output_file << h_atoms_id[i]<< " "
  //   << h_force_anglex[i] << " " << h_force_angley[i]  << " " << h_force_anglez[i]
  //   << std::endl;
  // }
  // output_file.close();

}

void CVFF::ComputeDihedralForce()
{
  thrust::fill(_device_data->_d_force_dihedral_x.begin(),
    _device_data->_d_force_dihedral_x.end(), 0.0f);
  thrust::fill(_device_data->_d_force_dihedral_y.begin(),
    _device_data->_d_force_dihedral_y.end(), 0.0f);
  thrust::fill(_device_data->_d_force_dihedral_z.begin(),
    _device_data->_d_force_dihedral_z.end(), 0.0f);

  //thrust::device_vector<int4> dihedral_list;
  auto atom_id_to_idx =
    LinkedCellLocator::GetInstance().GetLinkedCell()->_atom_id_to_idx;

  rbmd::Real h_energy_dihedral= 0.0;
  CHECK_RUNTIME(MEMSET(_d_total_edihedral, 0, sizeof(rbmd::Real)));

  auto num_dihedrals = *(_structure_info_data->_num_dihedrals);
  op::ComputeDihedralForceOp<device::DEVICE_GPU> dihedral_force_op;
  dihedral_force_op(_device_data->_d_box,num_dihedrals,
    thrust::raw_pointer_cast(atom_id_to_idx.data()),
    thrust::raw_pointer_cast(_device_data->_d_dihedral_coeffs_k.data()),
    thrust::raw_pointer_cast(_device_data->_d_dihedral_coeffs_sign.data()),
    thrust::raw_pointer_cast(_device_data->_d_dihedral_coeffs_multiplicity.data()),
    thrust::raw_pointer_cast(_device_data->_d_dihedral_type.data()),
    thrust::raw_pointer_cast(_device_data->_d_dihedral_id0.data()),
    thrust::raw_pointer_cast(_device_data->_d_dihedral_id1.data()),
    thrust::raw_pointer_cast(_device_data->_d_dihedral_id2.data()),
    thrust::raw_pointer_cast(_device_data->_d_dihedral_id3.data()),
    thrust::raw_pointer_cast(_device_data->_d_px.data()),
    thrust::raw_pointer_cast(_device_data->_d_py.data()),
    thrust::raw_pointer_cast(_device_data->_d_pz.data()),
    thrust::raw_pointer_cast(_device_data->_d_force_dihedral_x.data()),
    thrust::raw_pointer_cast(_device_data->_d_force_dihedral_y.data()),
    thrust::raw_pointer_cast(_device_data->_d_force_dihedral_z.data()),
    _d_total_edihedral);

  CHECK_RUNTIME(MEMCPY(&h_energy_dihedral,_d_total_edihedral , sizeof(rbmd::Real), D2H));

  // 打印累加后的总能量
  _ave_dihedral = h_energy_dihedral/num_dihedrals;

  std::cout << "test_current_step:" << test_current_step <<  " ,"
   << "average_dihedral_energy:" << _ave_dihedral << std::endl;

  //out
  std::ofstream outfile("ave_energy_dihedral.txt", std::ios::app);
  outfile << test_current_step << " " << _ave_dihedral << std::endl;
  outfile.close();


  //append force dihedral
  auto num_atoms = *(_structure_info_data->_num_atoms);
  std::vector<rbmd::Real> h_force_dihedralx(num_atoms);
  std::vector<rbmd::Real> h_force_dihedraly(num_atoms);
  std::vector<rbmd::Real> h_force_dihedralz(num_atoms);

  thrust::copy(_device_data->_d_force_dihedral_x.begin(),
    _device_data->_d_force_dihedral_x.end(), h_force_dihedralx.begin());
  thrust::copy(_device_data->_d_force_dihedral_y.begin(),
  _device_data->_d_force_dihedral_y.end(), h_force_dihedraly.begin());
  thrust::copy(_device_data->_d_force_dihedral_z.begin(),
  _device_data->_d_force_dihedral_z.end(), h_force_dihedralz.begin());

  thrust::host_vector<rbmd::Real> h_atoms_id = _device_data->_d_atoms_id;
  std::ofstream output_file("output_force_dihedral.txt");
  for (size_t i = 0; i < h_force_dihedralx.size(); ++i)
  {
    output_file << h_atoms_id[i]<< " "
    << h_force_dihedralx[i] << " " << h_force_dihedraly[i]  << " " << h_force_dihedralz[i]
    << std::endl;
  }
  output_file.close();

}

void CVFF::ReduceByKey()
{
  // 先对 temp_atom_ids 和 temp_forces_x 进行排序
  thrust::sort_by_key(_device_data->_d_temp_atom_ids.begin(),
    _device_data->_d_temp_atom_ids.end(), _device_data->_d_temp_forces_bondx.begin());
  thrust::sort_by_key(_device_data->_d_temp_atom_ids.begin(),
    _device_data->_d_temp_atom_ids.end(), _device_data->_d_temp_forces_bondy.begin());
  thrust::sort_by_key(_device_data->_d_temp_atom_ids.begin(),
    _device_data->_d_temp_atom_ids.end(), _device_data->_d_temp_forces_bondz.begin());

  auto num_atoms = *(_structure_info_data->_num_atoms);
  thrust::device_vector<rbmd::Id> reduce_atom_ids_x;
  thrust::device_vector<rbmd::Id> reduce_atom_ids_y;
  thrust::device_vector<rbmd::Id> reduce_atom_ids_z;
  reduce_atom_ids_x.resize(num_atoms);
  reduce_atom_ids_y.resize(num_atoms);
  reduce_atom_ids_z.resize(num_atoms);

    thrust::reduce_by_key(
      _device_data->_d_temp_atom_ids.begin(), _device_data->_d_temp_atom_ids.end(),
      _device_data->_d_temp_forces_bondx.begin(),   // 输入的 x 方向力
      reduce_atom_ids_x.begin(),   // 按键聚合后的原子id
      _device_data->_d_force_bond_x.begin()    // 聚合后的 x 方向力
    );
    thrust::reduce_by_key(
      _device_data->_d_temp_atom_ids.begin(), _device_data->_d_temp_atom_ids.end(),
      _device_data->_d_temp_forces_bondy.begin(),
      reduce_atom_ids_y.begin(),
      _device_data->_d_force_bond_y.begin()
    );
    thrust::reduce_by_key(
      _device_data->_d_temp_atom_ids.begin(), _device_data->_d_temp_atom_ids.end(),
      _device_data->_d_temp_forces_bondz.begin(),
      reduce_atom_ids_z.begin(),
      _device_data->_d_force_bond_z.begin()
     );
}




