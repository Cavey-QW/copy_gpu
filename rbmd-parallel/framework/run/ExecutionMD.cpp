#include "FieldName.h"
#include "ExecutionMD.h"
#include "run/worklet/MolecularWorklet.h"
#include "run/worklet/RunWorklet.h"
#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/ArrayHandleGroupVec.h>
#include "Executioner.h"
#include "ERFTable.h"
#include <vtkm/worklet/Keys.h>
#include <vtkm/cont/ArrayHandleConstant.h>
#include <vtkm/cont/DeviceAdapter.h>
#include <vtkm/cont/ArrayHandleSOA.h>
#include <vtkm/cont/EnvironmentTracker.h>

//#define NEIGHBOR_GHOST 
//#define INIT_NEIGHBOR_GHOST

struct ComputeNearForceVerlet_MPIWorklet : vtkm::worklet::WorkletMapField
{
  ComputeNearForceVerlet_MPIWorklet(const Real& cut_off)
    : _cut_off(cut_off)
  {
  }

  using ControlSignature = void(FieldIn atoms_id,
                                ExecObject locator,
                                ExecObject topology,
                                ExecObject force_function,
                                FieldIn group_j,
                                FieldIn num_j,
                                FieldIn position_group,
                                FieldInOut force_reduce,
                                FieldOut nearforce);
  using ExecutionSignature = void(_1, _2, _3, _4, _5, _6, _7, _8, _9);

  template<typename NeighbourGroupVecType, typename ForceReduce, typename PositionGroup>
  VTKM_EXEC void operator()(const Id atoms_id,
                            const ExecPointLocator& locator,
                            const ExecTopology& topology,
                            const ExecForceFunction& force_function,
                            const NeighbourGroupVecType& group_j,
                            const Id& num_j,
                            const PositionGroup& position_group,
                            ForceReduce& force_reduce,
                            Vec3f& nearforce) const
  {
    Vec3f LJ_force = { 0, 0, 0 };
    Vec3f ele_force = { 0, 0, 0 };

    if (atoms_id == -1)
      return;

    nearforce = { 0, 0, 0 };
    const auto& molecular_id_i = topology.GetMolecularId(atoms_id);
    const auto& pts_type_i = topology.GetAtomsType(atoms_id);
    auto eps_i = topology.GetEpsilon(pts_type_i);
    auto sigma_i = topology.GetSigma(pts_type_i);
    auto charge_pi = topology.GetCharge(atoms_id);

    auto function = [&](const Vec3f& p_i, const Vec3f& p_j, const Id& pts_id_j)
    {
      auto charge_pj = topology.GetCharge(pts_id_j);
      auto molecular_id_j = topology.GetMolecularId(pts_id_j);
      auto pts_type_j = topology.GetAtomsType(pts_id_j);
      auto eps_j = topology.GetEpsilon(pts_type_j);
      auto sigma_j = topology.GetSigma(pts_type_j);
      auto r_ij = p_j - p_i;

      auto f_ele_force = force_function.ComputeNearEnergyForce(p_i, p_j, charge_pi, charge_pj);
      ele_force += f_ele_force;

      force_reduce[pts_id_j] = f_ele_force;
      if (molecular_id_i == molecular_id_j)
        return;
      auto f_lj = force_function.ComputeLJForce(r_ij, eps_i, eps_j, sigma_i, sigma_j, _cut_off);
      LJ_force += f_lj;
      force_reduce[pts_id_j] += f_lj;
    };

    auto p_i = locator.GetPtsPosition(atoms_id);

    for (Id p = 0; p < num_j; p++)
    {
      auto idj = group_j[p];
      auto p_j = position_group[idj]; // locator.GetPtsPosition(idj) - coord_offset_j[p];
      function(p_i, p_j, idj);
    }
    nearforce = LJ_force + ele_force;
  }
  Real _cut_off;
};

struct ComputePoentialEnWorklet : vtkm::worklet::WorkletMapField
{
  ComputePoentialEnWorklet(const Real& cut_off)
    : _cut_off(cut_off)
  {
  }

  using ControlSignature = void(FieldIn atoms_id,
                                ExecObject locator,
                                ExecObject topology,
                                ExecObject force_function,
                                FieldIn group_j,
                                FieldIn num_j,
                                FieldIn coord_offset_j,
                                FieldIn position_group,
                                FieldInOut potential_en_group,
                                FieldOut Potential_ruduce,
                                FieldOut PotentialEnergy);
  using ExecutionSignature = void(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11);

  template<typename NeighbourGroupVecType,
           typename CoordOffsetj,
           typename ForceGroup,
           typename PotentialReduce,
           typename PositionGroup>
  VTKM_EXEC void operator()(const Id atoms_id,
                            const ExecPointLocator& locator,
                            const ExecTopology& topology,
                            const ExecForceFunction& force_function,
                            const NeighbourGroupVecType& group_j,
                            const Id& num_j,
                            const CoordOffsetj& coord_offset_j,
                            const PositionGroup& position_group,
                            ForceGroup& potential_en_group,
                            PotentialReduce& Potential_ruduce,
                            Real& PotentialEnergy) const
  {
    PotentialEnergy = 0;
    if (atoms_id == -1)
      return;
    Real PE_ij = 0;
    const auto& molecular_id_i = topology.GetMolecularId(atoms_id);
    const auto& pts_type_i = topology.GetAtomsType(atoms_id);
    auto eps_i = topology.GetEpsilon(pts_type_i);
    auto sigma_i = topology.GetSigma(pts_type_i);

    auto function = [&](const Vec3f& p_i, const Vec3f& p_j, const Id& pts_id_j)
    {
      auto pts_type_j = topology.GetAtomsType(pts_id_j);
      auto eps_j = topology.GetEpsilon(pts_type_j);
      auto sigma_j = topology.GetSigma(pts_type_j);
      auto r_ij = p_j - p_i;
      auto eps_ij = vtkm::Sqrt(eps_i * eps_j);
      auto sigma_ij = (sigma_i + sigma_j) / 2;

      // special lj part
      auto molecular_id_j = topology.GetMolecularId(pts_id_j);
      IdComponent force_factor_ij = (molecular_id_i == molecular_id_j) ? 0 : 1.0;

      auto PE_test =
        force_factor_ij * force_function.ComputePotentialEn(r_ij, eps_ij, sigma_ij, _cut_off);
      PE_ij += PE_test;

      Potential_ruduce[pts_id_j] = PE_test;
      potential_en_group[pts_id_j] = PE_test;
    };

    auto p_i = locator.GetPtsPosition(atoms_id);

    for (Id p = 0; p < num_j; p++)
    {
      auto idj = group_j[p];
      auto p_j = position_group[idj];
      //auto p_j = locator.GetPtsPosition(idj) - coord_offset_j[p];
      function(p_i, p_j, idj);
    }
    //potential_en_group[atoms_id] = PE_ij; 注意：这里要不要删除！
    PotentialEnergy = PE_ij;
  }
  Real _cut_off;
};

struct SetIndex : vtkm::worklet::WorkletMapField
{
  SetIndex(const Id& num)
    : _pos_num(num)
  {
  }
  using ControlSignature = void(FieldIn ids, FieldOut p);
  using ExecutionSignature = void(_1, _2);

  //template<typename CoordType>
  VTKM_EXEC void operator()(const Id& ids, Id& p) const { p = ids / _pos_num; }

  Id _pos_num;
};

struct ReduceWorklet : vtkm::worklet::WorkletReduceByKey
{
  using ControlSignature = void(KeysIn p, ValuesIn value, ReducedValuesOut reduce_value);

  using ExecutionSignature = _3(_2);
  using InputDomain = _1;
  template<typename ForceVecType>
  VTKM_EXEC typename ForceVecType::ComponentType operator()(const ForceVecType& value) const
  {
    typename ForceVecType::ComponentType sum = 0;
    for (vtkm::IdComponent index = 0; index < value.GetNumberOfComponents(); index++)
    {
      sum = sum + value[index];
    }
    return sum;
  }
};

struct SetValue : vtkm::worklet::WorkletMapField
{
  using ControlSignature = void(FieldInOut array_);
  using ExecutionSignature = void(_1);

  template<typename CoordType>
  VTKM_EXEC void operator()(CoordType& array_) const
  {
    array_ = 0;
  }
};

struct SetValue_New : vtkm::worklet::WorkletMapField
{
  SetValue_New(const Id& pos, const Id& num)
    : _pos_num(pos)
    , _pnumber(num)
  {
  }

  using ControlSignature = void(FieldIn ids, WholeArrayInOut array_);
  using ExecutionSignature = void(_1, _2);

  template<typename CoordType>
  VTKM_EXEC void operator()(const Id& id, CoordType& array_) const
  {
    for (int i = 0; i < _pnumber; ++i)
    {
      auto index = id + i * _pos_num;
      array_.Set(index, 0);
    }
  }

  Id _pos_num;
  Id _pnumber;
};

ExecutionMD::ExecutionMD(const Configuration& cfg)
  : Execution(cfg)
  //, rbmd_parallel(_para.GetParameter<std::array<double, 3>>(PARA_BOXMIN),
  //                _para.GetParameter<std::array<double, 3>>(PARA_BOXMAX),
  //                static_cast<double> (_para.GetParameter<Real>(PARA_CUTOFF)))
{
  num_cont = 0;
}

void ExecutionMD::Init()
{
  _position = _para.GetFieldAsArrayHandle<Vec3f>(field::position);
  InitPointLocator(); 
  SetForceFunction();
  SetTopology();

  if (_para.GetParameter<bool>(PARA_FAR_FORCE))
  {
    auto N = _position.GetNumberOfValues();
    _rhok_Re.AllocateAndFill(_RBE_P * N, 0.0);
    _rhok_Im.AllocateAndFill(_RBE_P * N, 0.0);

    vtkm::cont::ArrayHandleIndex indexArray(N * _RBE_P);
    vtkm::cont::Invoker{}(SetIndex{ N }, indexArray, _psamplekey);
  }
  _EleNearPairtimer_counting = 0.0;
}

void ExecutionMD::Execute() 
{
  _timer.Start();
  PreSolve();
  Solve();
  PostSolve();
}

void ExecutionMD::SetParameters()
{
  _para.SetParameter(PARA_ENSEMBLE, Get<std::string>("ensemble"));
  _para.SetParameter(PARA_TEMP_CTRL_TYPE, Get<std::string>("temp_ctrl_type"));
  _para.SetParameter(PARA_PRESS_CTRL_TYPE, Get<std::string>("press_ctrl_type"));
  _para.SetParameter(PARA_TIMESTEP, Get<Real>("timestep"));
  _para.SetParameter(PARA_NUM_STEPS, Get<Real>("num_steps"));
  _para.SetParameter(PARA_TEMPERATURE, GetVectorValue<Real>("temperature"));
  _para.SetParameter(PARA_PRESSURE, GetVectorValue<Real>("pressure"));
}

void ExecutionMD::InitialCondition()
{
  // 物理场初始化
  if (_init_way == "inbuild")
  {
    if (_para.GetParameter<bool>(PARA_FAR_FORCE))
    {
      _charge = _para.GetFieldAsArrayHandle<Real>(field::charge);
      auto n = _position.GetNumberOfValues();
      _charge.Allocate(n);
      _charge.Fill(-1.0, 0);
      _charge.Fill(1.0, n / 2);
      _topology.SetCharge(_charge);
    }
    else
    {
      _charge = _para.GetFieldAsArrayHandle<Real>(field::charge);
      auto n = _position.GetNumberOfValues();
      _charge.AllocateAndFill(n, 0);
      _topology.SetCharge(_charge);
    }
  }
  else if (_para.GetParameter<std::string>(PARA_FILE_TYPE) == "EAM")
  {
    _charge = _para.GetFieldAsArrayHandle<Real>(field::charge);
    auto n = _position.GetNumberOfValues();
    _charge.AllocateAndFill(n, 0);
    _topology.SetCharge(_charge);
  }
  else
  {
    _charge = _para.GetFieldAsArrayHandle<Real>(field::charge);
    _topology.SetCharge(_charge);
  }

  _velocity = _para.GetFieldAsArrayHandle<Vec3f>(field::velocity);
  _mass = _para.GetFieldAsArrayHandle<Real>(field::mass);

  _molecule_id = _para.GetFieldAsArrayHandle<Id>(field::molecule_id);
  _atoms_id = _para.GetFieldAsArrayHandle<Id>(field::atom_id);
}

void ExecutionMD::SetForceFunction() {}

void ExecutionMD::SetTopology() {}

void ExecutionMD::InitERF() 
{
  if (_para.GetParameter<bool>(PARA_FAR_FORCE)/*_use_erf == true*/)
  {     
     _static_table.SetTableIndex1(vtkm::cont::make_ArrayHandle(erf_index_leq2_v));
     _static_table.SetTableIndex2(vtkm::cont::make_ArrayHandle(erf_index_geq2_v));
     _static_table.SetTableRij1(vtkm::cont::make_ArrayHandle(erf_dis_leq2_v));
     _static_table.SetTableRij2(vtkm::cont::make_ArrayHandle(erf_dis_geq2_v));
     _static_table.SetTabledRij1(vtkm::cont::make_ArrayHandle(erf_dis_dif_leq2_v));
     _static_table.SetTabledRij2(vtkm::cont::make_ArrayHandle(erf_dis_dif_geq2_v));
     _static_table.SetTableFunctionRij1(vtkm::cont::make_ArrayHandle(erf_gnear_leq2_v));
     _static_table.SetTableFunctionRij2(vtkm::cont::make_ArrayHandle(erf_gnear_geq2_v));
     _static_table.SetTabledFunctionRij1(vtkm::cont::make_ArrayHandle(erf_gnear_der_leq2_v));
     _static_table.SetTabledFunctionRij2(vtkm::cont::make_ArrayHandle(erf_gnear_der_geq2_v));
  }
}

void ExecutionMD::UpdateVerletList() 
{
  auto cut_off = _para.GetParameter<Real>(PARA_CUTOFF);
  auto N = _position.GetNumberOfValues();
  auto verletlist_num = N * N;

  ArrayHandle<Id> id_verletlist;
  id_verletlist.Allocate(verletlist_num);

  std::vector<Id> offset_vec(N + 1);
  Id inc = 0;
  std::generate(
    offset_vec.begin(), offset_vec.end(), [&](void) -> Id { return (inc++) * N; });
  ArrayHandle<Id> temp_offset = vtkm::cont::make_ArrayHandle(offset_vec);

  auto id_verletlist_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(id_verletlist, temp_offset);

  ArrayHandle<Id> num_verletlist;
  ArrayHandle<Vec3f> offset_verletlist;
  offset_verletlist.Allocate(N * N);

  auto offset_verletlist_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(offset_verletlist, temp_offset);

  RunWorklet::ComputeNeighbours(
    cut_off, _atoms_id, _locator, id_verletlist_group, num_verletlist, offset_verletlist_group);

  _locator.SetVerletListInfo(num_verletlist, id_verletlist_group, offset_verletlist_group);
}

void ExecutionMD::ComputeCorrForce(vtkm::cont::ArrayHandle<Vec3f>& corr_force)
{
  auto rc = _para.GetParameter<Real>(PARA_CUTOFF);
  auto rs = _para.GetParameter<Real>(PARA_R_CORE);
  auto rho_system = _para.GetParameter<Real>(PARA_RHO);
  auto N = _position.GetNumberOfValues();
  auto box = _para.GetParameter<Vec3f>(PARA_BOX);

  vtkm::cont::ArrayHandle<vtkm::Id> num_verletlist;
  vtkm::cont::ArrayHandle<vtkm::Id> id_verletlist;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> offset_verletlist;

  Id rs_num = rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * rs * rs * rs) + 1;
  Id rc_num = rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * rc * rc * rc) + 1;
  Id rcs_num = rc_num - rs_num;
  Real random_num = _para.GetParameter<Real>(PARA_NEIGHBOR_SAMPLE_NUM);
  Real random_rate = random_num / (rcs_num);
  Id pice_num = std::ceil(1.0 / random_rate);
  Id cutoff_num = rs_num + random_num;
  Id verletlist_num = N * cutoff_num;

  num_verletlist.Allocate(N * 2);
  id_verletlist.Allocate(verletlist_num);
  offset_verletlist.Allocate(verletlist_num);

  std::vector<Id> stride_vec(N + 1);
  Id inc = 0;
  std::generate(stride_vec.begin(), stride_vec.end(), [&](void) -> Id { return (inc++) * cutoff_num; });
  vtkm::cont::ArrayHandle<vtkm::Id> stride_array = vtkm::cont::make_ArrayHandle(stride_vec);

  auto id_verletlist_group = vtkm::cont::make_ArrayHandleGroupVecVariable(id_verletlist, stride_array);
  auto offset_verletlist_group = vtkm::cont::make_ArrayHandleGroupVecVariable(offset_verletlist, stride_array);
  auto num_verletlist_group = vtkm::cont::make_ArrayHandleGroupVec<2>(num_verletlist);

  RunWorklet::ComputeRBLNeighboursOnce(rs_num,
                                          pice_num,
                                          _atoms_id,
                                          _locator,
                                          id_verletlist_group,
                                          num_verletlist_group,
                                          offset_verletlist_group);

  RunWorklet::NearForceRBLERF(rs_num,
                              pice_num,
                              _unit_factor._qqr2e,
                              box,
                              _atoms_id,
                              _locator,
                              _topology,
                              _force_function,
                              _static_table,
                              id_verletlist_group,
                              num_verletlist_group,
                              offset_verletlist_group,
                              corr_force);
}

std::vector<Vec2f> ExecutionMD::ComputeChargeStructureFactorRBE(Real& _Vlength, ArrayHandle<Vec3f>& _psample)
{
  std::vector<Vec2f> rhok;
  auto p_number = _psample.GetNumberOfValues();
  ArrayHandle<Real> density_real;
  ArrayHandle<Real> density_image;

  //auto sum = 0.0f;
  //vtkm::cont::Timer timer;
  //#pragma omp parallel for
  for (Id i = 0; i < p_number; i++)
  {
    auto kl = _psample.ReadPortal().Get(i);
    kl = 2 * vtkm::Pi() * kl / _Vlength;
    RunWorklet::ComputeChargeStructureFactorComponent(
      kl, _position, _charge, density_real, density_image);
    //timer.Start();
    Real value_Re =
      vtkm::cont::Algorithm::Reduce(density_real, vtkm::TypeTraits<Real>::ZeroInitialization());
    Real value_Im =
      vtkm::cont::Algorithm::Reduce(density_image, vtkm::TypeTraits<Real>::ZeroInitialization());
    //sum += timer.GetElapsedTime();
    rhok.push_back({ value_Re, value_Im });
  }
  //std::cout << "timer: " << sum << std::endl;
  return rhok;
}

void ExecutionMD::InitParameters() 
{
  _unit = _para.GetParameter<std::string>(PARA_UNIT);
  if (_unit == "REAL")
  {
    _unit_factor._kB = 1.9872067 * vtkm::Pow(10.0, -3);
    _unit_factor._fmt2v = 4.186 * vtkm::Pow(10.0, -4);
    _unit_factor._mvv2e = 1.0 / (4.186 * vtkm::Pow(10.0, -4));
    _unit_factor._qqr2e = 332.06371;
  }
  else if (_unit == "LJ")
  {
    _unit_factor._kB = 1.0;
    _unit_factor._fmt2v = 1.0;
    _unit_factor._mvv2e = 1.0;
    _unit_factor._qqr2e = 1.0;
  }
  else if (_unit == "METAL")
  {
    _unit_factor._kB = 8.617343e-5;
    _unit_factor._fmt2v = 1.0 / 1.0364269e-4;
    _unit_factor._mvv2e = 1.0364269e-4;
    _unit_factor._qqr2e = 14.399645;
  }
  _para.SetParameter(PARA_UNIT_FACTOR, _unit_factor);
}

void ExecutionMD::ComputeEwaldEleForce(IdComponent& Kmax, ArrayHandle<Vec3f>& Ewald_ele_force)
{
  auto Vlength = _para.GetParameter<Real>(PARA_VLENGTH);
  auto box = _para.GetParameter<Vec3f>(PARA_BOX);
  auto rhok = ComputeChargeStructureFactorEwald(box, Kmax);
  ArrayHandle<Vec2f> whole_rhok = vtkm::cont::make_ArrayHandle(rhok);
  RunWorklet::ComputeNewFarElectrostatics(
    Kmax, _atoms_id, whole_rhok, _force_function, _topology, _locator, Ewald_ele_force);
}

std::vector<Vec2f> ExecutionMD::ComputeChargeStructureFactorEwald( Vec3f& box, IdComponent& Kmax)
{
  std::vector<Vec2f> rhok;
  ArrayHandle<Real> density_real;
  ArrayHandle<Real> density_image;
  for (Id i = -Kmax; i <= Kmax; i++)
  {
    for (Id j = -Kmax; j <= Kmax; j++)
    {
      for (Id k = -Kmax; k <= Kmax; k++)
      {
        if (!(i == 0 && j == 0 && k == 0))
        {
          Vec3f K = { Real(i), Real(j), Real(k) };
          //
          K = { Real(2 * vtkm::Pi() * K[0] / box[0]),
                Real(2 * vtkm::Pi() * K[1] / box[1]),
                Real(2 * vtkm::Pi() * K[2] / box[2]) };
          RunWorklet::ComputeChargeStructureFactorComponent(
            K, _position, _charge, density_real, density_image);
          Real value_Re = vtkm::cont::Algorithm::Reduce(
            density_real, vtkm::TypeTraits<Real>::ZeroInitialization());
          Real value_Im = vtkm::cont::Algorithm::Reduce(
            density_image, vtkm::TypeTraits<Real>::ZeroInitialization());

          rhok.push_back({ value_Re, value_Im });
        }
      }
    }
  }
  return rhok;
}

void ExecutionMD::ComputeRBEEleForce(ArrayHandle<Vec3f>& psample,
                                  IdComponent& RBE_P,
                                  ArrayHandle<Vec3f>& RBE_ele_force)
{
 

  auto Vlength = _para.GetParameter<Real>(PARA_VLENGTH);
  auto box = _para.GetParameter<Vec3f>(PARA_BOX);

  ArrayHandle<Vec2f> new_whole_rhok;
  ComputeNewChargeStructureFactorRBE(box, psample, new_whole_rhok);
  
  RunWorklet::ComputeNewRBEForce(RBE_P,
                                    _atoms_id,
                                    psample,
                                    /*whole_rhok,*/ new_whole_rhok, _force_function,
                                    _topology,
                                    _locator,
                                    RBE_ele_force);
 

}


void ExecutionMD::ComputeNewChargeStructureFactorRBE(Vec3f& _box,
                                                     ArrayHandle<Vec3f>& _psample,
                                                     ArrayHandle<Vec2f>& new_rhok)
{

  auto N = _position.GetNumberOfValues();

  const Id p_number = _psample.GetNumberOfValues();



  RunWorklet::ComputePnumberChargeStructureFactor(_box,
                                                  p_number,
                                                  N,
                                                  _atoms_id,
                                                  _position,
                                                  _charge,
                                                  _psample,
                                                   _rhok_Re,
                                                   _rhok_Im);


  vtkm::cont::ArrayHandle<Id> psamplekey_out;
  vtkm::cont::ArrayHandle<Real> rhok_Re_reduce;
  vtkm::cont::ArrayHandle<Real> rhok_Im_reduce;


  vtkm::cont::Algorithm::ReduceByKey<Id, Real>(
    _psamplekey,
    _rhok_Re,
    psamplekey_out,
    rhok_Re_reduce,
    vtkm::Add());
  vtkm::cont::Algorithm::ReduceByKey<Id, Real>(
    _psamplekey,
    _rhok_Im,
    psamplekey_out,
    rhok_Im_reduce,
    vtkm::Add());



  RunWorklet::ChangePnumberChargeStructureFactor(rhok_Re_reduce, rhok_Im_reduce, new_rhok);
}

void ExecutionMD::ComputeRBLNearForce(ArrayHandle<Vec3f>& nearforce)
{
  auto N = _position.GetNumberOfValues();
  vtkm::cont::ArrayHandle<Vec3f> corr_force;
  corr_force.Allocate(N);

  auto rc = _para.GetParameter<Real>(PARA_CUTOFF);
  auto rs = _para.GetParameter<Real>(PARA_R_CORE);
  auto box = _para.GetParameter<Vec3f>(PARA_BOX);
  auto rho_system = _para.GetParameter<Real>(PARA_RHO);

  vtkm::cont::ArrayHandle<vtkm::Id> num_verletlist;
  vtkm::cont::ArrayHandle<vtkm::Id> id_verletlist;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> offset_verletlist;

  //0.05 is a temp parameter, to fit with the low density of system requirement
  //so as the RBL for LJ

  Real coeff_rcs = 1.0 + (0.05 / rho_system - 0.05);
  Id Id_coeff_rcs = vtkm::Round(coeff_rcs);
  Id rs_num = Id_coeff_rcs * rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * rs * rs * rs) + 1;
  Id rc_num = Id_coeff_rcs * rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * rc * rc * rc) + 1;
  Id rcs_num = rc_num - rs_num;
  Real random_num = _para.GetParameter<Real>(PARA_NEIGHBOR_SAMPLE_NUM);
  Real random_rate = random_num / (rcs_num);
  Id pice_num = std::ceil(1.0 / random_rate);
  //std::cout << "rs_num = " << rs_num << " rc_num = " << rc_num << " rcs_num = " << rcs_num
  //          << " random_num = " << random_num << " random_rate = " << random_rate << " pice_num = " << pice_num  << std::endl;
  //Id cutoff_num = rs_num + random_num;
  Id cutoff_num = rc_num + random_num;
  Id verletlist_num = N * cutoff_num;

  num_verletlist.Allocate(N * 2);
  id_verletlist.Allocate(verletlist_num);
  offset_verletlist.Allocate(verletlist_num);

  std::vector<Id> stride_vec(N + 1);
  Id inc = 0;
  std::generate(
    stride_vec.begin(), stride_vec.end(), [&](void) -> Id { return (inc++) * cutoff_num; });
  vtkm::cont::ArrayHandle<vtkm::Id> stride_array = vtkm::cont::make_ArrayHandle(stride_vec);

  auto id_verletlist_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(id_verletlist, stride_array);
  auto offset_verletlist_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(offset_verletlist, stride_array);
  auto num_verletlist_group = vtkm::cont::make_ArrayHandleGroupVec<2>(num_verletlist);

  RunWorklet::ComputeRBLNeighboursOnce(rs_num,
                                       pice_num,
                                       _atoms_id,
                                       _locator,
                                       id_verletlist_group,
                                       num_verletlist_group,
                                       offset_verletlist_group);

  auto special_offsets = _para.GetFieldAsArrayHandle<Id>(field::special_offsets);
  auto special_weights = _para.GetFieldAsArrayHandle<Real>(field::special_weights);
  auto specoal_ids = _para.GetFieldAsArrayHandle<Id>(field::special_ids);
  auto ids_group = vtkm::cont::make_ArrayHandleGroupVecVariable(specoal_ids, special_offsets);
  auto weight_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(special_weights, special_offsets);
  
  
  _EleNearPairtimer.Start();

  // 当前默认使用ERF 注意：这里要看NearForceRBLERFSpecialBonds和NearForceRBLERF
  //if (_use_erf == true)
  //{
   // RunWorklet::NearForceRBLERF(rs_num,
   //                                pice_num,
   //                                _unit_factor._qqr2e,
   //                                _atoms_id,
   //                                _locator,
   //                                _topology,
   //                                _force_function,
   //                                _static_table,
   //                                id_verletlist_group,
   //                                num_verletlist_group,
   //                                offset_verletlist_group,
   //                                corr_force);
  // 这里的判断待定！！
  if (_para.GetParameter<bool>(PARA_DIHEDRALS_FORCE))
  {
    RunWorklet::NearForceRBLERFSpecialBonds(rs_num,
                                            pice_num,
                                            _unit_factor._qqr2e,
                                             box,
                                            _atoms_id,
                                            _locator,
                                            _topology,
                                            _force_function,
                                            _static_table,
                                            id_verletlist_group,
                                            num_verletlist_group,
                                            offset_verletlist_group,
                                            ids_group,
                                            weight_group,
                                            corr_force);
  }
  else
  {
    RunWorklet::NearForceRBLERF(rs_num,
                                pice_num,
                                _unit_factor._qqr2e,
                                 box,
                                _atoms_id,
                                _locator,
                                _topology,
                                _force_function,
                                _static_table,
                                id_verletlist_group,
                                num_verletlist_group,
                                offset_verletlist_group,
                                corr_force);
  }

  //}
  //else
  //{
  //  RunWorklet::NearForceRBL(rs_num,
  //                              pice_num,
  //                              _unit_factor._qqr2e,
  //                              _atoms_id,
  //                              _locator,
  //                              _topology,
  //                              _force_function,
  //                              id_verletlist_group,
  //                              num_verletlist_group,
  //                              offset_verletlist_group,
  //                              corr_force);
  //}

  Vec3f corr_value = vtkm::cont::Algorithm::Reduce(corr_force, vtkm::TypeTraits<Vec3f>::ZeroInitialization()) / N;
  RunWorklet::SumRBLCorrForce(corr_value, corr_force, nearforce);

  //auto N = _position.GetNumberOfValues();
  //vtkm::cont::ArrayHandle<Vec3f> corr_force;
  //corr_force.Allocate(N);
  //ComputeCorrForce(corr_force);
  //Vec3f corr_value =
  //  vtkm::cont::Algorithm::Reduce(corr_force, vtkm::TypeTraits<Vec3f>::ZeroInitialization()) / N;
  //RunWorklet::SumRBLCorrForce(corr_value, corr_force, nearforce);
}

void ExecutionMD::ComputeRBLLJForce(ArrayHandle<Vec3f>& LJforce)
{
  auto N = _position.GetNumberOfValues();
  vtkm::cont::ArrayHandle<Vec3f> corr_ljforce;
  corr_ljforce.Allocate(N);
  auto box = _para.GetParameter<Vec3f>(PARA_BOX);

  auto rc = _para.GetParameter<Real>(PARA_CUTOFF);
  auto rs = _para.GetParameter<Real>(PARA_R_CORE);
  auto rho_system = _para.GetParameter<Real>(PARA_RHO);

  vtkm::cont::ArrayHandle<vtkm::Id> num_verletlist;
  vtkm::cont::ArrayHandle<vtkm::Id> id_verletlist;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> offset_verletlist;

  //0.05 is a temp parameter, to fit with the low density of system requirement
  //so as the RBL for LJ
  Real coeff_rcs = 1.0 + (0.05 / rho_system - 0.05);
  Id Id_coeff_rcs = vtkm::Round(coeff_rcs);
  Id rs_num = Id_coeff_rcs * rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * rs * rs * rs) + 1;
  Id rc_num = Id_coeff_rcs * rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * rc * rc * rc) + 1;
  Id rcs_num = rc_num - rs_num;
  Real random_num = _para.GetParameter<Real>(PARA_NEIGHBOR_SAMPLE_NUM);
  Real random_rate = random_num / (rcs_num);
  Id pice_num = std::ceil(1.0 / random_rate);
  //Id cutoff_num = rs_num + random_num;
  Id cutoff_num = rc_num + random_num;
  Id verletlist_num = N * cutoff_num;

  num_verletlist.Allocate(N * 2);
  id_verletlist.Allocate(verletlist_num);
  offset_verletlist.Allocate(verletlist_num);

  std::vector<Id> stride_vec(N + 1);
  Id inc = 0;
  std::generate(
    stride_vec.begin(), stride_vec.end(), [&](void) -> Id { return (inc++) * cutoff_num; });
  vtkm::cont::ArrayHandle<vtkm::Id> stride_array = vtkm::cont::make_ArrayHandle(stride_vec);

  auto id_verletlist_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(id_verletlist, stride_array);
  auto offset_verletlist_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(offset_verletlist, stride_array);
  auto num_verletlist_group = vtkm::cont::make_ArrayHandleGroupVec<2>(num_verletlist);

  RunWorklet::ComputeRBLNeighboursOnce(rs_num,
                                          pice_num,
                                          _atoms_id,
                                          _locator,
                                          id_verletlist_group,
                                          num_verletlist_group,
                                          offset_verletlist_group);

  RunWorklet::LJForceRBL(rs_num,
                               pice_num,
                                box,
                               _atoms_id,
                               _locator,
                               _topology,
                               _force_function,
                               id_verletlist_group,
                               num_verletlist_group,
                               offset_verletlist_group,
                               corr_ljforce);

  Vec3f corr_value =
    vtkm::cont::Algorithm::Reduce(corr_ljforce, vtkm::TypeTraits<Vec3f>::ZeroInitialization()) / N;
  RunWorklet::SumRBLCorrForce(corr_value, corr_ljforce, LJforce);
}

void ExecutionMD::ComputeVerletlistNearForce(ArrayHandle<Vec3f>& nearforce)
{
  auto cut_off = _para.GetParameter<Real>(PARA_CUTOFF);
  auto box = _para.GetParameter<Vec3f>(PARA_BOX);
  auto N = _position.GetNumberOfValues();
  auto rho_system = _para.GetParameter<Real>(PARA_RHO);
  auto max_j_num = rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * cut_off * cut_off * cut_off) + 1;
  auto verletlist_num = N * N;

  ArrayHandle<Id> id_verletlist;
  ArrayHandle<Vec3f> offset_verletlist;
  offset_verletlist.Allocate(verletlist_num);
  id_verletlist.Allocate(verletlist_num);

  std::vector<Id> temp_vec(N + 1);
  Id inc = 0;
  std::generate(temp_vec.begin(), temp_vec.end(), [&](void) -> Id { return (inc++) * N; });
  vtkm::cont::ArrayHandle<vtkm::Id> temp_offset = vtkm::cont::make_ArrayHandle(temp_vec);

  auto id_verletlist_group = vtkm::cont::make_ArrayHandleGroupVecVariable(id_verletlist, temp_offset);
  auto offset_verletlist_group = vtkm::cont::make_ArrayHandleGroupVecVariable(offset_verletlist, temp_offset);
  vtkm::cont::ArrayHandle<vtkm::Id> num_verletlist;

  RunWorklet::ComputeNeighbours( cut_off, _atoms_id, _locator, id_verletlist_group, num_verletlist, offset_verletlist_group);

  // 默认使用erf
  //if (_use_erf == true)
  //{
    //RunWorklet::NearForceVerletERF(cut_off,
    //                                  _atoms_id,
    //                                  _locator,
    //                                  _topology,
    //                                  _force_function,
    //                                  _static_table,
    //                                  id_verletlist_group,
    //                                  num_verletlist,
    //                                  offset_verletlist_group,
    //                                  nearforce);
  //}
  //else
  //{
    RunWorklet::NearForceVerlet(cut_off,
                                 box,
                                 _atoms_id,
                                 _locator,
                                 _topology,
                                 _force_function,
                                 id_verletlist_group,
                                 num_verletlist,
                                 offset_verletlist_group,
                                 nearforce);
  //}
}

void ExecutionMD::ComputeVerletlistLJForce(ArrayHandle<Vec3f>& ljforce)
{
  auto cut_off = _para.GetParameter<Real>(PARA_CUTOFF);
  auto box = _para.GetParameter<Vec3f>(PARA_BOX);
  auto N = _position.GetNumberOfValues();
  auto rho_system = _para.GetParameter<Real>(PARA_RHO);
  auto max_j_num = rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * cut_off * cut_off * cut_off) + 1;
  auto verletlist_num = N * max_j_num;

  ArrayHandle<Id> id_verletlist;
  ArrayHandle<Vec3f> offset_verletlist;
  offset_verletlist.Allocate(verletlist_num);
  id_verletlist.Allocate(verletlist_num);

  std::vector<Id> temp_vec(N + 1);
  Id inc = 0;
  std::generate(temp_vec.begin(), temp_vec.end(), [&](void) -> Id { return (inc++) * max_j_num; });
  vtkm::cont::ArrayHandle<vtkm::Id> temp_offset = vtkm::cont::make_ArrayHandle(temp_vec);

  vtkm::cont::ArrayHandle<vtkm::Id> num_verletlist;
  auto id_verletlist_group = vtkm::cont::make_ArrayHandleGroupVecVariable(id_verletlist, temp_offset);
  auto offset_verletlist_group = vtkm::cont::make_ArrayHandleGroupVecVariable(offset_verletlist, temp_offset);

  RunWorklet::ComputeNeighbours( cut_off, _atoms_id, _locator, id_verletlist_group, num_verletlist, offset_verletlist_group);

  RunWorklet::LJForceVerlet(cut_off,
                              box,
                            _atoms_id,
                            _locator,
                            _topology,
                            _force_function,
                            id_verletlist_group,
                            num_verletlist,
                            offset_verletlist_group,
                            ljforce);
}

void ExecutionMD::ComputeOriginalLJForce(ArrayHandle<Vec3f>& ljforce)
{
  auto cut_off = _para.GetParameter<Real>(PARA_CUTOFF);

  RunWorklet::LJForceWithPeriodicBC(
    cut_off, _atoms_id, _locator, _topology, _force_function, ljforce);
}

void ExecutionMD::ComputeRBLEAMForce(ArrayHandle<Vec3f>& force)
{
  ArrayHandle<Real> fp;
  auto Vlength = _para.GetParameter<Real>(PARA_VLENGTH);
  auto box = _para.GetParameter<Vec3f>(PARA_BOX);
  auto rhor_spline = _para.GetFieldAsArrayHandle<Vec7f>(field::rhor_spline);
  auto frho_spline = _para.GetFieldAsArrayHandle<Vec7f>(field::frho_spline);
  auto z2r_spline = _para.GetFieldAsArrayHandle<Vec7f>(field::z2r_spline);
  auto N = _position.GetNumberOfValues();

  auto rc = _para.GetParameter<Real>(PARA_CUTOFF);
  auto rs = _para.GetParameter<Real>(PARA_R_CORE);
  auto rho_system = _para.GetParameter<Real>(PARA_RHO);
  vtkm::cont::ArrayHandle<Vec3f> corr_force;
  corr_force.Allocate(N);

  vtkm::cont::ArrayHandle<vtkm::Id> num_verletlist;
  vtkm::cont::ArrayHandle<vtkm::Id> id_verletlist;
  vtkm::cont::ArrayHandle<vtkm::Vec3f> offset_verletlist;

  Id rs_num = rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * rs * rs * rs) + 1;
  Id rc_num = rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * rc * rc * rc) + 1;
  Id rcs_num = rc_num - rs_num;
  Real random_num = _para.GetParameter<Real>(PARA_NEIGHBOR_SAMPLE_NUM);
  Real random_rate = random_num / (rcs_num);
  Id pice_num = std::ceil(1.0 / random_rate);
  //Id cutoff_num = rs_num + random_num;
  Id cutoff_num = rc_num + random_num;
  Id verletlist_num = N * cutoff_num;

  num_verletlist.Allocate(N * 2);
  id_verletlist.Allocate(verletlist_num);
  offset_verletlist.Allocate(verletlist_num);

  std::vector<Id> stride_vec(N + 1);
  Id inc = 0;
  std::generate(
    stride_vec.begin(), stride_vec.end(), [&](void) -> Id { return (inc++) * cutoff_num; });
  vtkm::cont::ArrayHandle<vtkm::Id> stride_array = vtkm::cont::make_ArrayHandle(stride_vec);

  auto id_verletlist_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(id_verletlist, stride_array);
  auto offset_verletlist_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(offset_verletlist, stride_array);
  auto num_verletlist_group = vtkm::cont::make_ArrayHandleGroupVec<2>(num_verletlist);

  RunWorklet::ComputeRBLNeighboursOnce(rs_num,
                                          pice_num,
                                          _atoms_id,
                                          _locator,
                                          id_verletlist_group,
                                          num_verletlist_group,
                                          offset_verletlist_group);

  RunWorklet::EAMfp(rc,
                    box,
                    rs,
                    rs_num,
                    pice_num,
                    _atoms_id,
                    rhor_spline,
                    frho_spline,
                    _locator,
                    _force_function,
                    id_verletlist_group,
                    num_verletlist_group,
                    offset_verletlist_group,
                    fp);

 RunWorklet::EAMRBLForce(rc,
                         box,
                         rs,
                         rs_num,                            
                         pice_num,
                         _atoms_id,
                         rhor_spline,
                         z2r_spline,
                         fp,
                         _locator,
                         _force_function,
                         id_verletlist_group,
                         num_verletlist_group,
                         offset_verletlist_group,
                         corr_force);

  Vec3f corr_value =
    vtkm::cont::Algorithm::Reduce(corr_force, vtkm::TypeTraits<Vec3f>::ZeroInitialization()) / N;
 RunWorklet::SumRBLCorrForce(corr_value, corr_force, force);
}

void ExecutionMD::ComputeVerletlistEAMForce(ArrayHandle<Vec3f>& force)
{
 ArrayHandle<Real> fp;
 auto Vlength = _para.GetParameter<Real>(PARA_VLENGTH);
 auto box = _para.GetParameter<Vec3f>(PARA_BOX);
 auto rhor_spline = _para.GetFieldAsArrayHandle<Vec7f>(field::rhor_spline);
 auto frho_spline = _para.GetFieldAsArrayHandle<Vec7f>(field::frho_spline);
 auto z2r_spline = _para.GetFieldAsArrayHandle<Vec7f>(field::z2r_spline);

  auto cut_off = _para.GetParameter<Real>(PARA_CUTOFF);
  auto N = _position.GetNumberOfValues();
  auto rho_system = _para.GetParameter<Real>(PARA_RHO);
  auto max_j_num = rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * cut_off * cut_off * cut_off) + 1;
  auto verletlist_num = N * max_j_num;

  ArrayHandle<Id> id_verletlist;
  ArrayHandle<Vec3f> offset_verletlist;
  offset_verletlist.Allocate(verletlist_num);
  id_verletlist.Allocate(verletlist_num);

  std::vector<Id> temp_vec(N + 1);
  Id inc = 0;
  std::generate(temp_vec.begin(), temp_vec.end(), [&](void) -> Id { return (inc++) * max_j_num; });
  vtkm::cont::ArrayHandle<vtkm::Id> temp_offset = vtkm::cont::make_ArrayHandle(temp_vec);

  vtkm::cont::ArrayHandle<vtkm::Id> num_verletlist;
  auto id_verletlist_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(id_verletlist, temp_offset);
  auto offset_verletlist_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(offset_verletlist, temp_offset);

  RunWorklet::ComputeNeighbours(
    cut_off, _atoms_id, _locator, id_verletlist_group, num_verletlist, offset_verletlist_group);

  RunWorklet::EAMfpVerlet(cut_off,
                          box,
                          _atoms_id,
                          rhor_spline,
                          frho_spline,
                          _locator,
                          _force_function,
                          id_verletlist_group,
                          num_verletlist,
                          offset_verletlist_group,
                          fp);

  RunWorklet::EAMForceVerlet(cut_off,
                             box,
                             _atoms_id,
                              rhor_spline,
                              z2r_spline,
                              fp,
                             _locator,
                             _force_function,
                             id_verletlist_group,
                             num_verletlist,
                             offset_verletlist_group,
                             force);
}

void ExecutionMD::ComputeOriginalEAMForce(ArrayHandle<Vec3f>& force)
{
  auto cut_off = _para.GetParameter<Real>(EAM_PARA_CUTOFF);
  auto Vlength = _para.GetParameter<Real>(PARA_VLENGTH);
  auto box = _para.GetParameter<Vec3f>(PARA_BOX);

  auto rhor_spline = _para.GetFieldAsArrayHandle<Vec7f>(field::rhor_spline);
  auto frho_spline = _para.GetFieldAsArrayHandle<Vec7f>(field::frho_spline);
  auto z2r_spline = _para.GetFieldAsArrayHandle<Vec7f>(field::z2r_spline);

  ArrayHandle<Real> EAM_rho;
  ArrayHandle<Real> fp;

  //1:compute _EAM_rho   = density at each atom
  RunWorklet::EAM_rho(
    cut_off, box, _atoms_id, rhor_spline, _locator, _topology, _force_function, EAM_rho);

  // 2:compute fp    = derivative of embedding energy at each atom
  RunWorklet::EAM_fp(_atoms_id, EAM_rho, frho_spline, _locator, _topology, _force_function, fp);

  // 3:compute force  = EAM_force
  RunWorklet::EAM_force(cut_off,
                        box,
                        _atoms_id,
                        fp,
                        rhor_spline,
                        z2r_spline,
                        _locator,
                        _topology,
                        _force_function,
                        force);
}

void ExecutionMD::ComputeSpecialBondsLJForce(ArrayHandle<Vec3f>& ljforce) 
{
  auto cut_off = _para.GetParameter<Real>(PARA_CUTOFF);

  auto special_offsets = _para.GetFieldAsArrayHandle<Id>(field::special_offsets);
  auto special_weights = _para.GetFieldAsArrayHandle<Real>(field::special_weights);
  auto specoal_ids = _para.GetFieldAsArrayHandle<Id>(field::special_ids);
  auto ids_group = vtkm::cont::make_ArrayHandleGroupVecVariable(specoal_ids, special_offsets);
  auto weight_group = vtkm::cont::make_ArrayHandleGroupVecVariable(special_weights, special_offsets);

  RunWorklet::SpecicalBondsLJForce(
    cut_off, _atoms_id, _locator, _topology, _force_function, ids_group, weight_group, ljforce);
}

void ExecutionMD::InitPointLocator()
{
  vtkm::Vec<vtkm::Range, 3> range = _para.GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  vtkm::Vec3f left_bottom{ {
                             static_cast<vtkm::FloatDefault>(range[0].Min),
                           },
                           {
                             static_cast<vtkm::FloatDefault>(range[1].Min),
                           },
                           {
                             static_cast<vtkm::FloatDefault>(range[2].Min),
                           } };
  vtkm::Vec3f right_top{ {
                           static_cast<vtkm::FloatDefault>(range[0].Max),
                         },
                         {
                           static_cast<vtkm::FloatDefault>(range[1].Max),
                         },
                         {
                           static_cast<vtkm::FloatDefault>(range[2].Max),
                         } };
  _locator.SetRange(left_bottom, right_top);

  _locator.SetCutOff(_para.GetParameter<Real>(PARA_CUTOFF));

  _locator.SetRs(_para.GetParameter<Real>(PARA_R_CORE));

  _locator.SetPosition(_position);
}

void ExecutionMD::ComputePotentialEn_test(ArrayHandle<Real>& potential_en)
{
  auto cut_off = _para.GetParameter<Real>(PARA_CUTOFF);
  auto N = _position.GetNumberOfValues();
  auto rho_system = _para.GetParameter<Real>(PARA_RHO);
  auto max_j_num =
    N; //rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * cut_off * cut_off * cut_off) + 1;
  auto verletlist_num = N * max_j_num;

  ArrayHandle<Id> id_verletlist;
  ArrayHandle<Vec3f> offset_verletlist;
  offset_verletlist.AllocateAndFill(verletlist_num, { 0, 0, 0 });
  id_verletlist.AllocateAndFill(verletlist_num, 0);

  std::vector<Id> temp_vec(N + 1);
  Id inc = 0;
  std::generate(temp_vec.begin(), temp_vec.end(), [&](void) -> Id { return (inc++) * max_j_num; });
  vtkm::cont::ArrayHandle<vtkm::Id> temp_offset = vtkm::cont::make_ArrayHandle(temp_vec);

  vtkm::cont::ArrayHandle<vtkm::Id> num_verletlist;
  auto id_verletlist_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(id_verletlist, temp_offset);
  auto offset_verletlist_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(offset_verletlist, temp_offset);

  //并行测试代码：
  vtkm::cont::ArrayHandle<vtkm::Id> num_list;
  std::vector<Id> num_list_vec(N, 0);
  vtkm::cont::ArrayHandle<vtkm::Id> num_size_valid;
  std::vector<Id> num_size_valid_vec(N, 0);
  ArrayHandle<Vec3f> position_temp;
  position_temp.AllocateAndFill(N * N, { 0, 0, 0 });

  auto position_temp_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(position_temp, temp_offset);

  auto id_listgroup_write = id_verletlist_group.WritePortal();
  auto position_temp_group_write = position_temp_group.WritePortal();
  auto position_read = _position.ReadPortal();

  // 加上诡原子
  for (auto it = list_id_map_new.begin(); it != list_id_map_new.end(); ++it)
  {
    auto i = it->first;
    int count_num = 0;
    auto n_neighbor = list_id_map_new[i].neighborIds.size();
    auto n_ghost = list_id_map_new[i].ghostNeighborIds.size();

    num_size_valid_vec[i] = n_neighbor;
    num_list_vec[i] = n_neighbor + n_ghost;

    for (size_t j = 0; j < n_neighbor; j++)
    {
      auto id_neighbor = list_id_map_new[i].neighborIds[j];
      id_listgroup_write.Get(i)[j] = id_neighbor; // 索引 0 至 j-1：为实际原子Id；
      position_temp_group_write.Get(i)[id_neighbor] = position_read.Get(id_neighbor);
    }
    for (size_t k = n_neighbor; k < (n_neighbor + n_ghost); k++)
    {
      auto id_ghost = list_id_map_new[i].ghostNeighborIds[count_num];
      id_listgroup_write.Get(i)[k] = id_ghost; // 索引 j 至 (j+k)-1 为诡原子
      position_temp_group_write.Get(i)[id_ghost] =
        list_id_map_new[i].ghostNeighborPositions[count_num];
      count_num++;
    }
  }
  num_list = vtkm::cont::make_ArrayHandle(num_list_vec);
  num_size_valid = vtkm::cont::make_ArrayHandle(num_size_valid_vec);



  ArrayHandle<Real> potential_en_list;
  potential_en_list.AllocateAndFill(N * N, 0);
  std::vector<Id> potential_en_vec(N + 1);

  vtkm::cont::ArrayHandle<vtkm::Id> potential_en_offset = vtkm::cont::make_ArrayHandle(temp_vec);
  auto potential_en_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(potential_en_list, potential_en_offset);

  ArrayHandle<Real> potential_reduce_list;
  potential_reduce_list.AllocateAndFill(N * N, 0);
  auto potential_reduce_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(potential_reduce_list, potential_en_offset);

  vtkm::cont::Invoker{}(ComputePoentialEnWorklet{ cut_off },
                        _atoms_id,
                        _locator,
                        _topology,
                        _force_function,
                        id_verletlist_group,
                        num_list,
                        offset_verletlist_group,
                        position_temp_group,
                        potential_en_group,
                        potential_reduce_group,
                        potential_en);

  /*外面累加ljfoce*/
  ArrayHandle<Real> potential_redue_lj;
  potential_redue_lj.AllocateAndFill(N, 0);
  auto potential_redue_lj_write = potential_redue_lj.WritePortal();
  for (size_t i = 0; i < N; i++)
  {
    auto potential_reduce_group_read = potential_reduce_group.ReadPortal().Get(i);
    Real temp_potential = 0.0;
    for (size_t j = 0; j < N; j++)
    {
      temp_potential += potential_reduce_group_read[j];
    }
    //auto temp_ljpotential = potential_redue_read.Get(i);
    potential_redue_lj_write.Set(i, temp_potential /*+ temp_ljpotential*/);
  }

  vtkm::cont::ArrayCopy(potential_redue_lj, potential_en);
}

void ExecutionMD::ComputeVerletlistNearForceMPITest(ArrayHandle<Vec3f>& nearforce)
{
  auto cut_off = _para.GetParameter<Real>(PARA_CUTOFF);
  auto N = _position.GetNumberOfValues();
  auto rho_system = _para.GetParameter<Real>(PARA_RHO);
  auto max_j_num =
    rho_system * vtkm::Ceil(4.0 / 3.0 * vtkm::Pif() * cut_off * cut_off * cut_off) + 1;
  auto verletlist_num = N * N;

  ArrayHandle<Id> id_verletlist;
  ArrayHandle<Vec3f> offset_verletlist;
  offset_verletlist.Allocate(verletlist_num);
  id_verletlist.Allocate(verletlist_num);

  std::vector<Id> temp_vec(N + 1);
  Id inc = 0;
  std::generate(temp_vec.begin(), temp_vec.end(), [&](void) -> Id { return (inc++) * N; });
  vtkm::cont::ArrayHandle<vtkm::Id> temp_offset = vtkm::cont::make_ArrayHandle(temp_vec);

  auto id_verletlist_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(id_verletlist, temp_offset);
  auto offset_verletlist_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(offset_verletlist, temp_offset);
  vtkm::cont::ArrayHandle<vtkm::Id> num_verletlist;

  // 并行测试代码：
  vtkm::cont::ArrayHandle<vtkm::Id> num_list;
  std::vector<Id> num_list_vec(N, 0);
  vtkm::cont::ArrayHandle<vtkm::Id> num_size_valid;
  std::vector<Id> num_size_valid_vec(N, 0);
  ArrayHandle<Vec3f> position_temp;
  position_temp.AllocateAndFill(N * N, { 0, 0, 0 });
  ArrayHandle<Vec3f> force_list;
  force_list.AllocateAndFill(N * N, { 0, 0, 0 });

  std::vector<Id> force_vec(N + 1);
  Id force_inc = 0;
  std::generate(force_vec.begin(), force_vec.end(), [&](void) -> Id { return (force_inc++) * N; });
  vtkm::cont::ArrayHandle<vtkm::Id> force_offset = vtkm::cont::make_ArrayHandle(force_vec);

  auto force_group = vtkm::cont::make_ArrayHandleGroupVecVariable(force_list, force_offset);
  auto position_temp_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(position_temp, force_offset);

  auto id_listgroup_write = id_verletlist_group.WritePortal();
  auto position_temp_group_write = position_temp_group.WritePortal();
  auto position_read = _position.ReadPortal();

  // 加上诡原子
  for (auto it = list_id_map_new.begin(); it != list_id_map_new.end(); ++it)
  {
    auto i = it->first;
    int count_num = 0;
    auto n_neighbor = list_id_map_new[i].neighborIds.size();
    auto n_ghost = list_id_map_new[i].ghostNeighborIds.size();

    num_size_valid_vec[i] = n_neighbor;     //it->second[0].size();
    num_list_vec[i] = n_neighbor + n_ghost; // 总数 = 实际+诡原子

    for (size_t j = 0; j < n_neighbor; j++)
    {
      auto id_neighbor = list_id_map_new[i].neighborIds[j];
      id_listgroup_write.Get(i)[j] = id_neighbor; // 索引 0 至 j-1：为实际原子Id；
      position_temp_group_write.Get(i)[id_neighbor] = position_read.Get(id_neighbor);
    }
    for (size_t k = n_neighbor; k < (n_neighbor + n_ghost); k++)
    {
      auto id_ghost = list_id_map_new[i].ghostNeighborIds[count_num];
      id_listgroup_write.Get(i)[k] = id_ghost; // 索引 j 至 (j+k)-1 为诡原子
      position_temp_group_write.Get(i)[id_ghost] =
        list_id_map_new[i].ghostNeighborPositions[count_num];
      count_num++;
    }
  }
  num_list = vtkm::cont::make_ArrayHandle(num_list_vec);
  num_size_valid = vtkm::cont::make_ArrayHandle(num_size_valid_vec);

  ArrayHandle<Vec3f> force_reduce_list;
  force_reduce_list.AllocateAndFill(N * N, { 0, 0, 0 });
  auto force_reduce_group =
    vtkm::cont::make_ArrayHandleGroupVecVariable(force_reduce_list, force_offset);



  vtkm::cont::Invoker{}(ComputeNearForceVerlet_MPIWorklet{ cut_off },
                        _atoms_id,
                        _locator,
                        _topology,
                        _force_function,
                        id_verletlist_group,
                        num_list,
                        position_temp_group,
                        force_reduce_group,
                        nearforce);

  /*外面累加ljfoce*/
  ArrayHandle<Vec3f> force_redue_lj;
  force_redue_lj.AllocateAndFill(N, { 0.0, 0.0, 0.0 });
  auto force_redue_lj_write = force_redue_lj.WritePortal();
  for (size_t i = 0; i < N; i++)
  {
    auto force_reduce_group_read = force_reduce_group.ReadPortal().Get(i);
    Vec3f temp_force = { 0.0, 0.0, 0.0 };
    for (size_t j = 0; j < N; j++)
    {
      temp_force[0] += force_reduce_group_read[j][0];
      temp_force[1] += force_reduce_group_read[j][1];
      temp_force[2] += force_reduce_group_read[j][2];
    }
    //auto temp_ljforce = force_redue_read.Get(i);
    force_redue_lj_write.Set(i, temp_force /*+ temp_ljforce*/);
  }
  vtkm::cont::ArrayCopy(force_redue_lj, nearforce);

#ifdef NEIGHBOR_GHOST
  if (num_cont == 1 || num_cont == 100)
  {
    /*输出邻居粒子ID和诡原子ID 测试是否输出所有的邻居粒子*/
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (list_id_map_new.find(6) != list_id_map_new.end())
    {
      std::ofstream test_6;
      test_6.open(
        "NEIGHBOR_GHOST_6_" + std::to_string(num_cont) +
        ".txt"); // 因为有可能是不同进程的，这里的文件打开 和生成要在 if 里面完成，因为可能会被其他进程打开 导致无法填充数据
      auto n_neighbor = list_id_map_new[6].neighborIds.size();
      auto n_ghost = list_id_map_new[6].ghostNeighborIds.size();
      for (size_t i = 0; i < n_neighbor; i++)
      {
        auto id = list_id_map_new[6].neighborIds[i];
        test_6 << list_id_map_new[6].neighborIds[i] << "," << _position.ReadPortal().Get(id)[0]
               << "," << _position.ReadPortal().Get(id)[1] << ","
               << _position.ReadPortal().Get(id)[2] << std::endl;
      }

      for (size_t i = 0; i < n_ghost; i++)
      {
        test_6 << list_id_map_new[6].ghostNeighborIds[i] << ","
               << list_id_map_new[6].ghostNeighborPositions[i][0] << ","
               << list_id_map_new[6].ghostNeighborPositions[i][1] << ","
               << list_id_map_new[6].ghostNeighborPositions[i][2] << std::endl;
      }
    }

    if (list_id_map_new.find(1) != list_id_map_new.end())
    {
      std::ofstream test_1;
      test_1.open(
        "NEIGHBOR_GHOST_1_" + std::to_string(num_cont) +
        ".txt"); // 因为有可能是不同进程的，这里的文件打开 和生成要在 if 里面完成，因为可能会被其他进程打开 导致无法填充数据
      auto n_neighbor = list_id_map_new[1].neighborIds.size();
      auto n_ghost = list_id_map_new[1].ghostNeighborIds.size();
      for (size_t i = 0; i < n_neighbor; i++)
      {
        auto id = list_id_map_new[1].neighborIds[i];
        test_1 << list_id_map_new[1].neighborIds[i] << "," << _position.ReadPortal().Get(id)[0]
               << "," << _position.ReadPortal().Get(id)[1] << ","
               << _position.ReadPortal().Get(id)[2] << std::endl;
      }

      for (size_t i = 0; i < n_ghost; i++)
      {
        test_1 << list_id_map_new[1].ghostNeighborIds[i] << ","
               << list_id_map_new[1].ghostNeighborPositions[i][0] << ","
               << list_id_map_new[1].ghostNeighborPositions[i][1] << ","
               << list_id_map_new[1].ghostNeighborPositions[i][2] << std::endl;
      }
    }

    if (list_id_map_new.find(998) != list_id_map_new.end())
    {
      std::ofstream test_998;
      test_998.open(
        "NEIGHBOR_GHOST_998_" + std::to_string(num_cont) +
        ".txt"); // 因为有可能是不同进程的，这里的文件打开 和生成要在 if 里面完成，因为可能会被其他进程打开 导致无法填充数据
      auto n_neighbor = list_id_map_new[998].neighborIds.size();
      auto n_ghost = list_id_map_new[998].ghostNeighborIds.size();
      for (size_t i = 0; i < n_neighbor; i++)
      {
        auto id = list_id_map_new[998].neighborIds[i];
        test_998 << list_id_map_new[998].neighborIds[i] << "," << _position.ReadPortal().Get(id)[0]
                 << "," << _position.ReadPortal().Get(id)[1] << ","
                 << _position.ReadPortal().Get(id)[2] << std::endl;
      }

      for (size_t i = 0; i < n_ghost; i++)
      {
        test_998 << list_id_map_new[998].ghostNeighborIds[i] << ","
                 << list_id_map_new[998].ghostNeighborPositions[i][0] << ","
                 << list_id_map_new[998].ghostNeighborPositions[i][1] << ","
                 << list_id_map_new[998].ghostNeighborPositions[i][2] << std::endl;
      }
    }
  }
#endif

#ifdef INIT_NEIGHBOR_GHOST
  if (num_cont == 0)
  {
    /*输出邻居粒子ID和诡原子ID 测试是否输出所有的邻居粒子*/
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (list_id_map_new.find(6) != list_id_map_new.end())
    {
      std::ofstream test_6;
      test_6.open("INIT_NEIGHBOR_GHOST_6_" + std::to_string(num_cont) + ".txt");
      auto n_neighbor = list_id_map_new[6].neighborIds.size();
      auto n_ghost = list_id_map_new[6].ghostNeighborIds.size();
      for (size_t i = 0; i < n_neighbor; i++)
      {
        auto id = list_id_map_new[6].neighborIds[i];
        test_6 << list_id_map_new[6].neighborIds[i] << "," << _position.ReadPortal().Get(id)[0]
               << "," << _position.ReadPortal().Get(id)[1] << ","
               << _position.ReadPortal().Get(id)[2] << std::endl;
      }

      for (size_t i = 0; i < n_ghost; i++)
      {
        test_6 << list_id_map_new[6].ghostNeighborIds[i] << ","
               << list_id_map_new[6].ghostNeighborPositions[i][0] << ","
               << list_id_map_new[6].ghostNeighborPositions[i][1] << ","
               << list_id_map_new[6].ghostNeighborPositions[i][2] << std::endl;
      }
    }

    if (list_id_map_new.find(998) != list_id_map_new.end())
    {
      std::ofstream test_998;
      test_998.open("INIT_NEIGHBOR_GHOST_998_" + std::to_string(num_cont) + ".txt");
      auto n_neighbor = list_id_map_new[998].neighborIds.size();
      auto n_ghost = list_id_map_new[998].ghostNeighborIds.size();
      for (size_t i = 0; i < n_neighbor; i++)
      {
        auto id = list_id_map_new[998].neighborIds[i];
        test_998 << list_id_map_new[998].neighborIds[i] << "," << _position.ReadPortal().Get(id)[0]
                 << "," << _position.ReadPortal().Get(id)[1] << ","
                 << _position.ReadPortal().Get(id)[2] << std::endl;
      }

      for (size_t i = 0; i < n_ghost; i++)
      {
        test_998 << list_id_map_new[998].ghostNeighborIds[i] << ","
                 << list_id_map_new[998].ghostNeighborPositions[i][0] << ","
                 << list_id_map_new[998].ghostNeighborPositions[i][1] << ","
                 << list_id_map_new[998].ghostNeighborPositions[i][2] << std::endl;
      }
    }

    if (list_id_map_new.find(1) != list_id_map_new.end())
    {
      std::ofstream test_1;
      test_1.open("INIT_NEIGHBOR_GHOST_1_" + std::to_string(num_cont) + ".txt");
      auto n_neighbor = list_id_map_new[1].neighborIds.size();
      auto n_ghost = list_id_map_new[1].ghostNeighborIds.size();
      for (size_t i = 0; i < n_neighbor; i++)
      {
        auto id = list_id_map_new[1].neighborIds[i];
        test_1 << list_id_map_new[1].neighborIds[i] << "," << _position.ReadPortal().Get(id)[0]
               << "," << _position.ReadPortal().Get(id)[1] << ","
               << _position.ReadPortal().Get(id)[2] << std::endl;
      }

      for (size_t i = 0; i < n_ghost; i++)
      {
        test_1 << list_id_map_new[1].ghostNeighborIds[i] << ","
               << list_id_map_new[1].ghostNeighborPositions[i][0] << ","
               << list_id_map_new[1].ghostNeighborPositions[i][1] << ","
               << list_id_map_new[1].ghostNeighborPositions[i][2] << std::endl;
      }
    }
  }
#endif
  num_cont++; // 用来作统计的参数
}

void ExecutionMD::ComputeNewChargeStructureFactorRBE_MPItest(Vec3f& _box,
                                                     ArrayHandle<Vec3f>& _psample,
                                                     ArrayHandle<Vec2f>& new_rhok)
{
  auto N = _position.GetNumberOfValues();

  const Id p_number = _psample.GetNumberOfValues();

  RunWorklet::ComputePnumberChargeStructureFactor(_box, p_number, N, _atoms_id, _position, _charge, _psample, _rhok_Re, _rhok_Im);

  vtkm::cont::ArrayHandle<Id> psamplekey_out;
  vtkm::cont::ArrayHandle<Real> rhok_Re_reduce;
  vtkm::cont::ArrayHandle<Real> rhok_Im_reduce;


  vtkm::cont::Algorithm::ReduceByKey<Id, Real>(_psamplekey, _rhok_Re, psamplekey_out, rhok_Re_reduce, vtkm::Add());
  vtkm::cont::Algorithm::ReduceByKey<Id, Real>(_psamplekey, _rhok_Im, psamplekey_out, rhok_Im_reduce, vtkm::Add());

  std::vector<float> rhok_Re_reduce_vec;
  std::vector<float> rhok_Im_reduce_vec;

  auto read_rhok_Re_reduce = rhok_Re_reduce.ReadPortal();
  auto read_rhok_Im_reduce = rhok_Im_reduce.ReadPortal();
  for (size_t i = 0; i < rhok_Im_reduce.GetNumberOfValues(); i++)
  {
    rhok_Re_reduce_vec.push_back(read_rhok_Re_reduce.Get(i));
    rhok_Im_reduce_vec.push_back(read_rhok_Im_reduce.Get(i));
  }


  auto rhok_Re_reduce_temp = rbmd_parallel->GetGlobalRho(rhok_Re_reduce_vec); // 要更换接口
  auto rhok_Im_reduce_temp = rbmd_parallel->GetGlobalRho(rhok_Im_reduce_vec); // 要更换接口

  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(rhok_Re_reduce_temp), rhok_Re_reduce);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(rhok_Im_reduce_temp), rhok_Im_reduce);

  RunWorklet::ChangePnumberChargeStructureFactor(rhok_Re_reduce, rhok_Im_reduce, new_rhok);
}

void ExecutionMD::ComputeRBEEleForce_MPItest(ArrayHandle<Vec3f>& psample,
                                     IdComponent& RBE_P,
                                     ArrayHandle<Vec3f>& RBE_ele_force)
{


  auto Vlength = _para.GetParameter<Real>(PARA_VLENGTH);
  auto box = rbmd_parallel->GetBoxLength();
  ArrayHandle<Vec2f> new_whole_rhok;

  //int rank;
  //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //std::cout << "-------box-------" << box << "-------RBE_P-------" << RBE_P << std::endl;
  //if (rank == 0)
  //{
  //  auto _psample_read = psample.ReadPortal();
  //  for (size_t i = 0; i < psample.GetNumberOfValues(); i++)
  //  {
  //    std::cout << i << " ------psample_0--------- " << _psample_read.Get(i) << std::endl;
  //  }
  //}
  //else if (rank == 1)
  //{
  //  auto _psample_read = psample.ReadPortal();
  //  for (size_t i = 0; i < psample.GetNumberOfValues(); i++)
  //  {
  //    std::cout << i << " ------psample_1--------- " << _psample_read.Get(i) << std::endl;
  //  }
  //}  


  ComputeNewChargeStructureFactorRBE_MPItest(box, psample, new_whole_rhok);
  

  //if (rank == 0)
  //{
  //  auto _new_whole_rhok_read = new_whole_rhok.ReadPortal();
  //  for (size_t i = 0; i < new_whole_rhok.GetNumberOfValues(); i++)
  //  {
  //    std::cout << i << " ------new_whole_rhok_0--------- " << _new_whole_rhok_read.Get(i) << std::endl;
  //  }
  //}
  //else if (rank == 1)
  //{
  //  auto _new_whole_rhok_read = new_whole_rhok.ReadPortal();
  //  for (size_t i = 0; i < new_whole_rhok.GetNumberOfValues(); i++)
  //  {
  //    std::cout << i << " ------new_whole_rhok_1--------- " << _new_whole_rhok_read.Get(i) << std::endl;
  //  }
  //}  
  //
  //if (rank == 0)
  //{
  //  auto __atoms_id_read = _atoms_id.ReadPortal();
  //  for (size_t i = 0; i < _atoms_id.GetNumberOfValues(); i++)
  //  {
  //    std::cout << i << " ------_atoms_id_0--------- " << __atoms_id_read.Get(i)
  //              << std::endl;
  //  }
  //}
  //else if (rank == 1)
  //{
  //  auto __atoms_id_read = _atoms_id.ReadPortal();
  //  for (size_t i = 0; i < _atoms_id.GetNumberOfValues(); i++)
  //  {
  //    std::cout << i << " ------_atoms_id_1--------- " << __atoms_id_read.Get(i)
  //              << std::endl;
  //  }
  //}  

  //std::vector<Vec2f> rhok_new_vec;
  //auto num = new_whole_rhok.GetNumberOfValues();
  //auto read_rhok = new_whole_rhok.ReadPortal();
  //for (size_t i = 0; i < num; i++)
  //{
  //  rhok_new_vec.push_back(read_rhok.Get(i));
  //}

  RunWorklet::ComputeNewRBEForce(RBE_P,
                                 _atoms_id,
                                 psample,
                                 new_whole_rhok,
                                 _force_function,
                                 _topology,
                                 _locator,
                                 RBE_ele_force); 

  //if (rank == 0)
  //{
  //  auto _RBE_ele_force_read = RBE_ele_force.ReadPortal();
  //  for (size_t i = 0; i < RBE_ele_force.GetNumberOfValues(); i++)
  //  {
  //    std::cout << i << " ------RBE_ele_force_0--------- " << _RBE_ele_force_read.Get(i)
  //              << std::endl;
  //  }
  //}
  //else if (rank == 1)
  //{
  //  auto _RBE_ele_force_read = RBE_ele_force.ReadPortal();
  //  for (size_t i = 0; i < RBE_ele_force.GetNumberOfValues(); i++)
  //  {
  //    std::cout << i << " ------RBE_ele_force_1--------- " << _RBE_ele_force_read.Get(i)
  //              << std::endl;
  //  }
  //}  

  num_cont++;
}