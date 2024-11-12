﻿#include "ExecutionNVT.h"
#include "Application.h"
#include "FieldName.h"
#include "MeshFreeCondition.h"
#include "math/Math.h"
#include "math/Utiles.h"
#include "RBEPSample.h"
#include "run/worklet/MolecularWorklet.h"
#include "run/worklet/RunWorklet.h"
#include <cmath> // erfc(x)
#include <fstream>
#include <vtkm/Math.h>
#include <vtkm/Pair.h>
#include <vtkm/VectorAnalysis.h>
#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/ArrayHandleGroupVec.h>
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>
#include <vtkm/cont/ConvertNumComponentsToOffsets.h>
#include <vtkm/worklet/Keys.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/worklet/WorkletReduceByKey.h>
#include "output/worklet/OutPutWorklet.h"
//#define INIT_POS_FORCE 
//#define OUTPUT 
//#define POS_PRR_POS 
//#define POS_POST_POS

//#define INIT_FORCE
//#define FORCE
//#define PSAMPLE_RBEP

#define MIN(A, B) ((A) < (B) ? (A) : (B))
#define MAX(A, B) ((A) > (B) ? (A) : (B))

ExecutionNVT::ExecutionNVT(const Configuration& cfg)
  : ExecutionMD(cfg)
  , _executioner((_app.GetExecutioner()))
  //, cell_data(rbmd_parallel->GetCellData())

{
  ExecutionMD::SetParameters();
  InitParameters();
  _file.open("lj_potential_energy.csv");
}

void ExecutionNVT::Init()
{
  if (_para.GetParameter<std::string>(PARA_FILE_TYPE) == "EAM")
  {
    _potential_file.open(_para.GetParameter<std::string>(PARA_POTENTIAL_FILE));
    ReadPotentialFile(_potential_file);
    InitStyle();
  }

  ExecutionMD::Init();

  InitialCondition();
  //int rank;
  //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //if (rank == 0)
  //{
  //  auto _psample_read = _psample.ReadPortal();
  //  for (size_t i = 0; i < _psample.GetNumberOfValues(); i++)
  //  {
  //    std::cout << i << " ------MPILinkCell前：0--------- " << _psample_read.Get(i) << std::endl;
  //  }
  //}
  //else if (rank == 1)
  //{
  //  auto _psample_read = _psample.ReadPortal();
  //  for (size_t i = 0; i < _psample.GetNumberOfValues(); i++)
  //  {
  //    std::cout << i << " ------MPILinkCell前：1--------- " << _psample_read.Get(i) << std::endl;
  //  }
  //}  
  MPILinkCell();
  _all_force.AllocateAndFill(_position.GetNumberOfValues(), { 0, 0, 0 });
  ComputeForce(); // Presolve force

#ifdef INIT_FORCE
  /*输出循环力*/
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0)
  {
    std::ofstream ljforce_file0;
    ljforce_file0.open("INIT_force_MPI_" + std::to_string(num_cont) + "_0.txt");
    auto force_read = _all_force.ReadPortal();
    for (auto it = list_id_map_new.begin(); it != list_id_map_new.end(); ++it)
    {
      ljforce_file0 << it->first << "," << force_read.Get(it->first) << std::endl;
    }
  }
  else if (rank == 1)
  {
    std::ofstream ljforce_file1;
    ljforce_file1.open("INIT_force_MPI_" + std::to_string(num_cont) + "_1.txt");
    auto force_read = _all_force.ReadPortal();
    for (auto it = list_id_map_new.begin(); it != list_id_map_new.end(); ++it)
    {
      ljforce_file1 << it->first << "," << force_read.Get(it->first) << std::endl;
    }
  }
#endif
  //exit(1);
}

void ExecutionNVT::PreSolve()
{
  _locator.SetPosition(_position);
  if (_para.GetParameter<bool>(PARA_FAR_FORCE))
  {
    PreForce();
  }
  else
  {
    std::shared_ptr<Executioner>& executioner = _app.GetExecutioner();
    _dt = executioner->Dt();
  }
}

void ExecutionNVT::Solve()
{
  auto step=_app.GetExecutioner()->CurrentStep();
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // stage1:
  UpdateVelocity();

#ifdef OUTPUT
  /*输出速度文件*/
  if (42 < step && step < 51 /*step == 47 || step == 48 || step == 49 || step == 50*/)
  {
    if (rank == 0)
    {
      std::ofstream velocity_file0;
      velocity_file0.open("OUTPUT__velocity_MPI_1_" + std::to_string(step) + "_0.txt");
      auto read_velocity = _velocity.ReadPortal();
      for (auto it = list_id_map_new.begin(); it != list_id_map_new.end(); ++it)
      {
        velocity_file0 << it->first << "," << read_velocity.Get(it->first) << std::endl;
      }
    }
    else
    {
      std::ofstream velocity_file1;
      velocity_file1.open("OUTPUT__velocity_MPI_1_" + std::to_string(step) + "_1.txt");
      auto read_velocity = _velocity.ReadPortal();
      for (auto it = list_id_map_new.begin(); it != list_id_map_new.end(); ++it)
      {
        velocity_file1 << it->first << "," << read_velocity.Get(it->first) << std::endl;
      }
    }
  }
#endif

  // stage2:
  UpdatePosition();

#ifdef OUTPUT
  /*输出位置文件*/
  if (42 < step && step < 51 /*step == 47 || step == 48 || step == 49 || step == 50*/)
  {
    if (rank == 0)
    {
      std::ofstream position_file0;
      position_file0.open("OUTPUT__position_MPI_" + std::to_string(step) + "_0.txt");
      auto read_position = _position.ReadPortal();
      for (auto it = list_id_map_new.begin(); it != list_id_map_new.end(); ++it)
      {
        position_file0 << it->first << "," << read_position.Get(it->first) << std::endl;
      }
    }
    else
    {
      std::ofstream position_file1;
      position_file1.open("OUTPUT__position_MPI_" + std::to_string(step) + "_1.txt");
      auto read_position = _position.ReadPortal();
      for (auto it = list_id_map_new.begin(); it != list_id_map_new.end(); ++it)
      {
        position_file1 << it->first << "," << read_position.Get(it->first) << std::endl;
      }
    }
  }
#endif
  //New added
  if (_para.GetParameter<std::string>(PARA_FIX_SHAKE) == "true")
  {
    std::cout << "进入ConstraintA" << std::endl;
    ConstraintA();
  }

  PostCompute();

  // stage3:
  ComputeForce();

#ifdef OUTPUT
  /*输出循环力*/
  if (42 < step && step < 51 /*step == 47 || step == 48 || step == 49 || step == 50*/)
  {
    if (rank == 0)
    {
      std::ofstream ljforce_file0;
      ljforce_file0.open("OUTPUT__force_MPI_" + std::to_string(step) + "_0.txt");
      auto force_read = _all_force.ReadPortal();
      for (auto it = list_id_map_new.begin(); it != list_id_map_new.end(); ++it)
      {
        ljforce_file0 << it->first << "," << force_read.Get(it->first) << std::endl;
      }
    }
    else if (rank == 1)
    {
      std::ofstream ljforce_file1;
      ljforce_file1.open("OUTPUT__force_MPI_" + std::to_string(step) + "_1.txt");
      auto force_read = _all_force.ReadPortal();
      for (auto it = list_id_map_new.begin(); it != list_id_map_new.end(); ++it)
      {
        ljforce_file1 << it->first << "," << force_read.Get(it->first) << std::endl;
      }
    }
  }
#endif

  UpdateVelocity();

  //New added
  if (_para.GetParameter<std::string>(PARA_FIX_SHAKE) == "true")
  {
    std::cout << "进入ConstraintB" << std::endl;
    ConstraintB();
  }

  ComputeTempe();
  UpdateVelocityByTempConType();

#ifdef OUTPUT
  /*输出更新后的速度*/
  if (42 < step && step < 51 /*step == 47 || step == 48 || step == 49 || step == 50*/)
  {
    if (rank == 0)
    {
      std::ofstream velocity_file0;
      velocity_file0.open("OUTPUT__velocity_MPI_2_" + std::to_string(step) + "_0.txt");
      auto read_velocity = _velocity.ReadPortal();
      for (auto it = list_id_map_new.begin(); it != list_id_map_new.end(); ++it)
      {
        velocity_file0 << it->first << "," << read_velocity.Get(it->first) << std::endl;
      }
    }
    else
    {
      std::ofstream velocity_file1;
      velocity_file1.open("OUTPUT__velocity_MPI_2_" + std::to_string(step) + "_1.txt");
      auto read_velocity = _velocity.ReadPortal();
      for (auto it = list_id_map_new.begin(); it != list_id_map_new.end(); ++it)
      {
        velocity_file1 << it->first << "," << read_velocity.Get(it->first) << std::endl;
      }
    }
  }
#endif
}

void ExecutionNVT::PostSolve() 
{
   //2:self_potential_energy_avr
   ArrayHandle<Real> _self_energy;
   auto atoms_id = _para.GetFieldAsArrayHandle<Id>(field::atom_id);
   auto position = _para.GetFieldAsArrayHandle<Vec3f>(field::position);
   Id N = position.GetNumberOfValues();
   Real _far_ele_potential_energy_avr = 0.0;
   //auto charge = _para.GetFieldAsArrayHandle<Real>(field::charge);
   Real _self_potential_energy_avr;
   auto unit_factor = _para.GetParameter<UnitFactor>(PARA_UNIT_FACTOR);
 
   OutPut::ComputeSqCharge(_charge, _self_energy);
   //auto reduce__self_energy = vtkm::cont::Algorithm::Reduce(_self_energy, vtkm::TypeTraits<Real>::ZeroInitialization());
   auto self_potential_energy_total = -vtkm::Sqrt(_alpha / vtkm::Pi()) * vtkm::cont::Algorithm::Reduce(_self_energy, vtkm::TypeTraits<Real>::ZeroInitialization());
 
   auto self_energy_temp = rbmd_parallel->GetGlobalErengy(self_potential_energy_total);
   _self_potential_energy_avr = (self_energy_temp[0]/N) * unit_factor._qqr2e;
   //_self_potential_energy_avr = self_potential_energy_total / N;
   _self_potential_energy_avr = _self_potential_energy_avr * unit_factor._qqr2e;
 
   //3:_far_ele_potential_energy_avr
   Real Volume = vtkm::Pow(_Vlength, 3);
   Real Volume = _box[0] * _box[1] * _box[2];
   auto box = rbmd_parallel->GetBoxLength();
   ArrayHandle<Real> Density_Real;
   ArrayHandle<Real> Density_Image;
   Real far_ele_potential_energy_total = 0.0;
   for (Id i = -_Kmax; i <= _Kmax; i++)
   {
     for (Id j = -_Kmax; j <= _Kmax; j++)
     {
       for (Id k = -_Kmax; k <= _Kmax; k++)
       {
         if (!(i == 0 && j == 0 && k == 0))
         {
           Vec3f K = { Real(i), Real(j), Real(k) };
           //K = 2 * vtkm::Pi() * K / _Vlength;
           K = { Real(2 * vtkm::Pi() * K[0] / _box[0]),
                 Real(2 * vtkm::Pi() * K[1] / _box[1]),
                 Real(2 * vtkm::Pi() * K[2] / _box[2]) };
           Real Range_K = vtkm::Magnitude(K);
           OutPut::ComputeDensity(K, position, _charge, Density_Real, Density_Image); 
 
           Real Value_Re = vtkm::cont::Algorithm::Reduce(Density_Real, vtkm::TypeTraits<Real>::ZeroInitialization());
           Real Value_Im = vtkm::cont::Algorithm::Reduce(Density_Image, vtkm::TypeTraits<Real>::ZeroInitialization());
           auto energys_Re = rbmd_parallel->GetGlobalErengy(Value_Re);
           auto energys_Im = rbmd_parallel->GetGlobalErengy(Value_Im); 
 
          // Real Range_density2 = vtkm::Pow(Value_Re, 2) + vtkm::Pow(Value_Im, 2);
           Real Range_density2 = vtkm::Pow(energys_Re[0], 2) + vtkm::Pow(energys_Im[0], 2);
 
           far_ele_potential_energy_total += vtkm::Exp(-Range_K * Range_K / (4 * _alpha)) * Range_density2 / (Range_K * Range_K);
         }
       }
     }
   }
 
   Real Volume = box[0] * box[1] * box[2];
   far_ele_potential_energy_total = far_ele_potential_energy_total * (2 * vtkm::Pi() / Volume);
   auto energys = rbmd_parallel->GetGlobalErengy(far_ele_potential_energy_total);
 
   _far_ele_potential_energy_avr = energys[0]/N; //far_ele_potential_energy_total / N;
   _far_ele_potential_energy_avr = _far_ele_potential_energy_avr * unit_factor._qqr2e;
 
   _far_ele_potential_energy_avr = _far_ele_potential_energy_avr + _self_potential_energy_avr;
 
   //auto energys = rbmd_parallel->GetGlobalErengy(_far_ele_potential_energy_avr);
 
   try
   {
     if (rbmd_parallel->GetCurrentRank()==0)
     {
       _file << energys[0] << std::endl;
     }
 
     //if (!energys.empty())
     {
       _file << energys[0][0] << std::endl;
     }
   }
   catch (const std::exception& e)
   {
     _file.close();
     console::Error(e.what());
   }
 
   //  // 需要改为 之查看远程力的potential
   //auto position = _para.GetFieldAsArrayHandle<Vec3f>(field::position); 
   //
   //ArrayHandle<Real> lj_potential_energy; 
   //ComputePotentialEn_test(lj_potential_energy); 
   //auto _lj_potential_energy_total = vtkm::cont::Algorithm::Reduce(lj_potential_energy, vtkm::TypeTraits<Real>::ZeroInitialization()); 
   //auto energys = rbmd_parallel->GetGlobalErengy(_lj_potential_energy_total); 
  
   try
   {       
     auto energy_ave = energys[1];
     _file << energy_ave << std::endl;   
   }
   catch (const std::exception& e)
   {
     _file.close();
     console::Error(e.what());
   }
}

void ExecutionNVT::InitialCondition()
{
  ExecutionMD::InitialCondition();
  //_rho
  _Volume = _para.GetParameter<Real>(PARA_VOLUME);
  SetCenterTargetPositions();
  auto target_position = _para.GetFieldAsArrayHandle<Vec3f>(field::target_position);
  _para.SetParameter(PARA_RDF_RHO, target_position.GetNumberOfValues() / _Volume);
  rbmd_parallel = std::make_shared<RBMDParallelUtil>(_para.GetParameter<std::array<double, 3>>(PARA_BOXMIN),
                                       _para.GetParameter<std::array<double, 3>>(PARA_BOXMAX),
                                       static_cast<double>(_para.GetParameter<Real>(PARA_CUTOFF)));

  if (_para.GetParameter<bool>(PARA_FAR_FORCE))
  {
    auto p_number = _para.GetParameter<IdComponent>(PARA_COULOMB_SAMPLE_NUM);
    auto alpha = _para.GetParameter<Real>(PARA_ALPHA);
    rbmd_parallel->SetRbeConfigPnumber(p_number);
    rbmd_parallel->SetRbeConfigAlphaNumber(alpha);

    PreForce();
  }
  _nosehooverxi = 0.0;
}

void ExecutionNVT::ComputeForce()
{
  ComputeAllForce();
  TempConTypeForce();
}

void ExecutionNVT::ComputeAllForce()
{
  auto force_field = _para.GetParameter<std::string>(PARA_FORCE_FIELD_TYPE);
  if ("LJ/CUT" == force_field  )
  {
    //NearForceLJ();
    PreCompute();

    ComputeVerletlistNearForceMPITest(_all_force);
  }

  else if ("LJ/CUT/COUL/LONG" ==  force_field)
  {
    ////RunWorklet::SumFarNearForce(EleNewForce(), NearForce(), _all_force);
    PreCompute();
    ComputeRBEEleForce_MPItest(_psample, _RBE_P, _ele_new_force);
    ComputeVerletlistNearForceMPITest(_nearforce);
    RunWorklet::SumFarNearForce(_ele_new_force, _nearforce, _all_force);

#ifdef PSAMPLE_RBEP
    /*输出循环力*/
    //int rank;
    //MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
    {
       std::ofstream _psample_MPI_0;
       _psample_MPI_0.open("_psample_MPI_" + std::to_string(num_cont) + "_0.txt");
       auto _psample_read = _psample.ReadPortal();
       for (size_t i = 0; i < _psample.GetNumberOfValues(); i++)
       {
         _psample_MPI_0 << i << "," << _psample_read.Get(i) << std::endl;
       }      

       std::cout << num_cont<<"——rank=0 的RBE_P为：" << _RBE_P<<std::endl;

    }
    else if (rank == 1)
    {
       std::ofstream _psample_MPI_1;
       _psample_MPI_1.open("_psample_MPI_" + std::to_string(num_cont) + "_1.txt");
       auto _psample_read = _psample.ReadPortal();
       for (size_t i = 0; i < _psample.GetNumberOfValues(); i++)
       {
         _psample_MPI_1 << i << "," << _psample_read.Get(i) << std::endl;
       } 
       std::cout << num_cont << "——rank=1 的RBE_P为：" << _RBE_P << std::endl;
    }
#endif
    //exit(1);

#ifdef FORCE
    /*输出循环力*/
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
    {
       std::ofstream ljforce_file0;
       ljforce_file0.open("force_MPI_" + std::to_string(num_cont) + "_0.txt");
       auto force_read = _all_force.ReadPortal();
       for (auto it = list_id_map_new.begin(); it != list_id_map_new.end(); ++it)
       {
         ljforce_file0 << it->first << "," << force_read.Get(it->first) << std::endl;
       }
    }
    else if (rank == 1)
    {
       std::ofstream ljforce_file1;
       ljforce_file1.open("force_MPI_" + std::to_string(num_cont) + "_1.txt");
       auto force_read = _all_force.ReadPortal();
       for (auto it = list_id_map_new.begin(); it != list_id_map_new.end(); ++it)
       {
         ljforce_file1 << it->first << "," << force_read.Get(it->first) << std::endl;
       }
    }
#endif
  }

  else if ("CVFF" == force_field)
  {
    RunWorklet::SumFarNearForce(EleNewForce(), NearForce(), _all_force);

    Invoker{}(MolecularWorklet::AddForceWorklet{}, SpecialCoulForce(), _all_force);

    //all force with bondforce
    Invoker{}(MolecularWorklet::AddForceWorklet{}, BondForce(), _all_force);

    //all force with angleforce
    Invoker{}(MolecularWorklet::AddForceWorklet{}, AngleForce(), _all_force);

    if (_para.GetParameter<bool>(PARA_DIHEDRALS_FORCE) &&
        _para.GetParameter<bool>(PARA_FILE_DIHEDRALS))
    {
      //all force with dihedral+force
      Invoker{}(MolecularWorklet::AddForceWorklet{}, DihedralsForce(), _all_force);
    }
  }

  else if ("EAM" == force_field)
  {
    NearForceEAM();
  }

}

void ExecutionNVT::UpdateVelocity()
{
  try
  {
    vtkm::cont::ArrayCopy(_velocity, _old_velocity);

    RunWorklet::UpdateVelocity(_dt, _unit_factor._fmt2v, _all_force, _mass, _velocity);
  }
  catch (const std::exception& e)
  {
    std::cout << e.what() << std::endl;
  }
}

void ExecutionNVT::UpdatePosition()
{
  //auto n = _position.GetNumberOfValues();
  // store old position
  //_old_position.Allocate(n);
  //for (int i = 0; i < n; i++)
  //{
  //  _old_position.WritePortal().Set(i, _position.ReadPortal().Get(i));
  //}

  vtkm::cont::ArrayCopy(_position, _old_position);

  if (_para.GetParameter<std::string>(PARA_FIX_SHAKE) == "null" || _init_way == "inbuild")
  {
    auto&& position_flag = _para.GetFieldAsArrayHandle<Id3>(field::position_flag);
    RunWorklet::UpdatePositionFlag(_dt, _velocity, _locator, _position, position_flag);
  }
  else
  {
    RunWorklet::UpdatePosition(_dt, _velocity, _locator, _position);
  }

  _locator.SetPosition(_position);
  SetCenterTargetPositions();
}

void ExecutionNVT::UpdateVelocityByTempConType()
{
    // 
  _temp_con_type = _para.GetParameter<std::string>(PARA_TEMP_CTRL_TYPE);
  if (_temp_con_type == "NOSE_HOOVER")
  {
    //Nose Hoover
    //Maybe tauT = 20.0 * dt is a good choice for dt = 2e-3 from the test.
    //Because the temperature curve is the smoothest of all test simulations.
    //In fact, 5.0 and 10.0 are also optional.
    //As long as the coefficent is not too large, such as larger than 100 * dt.
    RunWorklet::UpdateVelocityNoseHoover(
      _dt, _unit_factor._fmt2v, _nosehooverxi, _all_force, _mass, _velocity);
    //Real tauT = 20.0 * _dt;
    Real tauT = vtkm::Pow(10.0, -1) * _dt;
    _nosehooverxi += 0.5 * _dt * (_tempT / _kbT - 1.0) / tauT;
  }
  else if (_temp_con_type == "TEMP_RESCALE")
  {
    Real coeff_rescale = vtkm::Sqrt(_kbT / _tempT);
    RunWorklet::UpdateVelocityRescale(coeff_rescale, _velocity);
  }
  else if (_temp_con_type == "BERENDSEN")
  {
    //
    //Velocity Rescale: Berendsen
    //Maybe dt_divide_taut = 0.05 is a good choice for dt = 2e-3. 0.005, 0.01, 0.1 is optional.
    //The selection of dt_divide_taut determines the temperature equilibrium time.
    
    //Real dt_divide_taut = 0.02; for PEO
    //Real dt_divide_taut = 0.1; // 
    auto dt_divide_taut = _dt / _Tdamp;
    Real coeff_Berendsen = vtkm::Sqrt(1.0 + dt_divide_taut * (_kbT / _tempT - 1.0));
    RunWorklet::UpdateVelocityRescale(coeff_Berendsen, _velocity);    
  }
}

void ExecutionNVT::SetCenterTargetPositions()
{
  if (_init_way == "inbuild")
  {
    auto num_pos = _position.GetNumberOfValues();
    vtkm::cont::ArrayHandle<vtkm::Vec3f> center_position_temp;
    center_position_temp.Allocate(num_pos / 2);
    auto&& write_prot_center = center_position_temp.WritePortal();
    auto&& read_prot_center = _position.ReadPortal();

    vtkm::cont::ArrayHandle<vtkm::Vec3f> target_position_temp;
    target_position_temp.Allocate(num_pos / 2);
    auto&& write_prot_target = target_position_temp.WritePortal();
    auto&& read_prot_target = _position.ReadPortal();

    for (int i = 0; i < num_pos / 2; i++)
    {
      write_prot_center.Set(i, read_prot_center.Get(i));
      write_prot_target.Set(i, read_prot_target.Get(i + (num_pos / 2)));
    }

    auto center_position = _para.GetFieldAsArrayHandle<Vec3f>(field::center_position);
    vtkm::cont::ArrayCopy(center_position_temp, center_position);

    auto target_position = _para.GetFieldAsArrayHandle<Vec3f>(field::target_position);
    vtkm::cont::ArrayCopy(target_position_temp, target_position);
  }
  else
  {
    /*auto atom_id_center = _para.GetFieldAsArrayHandle<Id>(field::atom_id_center);
    auto atom_id_target = _para.GetFieldAsArrayHandle<Id>(field::atom_id_target);
    auto center_position = _para.GetFieldAsArrayHandle<Vec3f>(field::center_position);
    auto target_position = _para.GetFieldAsArrayHandle<Vec3f>(field::target_position);

    Invoker{}(MolecularWorklet::GetPositionByTypeWorklet{}, atom_id_center, _position, center_position);
    Invoker{}(MolecularWorklet::GetPositionByTypeWorklet{}, atom_id_target, _position, target_position);*/
    //std::map<Id, ArrayHandle<Vec3f>> atom_pair_position;
    //std::map<Id, ArrayHandle<Id>> atom_pair_id;
    //auto atom_pair_type = _para.GetParameter<std::vector<int>>(PARA_ATOMS_PAIR_TYPE);
    auto rdf_id = _para.GetFieldAsArrayHandle<Id>(field::atom_pair_id);
    auto atoms_pair_type_offsets = _para.GetFieldAsArrayHandle<Id>(field::atoms_pair_type_offsets);
    vtkm::cont::ArrayHandle<Vec3f> atom_type_position;
    atom_type_position.Allocate(rdf_id.GetNumberOfValues());
    Invoker{}(MolecularWorklet::GetPositionByTypeWorklet{}, rdf_id, _position, atom_type_position);
    auto atom_pair_position = _para.GetFieldAsArrayHandle<Vec3f>(field::atom_pair_position);
    vtkm::cont::ArrayCopy(atom_type_position, atom_pair_position);
  }
 }

void ExecutionNVT::PreForce()
{
  _Vlength = _para.GetParameter<Real>(PARA_VLENGTH);
  _box = rbmd_parallel->GetBoxLength()/*_para.GetParameter<Vec3f>(PARA_BOX)*/;
  Vec3f sigma = { static_cast<Real>(vtkm::Sqrt(_alpha / 2.0) * _box[0] / vtkm::Pi()),
                  static_cast<Real>(vtkm::Sqrt(_alpha / 2.0) * _box[1] / vtkm::Pi()),
                  static_cast<Real>(vtkm::Sqrt(_alpha / 2.0) * _box[2] / vtkm::Pi()) };
  _dt = _executioner->Dt();
  // prepare for RBE force
  auto velocity_type = _para.GetParameter<std::string>(gtest::velocity_type);
  auto random = (velocity_type != "TEST") ? true : false;
  RBEPSAMPLE rbe_presolve_psample = { _alpha, _Vlength, _box, _RBE_P };
  rbe_presolve_psample._RBE_random = random;
  //_psample = rbe_presolve_psample.Fetch_P_Sample(Real(0.0), sigma);

  std::vector<vtkm::Vec3f> result;
  rbmd_parallel->GetRbeRadomM_vec3f(result);  

  //for (size_t i = 0; i < result.size(); i++)
  //{
  //  std::cout << i << "   " << result[i] << std::endl;
  //}

  MPI_Barrier(MPI_COMM_WORLD);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(result), _psample);

  //auto _psample_read = _psample.ReadPortal();
  //for (size_t i = 0; i < _psample.GetNumberOfValues(); i++)
  //{
  //  std::cout << i << " --------------- " << _psample_read.Get(i) << std::endl;
  //}

  // prepare for Langevin dynamics
  RBEPSAMPLE sample_presolve_1d;
  sample_presolve_1d._RBE_random = random;
  Real gaussianx = sample_presolve_1d.FetchSample_1D(Real(0.0), Real(1.0));
  Real gaussiany = sample_presolve_1d.FetchSample_1D(Real(0.0), Real(1.0));
  Real gaussianz = sample_presolve_1d.FetchSample_1D(Real(0.0), Real(1.0));
  _gaussian = { gaussianx, gaussiany, gaussianz };

  //for (size_t i = 0; i < _psample.GetNumberOfValues(); i++)
  //{
  //  std::cout << i << " ------_gaussian后：--------- " << _psample_read.Get(i) << std::endl;
  //}

  //exit(1);
}

vtkm::cont::ArrayHandle<Vec3f> ExecutionNVT::LJForce()
{
  _nearforce_type = _para.GetParameter<std::string>(PARA_NEIGHBOR_TYPE);
  if (_nearforce_type == "RBL")
  {
    ComputeRBLLJForce(_LJforce);
  }
  if (_nearforce_type == "VERLETLIST")
  {
    ComputeVerletlistLJForce(_LJforce);
  }
  else if (_nearforce_type == "ORIGINAL")
  {
    //ComputeOriginalLJForce(_LJforce);
    ComputeSpecialBondsLJForce(_LJforce);
  }
  return _LJforce;
}

vtkm::cont::ArrayHandle<Vec3f> ExecutionNVT::EleNearForce()
{
  if (_use_erf == true)
  {
    RunWorklet::ComputeNearElectrostaticsERF(
      _atoms_id, _static_table, _locator, _topology, _force_function, _ele_near_force);
  }
  else
  {
    RunWorklet::ComputeNearElectrostatics(
      _atoms_id, _locator, _topology, _force_function, _ele_near_force);
  }
  Invoker{}(MolecularWorklet::UnitRescaleWorklet{ _unit_factor._qqr2e }, _ele_near_force);
  return _ele_near_force;
}

vtkm::cont::ArrayHandle<Vec3f> ExecutionNVT::NearForce()
{
  _nearforce_type = _para.GetParameter<std::string>(PARA_NEIGHBOR_TYPE);

  if (_nearforce_type == "RBL")
  {
    ComputeRBLNearForce(_nearforce);
  }
  else if (_nearforce_type == "VERLETLIST")
  {
    ComputeVerletlistNearForce(_nearforce);
  }
  else if (_nearforce_type == "ORIGINAL")
  {
    RunWorklet::SumFarNearForce(EleNearForce(), LJForce(), _nearforce);
  }
  return _nearforce;
}

vtkm::cont::ArrayHandle<Vec3f> ExecutionNVT::NearForceLJ()
{
  _nearforce_type = _para.GetParameter<std::string>(PARA_NEIGHBOR_TYPE);
  if (_nearforce_type == "RBL")
  {
    ComputeRBLLJForce(_all_force);
  }
  else if (_nearforce_type == "VERLETLIST")
  {
    ComputeVerletlistLJForce(_all_force);
  }
  else if (_nearforce_type == "ORIGINAL")
  {
    ComputeOriginalLJForce(_all_force);
  }
  return _all_force;
}

vtkm::cont::ArrayHandle<Vec3f> ExecutionNVT::NearForceEAM()
{
  _nearforce_type = _para.GetParameter<std::string>(PARA_NEIGHBOR_TYPE);
  if (_nearforce_type == "RBL")
  {
    ComputeRBLEAMForce(_all_force);
  }
  else if (_nearforce_type == "VERLETLIST")
  {
    ComputeVerletlistEAMForce(_all_force);
  }
  else if (_nearforce_type == "ORIGINAL")
  {
    ComputeOriginalEAMForce(_all_force);
  }
  return _all_force;
}

vtkm::cont::ArrayHandle<Vec3f> ExecutionNVT::BondForce()
{
  // original
  auto bondlist = _para.GetFieldAsArrayHandle<Id>(field::bond_atom_id);
  auto bond_type = _para.GetFieldAsArrayHandle<Id>(field::bond_type);
  vtkm::IdComponent bondlist_num = bondlist.GetNumberOfValues();
  auto&& bondlist_group = vtkm::cont::make_ArrayHandleGroupVec<2>(bondlist);
  auto bond_num = bondlist_group.GetNumberOfValues();


  // bond_coeffs_k   bond_coeffs_equilibrium
  auto bond_coeffs_k = _para.GetFieldAsArrayHandle<Real>(field::bond_coeffs_k);
  auto bond_coeffs_equilibrium = _para.GetFieldAsArrayHandle<Real>(field::bond_coeffs_equilibrium);

  // forcebond
  vtkm::cont::ArrayHandle<Vec3f> forcebond;
  vtkm::cont::ArrayHandle<Real> bond_energy;
  auto&& forcebond_group = vtkm::cont::make_ArrayHandleGroupVec<2>(forcebond);
  Invoker{}(MolecularWorklet::ComputeBondHarmonicWorklet{ _box },
            bond_type,
            bondlist_group,
            bond_coeffs_k,
            bond_coeffs_equilibrium,
            _position,
            forcebond_group,
            bond_energy,
            _locator);

  auto bond_energy_avr =
    vtkm::cont::Algorithm::Reduce(bond_energy, vtkm::TypeTraits<Real>::ZeroInitialization()) /
    bond_num;
  _para.SetParameter(PARA_BOND_ENERGY, bond_energy_avr);

  vtkm::cont::ArrayHandle<Vec3f> reduce_force_bond;
  if (_para.GetParameter<bool>(PARA_DIHEDRALS_FORCE))
  {
    //reduce bond force
    auto atom_id = _para.GetFieldAsArrayHandle<Id>(field::atom_id);
    auto atom_id_number = atom_id.GetNumberOfValues();
    auto original_number = bondlist.GetNumberOfValues();

    std::vector<Id> new_bondlist(original_number);
    memcpy(&new_bondlist[0], bondlist.ReadPortal().GetArray(), original_number * sizeof(vtkm::Id));

    std::vector<vtkm::Id> append_list(atom_id_number);
    memcpy(&append_list[0], atom_id.ReadPortal().GetArray(), atom_id_number * sizeof(vtkm::Id));

    //append key
    new_bondlist.insert(new_bondlist.end(), append_list.begin(), append_list.end());

    //append force bond
    std::vector<vtkm::Vec3f> new_forcebond(original_number);
    memcpy(
      &new_forcebond[0], forcebond.ReadPortal().GetArray(), original_number * sizeof(vtkm::Vec3f));

    std::vector<Vec3f> value(atom_id_number, vtkm::Vec3f{ 0.0, 0.0, 0.0 });
    new_forcebond.insert(new_forcebond.end(), value.begin(), value.end());

    vtkm::worklet::Keys<vtkm::Id> keys_bond(vtkm::cont::make_ArrayHandle(new_bondlist));
    Invoker{}(MolecularWorklet::ReduceForceWorklet{},
              keys_bond,
              vtkm::cont::make_ArrayHandle(new_forcebond),
              reduce_force_bond);
  }
  else
  {
    vtkm::worklet::Keys<vtkm::Id> keys_bond(bondlist);
    Invoker{}(MolecularWorklet::ReduceForceWorklet{}, keys_bond, forcebond, reduce_force_bond);
  }
  return reduce_force_bond;
}

vtkm::cont::ArrayHandle<Vec3f> ExecutionNVT::AngleForce()
{
  //angle
  auto angle_list = _para.GetFieldAsArrayHandle<Id>(field::angle_atom_id);
  auto angle_type = _para.GetFieldAsArrayHandle<Id>(field::angle_type);
  vtkm::IdComponent anglelist_num = angle_list.GetNumberOfValues();
  auto&& anglelist_group = vtkm::cont::make_ArrayHandleGroupVec<3>(angle_list);
  auto angle_num = anglelist_group.GetNumberOfValues();

  // angle_coeffs_k   angle_coeffs_equilibrium
  auto angle_coeffs_k = _para.GetFieldAsArrayHandle<Real>(field::angle_coeffs_k);
  auto angle_coeffs_equilibrium = _para.GetFieldAsArrayHandle<Real>(field::angle_coeffs_equilibrium);

  // force_angle
  vtkm::cont::ArrayHandle<Vec3f> force_angle;
  vtkm::cont::ArrayHandle<Real> angle_energy;
  auto&& forceangle_group = vtkm::cont::make_ArrayHandleGroupVec<3>(force_angle);
  Invoker{}(MolecularWorklet::ComputeAngleHarmonicWorklet{ _box },
            angle_type,
            anglelist_group,
            angle_coeffs_k,
            angle_coeffs_equilibrium,
            _position,
            forceangle_group,
            angle_energy,
            _locator);

  auto angle_energy_avr =
    vtkm::cont::Algorithm::Reduce(angle_energy, vtkm::TypeTraits<Real>::ZeroInitialization()) /
    angle_num;
  _para.SetParameter(PARA_ANGLE_ENERGY, angle_energy_avr);

  vtkm::cont::ArrayHandle<Vec3f> reduce_force_angle;
  if (_para.GetParameter<bool>(PARA_DIHEDRALS_FORCE))
  {
    //reduce angle force
    auto atom_id = _para.GetFieldAsArrayHandle<Id>(field::atom_id);
    auto atom_id_number = atom_id.GetNumberOfValues();
    auto original_number = angle_list.GetNumberOfValues();

    std::vector<Id> new_anglelist(original_number);
    memcpy(
      &new_anglelist[0], angle_list.ReadPortal().GetArray(), original_number * sizeof(vtkm::Id));

    std::vector<vtkm::Id> append_list(atom_id_number);
    memcpy(&append_list[0], atom_id.ReadPortal().GetArray(), atom_id_number * sizeof(vtkm::Id));

    //append key
    new_anglelist.insert(new_anglelist.end(), append_list.begin(), append_list.end());

    //append force bond
    std::vector<vtkm::Vec3f> new_force_angle(original_number);
    memcpy(&new_force_angle[0],
           force_angle.ReadPortal().GetArray(),
           original_number * sizeof(vtkm::Vec3f));
    std::vector<Vec3f> value(atom_id_number, vtkm::Vec3f{ 0.0, 0.0, 0.0 });
    new_force_angle.insert(new_force_angle.end(), value.begin(), value.end());

    vtkm::worklet::Keys<vtkm::Id> keys_angle(vtkm::cont::make_ArrayHandle(new_anglelist));
    Invoker{}(MolecularWorklet::ReduceForceWorklet{},
              keys_angle,
              vtkm::cont::make_ArrayHandle(new_force_angle),
              reduce_force_angle);
  }
  else
  {
    vtkm::worklet::Keys<vtkm::Id> keys_angle(angle_list);
    Invoker{}(MolecularWorklet::ReduceForceWorklet{}, keys_angle, force_angle, reduce_force_angle);
  }

  return reduce_force_angle;
}

vtkm::cont::ArrayHandle<Vec3f> ExecutionNVT::DihedralsForce()
{
  //dihedrals
  auto dihedrals_list = _para.GetFieldAsArrayHandle<Id>(field::dihedrals_atom_id);
  auto dihedrals_type = _para.GetFieldAsArrayHandle<Id>(field::dihedrals_type);
  vtkm::IdComponent dihedralslist_num = dihedrals_list.GetNumberOfValues();
  auto&& dihedralslist_group = vtkm::cont::make_ArrayHandleGroupVec<4>(dihedrals_list);
  auto dihedrals_num = dihedralslist_group.GetNumberOfValues();

  // dihedrals_coeffs_k   dihedrals_coeffs_sign     dihedrals_coeffs_multiplicity
  auto dihedrals_coeffs_k = _para.GetFieldAsArrayHandle<Real>(field::dihedrals_coeffs_k);
  auto dihedrals_coeffs_sign =
    _para.GetFieldAsArrayHandle<vtkm::IdComponent>(field::dihedrals_coeffs_sign);
  auto dihedrals_coeffs_multiplicity =
    _para.GetFieldAsArrayHandle<vtkm::IdComponent>(field::dihedrals_coeffs_multiplicity);

  // force_dihedrals
  vtkm::cont::ArrayHandle<Vec3f> force_dihedrals;
  vtkm::cont::ArrayHandle<Real> dihedrals_energy;
  auto&& forcedihedrals_group = vtkm::cont::make_ArrayHandleGroupVec<4>(force_dihedrals);

  //auto a1 = dihedrals_type.GetNumberOfValues();
  //auto a2 = dihedralslist_group.GetNumberOfValues();
  //auto a3 = dihedrals_coeffs_k.GetNumberOfValues();
  //auto a4 = dihedrals_coeffs_sign.GetNumberOfValues();
  //auto a5 = dihedrals_coeffs_multiplicity.GetNumberOfValues();
  Invoker{}(MolecularWorklet::ComputeDihedralHarmonicWorklet{ _box },
            dihedrals_type,
            dihedralslist_group,
            dihedrals_coeffs_k,
            dihedrals_coeffs_sign,
            dihedrals_coeffs_multiplicity,
            _position,
            forcedihedrals_group,
            dihedrals_energy,
            _locator);

  auto dihedrals_energy_avr =
    vtkm::cont::Algorithm::Reduce(dihedrals_energy, vtkm::TypeTraits<Real>::ZeroInitialization()) /
    dihedrals_num;
  _para.SetParameter(PARA_DIHEDRAL_ENERGY, dihedrals_energy_avr);

  vtkm::cont::ArrayHandle<Vec3f> reduce_force_dihedrals;
  //reduce dihedrals force
  auto atom_id = _para.GetFieldAsArrayHandle<Id>(field::atom_id);
  auto atom_id_number = atom_id.GetNumberOfValues();
  auto original_number = dihedrals_list.GetNumberOfValues();

  std::vector<Id> new_dihedralslist(original_number);
  memcpy(&new_dihedralslist[0],
         dihedrals_list.ReadPortal().GetArray(),
         original_number * sizeof(vtkm::Id));

  std::vector<vtkm::Id> append_list(atom_id_number);
  memcpy(&append_list[0], atom_id.ReadPortal().GetArray(), atom_id_number * sizeof(vtkm::Id));

  //append key
  new_dihedralslist.insert(new_dihedralslist.end(), append_list.begin(), append_list.end());

  //append force bond
  std::vector<vtkm::Vec3f> new_force_dihedrals(original_number);
  memcpy(&new_force_dihedrals[0],
         force_dihedrals.ReadPortal().GetArray(),
         original_number * sizeof(vtkm::Vec3f));
  std::vector<Vec3f> value(atom_id_number, vtkm::Vec3f{ 0.0, 0.0, 0.0 });
  new_force_dihedrals.insert(new_force_dihedrals.end(), value.begin(), value.end());

  vtkm::worklet::Keys<vtkm::Id> keys_dihedrals(vtkm::cont::make_ArrayHandle(new_dihedralslist));
  Invoker{}(MolecularWorklet::ReduceForceWorklet{},
            keys_dihedrals,
            vtkm::cont::make_ArrayHandle(new_force_dihedrals),
            reduce_force_dihedrals);

  return reduce_force_dihedrals;
}

vtkm::cont::ArrayHandle<Vec3f> ExecutionNVT::SpecialCoulForce()
{
  auto vLength = _para.GetParameter<Real>(PARA_VLENGTH);
  auto box = _para.GetParameter<Vec3f>(PARA_BOX);
  auto source_array = _para.GetFieldAsArrayHandle<Id>(field::special_source_array);
  auto offsets_array = _para.GetFieldAsArrayHandle<Id>(field::special_offsets_array);
  auto groupVecArray = vtkm::cont::make_ArrayHandleGroupVecVariable(source_array, offsets_array);

  if (_para.GetParameter<bool>(PARA_DIHEDRALS_FORCE))
  {
    auto special_offsets = _para.GetFieldAsArrayHandle<Id>(field::special_offsets);
    auto special_weights = _para.GetFieldAsArrayHandle<Real>(field::special_weights);
    auto specoal_ids = _para.GetFieldAsArrayHandle<Id>(field::special_ids);
    auto ids_group = vtkm::cont::make_ArrayHandleGroupVecVariable(specoal_ids, special_offsets);
    auto weight_group =
      vtkm::cont::make_ArrayHandleGroupVecVariable(special_weights, special_offsets);

    RunWorklet::ComputeSpecialCoulGeneral(box,
                                             _atoms_id,
                                             groupVecArray,
                                             _force_function,
                                             _topology,
                                             _locator,
                                             ids_group,
                                             weight_group,
                                             _spec_coul_force);
  }
  else
  {  
    RunWorklet::ComputeSpecialCoul(box, _atoms_id, groupVecArray, _force_function, _topology, _locator, _spec_coul_force);
  }
  return _spec_coul_force;
}

vtkm::cont::ArrayHandle<Vec3f> ExecutionNVT::EleNewForce()
{
  _farforce_type = _para.GetParameter<std::string>(PARA_COULOMB_TYPE);
  if (_farforce_type == "RBE")
  {
    // New RBE force part
    ComputeRBEEleForce(_psample, _RBE_P, _ele_new_force);
  }
  else if (_farforce_type == "EWALD")
  {
    // New Ewald far part
    ComputeEwaldEleForce(_Kmax, _ele_new_force);
  }

  Invoker{}(MolecularWorklet::UnitRescaleWorklet{ _unit_factor._qqr2e }, _ele_new_force);
  return _ele_new_force;
}

void ExecutionNVT::TempConTypeForce()
{
  vtkm::cont::ArrayHandle<Real> mass;
  mass.AllocateAndFill(_all_force.GetNumberOfValues(), 1.0);
  _temp_con_type = _para.GetParameter<std::string>(PARA_TEMP_CTRL_TYPE);
  if (_temp_con_type == "LANGEVIN")
  {
    // Underdamped Langevin
    Real kBT = 1.0;
    Real gamma = 100.0;
    RunWorklet::UnderdampedLangevin(_gaussian, kBT, gamma, _dt, mass, _velocity, _all_force);
  }
}

void ExecutionNVT::ComputeTempe()
{
  auto n = _position.GetNumberOfValues();
  ArrayHandle<Real> sq_velocity;
  sq_velocity.Allocate(n);

  RunWorklet::ComputerKineticEnergy(_velocity, _mass, sq_velocity);
  Invoker{}(MolecularWorklet::UnitRescaleWorklet{ _unit_factor._mvv2e }, sq_velocity);

  _tempT_sum = vtkm::cont::Algorithm::Reduce(sq_velocity, vtkm::TypeTraits<Real>::ZeroInitialization());

  Real _global_tempT_sum;
  MPI_Allreduce(&_tempT_sum, &_global_tempT_sum, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
  std::cout << "______" << _global_tempT_sum << std::endl;

  //////////////////////////////////////////////////////////////////////
  //dof_shake is important!!!
  //It is used only for shake option, which may be modified in the future.
  //When shake is called, dof_shake = n(number of atoms of water); otherwise, dof_shake = 0
  // 3 is the extra dof
  ///////////////////////////////////////////////////////////////////////
  auto shake = _para.GetParameter<std::string>(PARA_FIX_SHAKE);
  Real temperature_kB = _unit_factor._kB;
  if (_init_way == "inbuild") 
  {
    _tempT = 0.5 * _global_tempT_sum / (3 * n / 2.0);
  }
  else if (shake == "null") // peo
  {
    _tempT = 0.5 * _global_tempT_sum / ((3 * n - 3) * temperature_kB / 2.0);
  }
  else if (shake == "true")
  {
    _tempT = 0.5 * _global_tempT_sum / ((3 * n - n - 3) * temperature_kB / 2.0);
  }
  else if (shake == "false")
  {
    _tempT = 0.5 * _global_tempT_sum / ((3 * n) * temperature_kB / 2.0);
  }
  _para.SetParameter(PARA_TEMPT_SUM, _global_tempT_sum);
  _para.SetParameter(PARA_TEMPT, _tempT);
}

void ExecutionNVT::SetForceFunction()
{
  InitERF();
  if (_para.GetParameter<bool>(PARA_FAR_FORCE))
  {
    auto cut_off = _para.GetParameter<Real>(PARA_CUTOFF);
    auto volume = _para.GetParameter<Real>(PARA_VOLUME);
    auto vlength = _para.GetParameter<Real>(PARA_VLENGTH);
    auto box = _para.GetParameter<Vec3f>(PARA_BOX);
    _force_function.SetParameters(cut_off, _alpha, volume, vlength, box, _Kmax);
  }

  if (_para.GetParameter<std::string>(PARA_FILE_TYPE) == "EAM")
  {
    auto rhomax = _para.GetParameter<Real>(EAM_PARA_RHOMAX);
    auto nrho = _para.GetParameter<Id>(EAM_PARA_NRHO);
    auto drho = _para.GetParameter<Real>(EAM_PARA_DRHO);
    auto nr = _para.GetParameter<Id>(EAM_PARA_NR);
    auto dr = _para.GetParameter<Real>(EAM_PARA_DR);
    _force_function.SetEAMParameters(rhomax, nrho, drho, nr, dr);
  }

}

void ExecutionNVT::SetTopology()
{
  auto pts_type = _para.GetFieldAsArrayHandle<Id>(field::pts_type);
  auto molecule_id = _para.GetFieldAsArrayHandle<Id>(field::molecule_id);
  auto epsilon = _para.GetFieldAsArrayHandle<Real>(field::epsilon);
  auto sigma = _para.GetFieldAsArrayHandle<Real>(field::sigma);

  _topology.SetAtomsType(pts_type);
  _topology.SetMolecularId(molecule_id);
  _topology.SetEpsAndSigma(epsilon, sigma);

  if (_init_way != "inbuild")
  {
    auto source_array = _para.GetFieldAsArrayHandle<Id>(field::special_source_array);
    auto offsets_array = _para.GetFieldAsArrayHandle<Id>(field::special_offsets_array);
    _topology.SetSourceAndOffsets(source_array, offsets_array);
  }
}

void ExecutionNVT::InitParameters()
{
  ExecutionMD::InitParameters();
  if (_para.GetParameter<bool>(PARA_FAR_FORCE))
  {
    _RBE_P = _para.GetParameter<IdComponent>(PARA_COULOMB_SAMPLE_NUM);
    _alpha = _para.GetParameter<Real>(PARA_ALPHA);
    _Kmax = _para.GetParameter<IdComponent>(PARA_KMAX); 
  }
  _kbT = _para.GetParameter<std::vector<Real>>(PARA_TEMPERATURE)[0]; 
  _Tdamp = _para.GetParameter<std::vector<Real>>(PARA_TEMPERATURE)[2]; 
  _para.SetParameter(PARA_TEMPT_SUM, Real{ 0.0 });
  _para.SetParameter(PARA_TEMPT, Real{ 0.0 });
  _init_way = _para.GetParameter<std::string>(PARA_INIT_WAY);
}

void ExecutionNVT::TimeIntegration() {}

void ExecutionNVT::ConstraintA()
{
  //vtkm::cont::Timer timer4ConstraintA;
  //timer4ConstraintA.Start();
  auto angle_list = _para.GetFieldAsArrayHandle<Id>(field::angle_atom_id);
  auto&& anglelist_group = vtkm::cont::make_ArrayHandleGroupVec<3>(angle_list);

  auto data_range = _para.GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  Vec<Vec2f, 3> range{
    { static_cast<Real>(data_range[0].Min), static_cast<Real>(data_range[0].Max) },
    { static_cast<Real>(data_range[1].Min), static_cast<Real>(data_range[1].Max) },
    { static_cast<Real>(data_range[2].Min), static_cast<Real>(data_range[2].Max) }
  };

  auto&& position_flag = _para.GetFieldAsArrayHandle<Id3>(field::position_flag);

  Invoker{}(
    MolecularWorklet::NewConstraintAWaterBondAngleWorklet{
      _box, _dt, _unit_factor._fmt2v, range },
    anglelist_group,
    _old_position,
    _old_velocity,
    _all_force,
    _mass,
    position_flag,
    _locator);

  vtkm::cont::ArrayCopy(_old_position, _position);
  vtkm::cont::ArrayCopy(_old_velocity, _velocity);
}

void ExecutionNVT::ConstraintB()
{
  //vtkm::cont::Timer timer4ConstraintB;
  //timer4ConstraintB.Start();
  auto angle_list = _para.GetFieldAsArrayHandle<Id>(field::angle_atom_id);
  auto&& anglelist_group = vtkm::cont::make_ArrayHandleGroupVec<3>(angle_list);

  auto data_range = _para.GetParameter<vtkm::Vec<vtkm::Range, 3>>(PARA_RANGE);
  Vec<Vec2f, 3> range{
    { static_cast<Real>(data_range[0].Min), static_cast<Real>(data_range[0].Max) },
    { static_cast<Real>(data_range[1].Min), static_cast<Real>(data_range[1].Max) },
    { static_cast<Real>(data_range[2].Min), static_cast<Real>(data_range[2].Max) }
  };

  Invoker{}(
    MolecularWorklet::NewConstraintBWaterBondAngleWorklet{
      _Vlength, _dt, _unit_factor._fmt2v, range },
    anglelist_group,
    _position,
    _old_velocity,
    _all_force,
    _mass,
    _locator);



  vtkm::cont::ArrayCopy(_old_velocity, _velocity);
}

void ExecutionNVT::ReadPotentialFile(std::ifstream& input_file)
{
  // 
  if (!input_file.is_open())
  {
    std::cerr << "Unable to open the file." << std::endl;
  }

  // 
  for (int i = 0; i < 2; ++i)
  {
    input_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  // third row
  input_file >> file.nrho >> file.drho >> file.nr >> file.dr >> file.cut_off;

  //
  file.frho.resize(file.nrho + 1);
  file.zr.resize(file.nr + 1);
  file.rhor.resize(file.nrho + 1);

  //  frho  array

  for (int i = 0; i < file.nrho; ++i)
  {
    input_file >> file.frho[i];
  }

  //zr  array

  for (int i = 0; i < file.nr; ++i)
  {
    input_file >> file.zr[i];
  }

  //  rhor array

  for (int i = 0; i < file.nrho; ++i)
  {
    input_file >> file.rhor[i];
  }

  input_file.close();
}

void ExecutionNVT::AllocateEAM() {}

void ExecutionNVT::file2array()
{
  Id i, j, k, m, n;
  Real sixth = 1.0 / 6.0;
  // auto ntypes = _header._num_atoms_type;

  Real rmax;
  dr = drho = rmax = rhomax = 0.0;

  dr = MAX(dr, file.dr);
  drho = MAX(drho, file.drho);
  rmax = MAX(rmax, (file.nr - 1) * file.dr);
  rhomax = MAX(rhomax, (file.nrho - 1) * file.drho);

  // set nr,nrho from cutoff and spacings
  // 0.5 is for round-off in divide

  nr = static_cast<int>(rmax / dr + 0.5);
  nrho = static_cast<int>(rhomax / drho + 0.5);

  // ------------------------------------------------------------------
  // setup frho arrays
  // ------------------------------------------------------------------
  frho.resize(nrho + 1);

  Real r, p, cof1, cof2, cof3, cof4;
  for (m = 1; m <= nrho; m++)
  {
    r = (m - 1) * drho;
    p = r / file.drho + 1.0;
    k = static_cast<int>(p);
    k = MIN(k, file.nrho - 2);
    k = MAX(k, 2);
    p -= k;
    p = MIN(p, 2.0);
    cof1 = -sixth * p * (p - 1.0) * (p - 2.0);
    cof2 = 0.5 * (p * p - 1.0) * (p - 2.0);
    cof3 = -0.5 * p * (p + 1.0) * (p - 2.0);
    cof4 = sixth * p * (p * p - 1.0);
    frho[m] = cof1 * file.frho[k - 1] + cof2 * file.frho[k] + cof3 * file.frho[k + 1] +
      cof4 * file.frho[k + 2];
  }

  // ------------------------------------------------------------------
  // setup rhor arrays
  // ------------------------------------------------------------------
  rhor.resize(nrho + 1);
  for (m = 1; m <= nr; m++)
  {
    r = (m - 1) * dr;
    p = r / file.dr + 1.0;
    k = static_cast<int>(p);
    k = MIN(k, file.nr - 2);
    k = MAX(k, 2);
    p -= k;
    p = MIN(p, 2.0);
    auto cof1 = -sixth * p * (p - 1.0) * (p - 2.0);
    auto cof2 = 0.5 * (p * p - 1.0) * (p - 2.0);
    auto cof3 = -0.5 * p * (p + 1.0) * (p - 2.0);
    auto cof4 = sixth * p * (p * p - 1.0);
    rhor[m] = cof1 * file.rhor[k - 1] + cof2 * file.rhor[k] + cof3 * file.rhor[k + 1] +
      cof4 * file.rhor[k + 2];
  }

  // ------------------------------------------------------------------
  // setup z2r arrays
  // ------------------------------------------------------------------
  z2r.resize(nr + 1);

  double zri;
  for (m = 1; m <= nr; m++)
  {
    r = (m - 1) * dr;

    p = r / file.dr + 1.0;
    k = static_cast<int>(p);
    k = MIN(k, file.nr - 2);
    k = MAX(k, 2);
    p -= k;
    p = MIN(p, 2.0);
    cof1 = -sixth * p * (p - 1.0) * (p - 2.0);
    cof2 = 0.5 * (p * p - 1.0) * (p - 2.0);
    cof3 = -0.5 * p * (p + 1.0) * (p - 2.0);
    cof4 = sixth * p * (p * p - 1.0);
    zri = cof1 * file.zr[k - 1] + cof2 * file.zr[k] + cof3 * file.zr[k + 1] + cof4 * file.zr[k + 2];

    z2r[m] = 27.2 * 0.529 * zri * zri;
  }
}

void ExecutionNVT::interpolate(Id n, Real delta, std::vector<Real>& f, std::vector<Vec7f>& spline)
{
  for (int m = 1; m <= n; m++)
  {
    spline[m][6] = f[m];
  }

  spline[1][5] = spline[2][6] -
    spline[1][6]; //f'(x) = (f(x + h) - f(x)) / h    [5] is the coefficient of the first derivative(energy expression)
  spline[2][5] = 0.5 * (spline[3][6] - spline[1][6]);
  spline[n - 1][5] = 0.5 * (spline[n][6] - spline[n - 2][6]);
  spline[n][5] = spline[n][6] - spline[n - 1][6];

  for (int m = 3; m <= n - 2; m++)
  {
    spline[m][5] =
      ((spline[m - 2][6] - spline[m + 2][6]) + 8.0 * (spline[m + 1][6] - spline[m - 1][6])) /
      12.0; //further sample points for a more accurate estimate
  }

  for (int m = 1; m <= n - 1; m++)
  {
    spline[m][4] = 3.0 * (spline[m + 1][6] - spline[m][6]) - 2.0 * spline[m][5] -
      spline[m + 1][5]; //[4] is the coefficient of the second derivative
    spline[m][3] = spline[m][5] + spline[m + 1][5] -
      2.0 * (spline[m + 1][6] - spline[m][6]); // [3] is the coefficient of the third derivative
  }

  spline[n][4] = 0.0;
  spline[n][3] = 0.0;       //The second and third derivative coefficients at the last sample point are zero,
                             // To make the interpolation curve smoother at both ends, the higher derivative coefficient at the boundary can be set to zero.
                             // This is because spline interpolation typically uses higher-order polynomial interpolation 
                              // at the inner sample points and lower-order polynomials at the boundaries to ensure smoothness.

  for (int m = 1; m <= n; m++)
  {
    spline[m][2] = spline[m][5] / delta;       //The coefficient of the second derivative(force expression)  
    spline[m][1] = 2.0 * spline[m][4] / delta; //The coefficient of the first derivative
    spline[m][0] = 3.0 * spline[m][3] / delta; //The coefficient of the zero derivative (i.e. the value of the function).
  }
}

void ExecutionNVT::array2spline()
{
  frho_spline.resize(nrho + 1);
  rhor_spline.resize(nrho + 1);
  z2r_spline.resize(nr + 1);

  interpolate(nrho, drho, frho, frho_spline);
  interpolate(nr, dr, rhor, rhor_spline);
  interpolate(nr, dr, z2r, z2r_spline);
}

void ExecutionNVT::SetEAM()
{
  _para.SetParameter(EAM_PARA_CUTOFF, file.cut_off);
  _para.SetParameter(EAM_PARA_RHOMAX, rhomax);
  _para.SetParameter(EAM_PARA_NRHO, nrho);
  _para.SetParameter(EAM_PARA_DRHO, drho);
  _para.SetParameter(EAM_PARA_NR, nr);
  _para.SetParameter(EAM_PARA_DR, dr);
  //

  _para.AddField(field::rhor_spline, ArrayHandle<Vec7f>{});
  auto rhor_spline_get = _para.GetFieldAsArrayHandle<Vec7f>(field::rhor_spline);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle<Vec7f>(rhor_spline), rhor_spline_get);

  _para.AddField(field::frho_spline, ArrayHandle<Vec7f>{});
  auto frho_spline_get = _para.GetFieldAsArrayHandle<Vec7f>(field::frho_spline);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle<Vec7f>(frho_spline), frho_spline_get);

  _para.AddField(field::z2r_spline, ArrayHandle<Vec7f>{});
  auto z2r_spline_get = _para.GetFieldAsArrayHandle<Vec7f>(field::z2r_spline);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle<Vec7f>(z2r_spline), z2r_spline_get);
}

void ExecutionNVT::InitStyle()
{
  AllocateEAM();
  file2array();
  array2spline();
  SetEAM();
}

void ExecutionNVT::MPILinkCell()
{
  // Add to LinkCell
  //auto space = GetParameter<Vec3f>(PARA_SPACE);
  //auto cut_off = GetParameter<Real>(PARA_CUTOFF);
  //RBMDParallelUtil rbmd_parallel(
  //  { static_cast<double>(space[0]), static_cast<double>(space[1]), static_cast<double>(space[2]) },
  //  static_cast<double>(cut_off));

  auto read_positong = _position.ReadPortal();
  auto read_velocity = _velocity.ReadPortal();
  auto read_charge = _charge.ReadPortal();

  auto num_position = _position.GetNumberOfValues();

  for (size_t i = 0; i < num_position; i++)
  {
          //rbmd_parallel->AddAtom(i,
          //                      static_cast<double>(read_positong.Get(i)[0]),
          //                      static_cast<double>(read_positong.Get(i)[1]),
          //                      static_cast<double>(read_positong.Get(i)[2]);

    rbmd_parallel->AddAtomAndVelocity(i,
                                     static_cast<double>(read_positong.Get(i)[0]),
                                     static_cast<double>(read_positong.Get(i)[1]),
                                     static_cast<double>(read_positong.Get(i)[2]),
                                     static_cast<double>(read_velocity.Get(i)[0]),
                                     static_cast<double>(read_velocity.Get(i)[1]),
                                     static_cast<double>(read_velocity.Get(i)[2]),
                                     static_cast<double>(read_charge.Get(i)));

    // 注意：初始时也应该将速度传入 后面preupdate才能获取初始速度；
  }
}

void ExecutionNVT::PreCompute()
{
  //int rank;
  //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //if (rank == 0)
  //{
  //  auto _psample_read = _psample.ReadPortal();
  //  for (size_t i = 0; i < _psample.GetNumberOfValues(); i++)
  //  {
  //    std::cout << i << " ------分区后PrepareUpdate前：0--------- "
  //              << _psample_read.Get(i) << std::endl;
  //  }
  //}
  //else if (rank == 1)
  //{
  //  auto _psample_read = _psample.ReadPortal();
  //  for (size_t i = 0; i < _psample.GetNumberOfValues(); i++)
  //  {
  //    std::cout << i << " ------分区后PrepareUpdate前：1--------- "
  //              << _psample_read.Get(i) << std::endl;
  //  }
  //}  
  int step = _app.GetExecutioner()->CurrentStep();
  //rbmd_parallel->PrepareUpdateForPositionChange(); // 这里是做分区预处理，是已经分好区了？
  rbmd_parallel->PrepareUpdateForPositionChangeTest(step);
 
  //if (rank == 0)
  //{
  //  auto _psample_read = _psample.ReadPortal();
  //  for (size_t i = 0; i < _psample.GetNumberOfValues(); i++)
  //  {
  //    std::cout << i << " ------分区后PrepareUpdateForPositionChange0--------- " << _psample_read.Get(i) << std::endl;
  //  }
  //}
  //else if (rank == 1)
  //{
  //  auto _psample_read = _psample.ReadPortal();
  //  for (size_t i = 0; i < _psample.GetNumberOfValues(); i++)
  //  {
  //    std::cout << i << " ------分区后PrepareUpdateForPositionChange1--------- " << _psample_read.Get(i) << std::endl;
  //  }
  //}    
  MPI_Barrier(MPI_COMM_WORLD);
  list_id_map_new = rbmd_parallel->BuildVerletListFull();
  auto& cell_data = rbmd_parallel->GetCellData();

#ifdef POS_PRR_POS
  /*输出位置*/
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (42 < step && step < 51)
  {
    if (rank == 0)
    {
      std::ofstream POS_PRR_POS_MPI_10;
      POS_PRR_POS_MPI_10.open("POS_PRR_POS_MPI_1_" + std::to_string(step) + "_0.txt");
      auto read_position = _position.ReadPortal();
      for (auto it = list_id_map_new.begin(); it != list_id_map_new.end(); ++it)
      {
         POS_PRR_POS_MPI_10 << it->first << "," << read_position.Get(it->first) << std::endl;
      }
    }
    else
    {
      std::ofstream POS_PRR_POS_MPI_11;
      POS_PRR_POS_MPI_11.open("POS_PRR_POS_MPI_1_" + std::to_string(step) + "_1.txt");
      auto read_position = _position.ReadPortal();
      for (auto it = list_id_map_new.begin(); it != list_id_map_new.end(); ++it)
      {
         POS_PRR_POS_MPI_11 << it->first << "," << read_position.Get(it->first) << std::endl;
      }
    }
  }
#endif

  // 这里是 为了让后面的数据返回给 并行块代码做的准备；提前找到对应Id原子的 cell 和 atom /*开始的计算方式*/
  auto num_cell = cell_data.size();
  auto N = _position.GetNumberOfValues();
  std::vector<Id> atom_id_mpi_temp(N, -1);
  std::vector<Vec3f> position_temp(N, { 0, 0, 0 });
  std::vector<Vec3f> velocity_temp(N, { 0, 0, 0 });
  std::vector<Vec3f> charge_temp(N, 0);


  int n_test1 = 0;
  int n_test2 = 0;
  int n_test3 = 0;

  cell_atom_map.clear();
  for (size_t i = 0; i < num_cell; i++)
  {
    if (cell_data[i].IsHaloCell())
    {
      n_test3++;
      continue;
    }
    for (size_t j = 0; j < cell_data[i].GetAtoms().size(); j++)
    {
      n_test1++;
      auto id = cell_data[i].GetAtoms()[j].GetID();
      if (list_id_map_new.find(id) != list_id_map_new.end())
      {
        position_temp[id] = cell_data[i].GetAtoms()[j].Position3D();
        velocity_temp[id] = cell_data[i].GetAtoms()[j].Velocity3D();
        charge_temp[id] = cell_data[i].GetAtoms()[j]._charge;

        atom_id_mpi_temp[id] = id;
        cell_atom_map.insert(std::make_pair(id, std::vector<unsigned long>{ i, j }));
        n_test2++;
      }
    }
  }
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(atom_id_mpi_temp), _atoms_id);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(position_temp), _position);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(velocity_temp), _velocity);
  vtkm::cont::ArrayCopy(vtkm::cont::make_ArrayHandle(charge_temp), _charge);


  _locator.SetPosition(_position);
  MPI_Barrier(MPI_COMM_WORLD);

  //int rank;
  //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //
  //if (rank == 0)
  //{
  //  std::ofstream position_file0;
  //  position_file0.open("POS_PRR_POS_MPI_2_test" + std::to_string(num_cont) + "_0.txt");
  //  auto read_position = _position.ReadPortal();
  //  for (auto it = list_id_map_new.begin(); it != list_id_map_new.end(); ++it)
  //  {
  //    position_file0 << it->first << "," << read_position.Get(it->first) << std::endl;
  //  }
  //}
  //else
  //{
  //  std::ofstream position_file1;
  //  position_file1.open("POS_PRR_POS_MPI_2_test" + std::to_string(num_cont) + "_1.txt");
  //  auto read_position = _position.ReadPortal();
  //  for (auto it = list_id_map_new.begin(); it != list_id_map_new.end(); ++it)
  //  {
  //    position_file1 << it->first << "," << read_position.Get(it->first) << std::endl;
  //  }
  //}

#ifdef POS_PRR_POS
  /*输出位置*/
  if (42 < step && step < 51)
  {
    if (rank == 0)
    {
      std::ofstream POS_PRR_POS_MPI_20;
      POS_PRR_POS_MPI_20.open("POS_PRR_POS_MPI_2_" + std::to_string(step) + "_0.txt");
      auto read_position = _position.ReadPortal();
      for (auto it = list_id_map_new.begin(); it != list_id_map_new.end(); ++it)
      {
        POS_PRR_POS_MPI_20 << it->first << "," << read_position.Get(it->first) << std::endl;
      }
    }
    else
    {
      std::ofstream POS_PRR_POS_MPI_21;
      POS_PRR_POS_MPI_21.open("POS_PRR_POS_MPI_2_" + std::to_string(step) + "_1.txt");
      auto read_position = _position.ReadPortal();
      for (auto it = list_id_map_new.begin(); it != list_id_map_new.end(); ++it)
      {
        POS_PRR_POS_MPI_21 << it->first << "," << read_position.Get(it->first) << std::endl;
      }
    }
  }
#endif
  //exit(1);
}

void ExecutionNVT::PostCompute()
{
#ifdef POS_POST_POS
  auto step = _app.GetExecutioner()->CurrentStep();
  /*输出位置*/
  if (42 < step && step < 51)
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
      std::ofstream POS_POST_POS_MPI10;
      POS_POST_POS_MPI10.open("POS_POST_POS_MPI_1_" + std::to_string(step) + "_0.txt");
      auto read_position = _position.ReadPortal();
      for (auto it = list_id_map_new.begin(); it != list_id_map_new.end(); ++it)
      {
        POS_POST_POS_MPI10 << it->first << "," << read_position.Get(it->first) << std::endl;
      }
    }
    else
    {
      std::ofstream POS_POST_POS_MPI11;
      POS_POST_POS_MPI11.open("POS_POST_POS_MPI_1_" + std::to_string(step) + "_1.txt");
      auto read_position = _position.ReadPortal();
      for (auto it = list_id_map_new.begin(); it != list_id_map_new.end(); ++it)
      {
        POS_POST_POS_MPI11 << it->first << "," << read_position.Get(it->first) << std::endl;
      }
    }
  }
#endif
  rbmd_parallel->PostUpdateForPositionChange();
  auto read_position = _position.ReadPortal();
  auto read_velocity = _velocity.ReadPortal();
  auto read_charge = _charge.ReadPortal();

  auto& cell_data = rbmd_parallel->GetCellData();

  for (auto it = list_id_map_new.begin(); it != list_id_map_new.end(); ++it)
  {
    auto id = it->first;
    auto cell_id = cell_atom_map[id][0];
    auto atom_id = cell_atom_map[id][1];
    auto& atoms = cell_data[cell_id].GetAtoms();
    auto position = read_position.Get(id);
    auto velocity = read_velocity.Get(id);
    auto charge = read_charge.Get(id);

    atoms[atom_id]._r[0] = static_cast<double>(position[0]);
    atoms[atom_id]._r[1] = static_cast<double>(position[1]);
    atoms[atom_id]._r[2] = static_cast<double>(position[2]);

    atoms[atom_id]._v[0] = static_cast<double>(velocity[0]);
    atoms[atom_id]._v[1] = static_cast<double>(velocity[1]);
    atoms[atom_id]._v[2] = static_cast<double>(velocity[2]);

    atoms[atom_id]._charge = static_cast<double>(charge);
  }
#ifdef POS_POST_POS
  /*输出位置*/
  if (42 < step && step < 51)
  {
    auto step = _app.GetExecutioner()->CurrentStep();
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
      std::ofstream POS_POST_POS_MPI20;
      POS_POST_POS_MPI20.open("POS_POST_POS_MPI_2_" + std::to_string(step) + "_0.txt");
      auto read_position = _position.ReadPortal();
      for (auto it = list_id_map_new.begin(); it != list_id_map_new.end(); ++it)
      {
        POS_POST_POS_MPI20 << it->first << "," << read_position.Get(it->first) << std::endl;
      }
    }
    else
    {
      std::ofstream POS_POST_POS_MPI21;
      POS_POST_POS_MPI21.open("POS_POST_POS_MPI_2_" + std::to_string(step) + "_1.txt");
      auto read_position = _position.ReadPortal();
      for (auto it = list_id_map_new.begin(); it != list_id_map_new.end(); ++it)
      {
        POS_POST_POS_MPI21 << it->first << "," << read_position.Get(it->first) << std::endl;
      }
    }
  }
#endif
  //rbmd_parallel->PostUpdateForPositionChange();
  list_id_map_new.clear();
  MPI_Barrier(MPI_COMM_WORLD);


}