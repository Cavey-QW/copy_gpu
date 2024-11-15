﻿//==================================================================================
//  RBMD 1.0 is developed for random batch molecular dynamics calculation.
//
//  Copyright(C) 2024 SHANGHAI JIAOTONG UNIVERSITY CHONGQING RESEARCH INSTITUTE
//
//  This program is free software : you can redistribute it and /or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.If not, see < https://www.gnu.org/licenses/>.
//
//  The post-processing data produced by VASPKIT may be used in any publications 
//  provided that its use is explicitly acknowledged. A suitable reference for VASPKIT is:
//  [1] Gao W, Zhao T, Guo Y, et al.RBMD: A molecular dynamics package enabling to simulate 
//  10 million all - atom particles in a single graphics processing unit[J].arXiv preprint arXiv : 2407.09315, 2024.
// 
//  Contact Email : [support_wz@sciai.com.cn]
//==================================================================================

#pragma once
#include <vtkm/cont/ExecutionObjectBase.h>
#include "topology/ExecTopology.h"
#include "FieldName.h"
#include "Types.h"
#include <vtkm/cont/ConvertNumComponentsToOffsets.h>
#include <vtkm/cont/ArrayHandleGroupVecVariable.h>


class ContTopology : vtkm::cont::ExecutionObjectBase
{
public:
  using GroupVecType = typename vtkm::cont::ArrayHandleGroupVecVariable<vtkm::cont::ArrayHandle<vtkm::Id>,
                                                                        vtkm::cont::ArrayHandle<vtkm::Id>>;

  void SetSourceAndOffsets(const vtkm::cont::ArrayHandle<vtkm::Id>& source_array,
                           const vtkm::cont::ArrayHandle<vtkm::Id>& offsets_array)
  {
    this->_group_vec_array = GroupVecType(source_array, offsets_array);
  }
  void SetMolecularId(const vtkm::cont::ArrayHandle<Id>& molecular_id)
  {
    this->_molecular_id = molecular_id;
  }

  void SetAtomsType(const vtkm::cont::ArrayHandle<Id>& atoms_type)
  {
    this->_atoms_type = atoms_type;
  }
  void SetCharge(const vtkm::cont::ArrayHandle<Real>& charge) 
  { 
      this->_charge = charge;
  }

  void SetEpsAndSigma(const vtkm::cont::ArrayHandle<Real>& epsilon,
                      const vtkm::cont::ArrayHandle<Real>& sigma)
  {
    this->_epsilon = epsilon;
    this->_sigma = sigma;
  }

  //void SetNeighbourIdAndNum(const vtkm::cont::ArrayHandle<vtkm::Id>& neighbour_j_id,
  //                         const vtkm::cont::ArrayHandle<vtkm::Id>& neighbour_j_num)
  //{
  //  this->_group_vec_neighbour = GroupVecType(neighbour_j_id, neighbour_j_num);
  //}

  VTKM_CONT
  ExecTopology PrepareForExecution(vtkm::cont::DeviceAdapterId device,
                                        vtkm::cont::Token& token) const;

private:

  GroupVecType _group_vec_array;
  //GroupVecType _group_vec_neighbour;
  vtkm::cont::ArrayHandle<vtkm::Id> _molecular_id;
  vtkm::cont::ArrayHandle<vtkm::Id> _atoms_type;
  vtkm::cont::ArrayHandle<Real> _charge;
  vtkm::cont::ArrayHandle<Real> _epsilon;
  vtkm::cont::ArrayHandle<Real> _sigma;
};