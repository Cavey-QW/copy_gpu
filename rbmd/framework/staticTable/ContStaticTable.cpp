﻿//==================================================================================
//  RBMD 2.1.0 is developed for random batch molecular dynamics calculation.
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

#include "staticTable/ContStaticTable.h"
#include <vtkm/cont/ArrayHandle.h>

VTKM_CONT ExecStaticTable ContStaticTable::PrepareForExecution(vtkm::cont::DeviceAdapterId device,
                                       vtkm::cont::Token& token) const
{
  return ExecStaticTable(this->_table_index_1.PrepareForInput(device, token),
                           this->_table_rij_1.PrepareForInput(device, token),
                           this->_table_drij_1.PrepareForInput(device, token),
                           this->_table_function_rij_1.PrepareForInput(device, token),
                           this->_table_dfunction_rij_1.PrepareForInput(device, token),
                           this->_table_index_2.PrepareForInput(device, token),
                           this->_table_rij_2.PrepareForInput(device, token),
                           this->_table_drij_2.PrepareForInput(device, token),
                           this->_table_function_rij_2.PrepareForInput(device, token),
                           this->_table_dfunction_rij_2.PrepareForInput(device, token)
                           );
}