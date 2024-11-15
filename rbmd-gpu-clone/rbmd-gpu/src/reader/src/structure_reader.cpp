#include "../include/structure_reader.h"

#include <hip/hip_runtime.h>

#include <sstream>

#include "../Utilities/string_util.h"
#include "model/md_data.h"
#include "data_manager.h"

#define HIP_CHECK(call)                                              \
  {                                                                  \
    hipError_t err = call;                                           \
    if (err != hipSuccess) {                                         \
      std::cerr << "HIP error: " << hipGetErrorString(err) << " at " \
                << __FILE__ << ":" << __LINE__ << std::endl;         \
      exit(err);                                                     \
    }                                                                \
  }
StructureReder::StructureReder(const std::string& filePath, MDData& data)
    : MmapReader(filePath), _md_data(data) {}

int StructureReder::Execute() {
  if (-1 == MmapReader::Execute()) {
    // log
    return -1;
  }

  if (-1 == ReadHeader()) {
    // log
    return -1;
  }

  if (-1 == ReadForceField()) {
    // log
    return -1;
  }

  AllocateDataSpace();
  if (-1 == ReadData()) {
    // log
    return -1;
  }

  return 0;
}

int StructureReder::ReadHeader() {
  try {
    auto& info = _md_data._structure_info_data;
    HIP_CHECK(MALLOCHOST(&(info->_num_atoms), sizeof(rbmd::Id)));
    HIP_CHECK(MALLOCHOST(&(info->_num_bonds), sizeof(rbmd::Id)));
    HIP_CHECK(MALLOCHOST(&(info->_num_angles), sizeof(rbmd::Id)));
    HIP_CHECK(MALLOCHOST(&(info->_num_dihedrals), sizeof(rbmd::Id)));
    HIP_CHECK(MALLOCHOST(&(info->_num_impropers), sizeof(rbmd::Id)));
    HIP_CHECK(MALLOCHOST(&(info->_num_atoms_type), sizeof(rbmd::Id)));
    HIP_CHECK(MALLOCHOST(&(info->_num_bounds_type), sizeof(rbmd::Id)));
    HIP_CHECK(MALLOCHOST(&(info->_num_angles_type), sizeof(rbmd::Id)));
    HIP_CHECK(MALLOCHOST(&(info->_num_dihedrals_type), sizeof(rbmd::Id)));
    HIP_CHECK(MALLOCHOST(&(info->_range), sizeof(rbmd::Range)));
    auto& box = _md_data._h_box;
    rbmd::Real coord_min[3];
    rbmd::Real coord_max[3];

    for (; _locate < _file_size; ++_locate) {
      if (_mapped_memory[_locate] == '\n') {
        auto line = std::string(_line_start, &_mapped_memory[_locate]);
        std::istringstream iss(line);
        if (rbmd::IsLegalLine(line)) {
          if (line.find("atoms") != std::string::npos) {
            iss >> *(info->_num_atoms);
            // std::cout << *(info->_num_atoms) << " atoms" << std::endl;
          } else if (line.find("bonds") != std::string::npos) {
            iss >> *(info->_num_bonds);
            // std::cout << *(info->_num_bonds) << " bonds" << std::endl;
          } else if (line.find("angles") != std::string::npos) {
            iss >> *(info->_num_angles);
            // std::cout << *(info->_num_angles) << " angles" << std::endl;
          } else if (line.find("dihedrals") != std::string::npos) {
            iss >> *(info->_num_dihedrals);
             //std::cout << *(info->_num_dihedrals) << " dihedrals" << std::endl;
          } else if (line.find("impropers") != std::string::npos) {
            iss >> *(info->_num_impropers);
            //std::cout << *(info->_num_impropers) << " impropers" << std::endl;
          } else if (line.find("atom types") != std::string::npos) {
            iss >> *(info->_num_atoms_type);
            //std::cout << *(info->_num_atoms_type) << " atom types" << std::endl;
          } else if (line.find("bond types") != std::string::npos) {
            iss >> *(info->_num_bounds_type);
            //std::cout << *(info->_num_bounds_type) << " bond types" << std::endl;
          } else if (line.find("angle types") != std::string::npos) {
            iss >> *(info->_num_angles_type);
            // std::cout << *(info->_num_angles_type) << " angle types" << std::endl;
          }
          else if (line.find("dihedral types") != std::string::npos) {
              iss >> *(info->_num_dihedrals_type);
              //std::cout << *(info->_num_angles_type) << " angle types" << std::endl;
          } else if (line.find("xlo xhi") != std::string::npos) {
            iss >> coord_min[0] >> coord_max[0];
            // std::cout << coord_min[0] << " " << coord_max[0] << "xlo xhi" <<
            // std::endl;
            (*info->_range)[0][0] = coord_min[0];
            (*info->_range)[0][1] = coord_max[0];
          } else if (line.find("ylo yhi") != std::string::npos) {
            iss >> coord_min[1] >> coord_max[1];
            // std::cout << coord_min[1] << " " << coord_max[1] << "ylo yhi" <<
            // std::endl;
            (*info->_range)[1][0] = coord_min[1];
            (*info->_range)[1][1] = coord_max[1];
          } else if (line.find("zlo zhi") != std::string::npos) {
            iss >> coord_min[2] >> coord_max[2];
            // std::cout << coord_min[2] << " " << coord_max[2] << "zlo zhi" <<
            // std::endl;
            (*info->_range)[2][0] = coord_min[2];
            (*info->_range)[2][1] = coord_max[2];
            bool pbc[3] = {1, 1, 1};
            box->Init(box->_type, coord_min, coord_max, pbc);
            _line_start = &_mapped_memory[_locate];
            break;
          }
        }

        _line_start = &_mapped_memory[_locate];
      }
    }
  } catch (const std::exception& e) {
    // log
    return -1;
  }

  return 0;
}

int StructureReder::ReadForceField() {
  try {
    auto& info = _md_data._structure_info_data;
    for (; _locate < _file_size; ++_locate) {
      if (_mapped_memory[_locate] == '\n') {
        auto line = std::string(_line_start, &_mapped_memory[_locate]);
        std::istringstream iss(line);

        if (rbmd::IsLegalLine(line)) {
          if (line.find("Masses") != std::string::npos) {
            // std::cout << "Masses" << std::endl;
            ReadMass(*(info->_num_atoms_type));
          } else if (line.find("Pair Coeffs") != std::string::npos) {
            // std::cout << "Pair Coeffs" << std::endl;
            ReadPairCoeffs(*(info->_num_atoms_type));
          } else if (line.find("Bond Coeffs") != std::string::npos) {
            ReadBondCoeffs(*(info->_num_bounds_type));
            // std::cout << "Bond Coeffs" << std::endl;
          } else if (line.find("Angle Coeffs") != std::string::npos) {
            ReadAngleCoeffs(*(info->_num_angles_type));
            // std::cout << "Angle Coeffs" << std::endl;
          } else if (line.find("Dihedral Coeffs") != std::string::npos) {
            ReadDihedralsCoeffs(*(info->_num_dihedrals_type));
          } else if (line.find("group") != std::string::npos) {
            break;
          } else if (line.find("Atoms") != std::string::npos) {
            break;
          }
        }
      }
    }
  } catch (const std::exception& e) {
    // log
    return -1;
  }

  return 0;
}

int StructureReder::ReadMass(const rbmd::Id& numAtomTypes) {
    auto force_style = DataManager::getInstance().getConfigData()->Get<std::string>( "type", "hyper_parameters", "force_field");
    if ("CVFF" == force_style) {
        try {
            auto force_filed =
                std::dynamic_pointer_cast<CVFFForceFieldData>(_md_data._force_field_data);
            auto& mass = force_filed->_h_mass;
            HIP_CHECK(MALLOCHOST(&mass, numAtomTypes * sizeof(rbmd::Real)));
            rbmd::Id atom_type;
            rbmd::Real value;

            _line_start = &_mapped_memory[_locate];
            for (auto num = 0; _locate < _file_size && num < numAtomTypes; ++_locate) {
                if (_mapped_memory[_locate] == '\n') {
                    auto line = std::string(_line_start, &_mapped_memory[_locate]);
                    std::istringstream iss(line);
                    if (rbmd::IsLegalLine(line)) {
                        iss >> atom_type >> value;
                        mass[atom_type - 1] = value;
                        // std::cout << atom_type << " " << force_filed->_h_mass[atom_type -
                        // 1] << std::endl;
                        ++num;
                    }

                    _line_start = &_mapped_memory[_locate];
                }
            }
            auto& mass_1 = force_filed->_h_mass;
            //std::cout << "mass[0]=" << force_filed->_h_mass[0] << ","<< "mass[1]=" << force_filed->_h_mass[1] << std::endl;
        }
        catch (const std::exception& e) {
            // log
            return -1;
        }
    }else
    {
        try {
            auto force_filed =
                std::dynamic_pointer_cast<LJForceFieldData>(_md_data._force_field_data);
            auto& mass = force_filed->_h_mass;
            HIP_CHECK(MALLOCHOST(&mass, numAtomTypes * sizeof(rbmd::Id)));
            rbmd::Id atom_type;
            rbmd::Real value;

            _line_start = &_mapped_memory[_locate];
            for (auto num = 0; _locate < _file_size && num < numAtomTypes; ++_locate) {
                if (_mapped_memory[_locate] == '\n') {
                    auto line = std::string(_line_start, &_mapped_memory[_locate]);
                    std::istringstream iss(line);
                    if (rbmd::IsLegalLine(line)) {
                        iss >> atom_type >> value;
                        mass[atom_type - 1] = value;
                        // std::cout << atom_type << " " << force_filed->_h_mass[atom_type -
                        // 1] << std::endl;
                        ++num;
                    }

                    _line_start = &_mapped_memory[_locate];
                }
            }
            auto& mass_1 = force_filed->_h_mass;
            /*std::cout << "mass[0]=" << force_filed->_h_mass[0] << ","
                << "mass[1]=" << force_filed->_h_mass[1] << std::endl;*/
        }
        catch (const std::exception& e) {
            // log
            return -1;
        }
    }
  

  return 0;
}

int StructureReder::ReadPairCoeffs(const rbmd::Id& numAtomTypes) {
    auto force_style = DataManager::getInstance().getConfigData()->Get<std::string>("type", "hyper_parameters", "force_field");
    if ("CVFF" == force_style) {
        try {
            auto force_filed =
                std::dynamic_pointer_cast<CVFFForceFieldData>(_md_data._force_field_data);
            auto& eps = force_filed->_h_eps;
            auto& sigma = force_filed->_h_sigma;
            HIP_CHECK(MALLOCHOST(&eps, numAtomTypes * sizeof(rbmd::Id)));
            HIP_CHECK(MALLOCHOST(&sigma, numAtomTypes * sizeof(rbmd::Id)));
            rbmd::Id atom_type;
            rbmd::Real eps_value;
            rbmd::Real sigma_value;

            _line_start = &_mapped_memory[_locate];
            for (auto num = 0; _locate < _file_size && num < numAtomTypes; ++_locate) {
                if (_mapped_memory[_locate] == '\n') {
                    auto line = std::string(_line_start, &_mapped_memory[_locate]);
                    std::istringstream iss(line);
                    if (rbmd::IsLegalLine(line)) {
                        iss >> atom_type >> eps_value >> sigma_value;
                        eps[atom_type - 1] = eps_value;
                        sigma[atom_type - 1] = sigma_value;
                        // std::cout << atom_type << " " << force_filed->_h_eps[atom_type - 1]
                        // << " " << force_filed->_h_sigma[atom_type - 1] << std::endl;
                        ++num;
                    }
                    _line_start = &_mapped_memory[_locate];
                }
            }
        }
        catch (const std::exception& e) {
            // log
            return -1;
        }
    }
    else {
        try {
            auto force_filed =
                std::dynamic_pointer_cast<LJForceFieldData>(_md_data._force_field_data);
            auto& eps = force_filed->_h_eps;
            auto& sigma = force_filed->_h_sigma;
            HIP_CHECK(MALLOCHOST(&eps, numAtomTypes * sizeof(rbmd::Id)));
            HIP_CHECK(MALLOCHOST(&sigma, numAtomTypes * sizeof(rbmd::Id)));
            rbmd::Id atom_type;
            rbmd::Real eps_value;
            rbmd::Real sigma_value;

            _line_start = &_mapped_memory[_locate];
            for (auto num = 0; _locate < _file_size && num < numAtomTypes; ++_locate) {
                if (_mapped_memory[_locate] == '\n') {
                    auto line = std::string(_line_start, &_mapped_memory[_locate]);
                    std::istringstream iss(line);
                    if (rbmd::IsLegalLine(line)) {
                        iss >> atom_type >> eps_value >> sigma_value;
                        eps[atom_type - 1] = eps_value;
                        sigma[atom_type - 1] = sigma_value;
                        // std::cout << atom_type << " " << force_filed->_h_eps[atom_type - 1]
                        // << " " << force_filed->_h_sigma[atom_type - 1] << std::endl;
                        ++num;
                    }
                    _line_start = &_mapped_memory[_locate];
                }
            }
        }
        catch (const std::exception& e) {
            // log
            return -1;
        }
    }
  

  return 0;
}

int StructureReder::ReadBondCoeffs(const rbmd::Id& numBondTypes) {
  try {
      auto force_filed = std::dynamic_pointer_cast<CVFFForceFieldData>(_md_data._force_field_data);
      auto& bond_coeffs_k = force_filed->_h_bond_coeffs_k;
      auto& bond_coeffs_equilibrium = force_filed->_h_bond_coeffs_equilibrium;
      HIP_CHECK(MALLOCHOST(&bond_coeffs_k, numBondTypes * sizeof(rbmd::Real)));
      HIP_CHECK(MALLOCHOST(&bond_coeffs_equilibrium, numBondTypes * sizeof(rbmd::Real)));
      rbmd::Id bound_type;
      rbmd::Real bond_coeffs_k_value;
      rbmd::Real equilibrium_value;
      
      _line_start = &_mapped_memory[_locate];
      for (auto num = 0; _locate < _file_size && num < numBondTypes; ++_locate)
      {
          if (_mapped_memory[_locate] == '\n')
          {
              auto line = std::string(_line_start,&_mapped_memory[_locate]); 
              std::istringstream iss(line); 
              if(rbmd::IsLegalLine(line))
              {
                  iss >> bound_type >> bond_coeffs_k_value >> equilibrium_value;
                  //std::cout << bound_type << " " <<bond_coeffs_k_value << " " << equilibrium_value << std::endl;
                  bond_coeffs_k[bound_type - 1] = bond_coeffs_k_value;
                  bond_coeffs_equilibrium[bound_type - 1] =equilibrium_value;
                  ++num;
               }
                  _line_start = &_mapped_memory[_locate];
           }
      }

  } catch (const std::exception& e) {
    // log
    return -1;
  }

  return 0;
}

int StructureReder::ReadAngleCoeffs(const rbmd::Id& numAngleTypes)
{
  try {
      auto force_filed = std::dynamic_pointer_cast<CVFFForceFieldData>(_md_data._force_field_data);
      auto& angle_coeffs_k = force_filed->_h_angle_coeffs_k;
      auto& angle_coeffs_equilibrium = force_filed->_h_angle_coeffs_equilibrium;
      HIP_CHECK(MALLOCHOST(&angle_coeffs_k, numAngleTypes * sizeof(rbmd::Real)));
      HIP_CHECK(MALLOCHOST(&angle_coeffs_equilibrium, numAngleTypes * sizeof(rbmd::Real)));
      rbmd::Id angle_type;
      rbmd::Real angle_coeffs_k_value;
      rbmd::Real equilibrium_value;
      
      _line_start = &_mapped_memory[_locate];
      for (auto num = 0; _locate < _file_size && num < numAngleTypes; ++_locate)
      {
         if (_mapped_memory[_locate] == '\n')
         {
             auto line = std::string(_line_start,&_mapped_memory[_locate]);
             std::istringstream iss(line); 
             if(rbmd::IsLegalLine(line))
             {
                 iss >> angle_type >> angle_coeffs_k_value >>equilibrium_value;
                 //std::cout << angle_type << " " <<angle_coeffs_k_value << " " << equilibrium_value << std::endl;
                 angle_coeffs_k[angle_type - 1] = angle_coeffs_k_value; 
                 angle_coeffs_equilibrium[angle_type - 1] = equilibrium_value;
                 ++num;
             }
             _line_start = &_mapped_memory[_locate];
         }
      }
  } catch (const std::exception& e) {
    // log
    return -1;
  }

  return 0;
}

int StructureReder::ReadDihedralsCoeffs(const rbmd::Id& numDihedralsTypes)
{
    try {
        auto force_filed = std::dynamic_pointer_cast<CVFFForceFieldData>(_md_data._force_field_data);
        auto& dihedral_coeffs_k = force_filed->_h_dihedral_coeffs_k;
        auto& dihedral_coeffs_sign = force_filed->_h_dihedral_coeffs_sign;
        auto& dihedral_coeffs_multiplicity = force_filed->_h_dihedral_coeffs_multiplicity;
        HIP_CHECK(MALLOCHOST(&dihedral_coeffs_k, numDihedralsTypes * sizeof(rbmd::Real)));
        HIP_CHECK(MALLOCHOST(&dihedral_coeffs_sign, numDihedralsTypes * sizeof(rbmd::Real)));
        HIP_CHECK(MALLOCHOST(&dihedral_coeffs_multiplicity, numDihedralsTypes * sizeof(rbmd::Real)));
        rbmd::Id dihedral_type;
        rbmd::Real dihedral_coeffs_k_value;
        rbmd::Real dihedral_coeffs_sign_value;
        rbmd::Real dihedral_coeffs_multiplicity_value;

        _line_start = &_mapped_memory[_locate];
        for (auto num = 0; _locate < _file_size && num < numDihedralsTypes; ++_locate)
        {
            if (_mapped_memory[_locate] == '\n')
            {
                auto line = std::string(_line_start, &_mapped_memory[_locate]); std::istringstream iss(line);
                if (rbmd::IsLegalLine(line))
                {
                    iss >> dihedral_type >> dihedral_coeffs_k_value >> dihedral_coeffs_sign_value >> dihedral_coeffs_multiplicity_value;
                    dihedral_coeffs_k[dihedral_type - 1] = dihedral_coeffs_k_value;
                    dihedral_coeffs_sign[dihedral_type - 1] = dihedral_coeffs_sign_value;
                    dihedral_coeffs_multiplicity[dihedral_type - 1] = dihedral_coeffs_multiplicity_value;
                    //std::cout << dihedral_type << " " << dihedral_coeffs_k_value << " " << dihedral_coeffs_sign_value << " " << dihedral_coeffs_multiplicity_value << std::endl;
                    ++num;
                }
                _line_start = &_mapped_memory[_locate];
            }
        }
    }
    catch (const std::exception& e) {
        // log
        return -1;
    }

    return 0;
}

