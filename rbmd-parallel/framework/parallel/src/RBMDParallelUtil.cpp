#include "RBMDParallelUtil.h"

#include <mpi.h>

#include <csignal>
#include <random>

#include "DomainDecompositionMPI.h"
#include "DomainServiceLocator.h"
#define PSAMPLE_RBEP
#include <fstream>
RBMDParallelUtil::RBMDParallelUtil(std::array<double, 3> min,std::array<double, 3> max, double cutoff){
  int current_rank, total_rank;

  MPI_Comm_size(MPI_COMM_WORLD, &total_rank);
  MPI_Comm_rank(MPI_COMM_WORLD, &current_rank);

  domdec_ = new DomainDecompositionMPI(current_rank, total_rank);
  DomainServiceLocator::Instance().SetDomain(new Domain(current_rank));
  auto domain = DomainServiceLocator::Instance().GetDomain();
  domain->SetGlobaMinMax(min,max);
  linked_cell_ = new LinkedCell();
  linked_cell_->SetCutoff(cutoff);
  double bBoxMin[3];
  double bBoxMax[3];
  domdec_->GetBoundingBoxMinMax(bBoxMin, bBoxMax);
  linked_cell_->rebuild(bBoxMin, bBoxMax);
  // 无法解决重复的问题
  // auto send_particles_together =  linked_cell_->rebuild(bBoxMin, bBoxMax);
  // domdec_->SetSendLeavingWithCopiesStrategy(send_particles_together);
  std::cout << "Bounding box in domian: " << current_rank << " : " << "["
            << bBoxMin[0] << ", " << bBoxMax[0] << "]" << " x " << "["
            << bBoxMin[1] << ", " << bBoxMax[1] << "]" << " x " << "["
            << bBoxMin[2] << ", " << bBoxMax[2] << "]" << std::endl;

  std::cout << "Halo Bounding box in domian: " << current_rank << " : " << "["
            << linked_cell_->_haloBoundingBoxMin[0] << ", " << linked_cell_->_haloBoundingBoxMax[0] << "]" << " x " << "["
            << linked_cell_->_haloBoundingBoxMin[1] << ", " << linked_cell_->_haloBoundingBoxMax[1] << "]" << " x " << "["
            << linked_cell_->_haloBoundingBoxMin[2] << ", " << linked_cell_->_haloBoundingBoxMax[2] << "]" << std::endl;

  domdec_->InitCommunicationPartners(cutoff, linked_cell_);
  this->computeOffsets();
  this->GetOffset();
}

RBMDParallelUtil::~RBMDParallelUtil() {
    delete domdec_;
    delete linked_cell_;
    domdec_ = nullptr;
    linked_cell_ = nullptr;
    // MPI_Finalize();
}

void RBMDParallelUtil::PrepareUpdateForPositionChange() {
    // TODO 暂时没有做力的交换
    linked_cell_->update();
    domdec_->BalanceAndExchange(1, false, linked_cell_);
}

void RBMDParallelUtil::PrepareUpdateForPositionChangeTest(int current_step) {
    step = current_step;
    //if (step > 42) {
    //    printAtomInfo(51, linked_cell_, "update之前");
    //    //printAtomInfo(51, linked_cell_, "update之前");
    //}
    linked_cell_->update();
    //if (step > 42) {
    //    printAtomInfo(51, linked_cell_, "update之后 Balance之前");
    //    //printAtomInfo(944, linked_cell_, "update之后 Balance之前");
    //}
    domdec_->BalanceAndExchange(1, false, linked_cell_);
    //if (step > 42)
    //{
    //    printAtomInfo(51, linked_cell_, "Balance之后");
    //    //printAtomInfo(944, linked_cell_, "Balance之后");
    //}
}

void RBMDParallelUtil::PostUpdateForPositionChange() {
    //if (step > 20) {
    //    printAtomInfo(500, linked_cell_, "Post之前");
    //    printAtomInfo(944, linked_cell_, "Post之前");
    //}
    linked_cell_->DeleteOuterParticles();
    //if (step > 20) {
    //    printAtomInfo(500, linked_cell_, "Post之后 ");
    //    printAtomInfo(944, linked_cell_, "Post之后 ");
    //}
}

void RBMDParallelUtil::AddAtom(int id, double x, double y, double z) {
    auto temp_atom = new Atom(id, x, y, z);
    double loc[3]{x, y, z};
    if (this->linked_cell_->isInBoundingBox(loc)) {
        this->linked_cell_->AddParticle(*temp_atom);
        // printf("process %d insert \n", myid);
    } else {
        delete temp_atom;
    }
}

void RBMDParallelUtil::AddAtomAndVelocity(int id,
                                          double x,
                                          double y,
                                          double z,
                                          double vx,
                                          double vy,
                                          double vz,
                                          double charge)
{
    auto temp_atom = new Atom(id, x, y, z, vx, vy, vz, charge);
    double loc[3]{ x, y, z };
    if (this->linked_cell_->isInBoundingBox(loc))
    {
        this->linked_cell_->AddParticle(*temp_atom);
        // printf("process %d insert \n", myid);
    }
    else
    {
        delete temp_atom;
    }
}

std::vector<Cell> &RBMDParallelUtil::GetCellData() {
    return this->linked_cell_->GetCells();
}

void RBMDParallelUtil::GetOffset() {
    int forwardNeighbourIndex = 0, backwardNeighbourIndex = 0;

    long maxNeighbourOffset = 0;
    long minNeighbourOffset = 0;

    std::array<long, 3> dims = {
        this->linked_cell_->GetCellPerDim(0),
        this->linked_cell_->GetCellPerDim(1),
        this->linked_cell_->GetCellPerDim(2)
    };

    std::array<long, 3> r{};
    for (r[2] = -1; r[2] <= 1; r[2]++) {
        for (r[1] = -1; r[1] <= 1; r[1]++) {
            for (r[0] = -1; r[0] <= 1; r[0]++) {
                long offset = ThreeToOneD(r, dims);

                if (offset > 0) {
                    forward_neighbour_offsets_[forwardNeighbourIndex] = offset;
                    ++forwardNeighbourIndex;
                    if (offset > maxNeighbourOffset) {
                        maxNeighbourOffset = offset;
                    }
                }
                if (offset < 0) {
                    backward_neighbour_offsets_[backwardNeighbourIndex] =
                            std::abs(offset);
                    ++backwardNeighbourIndex;
                    if (std::abs(offset) > minNeighbourOffset) {
                        minNeighbourOffset = std::abs(offset);
                    }
                }
            }
        }
    }
}

void RBMDParallelUtil::computeOffsets() {
    using std::make_pair;

    std::array<long, 3> dims = {
        this->linked_cell_->GetCellPerDim(0),
        this->linked_cell_->GetCellPerDim(1),
        this->linked_cell_->GetCellPerDim(2)
    };
    long int o = ThreeToOneD(0l, 0l, 0l, dims); // origin
    long int x = ThreeToOneD(1l, 0l, 0l, dims); // displacement to the right
    long int y = ThreeToOneD(0l, 1l, 0l, dims); // displacement ...
    long int z = ThreeToOneD(0l, 0l, 1l, dims);
    long int xy = ThreeToOneD(1l, 1l, 0l, dims);
    long int yz = ThreeToOneD(0l, 1l, 1l, dims);
    long int xz = ThreeToOneD(1l, 0l, 1l, dims);
    long int xyz = ThreeToOneD(1l, 1l, 1l, dims);

    int i = 0;
    // if incrementing along X, the following order will be more cache-efficient:
    _cellPairOffsets8Pack[i++] = std::make_pair(o, o);
    _cellPairOffsets8Pack[i++] = std::make_pair(o, y);
    _cellPairOffsets8Pack[i++] = std::make_pair(y, z);
    _cellPairOffsets8Pack[i++] = std::make_pair(o, z);
    _cellPairOffsets8Pack[i++] = std::make_pair(o, yz);

    _cellPairOffsets8Pack[i++] = std::make_pair(x, yz);
    _cellPairOffsets8Pack[i++] = std::make_pair(x, y);
    _cellPairOffsets8Pack[i++] = std::make_pair(x, z);
    _cellPairOffsets8Pack[i++] = std::make_pair(o, x);
    _cellPairOffsets8Pack[i++] = std::make_pair(o, xy);
    _cellPairOffsets8Pack[i++] = std::make_pair(xy, z);
    _cellPairOffsets8Pack[i++] = std::make_pair(y, xz);
    _cellPairOffsets8Pack[i++] = std::make_pair(o, xz);
    _cellPairOffsets8Pack[i++] = std::make_pair(o, xyz);

    i = 0;
    _cellOffsets8Pack[i++] = o;
    _cellOffsets8Pack[i++] = y;
    _cellOffsets8Pack[i++] = z;
    _cellOffsets8Pack[i++] = yz;

    _cellOffsets8Pack[i++] = x;
    _cellOffsets8Pack[i++] = xy;
    _cellOffsets8Pack[i++] = xz;
    _cellOffsets8Pack[i++] = xyz;
}

unsigned int RBMDParallelUtil::GetNumAtoms(const bool include_halo) {
    unsigned long N = 0;
    auto cells = this->GetCellData();
    unsigned long numCells = cells.size();

    for (unsigned long i = 0; i < numCells; ++i) {
        if (include_halo)
            N += cells.at(i).GetAtomsCount();
        else if (not cells.at(i).IsHaloCell()) {
            N += cells.at(i).GetAtomsCount();
        }
    }
    return N;
}

unsigned int RBMDParallelUtil::GetRbeGlobalRandomSeed() {
    unsigned int random_seed = 0u;
    auto domain = DomainServiceLocator::Instance().GetDomain();
    if (domdec_->GetCurrentRank() == 0) {
        random_seed = rd();
    }
    //    domdec_->CommunicatorBroadcast() 未实现uint 暂时先不用这个
    MPI_Bcast(&random_seed, 1, MPI_UNSIGNED, 0, domdec_->GetMPIComm());
    return random_seed;
}

void RBMDParallelUtil::SetRbeConfigPnumber(int rbe_p_number) {
    DomainServiceLocator::Instance().GetDomain()->SetGlobalRbePnumber(
        rbe_p_number);
}

int RBMDParallelUtil::GetRbeConfigPnumber() {
    return DomainServiceLocator::Instance().GetDomain()->GetGlobalRbePnumber();
}

std::vector<vtkm::Vec2f> RBMDParallelUtil::GetGlobalRho(const std::vector<vtkm::Vec2f>& local_rho)
{
    // TODO:后续改为返回引用吧
    int size = GetRbeConfigPnumber();
    int num_elements = local_rho.size() * 2;
    std::vector<float> local_elements(num_elements);
    std::vector<vtkm::Vec2f> global_rho(local_rho.size());
    std::vector<float> global_elements(num_elements);
    // 将local_data展开到local_elements
    for (size_t i = 0; i < local_rho.size(); ++i) {
      local_elements[i * 2] = local_rho[i][0];
      local_elements[i * 2 + 1] = local_rho[i][1];
    }
    MPI_Allreduce(local_elements.data(), global_elements.data(), num_elements, MPI_FLOAT, MPI_SUM, domdec_->GetMPIComm());
    // 将global_elements重新组装到global_rho
    for (size_t i = 0; i < global_rho.size(); ++i)
    {
      global_rho[i][0] = global_elements[i * 2];
      global_rho[i][1] = global_elements[i * 2 + 1];
    }
    return global_rho;
}

void RBMDParallelUtil::SetRbeConfigAlphaNumber(double rbe_alpha_number) {
    DomainServiceLocator::Instance().GetDomain()->SetGlobalRbeAlphaPnumber(
        rbe_alpha_number);
}

double RBMDParallelUtil::GetRbeConfigAlphaNumber() {
    return DomainServiceLocator::Instance()
            .GetDomain()
            ->GetGlobalRbeAlphaPnumber();
}

std::vector<std::vector<float> > RBMDParallelUtil::GetRbeRadomM() {
     std::vector<std::vector<float> > result;
    std::vector<float> temp_result;
    temp_result.reserve(GetRbeConfigPnumber());
    // 以随机数引擎作为种子的伪随机数生成器
    std::mt19937 gen(GetRbeGlobalRandomSeed());
    constexpr double pi = 3.141592653589793;
    // 定义正态分布
    for (int dim = 0; dim < 3; ++dim) {
        auto global_length =
                DomainServiceLocator::Instance().GetDomain()->GetGlobalLength(dim);
        std::normal_distribution<float> d(
            0, sqrt(GetRbeConfigAlphaNumber() * global_length * global_length /
                    (2 * pi * pi)));

        for (int i = 0; i < GetRbeConfigPnumber(); ++i) {
            temp_result.emplace_back(d(gen));
        }
        result.emplace_back(temp_result);
        temp_result.clear();
    }
    return result;
}

void  RBMDParallelUtil::GetRbeRadomM_vec3f(std::vector<vtkm::Vec3f>& result)
{
    // 以随机数引擎作为种子的伪随机数生成器
    std::mt19937 gen(GetRbeGlobalRandomSeed());
    constexpr double pi = 3.141592653589793;
    // 定义正态分布
    vtkm::Vec3f temp_result;

    auto global_length_x = DomainServiceLocator::Instance().GetDomain()->GetGlobalLength(0);
    std::normal_distribution<float> d_x(
      0, sqrt(GetRbeConfigAlphaNumber() * global_length_x * global_length_x / (2 * pi * pi)));

    auto global_length_y = DomainServiceLocator::Instance().GetDomain()->GetGlobalLength(1);
    std::normal_distribution<float> d_y(
      0, sqrt(GetRbeConfigAlphaNumber() * global_length_y * global_length_y / (2 * pi * pi)));

    
    auto global_length_z = DomainServiceLocator::Instance().GetDomain()->GetGlobalLength(2);
    std::normal_distribution<float> d_z(
      0, sqrt(GetRbeConfigAlphaNumber() * global_length_z * global_length_z / (2 * pi * pi)));

   for (int i = 0; i < GetRbeConfigPnumber(); ++i)
    {
        do
        {
            temp_result = { round(d_x(gen)), round(d_y(gen)), round(d_z(gen)) };
        } while (std::abs(temp_result[0]) == 0 && std::abs(temp_result[1]) == 0 &&
                 std::abs(temp_result[2]) == 0);

        result.push_back(temp_result);
    }
                
    //int rank;
    //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //
    //if (rank == 0)
    //{
    //  std::ofstream _psample_MPI_0;
    //  _psample_MPI_0.open("Parallel_psample_MPI_0.txt");
    //
    //  for (size_t i = 0; i < result.size(); i++)
    //  {
    //    _psample_MPI_0 << i << "," << result[i] << std::endl;
    //  }
    //}
    //else if (rank == 1)
    //{
    //  std::ofstream _psample_MPI_1;
    //  _psample_MPI_1.open("Parallel_psample_MPI_1.txt");
    //
    //  for (size_t i = 0; i < result.size(); i++)
    //  {
    //    _psample_MPI_1 << i << "," << result[i] << std::endl;
    //  }
    //}

    return ;
}

double calculateDistance(Atom &a, Atom &b) {
    return (std::pow(a.r(0) - b.r(0), 2) + std::pow(a.r(1) - b.r(1), 2) +
                     std::pow(a.r(2) - b.r(2), 2));
}

void RBMDParallelUtil::ProcessCurrentCell(
    Cell &current_cell,
    std::unordered_map<unsigned long, NeighborData> &verlet_list
) {
    if (!current_cell.IsInnerCell() && !current_cell.IsBoundaryCell()) {
        return;
    }
    auto current_cell_atoms = current_cell.GetAtoms();
    const auto cutoff = this->linked_cell_->GetCutoff() * this->linked_cell_->GetCutoff() ;
    for (int i = 0; i < current_cell_atoms.size(); ++i) {
        for (int j = i + 1; j < current_cell_atoms.size(); ++j) {
            if (calculateDistance(current_cell_atoms[i], current_cell_atoms[j]) <( cutoff-0.0001)) {
                unsigned long atom_i_id = current_cell_atoms[i].GetID();
                unsigned long atom_j_id = current_cell_atoms[j].GetID();
                verlet_list[atom_i_id].neighborIds.push_back(atom_j_id);
                verlet_list[atom_j_id].neighborIds.push_back(atom_i_id); // 对称添加
            }
        }
    }
}

void RBMDParallelUtil::ProcessCellPair(
    Cell &current_cell, Cell &neighbour_cell,
    std::unordered_map<unsigned long, NeighborData> &verlet_list
) {
    // TODO 优先度：高 后续改为方便能量输出的
    auto current_cell_atoms = current_cell.GetAtoms();
    auto neighbour_cell_atoms = neighbour_cell.GetAtoms();
    const auto cutoff = this->linked_cell_->GetCutoff() * this->linked_cell_->GetCutoff() ;
    for (int i = 0; i < current_cell_atoms.size(); ++i) {
        for (int j = 0; j < neighbour_cell_atoms.size(); ++j) {
            if (current_cell_atoms[i].GetID() == neighbour_cell_atoms[j].GetID())
                continue;
            if (calculateDistance(current_cell_atoms[i], neighbour_cell_atoms[j]) <
                ( cutoff-0.0001)) {
                unsigned long atom_i_id = current_cell_atoms[i].GetID();
                unsigned long atom_j_id = neighbour_cell_atoms[j].GetID();
                if (neighbour_cell.IsHaloCell()) {
                    verlet_list[atom_i_id].ghostNeighborIds.push_back(
                        neighbour_cell_atoms[j].GetID());
                    verlet_list[atom_i_id].ghostNeighborPositions.emplace_back(
                        neighbour_cell_atoms[j].Position3D());
                } else {
                    if (std::find(verlet_list[atom_i_id].neighborIds.begin(), verlet_list[atom_i_id].neighborIds.end(),
                                  atom_j_id) == verlet_list[atom_i_id].neighborIds.end()) {
                        verlet_list[atom_i_id].neighborIds.push_back(atom_j_id);
                    }
                    if (std::find(verlet_list[atom_j_id].neighborIds.begin(), verlet_list[atom_j_id].neighborIds.end(),
                                  atom_i_id) == verlet_list[atom_j_id].neighborIds.end()) {
                        verlet_list[atom_j_id].neighborIds.push_back(atom_i_id);
                    }
                }
            }
        }
    }
}

bool check_in_range(const unsigned long index, const unsigned start_index,
                    const unsigned long end_index) {
    if (index >= start_index && index < end_index) {
        return true;
    }
    return false;
}

//
//
// //TODO:how to fix bug?
std::unordered_map<unsigned long, NeighborData> &
RBMDParallelUtil::BuildVerletListFull() {
    result.clear();
    auto &cell_data = GetCellData();
    unsigned long start = 0;
    unsigned long end = cell_data.size();
    for (unsigned long cell_index = start; cell_index < end;
         ++cell_index) {
        // 当前遍历的Cell
        auto &current_cell = cell_data.at(cell_index);
        if (!current_cell.IsHaloCell()) {
            ProcessCurrentCell(current_cell, result);
            for (auto &forward: forward_neighbour_offsets_) {
                if (check_in_range(cell_index + forward, start, end)) {
                    auto &neighbour_cell = cell_data.at(cell_index + forward);
                    ProcessCellPair(current_cell, neighbour_cell, result);
                }
            }
            for (auto backward: backward_neighbour_offsets_) {
                if (check_in_range(cell_index - backward, start, end)) {
                    auto &neighbour_cell = cell_data.at(cell_index - backward);
                    ProcessCellPair(current_cell, neighbour_cell, result);
                }
            }
        }
    }
    return result;
}

std::unordered_map<unsigned long, NeighborData> &
RBMDParallelUtil::BuildVerletListC08(bool us_eighth_shell) {
    result.clear();

    const std::array<unsigned long, 3> strides = {2, 2, 2};
    std::array<unsigned long, 3> dims = {
        this->linked_cell_->GetCellPerDim(0),
        this->linked_cell_->GetCellPerDim(1),
        this->linked_cell_->GetCellPerDim(2)
    };

    std::array<unsigned long, 3> end{};
    for (int d = 0; d < 3; ++d) {
        end[d] = dims[d] - 1;
    }
    auto cell_data = GetCellData();

    for (unsigned long col = 0; col < 8; ++col) {
        std::array<unsigned long, 3> starts = OneToThreeD(col, strides);
        if (us_eighth_shell) {
            // if we are using eighth shell, we start at 1,1,1 instead of 0,0,0
            for (unsigned short i = 0; i < 3; ++i) {
                starts[i] += 1;
            }
        }

        const unsigned long start_x = starts[0], start_y = starts[1],
                start_z = starts[2];
        const unsigned long end_x = end[0], end_y = end[1], end_z = end[2];
        const unsigned long stride_x = strides[0], stride_y = strides[1],
                stride_z = strides[2];

        if (cell_data.empty()) {
            throw std::runtime_error("error in get cell data");
        }
        for (unsigned long z = start_z; z < end_z; z += stride_z) {
            for (unsigned long y = start_y; y < end_y; y += stride_y) {
                for (unsigned long x = start_x; x < end_x; x += stride_x) {
                    // 开始处理
                    unsigned long baseIndex = ThreeToOneD(x, y, z, dims);
                    using std::pair;
                    std::array<unsigned long, 3> threeDIndex =
                            OneToThreeD(baseIndex, dims);
                    for (int d = 0; d < 3; ++d) {
                        rbmd_assert(threeDIndex[d] != dims[d] - 1);
                    }
                    const int num_pairs = _cellPairOffsets8Pack.size();
                    for (int j = 0; j < num_pairs; ++j) {
                        std::pair<unsigned long, unsigned long> current_pair =
                                _cellPairOffsets8Pack[j];

                        unsigned offset1 = current_pair.first;
                        unsigned cellIndex1 = baseIndex + offset1;

                        unsigned offset2 = current_pair.second;
                        unsigned cellIndex2 = baseIndex + offset2;
                        if (cellIndex1 >= cell_data.size() ||
                            cellIndex2 >= cell_data.size() || cellIndex1 < 0 ||
                            cellIndex2 < 0) {
                            throw std::runtime_error("cell index error in c08");
                        }
                        auto cell1 = cell_data.at(cellIndex1);
                        auto cell2 = cell_data.at(cellIndex2);

                        if ((not us_eighth_shell) and cell1.IsHaloCell() and
                            cell2.IsHaloCell()) {
                            continue;
                        }

                        if (cellIndex1 == cellIndex2) {
                            ProcessCurrentCell(cell1, result);
                            // std::cout<< "ProcessCurrentCell" <<std::endl;
                        } else {
                            if (not us_eighth_shell) {
                                if (!cell1.IsHaloCell()) {
                                    // std::cout<< "ProcessCellPair" <<std::endl;
                                    ProcessCellPair(cell1, cell2, result);
                                } else {
                                    // std::cout<< "ProcessCellPair" <<std::endl;
                                    ProcessCellPair(cell2, cell1, result);
                                }
                                // } else {
                                //     // if we use eighthShell, we have to sum everything.
                                //     also we don't care about order of the cells!
                                //     ProcessCellPair(cell1, cell2, true);
                                // }
                            }
                        }
                    }
                }
            }
        }
    }
    return result;
}

vtkm::Vec2f RBMDParallelUtil::GetGlobalErengy(float locl_erengy)
{
    vtkm::Vec2f global_energy;
    float global_erengy_temp = 0;
    unsigned int local_atoms_number = this->GetNumAtoms(false);
    unsigned int global_atoms_number = 0;

    MPI_Allreduce(&locl_erengy, &global_erengy_temp, 1, MPI_FLOAT, MPI_SUM, domdec_->GetMPIComm());
    MPI_Allreduce(&local_atoms_number, &global_atoms_number, 1, MPI_UNSIGNED, MPI_SUM, domdec_->GetMPIComm());

    if (domdec_->GetCurrentRank() == 0) {
        // std::cout << "总的粒子数：" << global_atoms_number << std::endl;
        global_energy[0] = (global_erengy_temp);
        global_energy[1] = (global_erengy_temp / (static_cast<float>(global_atoms_number)));
    }
    return global_energy;
}

int RBMDParallelUtil::GetCurrentRank() const {
    int rank;
    MPI_Comm_rank(domdec_->GetMPIComm(), &rank);
    return rank;
}


vtkm::Vec3f  RBMDParallelUtil::GetBoxLength(){
  auto domain = DomainServiceLocator::Instance().GetDomain();
  double bBoxMin[3];
  double bBoxMax[3];
  domdec_->GetBoundingBoxMinMax(bBoxMin, bBoxMax);
  vtkm::Vec3f res = {static_cast<float>(bBoxMax[0]-bBoxMin[0]),static_cast<float>(bBoxMax[1]-bBoxMin[1]),static_cast<float>(bBoxMax[2]-bBoxMin[2])};
  return res;
}


std::vector<float> RBMDParallelUtil::GetGlobalRho(const std::vector<float>& local_rho)
{ // TODO:后续改为返回引用吧
  int size = GetRbeConfigPnumber();
  std::vector<float> global_rho(size);

  MPI_Allreduce(local_rho.data(),
                global_rho.data(),
                global_rho.size() /* size*/,
                MPI_FLOAT,
                MPI_SUM,
                domdec_->GetMPIComm());
  return global_rho;
}
