#pragma once
#include <array>
#include "DomainDecompositionMPI.h"
#include "Cell.h"
#include <random>
#include <vtkm/Types.h>

struct NeighborData {
    std::vector<unsigned long> neighborIds; // 邻居id数组
    std::vector<unsigned long> ghostNeighborIds; // 鬼邻居id数组
    std::vector<vtkm::Vec3f> ghostNeighborPositions; // 鬼邻居位置数组
};

static int step;

class RBMDParallelUtil {
private:
#pragma region 变量
    DomainDecompositionMPI *domdec_ = nullptr;
    LinkedCell *linked_cell_ = nullptr;
    std::random_device rd;
    std::unordered_map<unsigned long, NeighborData> result = std::unordered_map<unsigned long, NeighborData>{};
#pragma endregion
 
    void ProcessCurrentCell(Cell &current_cell,
                            std::unordered_map<unsigned long, NeighborData> &verlet_list);

    void ProcessCellPair(Cell &current_cell,
                         Cell &neighbour_cell,
                         std::unordered_map<unsigned long, NeighborData> &verlet_list);

    std::array<std::pair<unsigned long, unsigned long>, 14> _cellPairOffsets8Pack{};
    std::array<unsigned long, 8> _cellOffsets8Pack{};
    std::array<unsigned long, 13> forward_neighbour_offsets_{};
    std::array<unsigned long, 13> backward_neighbour_offsets_{};

    void computeOffsets();

    void GetOffset();

public:

    static void printAtomInfo(unsigned long target_id, LinkedCell *linked_cell,
                              std::string msg) {
        const auto &cells = linked_cell->GetCells();
        for (auto cell: cells) {
            const auto &atoms = cell.GetAtoms();
            for (const auto &atom: atoms) {
                if (atom._id == target_id) {
                    std::cout << "当前第  " << step - 1 << "步" << std::endl;
                    std::cout << msg << "-------------------";
                    std::cout << "  当前原子ID： " << atom._id << "  原子位置  " << "["
                            << atom._r[0] << "," << atom._r[1] << "," << atom._r[2] << "]"
                              << "当前Cell ID:" << cell.GetID() << "是否鬼原子："
                              << (cell.IsHaloCell() ? "是" : "否") << std::endl;
                }
            }
        }
    }

    /*!
       * 构造函数 完成了Domian的划分，数据结构的初始化。构造之前需要确保已经进行了MPI_Init.非并行程序请勿调用。
       * @param x_y_z_length 模拟盒子的大小xyz
       * @param cutoff cutoff
       */
    RBMDParallelUtil(std::array<double, 3> min,std::array<double, 3> max, double cutoff);

    ~RBMDParallelUtil();

    /*!
       * 数据读取时使用，写入原子数据
       * @param id 原子ID
       * @param x 坐标X
       * @param y 坐标Y
       * @param z 坐标Z
       */
    void AddAtom(int id, double x, double y, double z);

    void AddAtomAndVelocity(int id,
                            double x,
                            double y,
                            double z,
                            double vx,
                            double vy,
                            double vz,
                            double charge);

    /*!
       * 返回Cell
       * Cell内包含属于该Cell的粒子
       * 访问邻居Cell:
       * for (auto& neighbourOffset : this->forward_neighbour_offsets_) {
       * 	    auto&  neighbourCell = this->_cells->at(cellIndex + neighbourOffset)
       * 	};
       * 	通常每个Cell有13个前向邻居，13个后向邻居。
       * 	对于InnerCell只需计算与前向邻居的相互作用，对于boundy需要计算前向邻居与后向邻居的相互作用。
       * 	计算atom i，j相互作用时:
       * 	i.f += f
       * 	if j is not in halo cell:
       * 	    j -= f
       * 	i不会来自halo cell，i取自当前cell
       * @return
       */
    std::vector<Cell> &GetCellData();


    /*!
       * 计算力并更新位置之前调用该函数实现并行需要的预更新。  NOTE：访问Cell之前必须调用这个，否则无法访问HaloCell
       */
    void PrepareUpdateForPositionChange();

    void PrepareUpdateForPositionChangeTest(int current_step);


    /*!
       * 更新完成位置后调用该函数实现并行需要的后更新   NOTE：调用这个会清除HaloCell
       */
    void PostUpdateForPositionChange();


    unsigned int GetNumAtoms(bool include_halo = false); //TODO 修改 获取全部原子

    /*！
       * 用于RBE的m生成的随机数种子，可以保证每个进程使用同一个种子生成同样的随机数
       */
    unsigned int GetRbeGlobalRandomSeed();


    /*！
       * 从配置文件中读取并设置RBE的p的值
       */
    void SetRbeConfigPnumber(int rbe_p_number);

    /*!
       * 获取配置的RBE的p的值
       * @return 用于RBE计算的p值
       */
    int GetRbeConfigPnumber();


    /*！
       * 从配置文件中读取并设置RBE的alpha的值
       */
    void SetRbeConfigAlphaNumber(double rbe_alpha_number);

    /*!
       * 获取配置的RBE的alpha的值
       * @return 用于RBE计算的alpha值
       */
    double GetRbeConfigAlphaNumber();

    //TODO 暂时先这样，后面考虑domian直接同步全部？ 但是这个是在计算里面同步的，好像不是很合适。

    /*!
     *
     * @param local_rho 局部（当前进程运算结果）rho数组
     * @return RBE中全局的rho数组
     */
    std::vector<vtkm::Vec2f> GetGlobalRho(const std::vector<vtkm::Vec2f>& local_rho);


    /*!
     * 生成随机数组M，不同进程调用得到相同的结果，每次调用得到不同的随机数组M
     * @return
     */
    std::vector<std::vector<float>> GetRbeRadomM();
    void GetRbeRadomM_vec3f(std::vector<vtkm::Vec3f>& result);


    // 二维数组：数组1 原子id=1的邻居   一维数组：邻居的长度   //TODO 划分之后有可能ID=1的就不在当前domain了

    std::unordered_map<unsigned long, NeighborData> &BuildVerletListFull();

    std::unordered_map<unsigned long, NeighborData> &BuildVerletListC08(bool us_eighth_shell = false);


    /*!
     * 统计能量输出
     * @param locl_erengy 局部能量，就是当前进程运算的能量。
     * @return Erengy能量数组，如果当前进程不是主进程返回空的vector，如果是主进程返回size为2的vector，第一个是总能量，第二个是平均的能量。auto erengys = GetGlobalErengy(local_erengr) ---> if(!energys.empty()) std::cout
     */
    //std::vector<vtkm::Vec2f> GetGlobalErengy(float locl_erengy);
    vtkm::Vec2f GetGlobalErengy(float locl_erengy);

    int GetCurrentRank() const;

	vtkm::Vec3f  GetBoxLength();


    std::vector<float> GetGlobalRho(const std::vector<float>& local_rho);
};
