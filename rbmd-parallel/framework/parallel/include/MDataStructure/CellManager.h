#pragma once

#include <array>
#include <vector>

//! @brief
//! 一个线程安全的单例实现(C++11之后)，负责管理Cell的一些参数，后续需要考虑如何在Device端实现，暂时使用该类，后续抽象基类并继承Host和Device
class CellManager {
private:
    CellManager() {

    }

    ~CellManager()= default;

    CellManager(const CellManager &);

    CellManager &operator=(const CellManager &);

#pragma region 私有成员变量

    //! @brief 是否完成了初始化
    bool _initlized = false;

    //! @brief 每个维度的cell数量
    std::array<unsigned long, 3> _num_cells_per_dimension{};

    //! @brief
    std::array<double, 3> _bounding_box_min{}, _bounding_box_max{};

    //! @brief
    //! halo层的盒子最小（左下前方），最大坐标（右上后方）。Halo层是为了并行需要而添加的，层内包含ghost原子。
    std::array<double, 3> _halo_bounding_box_min{}, _halo_bounding_box_max{};

    //! @brief cell的长度，通常是cutoff的数值
    std::array<double, 3> _cell_length{};

    std::array<unsigned int, 3> _halo_width_in_num_cells{};

    //! @brief 与Cell数组长度相同，如果是HaloCell则为true
    std::vector<bool> _halo_cell_flags;

    enum class CellType {
        HALO = 0,
        BOUNDARY = 1,
        INNER = 2,
        INNERMOST = 3
    };  // TODO 理清楚并详细注释
#pragma endregion

#pragma region 私有函数

    /*!
     * @brief 将一维索引转换为三维索引          //TODO memory per?
     * 内存中数据应该是行主序存储的！
     * @param cell_index 一维cell索引
     * @return xyz的索引数组
     */
    std::array<unsigned long, 3> MapLinearTo3DIndex(
            unsigned long cell_index) const;

    /*!
     * @brief 判断Cell属于何种类型
     * @param cell_type 需要判断的Cell类型
     * @param cell_index Cell的索引
     * @return 如果是该类型则返回true,否则返回false
     */
    bool IsCellOfType(CellType cell_type, unsigned long cell_index) const;

    double GetCellBorder(unsigned long cell_index, int dimension) const;

#pragma endregion

public:
    /*!
     * @brief 单例只能通过该方法获取实例
     * @return 全局唯一的CellManager对象，值得注意的是C++之前本代码是非线程安全的
     */
    static CellManager &Instance();

    /*!
     * @brief 判断给定index的Cell是否是HaloCell
     * @param cell_index 需要判定的Cell的索引
     * @return 如果该Cell是HaloCell则返回true否则返回false
     */
    bool IsHaloCell(unsigned long cell_index) const;

    /*!
     * @brief 判断给定index的Cell是否是HaloCell
     * @param cell_index 需要判定的Cell的索引
     * @return 如果该Cell是BoundaryCell则返回true否则返回false
     */
    bool IsBoundaryCell(unsigned long cell_index) const;

    /*!
     * @brief 判断给定index的Cell是否是InnerCell
     * @param cell_index 需要判定的Cell的索引
     * @return 如果该Cell是InnerCell则返回true否则返回false
     */
    bool IsInnerCell(unsigned long cell_index) const;

    /*!
     * @brief 判断给定index的Cell是否是InnerMostCell
     * @param cell_index 需要判定的Cell的索引
     * @return 如果该Cell是InnerMostCell则返回true否则返回false
     */
    bool IsInnerMostCell(unsigned long cellIndex) const;

    double GetCellBoundingBoxMin(unsigned long cell_index, int dimension) const;

    double GetCellBoundingBoxMax(unsigned long cell_index, int dim) const;

    void init(int cellsPerDim[3],
            double haloBoxMin[3], double haloBoxMax[3],
            double boxMin[3], double boxMax[3],
            double cellLength[3], int haloWidthInNumCells[3]);

};
