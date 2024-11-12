#pragma once

#include <vector>
#include "Atom.h"
#include "SingleCellIterator.hpp"

class Cell {
public:
    Cell();

    ~Cell();

    //! @brief Cell的id,与在数组中的索引相同
    unsigned long _id{};

    /*!
     * @brief 获取Cell的ID
     * @return ID
     */
    unsigned long GetID() const;

    /*!
     * @brief 设置cell的ID
     * @param id ID的值
     */
    void SetID(unsigned long id);

    /*!
     * @brief 清空当前cell
     */
    void EmptyAtoms();

    /*!
     * @brief 添加单个原子到Cell中
     * @param atom 原子对象
     * @return 是否插入成功 -----TODO考虑要不要加上判断是否重复插入，遍历一次很耗时
     */
    bool AddAtomToCell(Atom &atom); //TODO

    /*!
     * 获取当前cell原子总数
     * @return 当前cell的原子数量
     */
    unsigned long GetAtomsCount() const;

    /*!
     * 判断当前Cell是否为空
     * @return true表示当前Cell为空，false表示当前Cell非空
     */
    bool IsEmpty() const;

    /*!
     * 获取当前Cell的左下前角坐标
     * @param dimension 维度
     * @return 当前Cell的左下前角坐标
     */
    double GetBoundingBoxMin(int dimension) const;

    double GetBoundingBoxMax(int dimension) const;


    bool IsHaloCell() const;

    bool IsBoundaryCell() const;

    bool IsInnerCell() const;

    bool IsInnerMostCell() const;

    void AssignCellToHaloRegion();

    bool TestPointInCell(const double point[3]) const;

    /*!
     * 判断原子是否位于当前Cell的box区域
     * @param atom 需要判断的原子对象
     * @return true表示在当前Cell的区域，false表示不在当前Cell的区域
     */
    bool InBox(const Atom &atom) const;

    bool IsNotEmpty() const;

    SingleCellIterator<Cell> iterator();

    void moleculesAtNew(size_t i, Atom *&multipurposePointer) {
        multipurposePointer = &_atoms.at(i);
    }

    bool deleteMoleculeByIndex(size_t index);

    void increaseMoleculeStorage(size_t numExtraMols);


    void PreUpdateLeavingMolecules();
    void UpdateLeavingMolecules(Cell &otherCell);

    void PostUpdateLeavingMolecules();

    //! @brief 模拟过程中离开当前Cell的粒子
    std::vector<Atom> _leaveing_atoms;

    std::vector<Atom>& GetAtoms();

#pragma region 暂时不用的 后续动态负载均衡要
    /*
     *  // setBoxMin setBoxMax 暂时用不到
     * void assignCellToBoundaryRegion() { mardyn_assert(isBoundaryCell()); }
	void assignCellToInnerRegion() { mardyn_assert(isInnerCell()); }
	void assignCellToInnerMostAndInnerRegion() { mardyn_assert(isInnerMostCell() and isInnerCell()); }

	void skipCellFromHaloRegion() { mardyn_assert(not isHaloCell()); }
	void skipCellFromBoundaryRegion() { mardyn_assert(not isBoundaryCell()); }
	void skipCellFromInnerRegion() { mardyn_assert(not isInnerCell()); }
	void skipCellFromInnerMostRegion() { mardyn_assert(not isInnerMostCell()); }
     moleculesAtNew rmm使用 降低内存消耗
     * */
#pragma endregion

private:
    //! @brief 每个cell中有若干粒子
    std::vector<Atom> _atoms{};
};
