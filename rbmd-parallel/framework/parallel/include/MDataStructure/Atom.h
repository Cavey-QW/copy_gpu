#pragma once

#include <array>
#include <vtkm/Types.h>
//#ifdef USE_DOUBLE
//typedef double Real_mpi;
//#else
//typedef float Real_mpi;
//#endif
//TODO 暂时未兼容
typedef double Real_mpi;



class Atom {
public:
#pragma region 构造函数
  explicit Atom(unsigned long id = 0,
          double rx = 0., double ry = 0., double rz = 0.,
          double vx = 0., double vy = 0., double vz = 0.,double charge = 0.
 );
#pragma endregion

#pragma region 变量
    //! @brief 位置坐标
    Real_mpi _r[3];
    //! @brief 力
    Real_mpi _F[3];
    //! @brief 速度
    Real_mpi _v[3];
    //! @brief 质量
    Real_mpi _m;
    Real_mpi _L[3];  /**< angular momentum */
    Real_mpi _Vi[3]; /** Virial tensor **/
    //TODO 其他变量
    //! @brief 唯一ID
    unsigned long _id;
    Real_mpi _M[3];  /**< torsional moment */
    Real_mpi _I[3]{0.,0.,0.},_invI[3]{0.,0.,0.};  // moment of inertia for principal axes and it's inverse

    Real_mpi _charge;

    double M(unsigned short d) const  { return _M[d]; }

    /** get the virial **/
    double Vi(unsigned short d) const  { return _Vi[d];}
#pragma endregion

public:
#pragma region 变量操作

    /*!
     * @brief 获取当前原子ID
     * @return 原子ID
     */
    unsigned long GetID() const;

    /*!
     * @brief 设置原子ID
     * @param id ID值
     */
    void SetID(unsigned long id);

    /*!
     * @brief 获取原子某个维度的坐标值
     * @param dimension 维度，通常是xyz三个维度，x为0,y为1,z为2
     * @return 该维度的原子坐标，使用浮点数表示
     */
    Real_mpi r(unsigned short dimension) const;

    vtkm::Vec3f Position3D() const;

    vtkm::Vec3f Velocity3D() const;
    /*!
     * @brief 设置该原子某个维度的坐标
     * @param dimension 维度
     * @param r 需要设置的坐标值
     */
    void Setr(unsigned short dimension, Real_mpi r);

    /*!
     * @brief 获取原子某个维度的速度
     * @param dimension 维度
     * @return 浮点数表示的速度
     */
    Real_mpi v(unsigned short dimension) const;

    /*!
     * @brief 设置原子某个维度的速度
     * @param dimension 维度
     * @param v 需要设置的数值
     */
    void Setv(unsigned short dimension, Real_mpi v);

    /*!
     * 判断原子是否位于某个box区域内
     * @param bounding_box_min box的左下前角坐标
     * @param bounding_box_max box的右上后角坐标
     * @return true表示位于该box区域内，否则返回false
     */
    bool InBox(const double bounding_box_min[3], const double bounding_box_max[3]) const;

    void setF(unsigned short d, double F) ;
    /** get coordinate of current force onto molecule */
    double F(unsigned short d) const ;

    std::array<double, 3> F_arr() ;


    std::array<double, 3> M_arr() {
        std::array<double, 3> m_ret{M(0), M(1), M(2)};
        return m_ret;
    }
    virtual std::array<double, 3> Vi_arr() {
        std::array<double, 3> vi_ret{Vi(0), Vi(1), Vi(2)};
        return vi_ret;
    }

    void Fadd(const double a[])  { for(unsigned short d=0;d<3;++d) _F[d]+=a[d]; }
    void Madd(const double a[])  { for(unsigned short d=0;d<3;++d) _M[d]+=a[d]; }
    void Viadd(const double a[])  { for(unsigned short d=0;d<3;++d) _Vi[d]+=a[d]; }
    void move(int d, double dr)  { _r[d] += dr; }
#pragma endregion



};

