#pragma once
#include <array>

/*
 * 三维的数据按行主序存储
 */

template <typename T>
T ThreeToOneD(T x, T y, T z, const std::array<T, 3> dims) {
  return (z * dims[1] + y) * dims[0] + x;
}

template <typename T>
T ThreeToOneD(const std::array<T, 3>& pos, const std::array<T, 3> dims) {
  return (pos[2] * dims[1] + pos[1]) * dims[0] + pos[0];
}

template <typename T>
std::array<T, 3> OneToThreeD(T ind, const std::array<T, 3> dims) {
  std::array<T, 3> pos;
  pos[2] = ind / (dims[0] * dims[1]);
  pos[1] = (ind - pos[2] * dims[0] * dims[1]) / dims[0];
  pos[0] = ind - dims[0] * (pos[1] + dims[1] * pos[2]);
  return pos;
}