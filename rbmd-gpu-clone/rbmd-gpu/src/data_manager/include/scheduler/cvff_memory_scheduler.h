#pragma once
#include "memory_scheduler.h"

class CVFFMemoryScheduler : public MemoryScheduler {
 public:
  CVFFMemoryScheduler(){};
  virtual ~CVFFMemoryScheduler() = default;
  /**
   * @brief async memeory host to device
   * @return error code
   */
  bool asyncMemoryH2D() override;

  /**
   * @brief async memory device to host
   * @return error code
   */
  bool asyncMemoryD2H() override;
};
