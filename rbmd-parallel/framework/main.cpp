#include "MDApplication.h"
#include <iostream>
#include <string>

#include <mpi/mpi.h>

const std::string RBMD_VERSION = "1.2.0";
int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  std::shared_ptr<Application> app = std::make_shared<MDApplication>(argc, argv);
  app->Run();
  MPI_Finalize();
  return 0;
};