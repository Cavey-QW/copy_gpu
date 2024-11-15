#include "../include/mmap_reader.h"

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

MmapReader::MmapReader(const std::string& filePath)
    : BaseReader(filePath),
      _file_size(0),
      _mapped_memory(nullptr),
      _line_start(nullptr),
      _locate(0) {}

MmapReader::~MmapReader() {
  if (-1 == munmap(_mapped_memory, _file_size)) {
    // log
  }
}

int MmapReader::Execute() {
  auto fd = ::open(_file_path.c_str(), O_RDONLY);
  if (-1 == fd) {
    // log
    return -1;
  }

  // get file size
  struct stat file_stat;
  if (-1 == fstat(fd, &file_stat)) {
    // log
    close(fd);
    return -1;
  }
  _file_size = file_stat.st_size;

  // mmap
  _mapped_memory = static_cast<char*>(
      mmap(nullptr, _file_size, PROT_READ, MAP_PRIVATE, fd, 0));
  if (MAP_FAILED == _mapped_memory) {
    // log
    close(fd);
    return -1;
  }

  _line_start = _mapped_memory;

  // close fd
  close(fd);

  return 0;
}
