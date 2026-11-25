/** @file native_files.hpp Binary snapshots, SHA-256 and exclusive output
 * creation. */
#pragma once
#include "third_party/picosha2.h"
#include <algorithm>
#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>
#include <vector>
#ifdef _WIN32
#include <io.h>
#include <windows.h>
#else
#include <unistd.h>
#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif
#endif
namespace swordfish {
namespace fs = std::filesystem;
inline std::string read_file(const fs::path &path) {
  std::ifstream in(path, std::ios::binary);
  if (!in)
    throw std::runtime_error("Cannot read file: " + path.string());
  std::string result;
  char buffer[16384];
  while (in) {
    in.read(buffer, sizeof buffer);
    result.append(buffer, size_t(in.gcount()));
    if (result.size() > 128 * 1024 * 1024)
      throw std::invalid_argument("Input exceeds 128 MiB: " + path.string());
  }
  if (!in.eof())
    throw std::runtime_error("File read failure: " + path.string());
  return result;
}
inline std::string sha256(const std::string &bytes) {
  return picosha2::hash256_hex_string(bytes);
}
inline void require_new_path(const fs::path &path) {
  if (fs::exists(path) || fs::is_symlink(fs::symlink_status(path)))
    throw std::invalid_argument("Refusing to replace " + path.string());
}
inline void exclusive_file(const fs::path &path, const std::string &bytes) {
#ifdef _WIN32
  int fd = _wopen(path.c_str(), _O_CREAT | _O_EXCL | _O_WRONLY | _O_BINARY,
                  _S_IREAD | _S_IWRITE);
#else
  int fd = ::open(path.c_str(), O_CREAT | O_EXCL | O_WRONLY, 0666);
#endif
  if (fd < 0)
    throw std::runtime_error("Cannot create " + path.string() + ": " +
                             std::strerror(errno));
  bool ok = true;
  size_t offset = 0;
  while (offset < bytes.size()) {
#ifdef _WIN32
    int written =
        _write(fd, bytes.data() + offset,
               unsigned(std::min<size_t>(bytes.size() - offset, 1 << 20)));
#else
    auto written = ::write(fd, bytes.data() + offset, bytes.size() - offset);
#endif
    if (written < 0 && errno == EINTR)
      continue;
    if (written <= 0) {
      ok = false;
      break;
    }
    offset += size_t(written);
  }
#ifdef _WIN32
  if (_commit(fd) != 0)
    ok = false;
  if (_close(fd) != 0)
    ok = false;
#else
  if (::fsync(fd) != 0)
    ok = false;
  if (::close(fd) != 0)
    ok = false;
#endif
  if (!ok) {
    std::error_code ec;
    fs::remove(path, ec);
    throw std::runtime_error("Cannot write " + path.string());
  }
}
inline fs::path executable_path() {
#ifdef _WIN32
  std::vector<wchar_t> buffer(32768);
  auto n = GetModuleFileNameW(nullptr, buffer.data(), DWORD(buffer.size()));
  if (n && n < buffer.size())
    return fs::path(std::wstring(buffer.data(), n));
#elif defined(__APPLE__)
  uint32_t size = 0;
  _NSGetExecutablePath(nullptr, &size);
  std::vector<char> buffer(size);
  if (_NSGetExecutablePath(buffer.data(), &size) == 0)
    return fs::canonical(buffer.data());
#else
  std::error_code ec;
  auto path = fs::read_symlink("/proc/self/exe", ec);
  if (!ec)
    return path;
#endif
  throw std::runtime_error("Cannot locate running executable for provenance");
}
} // namespace swordfish
