#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>

int main(int argc, char *argv[]) {
  int n = 10;
  double freq = 1.0;

  for (int i = 1; i < argc; ++i) {
    if (std::strcmp(argv[i], "--n") == 0 && i + 1 < argc) {
      n = std::atoi(argv[++i]);
    } else if (std::strcmp(argv[i], "--freq") == 0 && i + 1 < argc) {
      freq = std::atof(argv[++i]);
    }
  }

  for (int i = 0; i < n; ++i) {
    double t = i * 0.1;
    std::cout << std::sin(freq * t) << "\n";
  }

  return 0;
}
