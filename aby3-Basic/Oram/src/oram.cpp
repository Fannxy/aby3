#include "../include/oram.h"

void get_permutation(size_t len, std::vector<size_t> &permutation) {
  // generate the permutation for a len-sized array and save it into the
  // permutation.
  permutation.resize(len);
  for (size_t i = 0; i < len; i++) {
    permutation[i] = i;
  }
  std::random_shuffle(permutation.begin(), permutation.end());
  return;
}