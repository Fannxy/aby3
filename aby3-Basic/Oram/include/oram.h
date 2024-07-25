#ifndef ORAM_H
#define ORAM_H

#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "test.h"

// #define DEBUG
#define TEST

template <typename T>
class StashElement {
 public:
  T data;
  size_t logicalIndex;
  StashElement() {}
  StashElement(T data, size_t logicalIndex)
      : data(data), logicalIndex(logicalIndex) {}
};

template <typename T>
class PackedIndex {
 public:
  size_t pack;
  std::vector<T> packedIndices;
  T logicalIndex;

  PackedIndex() {}

  PackedIndex(size_t pack, T logicalIndex)
      : pack(pack), logicalIndex(logicalIndex) {
    packedIndices.resize(pack);
  }

  void set_packed_indices(std::vector<T> &indices) {
    this->packedIndices = indices;
  }
};

template <typename T>
struct is_packed_index : std::false_type {};

template <typename T>
struct is_packed_index<PackedIndex<T>> : std::true_type {};

template <typename T>
inline void print_vector_impl(std::vector<T> &data, std::false_type) {
  size_t len = data.size();
  for (size_t i = 0; i < len; i++) {
    std::cout << data[i] << " ";
  }
  std::cout << std::endl;
}

template <typename T>
inline void print_vector_impl(std::vector<PackedIndex<T>> &data,
                              std::true_type) {
  size_t len = data.size();
  for (size_t i = 0; i < len; i++) {
    std::cout << "index: " << data[i].logicalIndex << " ";
    std::cout << "data: ";
    for (size_t j = 0; j < data[i].pack; j++) {
      std::cout << data[i].packedIndices[j] << " ";
    }
  }
  std::cout << std::endl;
}

template <typename T>
inline void print_vector(std::vector<T> &data) {
  print_vector_impl(data, is_packed_index<T>{});
}

void get_permutation(size_t len, std::vector<size_t> &permutation);

template <typename T, typename D>
void permutate(std::vector<T> &permutation, std::vector<D> &data) {
  std::vector<D> temp(data.size());
  for (size_t i = 0; i < permutation.size(); i++) {
    temp[permutation[i]] = data[i];
  }
  data = temp;
}

template <typename T>  // index map data type T.
class PosMap {
 public:
  size_t n;                    // data size;
  size_t pack;                 // reduce size;
  size_t S;                    // stash size;
  size_t t;                    // current stash size;
  std::vector<T> permutation;  // permutation, size_t
  std::vector<PackedIndex<T>> packed_index;

  bool
      linear;  // if linear, then we can use a simple array to store the posmap.
  size_t map_len;
  std::vector<bool> usage_map;
  PosMap<T>
      *subPosMap;  // if not linear, then we need to construct a subPosMap.
  std::vector<PackedIndex<T>> stash;  // stash, init to empty vector.

  PosMap() {}

  PosMap(size_t n, size_t pack, size_t S) : n(n), pack(pack), S(S), t(0) {
    map_len = n / pack;
    if (map_len < S)
      linear = true;
    else
      linear = false;
  }

  PosMap(size_t n, size_t pack, size_t S, std::vector<T> permutation)
      : n(n), pack(pack), S(S), t(0), permutation(permutation) {
    size_t map_len = n / pack;
    if (map_len < S) {  // now construct a linear map.
      linear = true;
      usage_map.resize(n, false);  // init to false.
    } else {                       // resursively construct a posmap.
      linear = false;
      // // construct the packed posMap.
      for (size_t i = 0; i < map_len; i++) {
        std::vector<T> tmp;
        for (size_t j = 0; j < pack; j++) {
          tmp.push_back(permutation[i * pack + j]);
        }
        // packed_permutation.push_back(tmp);
        PackedIndex<T> packedIndex(pack, i);
        packedIndex.set_packed_indices(tmp);
        packed_index.push_back(packedIndex);
      }
      // shuffle the packed posMap.
      std::vector<size_t> sub_permutation;
      get_permutation(map_len, sub_permutation);
      permutate(sub_permutation, packed_index);
      
#ifdef DEBUG
      print_vector(sub_permutation);
      print_vector(packed_index);
#endif

      // free the space for permutation
      std::vector<T>().swap(permutation);

      // construct the posMap.
      this->subPosMap = new PosMap<T>(map_len, pack, S, sub_permutation);
    }
  }

  size_t access(size_t index, bool fake) {  // input the logical index.

#ifdef DEBUG
    std::cerr << "index = " << index << std::endl;
#endif
    size_t physical_index = -1;

    if (linear) {  // for linear oram.

#ifdef DEBUG
      std::cerr << "linear" << std::endl;
#endif
      bool done = false;
      for (int i = 0; i < this->n; i++) {
        bool s1 = !fake && (i == index);
        // bool s2 = fake && (!this->usage_map[i]) && (!done);
        bool s2 = fake && (!done);
        if (s1 || s2) {
          this->usage_map[i] = true;
          done = true;
          physical_index = this->permutation[i];
        }
      }
    } else {  // for recursive oram.
#ifdef DEBUG
      std::cerr << "recursive" << std::endl;
#endif
      // first look into the stash.
      bool found = false;
      size_t h = index / pack;
      size_t l = index % pack;

#ifdef DEBUG
      std::cerr << "h = " << h << std::endl;
      std::cerr << "l = " << l << std::endl;
#endif

      for (int i = 0; i < this->t; i++) {
        if (this->stash[i].logicalIndex == h) {
          found = true;
          physical_index = this->stash[i].packedIndices[l];
        }
      }

#ifdef DEBUG
      std::cerr << "found in stash -  " << found << std::endl;
      std::cerr << "l = " << l << std::endl;
      std::cerr << "physical_index in stash = " << physical_index << std::endl;
#endif
      // select the (maybe) random physical index from the next-level index map.
      size_t tmp_physical_index = this->subPosMap->access(h, fake || found);

      // update the stash.
      this->stash.push_back(this->packed_index[tmp_physical_index]);
      this->t++;
      if (fake || !found) {  // the hit happened beforehead, aggregassively
                             // return the first element in the stash.
        PackedIndex<T> tmp = this->stash[this->t - 1];
        physical_index = tmp.packedIndices[l];
      }
    }
#ifdef DEBUG
    std::cerr << "physical_index = " << physical_index << std::endl;
#endif
    return physical_index;
  }
};

template <typename T>
class SqrtOram {
 public:
  int n;                               // data size;
  int S;                               // stash size;
  int pack;                            // reduce size;
  int t;                               // current stash size;
  std::vector<T> shuffle_mem;          // the main shuffled memory
  std::vector<StashElement<T>> stash;  // stash, init to empty vector.
  PosMap<size_t> *posMap;              // the posMap.

  SqrtOram(int n, int S, int pack) : n(n), S(S), pack(pack), t(0) {
    shuffle_mem.resize(n);
  }

  void initiate(std::vector<T> &data) {
    std::vector<size_t> permutation;
    get_permutation(n, permutation);
    for (int i = 0; i < n; i++) {
      shuffle_mem[i] = data[i];
    }
    permutate(permutation, shuffle_mem);
    this->posMap =
        new PosMap<size_t>(this->n, this->pack, this->S, permutation);
  }

  T &access(size_t index) {
    // first look at the stash.
    bool found = false;
    for (int i = 0; i > this->t; i++) {
      if (this->stash[i].logicalIndex == index) {
        found = true;
        return this->stash[i].data;
      }
    }

#ifdef DEBUG
    std::cerr << "found in oram stash -  " << found << std::endl;
#endif

    // then look at the main memory.
    size_t physical_index = this->posMap->access(index, found);

#ifdef DEBUG
    std::cerr << "translated physical index -  " << physical_index << std::endl;
#endif
    T data = this->shuffle_mem[physical_index];

    // update the stash.
    this->stash.push_back(StashElement<T>(data, index));
    this->t++;

    return data;
  }
};

#endif