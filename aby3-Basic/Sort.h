#include "Basics.h"
#include "Shuffle.h"

#include <aby3/Circuit/CircuitLibrary.h>
#include "../aby3-RTR/debug.h"

int bc_sort_different(std::vector<aby3::sbMatrix> &data, std::vector<size_t> &lows, std::vector<size_t> &highs, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime, size_t max_size);

int bc_sort_different(aby3::si64Matrix &data, std::vector<size_t> &lows, std::vector<size_t> &highs, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime, size_t max_size);

int quick_sort_different(std::vector<aby3::sbMatrix> &data, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime, size_t min_size);

int quick_sort(std::vector<aby3::sbMatrix> &data, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime, size_t min_size);

int quick_sort_different(aby3::si64Matrix& data, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime, size_t min_size);

int quick_sort(aby3::si64Matrix& data, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime, size_t min_size);

int odd_even_merge(aby3::sbMatrix& data1, aby3::sbMatrix& data2, aby3::sbMatrix& res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

int odd_even_merge(aby3::si64Matrix& data1, aby3::si64Matrix& data2, aby3::si64Matrix& res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

int odd_even_multi_merge(std::vector<aby3::sbMatrix> &data, aby3::sbMatrix& sorted_res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

int odd_even_multi_merge(std::vector<aby3::si64Matrix> &data, aby3::si64Matrix& sorted_res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

/*
data1: k * m1.
data2: k * m2.
sorted_res: k * (m1 + m2).
*/
int high_dimensional_odd_even_merge(std::vector<aby3::sbMatrix> &data1, std::vector<aby3::sbMatrix>& data2, std::vector<aby3::sbMatrix>& sorted_res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);

/*
    data - (k * {m_1, m_2, ..., m_n})
    sorted_res - (k, \sum{m_i}).
    This function is used to merge the data in the k-dimensional space, accepting k sorted sub-arrays with size {m_1, ..., m_n}, and outputting k sorted array with size \sum{m_i}.
*/
int high_dimensional_odd_even_multi_merge(std::vector<std::vector<aby3::sbMatrix>> &data, std::vector<aby3::sbMatrix>& sorted_res, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);