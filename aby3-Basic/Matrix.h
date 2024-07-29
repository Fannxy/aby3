#include "Basics.h"

#include <aby3/Circuit/CircuitLibrary.h>
#include "../aby3-RTR/debug.h"

int fixed_matrix_mult(aby3::sf64Matrix<aby3::D8>& A, aby3::sf64Matrix<aby3::D8>& B, aby3::sf64Matrix<aby3::D8>& C, int pIdx, aby3::Sh3Encryptor& enc, aby3::Sh3Evaluator& eval, aby3::Sh3Runtime& runtime);