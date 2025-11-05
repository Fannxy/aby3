#pragma once
#include <aby3/sh3/Sh3Encryptor.h>
#include <aby3/sh3/Sh3Evaluator.h>
#include <aby3/sh3/Sh3FixedPoint.h>
#include <aby3/sh3/Sh3Runtime.h>
#include <aby3/sh3/Sh3Types.h>
#include <cryptoTools/Common/CLP.h>
#include <cryptoTools/Network/IOService.h>

#include "../aby3-RTR/debug.h"

int oblivious_idx_select_test(oc::CLP &cmd);

int genperm_test(oc::CLP &cmd);

int feddb_shuffle_test(oc::CLP &cmd);

int persist_test(oc::CLP &cmd);

int index_agg_test(oc::CLP &cmd);

int group_by_common_test(oc::CLP &cmd);

int group_by_test(oc::CLP &cmd);

int odd_even_merge_sort_test(oc::CLP &cmd);

int join_test(oc::CLP &cmd);