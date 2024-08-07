#pragma once
#include <aby3/sh3/Sh3Encryptor.h>
#include <aby3/sh3/Sh3Evaluator.h>
#include <aby3/sh3/Sh3FixedPoint.h>
#include <aby3/sh3/Sh3Runtime.h>
#include <aby3/sh3/Sh3Types.h>
#include <cryptoTools/Network/IOService.h>
#include <mpi.h>
#include <iomanip>

#include "./Pair_then_Reduce/include/datatype.h"
#include "./Pair_then_Reduce/include/tasks.h"
#include "BuildingBlocks.h"
#include "debug.h"

// aby3 subtask definition, currently we do not mainain this class.
template <typename NUMX, typename NUMY, typename NUMT, typename NUMR>
class SubABY3Task : public SubTask<NUMX, NUMY, NUMT, NUMR>{
public:
    // aby3 info.
    int pIdx;
    aby3::Sh3Encryptor* enc;
    aby3::Sh3Runtime* runtime;
    aby3::Sh3Evaluator* eval;

    size_t m_start = 0;
    size_t m_end = 0;

    using SubTask<NUMX, NUMY, NUMT, NUMR>::SubTask;

    SubABY3Task(const size_t optimal_block, const int task_id, const int pIdx,
            aby3::Sh3Encryptor& enc, aby3::Sh3Runtime& runtime,
            aby3::Sh3Evaluator& eval)
        : pIdx(pIdx),
            enc(&enc),
            runtime(&runtime),
            eval(&eval),
            SubTask<NUMX, NUMY, NUMT, NUMR>(optimal_block, task_id) {
    }

    size_t get_partial_m_lens(){
        this->m_start = this->table_start / this->n;
        this->m_end = (this->table_end - 1) / this->n;
        size_t partial_len = this->m_end - this->m_start + 1;
        if(this->m_end >= this->m - 1){
            return partial_len;
        }
        else{
            partial_len += this->lookahead;
        }
        return partial_len;
    }
};


// mpi function definition.
template <typename NUMX, typename NUMY, typename NUMT, typename NUMR,
          template <typename, typename, typename, typename> class TASK>
class ABY3MPITask : public MPIPTRTask<NUMX, NUMY, NUMT, NUMR, TASK> {
 public:
  // aby3 info
  int pIdx;
  aby3::Sh3Encryptor& enc;
  aby3::Sh3Runtime& runtime;
  aby3::Sh3Evaluator& eval;

  // setup all the aby3 environment variables, pIdx and rank.
  using MPIPTRTask<NUMX, NUMY, NUMT, NUMR, TASK>::MPIPTRTask;
  ABY3MPITask(int tasks, size_t optimal_block, const int pIdx,
                 aby3::Sh3Encryptor& enc, aby3::Sh3Runtime& runtime,
                 aby3::Sh3Evaluator& eval)
      : pIdx(pIdx),
        enc(enc),
        runtime(runtime),
        eval(eval),
        MPIPTRTask<NUMX, NUMY, NUMT, NUMR, TASK>(tasks, optimal_block) {}

  // override sub_task create.
  void create_sub_task(size_t optimal_block, int task_id, size_t table_start,
                       size_t table_end) override {
    auto subTaskPtr(
        new TASK<NUMX, NUMY, NUMT, NUMR>(optimal_block, task_id, this->pIdx,
                                         this->enc, this->runtime, this->eval));
    this->subTask.reset(subTaskPtr);
    this->subTask->initial_value = this->default_value;
    this->subTask->circuit_construct(this->shapeX, this->shapeY, table_start,
                                     table_end);
    for (int j = 0; j < this->n; j++)
      this->subTask->res[j] = this->default_value;
  }
};

// mpi function definition.
template <typename NUMX, typename NUMY, typename NUMT, typename NUMR,
          template <typename, typename, typename, typename> class TASK>
class ABY3MPIPairOnlyTask : public PairOnlyMPIPTRTask<NUMX, NUMY, NUMT, NUMR, TASK> {
 public:
  // aby3 info
  int pIdx;
  aby3::Sh3Encryptor& enc;
  aby3::Sh3Runtime& runtime;
  aby3::Sh3Evaluator& eval;

  // setup all the aby3 environment variables, pIdx and rank.
  using PairOnlyMPIPTRTask<NUMX, NUMY, NUMT, NUMR, TASK>::PairOnlyMPIPTRTask;
  ABY3MPIPairOnlyTask(int tasks, size_t optimal_block, const int pIdx,
                 aby3::Sh3Encryptor& enc, aby3::Sh3Runtime& runtime,
                 aby3::Sh3Evaluator& eval)
      : pIdx(pIdx),
        enc(enc),
        runtime(runtime),
        eval(eval),
        PairOnlyMPIPTRTask<NUMX, NUMY, NUMT, NUMR, TASK>(tasks, optimal_block) {}

  // override sub_task create.
  void create_sub_task(size_t optimal_block, int task_id, size_t table_start,
                       size_t table_end) override {
    auto subTaskPtr(
        new TASK<NUMX, NUMY, NUMT, NUMR>(optimal_block, task_id, this->pIdx,
                                         this->enc, this->runtime, this->eval));
    this->subTask.reset(subTaskPtr);
    this->subTask->initial_value = this->default_value;
    this->subTask->circuit_construct(this->shapeX, this->shapeY, table_start,
                                     table_end);
    for (int j = 0; j < this->n; j++)
      this->subTask->res[j] = this->default_value;
  }
};