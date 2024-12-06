#include "Test.h"

#include <chrono>
#include <random>
#include <thread>

#include "../aby3-Basic/Basics.h"
#include "../aby3-ML/Regression.h"
#include "../aby3-ML/LinearModelGen.h"
#include "../aby3-ML/aby3ML.h"
#include "../aby3-ML/PlainML.h"
#include "../aby3-RTR/BuildingBlocks.h"
#include "../aby3-Basic/Matrix.h"


using namespace oc;
using namespace aby3;

int logisticFunc(sf64Matrix<D8>& X, sf64Matrix<D8>& Y, Sh3Piecewise& mLogistic, int pIdx, Sh3Encryptor& enc, Sh3Evaluator& eval, Sh3Runtime& runtime) {
	if (X.rows() != Y.rows() || X.cols() != Y.cols()) {
		THROW_RUNTIME_ERROR("Matrix dimensions do not match.");
	}

	mLogistic.eval<D8>(runtime, X, Y, eval);
	
	return 0;
}

int logistic_regression_test(oc::CLP& cmd) {

    SPLITTED_TEST_INIT

    auto N = cmd.getOr<int>("dataSize", 10000);
    auto dim = cmd.getOr<int>("D", 16);
    
    PRNG prng(toBlock(1));

    RegressionParam params;
    params.mLearningRate = 1.0 / (1 << 3);

    const Decimal D = D8;
    sf64Matrix<D> X(N, dim), Y(N, 1), w(dim, 1);

	auto& mCastX = (si64Matrix&) X;
    auto& mCastY = (si64Matrix&) Y;
	auto& mCastw = (si64Matrix&) w;
    for(i64 i=0; i<mCastX.size(); i++) {
		mCastX.mShares[0](i) = 0;
		mCastX.mShares[1](i) = 0;
    }
	for(i64 i=0; i<mCastY.size(); i++) {
		mCastY.mShares[0](i) = 0;
		mCastY.mShares[1](i) = 0;
	}
	for(i64 i=0; i<mCastw.size(); i++) {
		mCastw.mShares[0](i) = 0;
		mCastw.mShares[1](i) = 0;
	}

	Sh3Piecewise mLogistic;
	mLogistic.mThresholds.resize(2);
	mLogistic.mThresholds[0] = -0.5;
	mLogistic.mThresholds[1] = 0.5;
	mLogistic.mCoefficients.resize(3);
	mLogistic.mCoefficients[1].resize(2);
	mLogistic.mCoefficients[1][0] = 0.5;
	mLogistic.mCoefficients[1][1] = 1;
	mLogistic.mCoefficients[2].resize(1);
	mLogistic.mCoefficients[2][0] = 1;

    u64 aB = std::log2(1 / (params.mLearningRate / N));

    sf64Matrix<D> xw(X.rows(), 1);
	fixed_matrix_mult(X, w, xw, 0, enc, eval, runtime);

	sf64Matrix<D> fxw(X.rows(), 1);
	logisticFunc(xw, fxw, mLogistic, 0, enc, eval, runtime);

	sf64Matrix<D> error = fxw - Y;

	X.transposeInPlace();

	sf64Matrix<D> update(X.rows(), 1);
	fixed_matrix_mult_shift(X, error, update, aB, 0, enc, eval, runtime);

	w = w - update;

    return 0;
}

int splitted_logistic_regression_test(oc::CLP& cmd) {

	int phase = -1;
    if(cmd.isSet("LogReg-0")) {
        phase = 0;
    } else if(cmd.isSet("LogReg-1")) {
        phase = 1;
    } else if(cmd.isSet("LogReg-2")) {
        phase = 2;
    } else {
        throw std::runtime_error(LOCATION);
    }

    SPLITTED_TEST_INIT

    auto N = cmd.getOr<int>("dataSize", 10000);
    auto dim = cmd.getOr<int>("D", 16);
    
    PRNG prng(toBlock(1));

    RegressionParam params;
    params.mLearningRate = 1.0 / (1 << 3);

    const Decimal D = D8;
    sf64Matrix<D> X(N, dim), Y(N, 1), w(dim, 1);

	auto& mCastX = (si64Matrix&) X;
    auto& mCastY = (si64Matrix&) Y;
	auto& mCastw = (si64Matrix&) w;
    for(i64 i=0; i<mCastX.size(); i++) {
		mCastX.mShares[0](i) = 0;
		mCastX.mShares[1](i) = 0;
    }
	for(i64 i=0; i<mCastY.size(); i++) {
		mCastY.mShares[0](i) = 0;
		mCastY.mShares[1](i) = 0;
	}
	for(i64 i=0; i<mCastw.size(); i++) {
		mCastw.mShares[0](i) = 0;
		mCastw.mShares[1](i) = 0;
	}

	Sh3Piecewise mLogistic;
	mLogistic.mThresholds.resize(2);
	mLogistic.mThresholds[0] = -0.5;
	mLogistic.mThresholds[1] = 0.5;
	mLogistic.mCoefficients.resize(3);
	mLogistic.mCoefficients[1].resize(2);
	mLogistic.mCoefficients[1][0] = 0.5;
	mLogistic.mCoefficients[1][1] = 1;
	mLogistic.mCoefficients[2].resize(1);
	mLogistic.mCoefficients[2][0] = 1;

    u64 aB = std::log2(1 / (params.mLearningRate / N));

    sf64Matrix<D> xw(X.rows(), 1);
    auto &mCastxw = (si64Matrix&) xw;
    for(i64 i=0; i<mCastxw.size(); i++) {
        mCastxw.mShares[0](i) = 0;
        mCastxw.mShares[1](i) = 0;
    }
    if(phase == 0) {
        fixed_matrix_mult(X, w, xw, 0, enc, eval, runtime);
    }

	sf64Matrix<D> fxw(X.rows(), 1);
    auto &mCastfxw = (si64Matrix&) fxw;
    for(i64 i=0; i<mCastfxw.size(); i++) {
        mCastfxw.mShares[0](i) = 0;
        mCastfxw.mShares[1](i) = 0;
    }
    if(phase == 1) {
        logisticFunc(xw, fxw, mLogistic, 0, enc, eval, runtime);
    }

    if(phase == 2) {
        sf64Matrix<D> error = fxw - Y;

        X.transposeInPlace();

        sf64Matrix<D> update(X.rows(), 1);
        fixed_matrix_mult_shift(X, error, update, aB, 0, enc, eval, runtime);

        w = w - update;
    }

    return 0;
}