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

    BASIC_TEST_INIT

    auto N = cmd.getOr<int>("dataSize", 10000);
    auto dim = cmd.getOr<int>("D", 1000);
    auto B = cmd.getOr<int>("B", 128);
    auto IT = cmd.getOr<int>("I", 10000);
    
    PRNG prng(toBlock(1));

    RegressionParam params;
    params.mBatchSize = B;
    params.mIterations = IT;
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

    u64 aB = std::log2(1 / (params.mLearningRate / params.mBatchSize));

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

	std::vector<int> role;
	role = cmd.getMany<int>("role");

	std::vector<int> rank;
	rank = cmd.getMany<int>("rank");
	if(rank.size() != role.size()) {
		throw std::runtime_error(LOCATION);
	}

	std::vector<std::string> p0_ip, p1_ip;
	p0_ip = cmd.getMany<std::string>("p0_ip");
	p1_ip = cmd.getMany<std::string>("p1_ip");
	if(p0_ip.size() != role.size() || p1_ip.size() != role.size()) {
		throw std::runtime_error(LOCATION);
	}

	std::vector<IOService> ios(role.size());
	std::vector<Sh3Encryptor> enc(role.size());
	std::vector<Sh3Evaluator> eval(role.size());
	std::vector<Sh3Runtime> runtime(role.size());
	for(size_t i=0; i<role.size(); i++) {
		splitted_setup((u64)role[i], rank[i], ios[i], enc[i], eval[i], runtime[i], p0_ip[i], p1_ip[i]);
	}

	std::vector<u64> dataSize;
	dataSize = cmd.getMany<u64>("dataSize");
	if(dataSize.size() != role.size()) {
		throw std::runtime_error(LOCATION);
	}

	u64 N = 0;
	for(size_t i=0; i<role.size(); i++) {
		N += dataSize[i];
	}
    auto dim = cmd.getOr<int>("D", 1000);
    auto B = cmd.getOr<int>("B", 128);
    auto IT = cmd.getOr<int>("I", 1000);
    
    PRNG prng(toBlock(1));

    RegressionParam params;
    params.mBatchSize = B;
    params.mIterations = IT;
    params.mLearningRate = 1.0 / (1 << 3);

    const Decimal D = D8;
    sf64Matrix<D> X(N, dim), Y(N, 1), w(dim, 1);

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

	u64 aB = std::log2(1 / (params.mLearningRate / params.mBatchSize));

    std::vector<std::thread> threads;

    sf64Matrix<D> xw(N, 1);
    if(phase == 0) {
        auto& mCastX = (si64Matrix&) X;
        auto& mCastY = (si64Matrix&) Y;
        auto& mCastw = (si64Matrix&) w;
        for(size_t i=0; i<mCastX.size(); i++) {
            mCastX.mShares[0](i) = 0;
            mCastX.mShares[1](i) = 0;
        }
        for(size_t i=0; i<mCastY.size(); i++) {
            mCastY.mShares[0](i) = 0;
            mCastY.mShares[1](i) = 0;
        }
        for(size_t i=0; i<mCastw.size(); i++) {
            mCastw.mShares[0](i) = 0;
            mCastw.mShares[1](i) = 0;
        }
        threads.clear();
        std::vector<sf64Matrix<D>> subX(role.size());
        std::vector<sf64Matrix<D>> subxw(role.size());
        int offset = 0;
        for(int i=0; i<role.size(); i++) {
            subX[i].resize(dataSize[i], dim);
            for(size_t j=0; j<subX[i].rows(); j++) {
                for(size_t k=0; k<subX[i].cols(); k++) {
                    subX[i](j, k, X(offset + j, k));
                }
            }
            subxw[i].resize(dataSize[i], 1);
            threads.push_back(std::thread(fixed_matrix_mult, std::ref(subX[i]), std::ref(w), std::ref(subxw[i]), 0, std::ref(enc[i]), std::ref(eval[i]), std::ref(runtime[i])));
            offset += dataSize[i];
        }
        for(auto& th : threads) {
            th.join();
        }
        offset = 0;
        for(int i=0; i<role.size(); i++) {
            for(size_t j=0; j<subxw[i].rows(); j++) {
                xw(offset + j, 0, subxw[i](j, 0));
            }
            offset += dataSize[i];
        }
    }

	sf64Matrix<D> fxw(X.rows(), 1);
    if(phase == 1) {
        auto& mCastxw = (si64Matrix&) xw;
        for(size_t i=0; i<mCastxw.size(); i++) {
            mCastxw.mShares[0](i) = 0;
            mCastxw.mShares[1](i) = 0;
        }
        threads.clear();
        std::vector<sf64Matrix<D>> subxw(role.size());
        std::vector<sf64Matrix<D>> subfxw(role.size());
        int offset = 0;
        for(int i=0; i<role.size(); i++) {
            subxw[i].resize(dataSize[i], 1);
            for(size_t j=0; j<subxw[i].rows(); j++) {
                subxw[i](j, 0, xw(offset + j, 0));
            }
            subfxw[i].resize(dataSize[i], 1);
            threads.push_back(std::thread(logisticFunc, std::ref(subxw[i]), std::ref(subfxw[i]), std::ref(mLogistic), 0, std::ref(enc[i]), std::ref(eval[i]), std::ref(runtime[i])));
            offset += dataSize[i];
        }
        for(auto& th : threads) {
            th.join();
        }
        offset = 0;
        for(int i=0; i<role.size(); i++) {
            for(size_t j=0; j<subfxw[i].rows(); j++) {
                fxw(offset + j, 0, subfxw[i](j, 0));
            }
            offset += dataSize[i];
        }
    }

    if(phase == 2) {
        auto& mCastX = (si64Matrix&) X;
        auto& mCastY = (si64Matrix&) Y;
        auto& mCastfxw = (si64Matrix&) fxw;
        for(size_t i=0; i<mCastX.size(); i++) {
            mCastX.mShares[0](i) = 0;
            mCastX.mShares[1](i) = 0;
        }
        for(size_t i=0; i<mCastY.size(); i++) {
            mCastY.mShares[0](i) = 0;
            mCastY.mShares[1](i) = 0;
        }
        for(size_t i=0; i<mCastfxw.size(); i++) {
            mCastfxw.mShares[0](i) = 0;
            mCastfxw.mShares[1](i) = 0;
        }

        sf64Matrix<D> error = fxw - Y;

        X.transposeInPlace();

        threads.clear();
        std::vector<sf64Matrix<D>> subX(role.size());
        std::vector<sf64Matrix<D>> suberror(role.size());
        std::vector<sf64Matrix<D>> subupdate(role.size());
        int offset = 0;
        for(int i=0; i<role.size(); i++) {
            subX[i].resize(dataSize[i], dim);
            for(size_t j=0; j<subX[i].rows(); j++) {
                for(size_t k=0; k<subX[i].cols(); k++) {
                    subX[i](j, k, X(offset + j, k));
                }
            }
            suberror[i].resize(dataSize[i], 1);
            for(size_t j=0; j<suberror[i].rows(); j++) {
                suberror[i](j, 0, error(offset + j, 0));
            }
            subupdate[i].resize(dim, 1);
            threads.push_back(std::thread(fixed_matrix_mult_shift, std::ref(subX[i]), std::ref(suberror[i]), std::ref(subupdate[i]), aB, 0, std::ref(enc[i]), std::ref(eval[i]), std::ref(runtime[i])));
            offset += dataSize[i];
        }
        for(auto& th : threads) {
            th.join();
        }

        sf64Matrix<D> update(dim, 1);
        for(int i=0; i<role.size(); i++) {
            update = update + subupdate[i];
        }

        w = w - update;
    }

    return 0;
}