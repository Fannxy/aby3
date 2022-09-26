#pragma once
#include <cryptoTools/Circuit/BetaCircuit.h>
#include <cryptoTools/Common/Defines.h>
#include <unordered_map>
#include "aby3/Common/Defines.h"
#include <cryptoTools/Circuit/BetaLibrary.h>

namespace aby3
{
    class CircuitLibrary : public oc::BetaLibrary
    {
    public:
		using BetaCircuit = oc::BetaCircuit;
		using BetaBundle = oc::BetaBundle;


		BetaCircuit* int_Sh3Piecewise_helper(u64 aSize, u64 numThesholds);

        BetaCircuit* convert_arith_to_bin(u64 n, u64 bits);

		static void int_Sh3Piecewise_build_do(
			BetaCircuit& cd,
			span<const BetaBundle> aa,
			const BetaBundle & b,
			span<const BetaBundle> c);

		static void Preproc_build(BetaCircuit& cd, u64 dec);
		static void argMax_build(BetaCircuit& cd, u64 dec, u64 numArgs);
		BetaCircuit* gt_build_do();

		BetaCircuit* int_comp_helper(u64 aSize);
        static void int_comp_build_do(
			BetaCircuit& cd,
			const BetaBundle & aa,
			const BetaBundle & b,
			const BetaBundle & cc);

		BetaCircuit* bits_nor_helper(u64 size);
		static void bits_nor_build_do(BetaCircuit& cd,
			const BetaBundle & a,
			const BetaBundle & b,
			const BetaBundle & c);

    };

}