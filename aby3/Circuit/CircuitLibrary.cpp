#include "CircuitLibrary.h"

#include <cryptoTools/Circuit/BetaLibrary.h>
#include <cryptoTools/Crypto/RandomOracle.h>

namespace aby3
{

	namespace {

		template<typename THead>
		void _hash(oc::RandomOracle& ro, THead h)
		{
			ro.Update(h);
			//auto hh = std::hash<THead>()(h);
			//return hh;
		}
		template<typename THead, typename... TTail>
		void _hash(oc::RandomOracle& ro, THead h, TTail... tail)
		{
			_hash(ro, h);
			_hash<TTail...>(ro, tail...);
		}

		template<typename... TTail>
		size_t hash(TTail... tail)
		{
			oc::RandomOracle ro(sizeof(size_t));

			_hash<TTail...>(ro, tail...);
			size_t ret;
			ro.Final(ret);
			return ret;
		}

	}

	oc::BetaCircuit * CircuitLibrary::int_Sh3Piecewise_helper(u64 size, u64 numThesholds)
	{
		size_t key = hash(__FUNCTION__, size, numThesholds);

		auto iter = mCirMap.find(key);

		if (iter == mCirMap.end())
		{
			auto* cd = new BetaCircuit;

			std::vector<BetaBundle> aa(numThesholds);
			for (auto& a : aa)
			{
				a.mWires.resize(size);
				cd->addInputBundle(a);
			}

			BetaBundle b(size);
			cd->addInputBundle(b);


			std::vector<BetaBundle> cc(numThesholds + 1);
			for (auto& c : cc)
			{
				c.mWires.resize(1);
				cd->addOutputBundle(c);
			}

			int_Sh3Piecewise_build_do(*cd,
				aa,
				b,
				cc);

			//cd->levelByAndDepth();

			iter = mCirMap.insert(std::make_pair(key, cd)).first;
		}

		return iter->second;
	}

	void CircuitLibrary::int_Sh3Piecewise_build_do(
		BetaCircuit & cd,
		span<const BetaBundle> aa,
		const BetaBundle & b,
		span<const BetaBundle> cc)
	{

		std::vector<BetaBundle>
			temps(aa.size()),
			thresholds(aa.size());

		for (u64 t = 0; t < thresholds.size(); ++t)
		{
			thresholds[t].mWires.resize(1);
			temps[t].mWires.resize(b.mWires.size() * 2);
			cd.addTempWireBundle(temps[t]);
		}

		for (u64 t = 1; t < thresholds.size() - 1; ++t)
			cd.addTempWireBundle(thresholds[t]);

		// the first region bit is just the thrshold. 
		thresholds[0] = cc[0];

		// compute all the signs
		for (u64 t = 0; t < thresholds.size(); ++t)
		{
			if (aa[t].size() != b.size())
				throw RTE_LOC;

			//int_int_add_msb
			extractBit_build(cd, aa[t], b, thresholds[t], temps[t], b.size() - 1, 
				IntType::TwosComplement, 
				AdderType::Addition,
				Optimized::Depth);
			//int_int_add_msb_build_do(cd, aa[t], b, thresholds[t], temps[t]);

			//std::cout << t << " @ " << thresholds[t].mWires[0] << std::endl;
			//cd.addPrint(std::to_string(t) + "  ");
			//cd.addPrint(thresholds[t]);
			//cd.addPrint("\n");
		}


		// Take the and of the signs
		for (u64 t = 1; t < thresholds.size(); ++t)
		{
			cd.addGate(
				thresholds[t - 1].mWires[0],
				thresholds[t].mWires[0],
				oc::GateType::na_And,
				cc[t].mWires[0]);
		}

		// mark the last region as the inverse of the last threshold bit.
		cd.addInvert(
			thresholds.back().mWires[0],
			cc[thresholds.size()].mWires[0]);
	}

	void CircuitLibrary::Preproc_build(BetaCircuit& cd, u64 dec)
	{

		const auto size = sizeof(i64) * 8;
		BetaBundle a(size);
		BetaBundle b0(size);
		BetaBundle c0(size);
		BetaBundle b1(size - dec);
		BetaBundle c1(size - dec);

		cd.addInputBundle(a);
		//std::cout << "a " << a.mWires[0] << "-" << a.mWires.back() << std::endl;
		cd.addInputBundle(b0);
		cd.addInputBundle(b1);
		cd.addOutputBundle(c0);
		cd.addOutputBundle(c1);

		BetaBundle t(3);
		cd.addTempWireBundle(t);

		//lib.int_int_add_build_so(cd, a, b0, c0, t);
		add_build(cd, a, b0, c0, t, IntType::TwosComplement, Optimized::Size);

		a.mWires.erase(a.mWires.begin(), a.mWires.begin() + dec);

		//lib.int_int_add_build_so(cd, a, b1, c1, t);
		add_build(cd, a, b1, c1, t, IntType::TwosComplement, Optimized::Size);

		if (cd.mNonlinearGateCount != 2 * (size - 1) - dec)
			throw std::runtime_error(LOCATION);

	}

	void CircuitLibrary::argMax_build(BetaCircuit& cd, u64 dec, u64 numArgs)
	{
		auto size = 64 - dec;
		std::vector<BetaBundle> a0(numArgs), a1(numArgs);

		std::vector<BetaBundle> a(numArgs);

		//std::vector<std::vector<BetaBundle>> select(log2ceil(numArgs));

		for (u64 i = 0; i < numArgs; ++i)
		{
			a[i].mWires.resize(size);
			a0[i].mWires.resize(size);
			a1[i].mWires.resize(size);

			cd.addInputBundle(a0[i]);
			cd.addInputBundle(a1[i]);
		}


		BetaBundle argMax(oc::log2ceil(numArgs));
		cd.addOutputBundle(argMax);

		// set argMax to zero initially...
		for (auto& w : argMax.mWires) cd.addConst(w, 0);


		BetaBundle t(3);
		cd.addTempWireBundle(t);
		for (u64 i = 0; i < numArgs; ++i)
		{
			cd.addTempWireBundle(a[i]);
			//lib.int_int_add_build_so(cd, a0[i], a1[i], a[i], t);
			add_build(cd, a0[i], a1[i], a[i], t, IntType::TwosComplement, Optimized::Size);

			cd.addPrint("a[" + std::to_string(i) + "] = ");
			cd.addPrint(a[i]);
			cd.addPrint("\n");
		}

		// maxPointer will equal 1 if the the rhs if greater.
		BetaBundle maxPointer(1);
		maxPointer[0] = argMax[0];
		auto& max = a[0];

		//lib.int_int_lt_build(cd, max, a[1], maxPointer);
		lessThan_build(cd, max, a[1], maxPointer,
			IntType::TwosComplement,
			Optimized::Depth);

		for (u64 i = 2; i < numArgs; ++i)
		{
			// currently, max = max(a[0],...,a[i-2]). Now make 
			//		max = max(a[0],...,a[i-1]) 
			//      max = maxPointer ? a[i - 1] : max;
			//lib.int_int_multiplex_build(cd, a[i - 1], max, maxPointer, max, t);
			multiplex_build(cd, a[i - 1], max, maxPointer, max, t);

			cd.addPrint("max({0, ...," + std::to_string(i - 1) + " }) = ");
			cd.addPrint(max);
			cd.addPrint(" @ ");
			cd.addPrint(argMax);
			cd.addPrint(" * ");
			cd.addPrint(maxPointer);
			cd.addPrint("\n");

			// now compute which is greater (a[i], max) and store the result in maxPointer
			cd.addTempWireBundle(maxPointer);
			//lib.int_int_lt_build(cd, max, a[i], maxPointer);
			lessThan_build(cd, max, a[i], maxPointer, 
				IntType::TwosComplement, 
				Optimized::Depth);

			// construct a const wire bundle that encodes the index i.
			BetaBundle idx(oc::log2ceil(numArgs));
			cd.addConstBundle(idx, oc::BitVector((u8*)& i, oc::log2ceil(numArgs)));

			// argMax = (max < a[i]) ? i : argMax ;
			// argMax =   maxPointer ? i : argMax ;
			multiplex_build(cd, idx, argMax, maxPointer, argMax, t);
		}

		cd.addPrint("max({0, ...," + std::to_string(numArgs) + " }) = _____ @ ");
		cd.addPrint(argMax);
		cd.addPrint("\n");

	}

    oc::BetaCircuit* CircuitLibrary::convert_arith_to_bin(u64 n, u64 bits)
    {
        //auto key = "convert_arith_to_bin" + std::to_string(n) + "_" + std::to_string(bits);
		auto key = hash(__FUNCTION__, n, bits);

        auto iter = mCirMap.find(key);

        if (iter == mCirMap.end())
        {
            auto* cd = new BetaCircuit;

            std::array<oc::BetaBundle, 2> inputs, inSubs;
            oc::BetaBundle outputs, outSub;

            inputs[0].mWires.resize(bits * n);
            inputs[1].mWires.resize(bits * n);
            outputs.mWires.resize(bits * n);
            cd->addInputBundle(inputs[0]);
            cd->addInputBundle(inputs[1]);
            cd->addOutputBundle(outputs);

            BetaBundle temp(bits * 2);
            auto iter0 = inputs[0].mWires.begin();
            auto iter1 = inputs[1].mWires.begin();
            auto iter2 = outputs.mWires.begin();
            for (u64 i = 0; i < inputs.size(); ++i)
            {
                inSubs[0].mWires.clear();
                inSubs[1].mWires.clear();
                outSub.mWires.clear();

                inSubs[0].mWires.insert(inSubs[0].mWires.end(), iter0, iter0 + bits);
                inSubs[1].mWires.insert(inSubs[1].mWires.end(), iter1, iter1 + bits);
                outSub.mWires.insert(   outSub.mWires.end(),    iter2, iter2 + bits);

                iter0 += bits;
                iter1 += bits;
                iter2 += bits;

                cd->addTempWireBundle(temp);
                //int_int_add_build_do(*cd, inSubs[0], inSubs[1], outSub, temp);
				add_build(*cd, inSubs[0], inSubs[1], outSub, temp,
					IntType::TwosComplement,
					Optimized::Depth);
            }

            iter = mCirMap.insert(std::make_pair(key, cd)).first;
        }

        return iter->second;
    }

	oc::BetaCircuit* CircuitLibrary::gt_build_do(){
		// auto key = "gt_circuit";
		auto key = hash(__FUNCTION__);
		auto iter = mCirMap.find(key);
		const auto unit_size = sizeof(i64) * 8;

		if (iter == mCirMap.end())
        {
			auto* cd = new BetaCircuit;
			BetaBundle ba(unit_size);
			BetaBundle bb(unit_size);
			BetaBundle msb(1);
			BetaBundle temps(bb.mWires.size() * 2);

			// temps.mWires.resize(bb.mWires.size() * 2);
			std::cout << "mWires size ba: " << ba.mWires.size() << std::endl; // 64
			std::cout << "mWires size bb: " << bb.mWires.size() << std::endl; // 64
			std::cout << "output size: " << msb.mWires.size() << std::endl; // 1
			std::cout << "tmp size: " << temps.mWires.size() << std::endl; // 128

			cd->addInputBundle(ba);
			cd->addInputBundle(bb);
			cd->addOutputBundle(msb);
			cd->addTempWireBundle(temps);

			std::cout << "succeed here" << std::endl;
			// oc::BetaLibrary::int_int_add_msb_build_do(*cd, ba, bb, msb, temps);
			extractBit_build(*cd, ba, bb, msb, temps, bb.size() - 1, 
				IntType::TwosComplement, 
				AdderType::Addition,
				Optimized::Depth);
			std::cout << " after build do " << std::endl;
			iter = mCirMap.insert(std::make_pair(key, cd)).first;
		}
		return iter->second;
	}


	oc::BetaCircuit * CircuitLibrary::int_comp_helper(u64 size)
	{
		// auto key = "int_comp_helper" + std::to_string(size);
		auto key = hash(__FUNCTION__);
		auto iter = mCirMap.find(key);

		if (iter == mCirMap.end())
		{
			auto* cd = new BetaCircuit;

			BetaBundle aa(size);
			cd->addInputBundle(aa);
			BetaBundle b(size);
			cd->addInputBundle(b);

			BetaBundle cc(1);
			cd->addOutputBundle(cc);

			int_comp_build_do(*cd,
				aa,
				b,
				cc);

			iter = mCirMap.insert(std::make_pair(key, cd)).first;
		}

		return iter->second;
	}

	void CircuitLibrary::int_comp_build_do(
		BetaCircuit & cd,
		const BetaBundle & aa,
		const BetaBundle & b,
		const BetaBundle & cc)
	{

		BetaBundle temps(aa.mWires.size()*2), thresholds(1);
		cd.addTempWireBundle(temps);
		thresholds = cc;
		// oc::BetaLibrary::int_int_add_msb_build_do(cd, aa, b, thresholds, temps);
		extractBit_build(cd, aa, b, thresholds, temps, aa.size() - 1, 
			IntType::TwosComplement, 
			AdderType::Addition,
			Optimized::Depth);
	}

	oc::BetaCircuit * CircuitLibrary::bits_nor_helper(u64 size){
		// auto key = "bit_or_helper" + std::to_string(size);
		auto key = hash(__FUNCTION__);
		auto iter = mCirMap.find(key);
		if (iter == mCirMap.end()){
			auto* cd = new BetaCircuit;

			BetaBundle a(size);
			cd->addInputBundle(a);
			BetaBundle b(size);
			cd->addInputBundle(b);
			BetaBundle c(size);
			cd->addOutputBundle(c);
			bits_nor_build_do(*cd,
				a,
				b,
				c);

			iter = mCirMap.insert(std::make_pair(key, cd)).first;
		}
		return iter->second;
	}

	void CircuitLibrary::bits_nor_build_do(
		BetaCircuit & cd,
		const BetaBundle & a,
		const BetaBundle & b,
		const BetaBundle & c)
	{
		for(u64 i=0; i < a.mWires.size(); i++){
			cd.addGate(a.mWires[i], b.mWires[i], oc::GateType::Nor, c.mWires[i]);
		}
	}
}