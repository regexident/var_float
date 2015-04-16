#include <iostream>

#define CATCH_CONFIG_MAIN

#include <catch.hpp>

#include "var_float.h"

namespace var_float_c_tests {

	const vf_profile_t profile_8 = vf_profile_default(8);
	const vf_profile_t profile_16 = vf_profile_default(16);
	const vf_profile_t profile_24 = vf_profile_default(24);
	const vf_profile_t profile_32 = vf_profile_default(32);
	const vf_profile_t profile_40 = vf_profile_default(40);
	const vf_profile_t profile_48 = vf_profile_default(48);
	const vf_profile_t profile_56 = vf_profile_default(56);
	const vf_profile_t profile_64 = vf_profile_default(64);
	
	const vf_base_t pos_value_32 = 0b01000010001010000000000000000000; // 42.0
	const vf_base_t pos_value_64 = 0b0100000001000101000000000000000000000000000000000000000000000000; // 42.0
	
	const vf_base_t neg_value_32 = 0b11000010001010000000000000000000; // -42.0
	const vf_base_t neg_value_64 = 0b1100000001000101000000000000000000000000000000000000000000000000; // -42.0

	const vf_base_t pos_nan_32 = vf_nan(profile_32);
	const vf_base_t neg_nan_32 = vf_negate(pos_nan_32, profile_32);
	const vf_base_t pos_inf_32 = vf_inf(profile_32);
	const vf_base_t neg_inf_32 = vf_negate(pos_inf_32, profile_32);
	const vf_base_t pos_subnormal_32 = vf_base_from_float(1.0 / FLT_MAX, profile_32, nullptr);
	const vf_base_t neg_subnormal_32 = vf_negate(pos_subnormal_32, profile_32);
	
	TEST_CASE("vf_profile_default(8)", "[vf_profile]") {
		SECTION("8 bit") {
			REQUIRE(profile_8.bias == 3);
			REQUIRE(profile_8.bits == 8);
			REQUIRE(profile_8.sgn_bits == 1);
			REQUIRE(profile_8.exp_bits == 3);
			REQUIRE(profile_8.mnt_bits == 4);
			REQUIRE(profile_8.sgn_mask == 0b10000000);
			REQUIRE(profile_8.exp_mask == 0b01110000);
			REQUIRE(profile_8.mnt_mask == 0b00001111);
		}
	
		SECTION("16 bit") {
			REQUIRE(profile_16.bias == 15);
			REQUIRE(profile_16.bits == 16);
			REQUIRE(profile_16.sgn_bits == 1);
			REQUIRE(profile_16.exp_bits == 5);
			REQUIRE(profile_16.mnt_bits == 10);
			REQUIRE(profile_16.sgn_mask == 0b1000000000000000);
			REQUIRE(profile_16.exp_mask == 0b0111110000000000);
			REQUIRE(profile_16.mnt_mask == 0b0000001111111111);
		}

		SECTION("32 bit") {
			REQUIRE(profile_32.bias == 127);
			REQUIRE(profile_32.bits == 32);
			REQUIRE(profile_32.sgn_bits == 1);
			REQUIRE(profile_32.exp_bits == 8);
			REQUIRE(profile_32.mnt_bits == 23);
			REQUIRE(profile_32.sgn_mask == 0b10000000000000000000000000000000);
			REQUIRE(profile_32.exp_mask == 0b01111111100000000000000000000000);
			REQUIRE(profile_32.mnt_mask == 0b00000000011111111111111111111111);
		}

		SECTION("64 bit") {
			REQUIRE(profile_64.bias == 1023);
			REQUIRE(profile_64.bits == 64);
			REQUIRE(profile_64.sgn_bits == 1);
			REQUIRE(profile_64.exp_bits == 11);
			REQUIRE(profile_64.mnt_bits == 52);
			REQUIRE(profile_64.sgn_mask == 0b1000000000000000000000000000000000000000000000000000000000000000);
			REQUIRE(profile_64.exp_mask == 0b0111111111110000000000000000000000000000000000000000000000000000);
			REQUIRE(profile_64.mnt_mask == 0b0000000000001111111111111111111111111111111111111111111111111111);
		}
	}
	
	TEST_CASE("vf_base_from_float", "[float]") {
		const float flt = 42.0f;
		vf_status_t status;
		uint32_t base = (uint32_t) vf_base_from_float(flt, profile_32, &status);
		
		REQUIRE(flt == *(float *) &base);
		REQUIRE(status == VF_STATUS_OK);
	}

	TEST_CASE("vf_base_from_double", "[double]") {
		const double dbl = 42.0;
		vf_status_t status;
		uint64_t base = (uint64_t) vf_base_from_double(dbl, profile_64, &status);
		
		REQUIRE(dbl == *(double *) &base);
		REQUIRE(status == VF_STATUS_OK);
	}

	TEST_CASE("vf_base_to_float", "[float]") {
		const vf_base_t base = pos_value_32; // 42.0
		vf_status_t status;
		float flt = vf_base_to_float(base, profile_32, &status);

		REQUIRE(flt == *(float *)&base);
		REQUIRE(status == VF_STATUS_OK);
	}

	TEST_CASE("vf_base_to_double", "[float]") {
		const vf_base_t base = pos_value_64; // 42.0
		vf_status_t status;
		double dbl = vf_base_to_double(base, profile_64, &status);
		
		REQUIRE(dbl == *(double *)&base);
		REQUIRE(status == VF_STATUS_OK);
	}

	TEST_CASE("vf_cast", "[cast]") {
		const uint64_t pos_overflow_64 = (uint64_t) vf_base_from_double(DBL_MAX, profile_64, nullptr);
		const uint64_t pos_underflow_64 = (uint64_t) vf_base_from_double(DBL_MIN, profile_64, nullptr);
		
		const uint64_t neg_overflow_64 = vf_negate(pos_overflow_64, profile_64);
		const uint64_t neg_underflow_64 = vf_negate(pos_underflow_64, profile_64);
		
		const vf_base_t pos_inf_64 = vf_inf(profile_64);
		const vf_base_t pos_inf_32 = vf_inf(profile_32);
		
		const vf_base_t neg_inf_64 = vf_negate(pos_inf_64, profile_64);
		const vf_base_t neg_inf_32 = vf_negate(pos_inf_32, profile_32);
		
		const vf_base_t pos_zero_32 = 0;
		const vf_base_t neg_zero_32 = vf_negate(pos_zero_32, profile_32);
		
		vf_status_t status;
		
		SECTION("Values that fit into target profile convert just fine") {
			REQUIRE(pos_value_32 == vf_cast(pos_value_64, profile_64, profile_32, &status));
			REQUIRE(status == VF_STATUS_OK);
			
			REQUIRE(pos_value_64 == vf_cast(pos_value_32, profile_32, profile_64, &status));
			REQUIRE(status == VF_STATUS_OK);
			
			REQUIRE(pos_value_64 == vf_cast(pos_value_64, profile_64, profile_64, &status));
			REQUIRE(status == VF_STATUS_OK);
			
			REQUIRE(pos_value_32 == vf_cast(pos_value_32, profile_32, profile_32, &status));
			REQUIRE(status == VF_STATUS_OK);
		}
		
		SECTION("Values that are too large for target profile overflow and become infinity") {
			REQUIRE(vf_cast(pos_overflow_64, profile_64, profile_32, &status) == pos_inf_32);
			REQUIRE(status == VF_STATUS_OVERFLOW);
			
			REQUIRE(vf_cast(vf_negate(pos_overflow_64, profile_64), profile_64, profile_32, &status) == neg_inf_32);
			REQUIRE(status == VF_STATUS_OVERFLOW);
		}
		
		SECTION("Values that are too small for target profile underflow and become zero") {
			REQUIRE(vf_cast(pos_underflow_64, profile_64, profile_32, &status) == pos_zero_32);
			REQUIRE(status == VF_STATUS_UNDERFLOW);
			
			REQUIRE(vf_cast(vf_negate(pos_underflow_64, profile_64), profile_64, profile_32, &status) == neg_zero_32);
			REQUIRE(status == VF_STATUS_UNDERFLOW);
		}
		
		SECTION("Infinity remains infinity") {
			REQUIRE(vf_cast(pos_inf_64, profile_64, profile_32, &status) == pos_inf_32);
			REQUIRE(status == VF_STATUS_OK);
			
			REQUIRE(vf_cast(neg_inf_64, profile_64, profile_32, &status) == neg_inf_32);
			REQUIRE(status == VF_STATUS_OK);
		
			REQUIRE(vf_cast(pos_inf_32, profile_32, profile_64, &status) == pos_inf_64);
			REQUIRE(status == VF_STATUS_OK);
			
			REQUIRE(vf_cast(neg_inf_32, profile_32, profile_64, &status) == neg_inf_64);
			REQUIRE(status == VF_STATUS_OK);
		}
		
		SECTION("Casting from float to double mirrors system behaviour") {
			size_t count = 1000000;
			srand(time(NULL));
			for (size_t i = 0; i < count; i++) {
				vf_base_t before_vf = rand();
				float before = vf_base_to_float(before_vf, profile_32, nullptr);
				vf_base_t after_vf = vf_cast(before_vf, profile_32, profile_64, nullptr);
				double after = (double)before;
				vf_base_t after_os = vf_base_from_double(after, profile_64, nullptr);
				if (after_vf != after_os) {
					REQUIRE(after_vf == after_os);
				}
			}
		}
		
		SECTION("Casting from double to float mirrors system behaviour") {
			size_t count = 1000000;
			srand(time(NULL));
			for (size_t i = 0; i < count; i++) {
				vf_base_t before_vf = ((vf_base_t)rand() << 32) | rand();
				double before = vf_base_to_double(before_vf, profile_64, nullptr);
				vf_base_t after_vf = vf_cast(before_vf, profile_64, profile_32, nullptr);
				float after = (float)before;
				vf_base_t after_os = vf_base_from_float(after, profile_32, nullptr);
				if (after_vf != after_os) {
					REQUIRE(after_vf == after_os);
				}
			}
		}
	}
	
	TEST_CASE("vf_rounding_error", "[vf_cast]") {
		double input = 1.234567890123456789;
		
		double error_8  = fabs(vf_rounding_error(input, profile_8));
		double error_16 = fabs(vf_rounding_error(input, profile_16));
		double error_24 = fabs(vf_rounding_error(input, profile_24));
		double error_32 = fabs(vf_rounding_error(input, profile_32));
		double error_40 = fabs(vf_rounding_error(input, profile_40));
		double error_48 = fabs(vf_rounding_error(input, profile_48));
		double error_56 = fabs(vf_rounding_error(input, profile_56));
		double error_64 = fabs(vf_rounding_error(input, profile_64));

		REQUIRE(vf_base_from_double(error_8,  profile_8,  nullptr) == 0b00000001);
		REQUIRE(vf_base_from_double(error_16, profile_16, nullptr) == 0b0000101001010010);
		REQUIRE(vf_base_from_double(error_24, profile_24, nullptr) == 0b001011010110111101011101);
		REQUIRE(vf_base_from_double(error_32, profile_32, nullptr) == 0b00110010001000110001011001111111);
		REQUIRE(vf_base_from_double(error_40, profile_40, nullptr) == 0b0010111101000101100111111011000000000000);
		REQUIRE(vf_base_from_double(error_48, profile_48, nullptr) == 0b001101011110011111101100000000000000000000000000);
		REQUIRE(vf_base_from_double(error_56, profile_56, nullptr) == 0b00111001101010000000000000000000000000000000000000000000);
		REQUIRE(vf_base_from_double(error_64, profile_56, nullptr) == 0b0000000000000000000000000000000000000000000000000000000000000000);
	}
	
	TEST_CASE("vf_optimal_profile", "[vf_cast]") {
		double input = 1.234567890123456789;
		const vf_profile_t profiles[] {
			profile_8,
			profile_16,
			profile_24,
			profile_32,
			profile_40,
			profile_48,
			profile_56
		};		
		REQUIRE(vf_optimal_profile(input, 0.1, (vf_profile_t *)profiles, 7).bits == 8);
		REQUIRE(vf_optimal_profile(input, 0.001, (vf_profile_t *)profiles, 7).bits == 16);
		REQUIRE(vf_optimal_profile(input, 0.00001, (vf_profile_t *)profiles, 7).bits == 24);
		REQUIRE(vf_optimal_profile(input, 0.00000001, (vf_profile_t *)profiles, 7).bits == 32);
		REQUIRE(vf_optimal_profile(input, 0.000000001, (vf_profile_t *)profiles, 7).bits == 40);
		REQUIRE(vf_optimal_profile(input, 0.00000000001, (vf_profile_t *)profiles, 7).bits == 48);
		REQUIRE(vf_optimal_profile(input, 0.00000000000001, (vf_profile_t *)profiles, 7).bits == 56);
	}
	
	TEST_CASE("vf_optimize", "[vf_cast]") {
		double input = 1.234567890123456789;
		const vf_profile_t profiles[] {
			profile_8,
			profile_16,
			profile_24,
			profile_32,
			profile_40,
			profile_48,
			profile_56
		};
		vf_profile_t profile_out;
		REQUIRE(vf_optimize(input, 0.1, (vf_profile_t *)profiles, 7, &profile_out) == 0b00110100);
		REQUIRE(vf_optimize(input, 0.001, (vf_profile_t *)profiles, 7, &profile_out) == 0b0011110011110000);
		REQUIRE(vf_optimize(input, 0.00001, (vf_profile_t *)profiles, 7, &profile_out) == 0b001111110011110000001101);
		REQUIRE(vf_optimize(input, 0.00000001, (vf_profile_t *)profiles, 7, &profile_out) == 0b00111111100111100000011001010010);
		REQUIRE(vf_optimize(input, 0.000000001, (vf_profile_t *)profiles, 7, &profile_out) == 0b0011111110011110000001100101001000010100);
		REQUIRE(vf_optimize(input, 0.00000000001, (vf_profile_t *)profiles, 7, &profile_out) == 0b001111111100111100000011001010010000101000110001);
		REQUIRE(vf_optimize(input, 0.00000000000001, (vf_profile_t *)profiles, 7, &profile_out) == 0b00111111111001111000000110010100100001010001100010110100);
	}
	
	TEST_CASE("vf_min", "[vf_util]") {
		REQUIRE(vf_min(profile_32) == vf_base_from_float(FLT_MIN, profile_32, nullptr));
		REQUIRE(vf_min(profile_64) == vf_base_from_double(DBL_MIN, profile_64, nullptr));
	}
	
	TEST_CASE("vf_max", "[vf_util]") {
		REQUIRE(vf_max(profile_32) == vf_base_from_float(FLT_MAX, profile_32, nullptr));
		REQUIRE(vf_max(profile_64) == vf_base_from_double(DBL_MAX, profile_64, nullptr));
	}
	
	TEST_CASE("vf_lowest", "[vf_util]") {
		REQUIRE(vf_lowest(profile_32) == vf_base_from_float(-FLT_MAX, profile_32, nullptr));
		REQUIRE(vf_lowest(profile_64) == vf_base_from_double(-DBL_MAX, profile_64, nullptr));
	}
	
//	TEST_CASE("vf_epsilon", "[vf_util]") {
//		REQUIRE(vf_epsilon(profile_32) == vf_base_from_float(FLT_EPSILON, profile_32, nullptr));
//		REQUIRE(vf_epsilon(profile_64) == vf_base_from_double(DBL_EPSILON, profile_64, nullptr));
//	}
	
//	TEST_CASE("vf_round_error", "[vf_util]") {
//		REQUIRE(vf_round_error(profile_32) == vf_base_from_float(0.5f, profile_32, nullptr));
//		REQUIRE(vf_round_error(profile_64) == vf_base_from_double(0.5, profile_64, nullptr));
//	}
	
	TEST_CASE("vf_qnan", "[vf_util]") {
		REQUIRE(vf_qnan(profile_32) == 0b01111111110000000000000000000000);
		REQUIRE(vf_qnan(profile_64) == 0b0111111111111000000000000000000000000000000000000000000000000000);
		REQUIRE(vf_qnan(profile_32) == vf_base_from_float((float)NAN, profile_32, nullptr));
		REQUIRE(vf_qnan(profile_64) == vf_base_from_double((double)NAN, profile_64, nullptr));
	}
	
	TEST_CASE("vf_snan", "[vf_util]") {
		REQUIRE(vf_snan(profile_32) == 0b01111111101000000000000000000000);
		REQUIRE(vf_snan(profile_64) == 0b0111111111110100000000000000000000000000000000000000000000000000);
	}
	
	TEST_CASE("vf_denorm_min", "[vf_util]") {
		REQUIRE(vf_denorm_min(profile_32) == 1);
		REQUIRE(vf_denorm_min(profile_64) == 1);
	}
	
	TEST_CASE("vf_digits", "[vf_util]") {
		REQUIRE(vf_digits(profile_32) == 24);
		REQUIRE(vf_digits(profile_64) == 53);
	}
	
	TEST_CASE("vf_digits10", "[vf_util]") {
		REQUIRE(vf_digits10(profile_32) == (int)(24 * 0.301));
		REQUIRE(vf_digits10(profile_64) == (int)(53 * 0.301));
	}
	
	TEST_CASE("vf_max_digits10", "[vf_util]") {
		REQUIRE(vf_max_digits10(profile_32) == (int)(24 * 0.301) / 2);
		REQUIRE(vf_max_digits10(profile_64) == (int)(53 * 0.301) / 2);
	}
	
	TEST_CASE("vf_min_exponent", "[vf_util]") {
		REQUIRE(vf_min_exponent(profile_32) == -profile_32.bias);
		REQUIRE(vf_min_exponent(profile_64) == -profile_64.bias);
	}
	
	TEST_CASE("vf_min_exponent10", "[vf_util]") {
		REQUIRE(vf_min_exponent10(profile_32) == (int)(vf_min_exponent(profile_32) * 0.301));
		REQUIRE(vf_min_exponent10(profile_64) == (int)(vf_min_exponent(profile_64) * 0.301));
	}
	
	TEST_CASE("vf_max_exponent", "[vf_util]") {
		REQUIRE(vf_max_exponent(profile_32) == profile_32.bias);
		REQUIRE(vf_max_exponent(profile_64) == profile_64.bias);
	}
	
	TEST_CASE("vf_max_exponent10", "[vf_util]") {
		REQUIRE(vf_max_exponent10(profile_32) == (int)(vf_max_exponent(profile_32) * 0.301));
		REQUIRE(vf_max_exponent10(profile_64) == (int)(vf_max_exponent(profile_64) * 0.301));
	}
	
	TEST_CASE("vf_nan", "[vf_util]") {
		REQUIRE(vf_nan(profile_32) == vf_base_from_float((float)NAN, profile_32, nullptr));
		REQUIRE(vf_nan(profile_64) == vf_base_from_double((double)NAN, profile_64, nullptr));
	}
	
	TEST_CASE("vf_inf", "[vf_util]") {
		REQUIRE(vf_inf(profile_32) == vf_base_from_float((float)INFINITY, profile_32, nullptr));
		REQUIRE(vf_inf(profile_64) == vf_base_from_double((double)INFINITY, profile_64, nullptr));
	}
	
	TEST_CASE("vf_abs", "[vf_util]") {
		REQUIRE(vf_abs(neg_value_32, profile_32) == pos_value_32);
		REQUIRE(vf_abs(neg_value_64, profile_64) == pos_value_64);
	}
	
	TEST_CASE("vf_compare", "[vf_util]") {
		REQUIRE(vf_compare(pos_value_32, neg_value_32, profile_32) == -1);
		REQUIRE(vf_compare(neg_value_32, pos_value_32, profile_32) == 1);
		REQUIRE(vf_compare(neg_value_32, neg_value_32, profile_32) == 0);
		REQUIRE(vf_compare(pos_value_32, pos_value_32, profile_32) == 0);
		REQUIRE(vf_compare(pos_inf_32, pos_value_32, profile_32) == -1);
		REQUIRE(vf_compare(pos_value_32, pos_inf_32, profile_32) == 1);
	}
	
	TEST_CASE("vf_isnormal", "[vf_util]") {
		REQUIRE(vf_isnormal(pos_value_32, profile_32) == 1);
		REQUIRE(vf_isnormal(neg_value_32, profile_32) == 1);
		REQUIRE(vf_isnormal(pos_nan_32, profile_32) == 0);
		REQUIRE(vf_isnormal(neg_nan_32, profile_32) == 0);
		REQUIRE(vf_isnormal(pos_inf_32, profile_32) == 0);
		REQUIRE(vf_isnormal(neg_inf_32, profile_32) == 0);
		REQUIRE(vf_isnormal(pos_subnormal_32, profile_32) == 0);
		REQUIRE(vf_isnormal(neg_subnormal_32, profile_32) == 0);
	}
	
	TEST_CASE("vf_isfinite", "[vf_util]") {
		REQUIRE(vf_isfinite(pos_value_32, profile_32) == 1);
		REQUIRE(vf_isfinite(neg_value_32, profile_32) == 1);
		REQUIRE(vf_isfinite(pos_nan_32, profile_32) == 0);
		REQUIRE(vf_isfinite(neg_nan_32, profile_32) == 0);
		REQUIRE(vf_isfinite(pos_inf_32, profile_32) == 0);
		REQUIRE(vf_isfinite(neg_inf_32, profile_32) == 0);
		REQUIRE(vf_isfinite(pos_subnormal_32, profile_32) == 1);
		REQUIRE(vf_isfinite(neg_subnormal_32, profile_32) == 1);
	}
	
	TEST_CASE("vf_isnan", "[vf_util]") {
		REQUIRE(vf_isnan(pos_value_32, profile_32) == 0);
		REQUIRE(vf_isnan(neg_value_32, profile_32) == 0);
		REQUIRE(vf_isnan(pos_nan_32, profile_32) == 1);
		REQUIRE(vf_isnan(neg_nan_32, profile_32) == 1);
		REQUIRE(vf_isnan(pos_inf_32, profile_32) == 0);
		REQUIRE(vf_isnan(neg_inf_32, profile_32) == 0);
		REQUIRE(vf_isnan(pos_subnormal_32, profile_32) == 0);
		REQUIRE(vf_isnan(neg_subnormal_32, profile_32) == 0);
	}
	
	TEST_CASE("vf_isinf", "[vf_util]") {
		REQUIRE(vf_isinf(pos_value_32, profile_32) == 0);
		REQUIRE(vf_isinf(neg_value_32, profile_32) == 0);
		REQUIRE(vf_isinf(pos_nan_32, profile_32) == 0);
		REQUIRE(vf_isinf(neg_nan_32, profile_32) == 0);
		REQUIRE(vf_isinf(pos_inf_32, profile_32) == 1);
		REQUIRE(vf_isinf(neg_inf_32, profile_32) == 1);
		REQUIRE(vf_isinf(pos_subnormal_32, profile_32) == 0);
		REQUIRE(vf_isinf(neg_subnormal_32, profile_32) == 0);
	}

}
