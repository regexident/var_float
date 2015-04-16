#include <stdint.h>
#include <limits.h>
#include <stddef.h>
#include <math.h>
#include <assert.h>

#include "var_float.h"

#define VF_ROUND_TIES_TO_EVEN 1

vf_profile_t vf_profile_default(const uint16_t bits) {
	assert(bits >= 8);
	static const uint16_t exponents[] = {3, 5, 7, 8, 8, 9, 10, 11};
	uint16_t exponent;
	if (bits <= 64) {
		exponent = exponents[(size_t)round((float)bits / 8) - 1];
	} else {
		// Formula for IEEE 754-1985 interchange formats:
		exponent = (uint16_t)(round(4 * (log(bits) / log(2))) - 13);
	}
	const vf_base_t bias = ((vf_base_t)1 << (exponent - 1)) - 1;
	return vf_profile_custom(bits, exponent, bias);
}

vf_profile_t vf_profile_custom(const uint16_t bits, const uint16_t exp, const vf_base_t bias) {
	assert(bits >= 8);
	assert(exp > 0);
	assert(exp < (bits - 2));
	assert(bias > 0);
	assert(bias < ((vf_base_t)1 << exp) - 1);
	return (vf_profile_t){
		.bias = bias,
		.bits = bits,
		.mask = (~(vf_base_t)0) >> ((sizeof(vf_base_t) * CHAR_BIT) - bits),
		.sgn_bits = 1,
		.exp_bits = exp,
		.mnt_bits = (uint16_t)(bits - 1 - exp),
		.sgn_mask = (vf_base_t)1 << (bits - 1),
		.exp_mask = (((vf_base_t)1 << (bits - 1)) - 1) & ~(((vf_base_t)1 << (bits - 1 - exp)) - 1),
		.mnt_mask = ((vf_base_t)1 << (bits - 1 - exp)) - 1
	};
}

int vf_profile_equal(vf_profile_t lhs, vf_profile_t rhs) {
	int bias_equal = lhs.bias == rhs.bias;
	int bits_equal = lhs.bits == rhs.bits;
	int sgn_bits_equal = lhs.sgn_bits == rhs.sgn_bits;
	int exp_bits_equal = lhs.exp_bits == rhs.exp_bits;
	int mnt_bits_equal = lhs.mnt_bits == rhs.mnt_bits;
	int sgn_mask_equal = lhs.sgn_mask == rhs.sgn_mask;
	int exp_mask_equal = lhs.exp_mask == rhs.exp_mask;
	int mnt_mask_equal = lhs.mnt_mask == rhs.mnt_mask;
	return bias_equal & bits_equal & sgn_bits_equal & exp_bits_equal & mnt_bits_equal & sgn_mask_equal & exp_mask_equal & mnt_mask_equal;
}

vf_status_t vf_validate(const vf_profile_t profile) {
	const int bases_valid = (sizeof(vf_base_t) * CHAR_BIT) == (sizeof(vf_float_base_t) * CHAR_BIT);
	const int bits_valid = (profile.bits > 0) && (profile.bits <= sizeof(vf_base_t) * CHAR_BIT);
	const int exponent_valid = (profile.exp_bits > 0) && (1 + profile.exp_bits + 1 < profile.bits);
	const int bias_valid = (profile.bias > 0) && (profile.bias < ((vf_base_t) 1 << profile.exp_bits) - 1);
	return (bases_valid && bits_valid && exponent_valid && bias_valid) ? VF_STATUS_OK : VF_STATUS_INVALID_PROFILE;
}

void vf_set_status(vf_status_t *status_ptr, vf_status_t status) {
	if (status_ptr != NULL) {
		*status_ptr = status;
	}
}

vf_base_t vf_deflate_with_overflow(vf_base_t input, vf_profile_t profile_in, vf_profile_t profile_out, vf_status_t *status) {
	const vf_base_t mnt_shift = profile_in.mnt_bits - profile_out.mnt_bits;
	const vf_base_t bits_delta = profile_in.bits - profile_out.bits;
	
	vf_base_t exp_in;
	vf_base_t mnt_in;
	vf_base_t sgn_out;
	
	sgn_out = (input & profile_in.sgn_mask) >> bits_delta;
	exp_in = input & profile_in.exp_mask;
	
	const vf_base_t overflow_threshold = (profile_in.sgn_mask >> 1) | ((((vf_base_t)1 << (profile_out.exp_bits - 1)) - 1) << profile_in.mnt_bits);
	assert(exp_in >= overflow_threshold);
	
	vf_set_status(status, VF_STATUS_OK);
	
	if (exp_in == profile_in.exp_mask) {
		// Inf or NaN
		mnt_in = input & profile_in.mnt_mask;
		if (mnt_in != 0) {
			vf_base_t exp_out = profile_out.exp_mask;
			vf_base_t mnt_out = mnt_in >> mnt_shift;
			// Add QNaN bit...
			mnt_out |= ((vf_base_t)1 << (profile_out.mnt_bits - 1));
			const vf_base_t output = (sgn_out + exp_out + mnt_out) & profile_out.mask;
			return output;
		}
	} else {
		// overflow to sgned inf
		vf_set_status(status, VF_STATUS_OVERFLOW);
	}
	vf_base_t output = (sgn_out + profile_out.exp_mask) & profile_out.mask;
	return output;
}

vf_base_t vf_deflate_with_underflow(vf_base_t input, vf_profile_t profile_in, vf_profile_t profile_out, vf_status_t *status) {
	vf_base_t exp_in;
	vf_base_t mnt_in;
	
	vf_base_t sgn_out;
	vf_base_t mnt_out;
	
	const vf_base_t bias_delta = profile_in.bias - profile_out.bias;
	const vf_base_t mnt_shift = profile_in.mnt_bits - profile_out.mnt_bits;
	const vf_base_t bits_delta = profile_in.bits - profile_out.bits;
	
	const vf_base_t rounding_mask = ((vf_base_t)1 << (profile_in.mnt_bits - profile_out.mnt_bits + 1)) - 1;
	
	sgn_out = (input & profile_in.sgn_mask) >> bits_delta;
	exp_in = input & profile_in.exp_mask;
	
	const vf_base_t underflow_threshold = (((vf_base_t)0x1 << (profile_in.exp_bits - profile_out.exp_bits)) - 1) << (profile_in.mnt_bits + profile_out.exp_bits - 1);
	assert(exp_in <= underflow_threshold);
	
	vf_set_status(status, VF_STATUS_OK);
	
	/*
	 * Signed zeros, subnormal floats, and floats with small
	 * exponents all convert to sgned zero halfs.
	 */
	vf_base_t threshold = ((profile_in.bias - profile_out.bias- profile_out.mnt_bits) << profile_in.mnt_bits);
	if (exp_in < threshold) {
		if ((input & (input & ~(profile_in.sgn_mask))) != 0) {
			vf_set_status(status, VF_STATUS_UNDERFLOW);
		}
		const vf_base_t output = (sgn_out) & profile_out.mask;
		return output;
	}
	// Make the subnormal sgnificand
	exp_in >>= profile_in.mnt_bits;
	mnt_in = ((vf_base_t)1 << (profile_in.mnt_bits)) + (input & profile_in.mnt_mask);
	// If it's not exactly represented, it underflowed
	if ((mnt_in & ((vf_base_t)1 << ((profile_in.bias - profile_out.bias) + (profile_in.mnt_bits - profile_out.mnt_bits) + 1 - exp_in)) - 1) != 0) {
		vf_set_status(status, VF_STATUS_UNDERFLOW);
	}
	mnt_in >>= (bias_delta + 1) - exp_in;
	// Handle rounding by adding 1 to the bit beyond half precision
#if VF_ROUND_TIES_TO_EVEN
	/*
	 * If the last bit in the half sgnificand is 0 (already even), and
	 * the remaining bit pattern is 1000...0, then we do not add one
	 * to the bit after the half sgnificand.  In all other cases, we do.
	 */
	const vf_base_t round_bit = profile_in.sgn_mask >> (1 + profile_in.exp_bits + profile_out.mnt_bits);
	if ((mnt_in & rounding_mask) != round_bit) {
		mnt_in += round_bit;
	}
#else
	mnt_in += round_bit;
#endif
	mnt_out = mnt_in >> mnt_shift;
	/*
	 * If the rounding causes a bit to spill into exp_out, it will
	 * increment exp_out from zero to one and mnt_out will be zero.
	 * This is the correct result.
	 */
	const vf_base_t output = (sgn_out + mnt_out) & profile_out.mask;
	return output;
}

vf_base_t vf_deflate_without_over_or_underflow(vf_base_t input, vf_profile_t profile_in, vf_profile_t profile_out, vf_status_t *status) {
	// casting between different profiles of same bit length not supported yet. Patches welcome.
	assert(profile_in.bits > profile_out.bits);
	assert(profile_in.exp_bits >= profile_out.exp_bits);
	assert(profile_in.mnt_bits >= profile_out.mnt_bits);
	
	vf_base_t exp_in;
	vf_base_t mnt_in;
	
	vf_base_t sgn_out;
	vf_base_t exp_out;
	vf_base_t mnt_out;
	
	const vf_base_t mnt_shift = profile_in.mnt_bits - profile_out.mnt_bits;
	const vf_base_t bits_delta = profile_in.bits - profile_out.bits;
	const vf_base_t underflow_threshold = (((vf_base_t)0x1 << (profile_in.exp_bits - profile_out.exp_bits)) - 1) << (profile_in.mnt_bits + profile_out.exp_bits - 1);
	
	vf_set_status(status, VF_STATUS_OK);
	
	sgn_out = (input & profile_in.sgn_mask) >> bits_delta;
	exp_in = input & profile_in.exp_mask;
	
	exp_out = ((exp_in - underflow_threshold) >> mnt_shift) & profile_out.exp_mask;
	mnt_in = input & profile_in.mnt_mask;
	// Handle rounding by adding 1 to the bit beyond half precision
	const vf_base_t round_bit = profile_in.sgn_mask >> (1 + profile_in.exp_bits + profile_out.mnt_bits);
#if VF_ROUND_TIES_TO_EVEN
	/*
	 * If the last bit in the half sgnificand is 0 (already even), and
	 * the remaining bit pattern is 1000...0, then we do not add one
	 * to the bit after the half sgnificand.  In all other cases, we do.
	 */
	const vf_base_t rounding_mask = ((vf_base_t)1 << (profile_in.mnt_bits - profile_out.mnt_bits + 1)) - 1;
	if ((mnt_in & rounding_mask) != round_bit) {
		mnt_in += round_bit;
	}
#else
	mnt_in += round_bit;
#endif
	mnt_out = mnt_in >> mnt_shift;
	/*
	 * If the rounding causes a bit to spill into exp_out, it will
	 * increment exp_out by one and mnt_out will be zero.  This is the
	 * correct result.  exp_out may increment to bias, at greatest, in
	 * which case the result overflows to a sgned inf.
	 */
	if (((mnt_out + exp_out) & profile_out.mask) == profile_out.exp_mask) {
		vf_set_status(status, VF_STATUS_OVERFLOW);
	}
	const vf_base_t output = (sgn_out + exp_out + mnt_out) & profile_out.mask;
	return output;
}

vf_base_t vf_deflate(vf_base_t input, vf_profile_t profile_in, vf_profile_t profile_out, vf_status_t *status) {
	// casting between different profiles of same bit length not supported yet. Patches welcome.
	assert(profile_in.bits > profile_out.bits);
	assert(profile_in.exp_bits >= profile_out.exp_bits);
	assert(profile_in.mnt_bits >= profile_out.mnt_bits);
	
	vf_base_t exp_in = input & profile_in.exp_mask;
	
	// Exponent overflow/NaN converts to sgned inf/NaN
	const vf_base_t overflow_threshold = (profile_in.sgn_mask >> 1) | ((((vf_base_t)1 << (profile_out.exp_bits - 1)) - 1) << profile_in.mnt_bits);
	if (exp_in >= overflow_threshold) {
		const vf_base_t output = vf_deflate_with_overflow(input, profile_in, profile_out, status);
		return output;
	}
	
	// Exponent underflow converts to a subnormal half or sgned zero
	const vf_base_t underflow_threshold = (((vf_base_t)0x1 << (profile_in.exp_bits - profile_out.exp_bits)) - 1) << (profile_in.mnt_bits + profile_out.exp_bits - 1);
	if (exp_in <= underflow_threshold) {
		const vf_base_t output = vf_deflate_with_underflow(input, profile_in, profile_out, status);
		return output;
	}
	
	// Regular case with no overflow or underflow
	const vf_base_t output = vf_deflate_without_over_or_underflow(input, profile_in, profile_out, status);
	return output;
}

vf_base_t vf_inflate(vf_base_t input, vf_profile_t profile_in, vf_profile_t profile_out) {
	// casting between different profiles of same bit length not supported yet. Patches welcome.
	assert(profile_in.bits < profile_out.bits);
	assert(profile_in.exp_bits <= profile_out.exp_bits);
	assert(profile_in.mnt_bits <= profile_out.mnt_bits);
	
	vf_base_t exp_in;
	vf_base_t mnt_in;
	
	vf_base_t sgn_out;
	vf_base_t exp_out;
	vf_base_t mnt_out;
	
	const vf_base_t bias_delta = profile_out.bias - profile_in.bias;
	const size_t mnt_shift = profile_out.mnt_bits - profile_in.mnt_bits;
	const size_t bits_delta = profile_out.bits - profile_in.bits;
	
	exp_in = input & profile_in.exp_mask;
	sgn_out = (input & profile_in.sgn_mask) << bits_delta;
	if (exp_in == 0) { // 0 or subnormal
		mnt_in = (input & profile_in.mnt_mask);
		// Signed zero
		if (mnt_in == 0) {
			const vf_base_t output = (sgn_out) & profile_out.mask;
			return output;
		}
		// Subnormal
		mnt_in <<= 1;
		const vf_base_t mask = ((vf_base_t)1 << profile_in.mnt_bits);
		while ((mnt_in & mask) == 0) {
			mnt_in <<= 1;
			exp_in++;
		}
		exp_out = (((bias_delta - exp_in)) << profile_out.mnt_bits) & profile_out.exp_mask;
		mnt_out = (mnt_in & profile_in.mnt_mask) << mnt_shift;
		const vf_base_t output = (sgn_out + exp_out + mnt_out) & profile_out.mask;
		return output;
	} else if (exp_in == profile_in.exp_mask) { // inf or NaN
		// All-ones exponent and a copy of the sgnificand
		exp_out = profile_out.exp_mask;
		mnt_out = (input & profile_in.mnt_mask) << mnt_shift;
		// QNaN
		if (mnt_out != 0) {
			mnt_out |= ((vf_base_t)1 << (profile_out.mnt_bits - 1));
		}
		const vf_base_t output = (sgn_out + exp_out + mnt_out) & profile_out.mask;
		return output;
	} else { // normalized
		// Just need to adjust the exponent and shift
		exp_out = ((input & profile_in.exp_mask) + (bias_delta << profile_in.mnt_bits)) << mnt_shift;
		mnt_out = (input & profile_in.mnt_mask) << mnt_shift;
		const vf_base_t output = (sgn_out + exp_out + mnt_out) & profile_out.mask;
		return output;
	}
}

vf_base_t vf_translate(vf_base_t input, vf_profile_t profile_in, vf_profile_t profile_out, vf_status_t *status) {
	// casting between different profiles of same bit length not supported yet. Patches welcome.
	assert(0);
	vf_set_status(status, VF_STATUS_UNSUPPORTED);
	return input;
}

vf_base_t vf_pass_thru(vf_base_t input, vf_profile_t profile_in, vf_profile_t profile_out, vf_status_t *status) {
	vf_set_status(status, VF_STATUS_OK);
	return input;
}

vf_base_t vf_cast(vf_base_t input, vf_profile_t profile_in, vf_profile_t profile_out, vf_status_t *status) {
	vf_base_t output;
	if (profile_in.bits < profile_out.bits) {
		vf_set_status(status, VF_STATUS_OK);
		output = vf_inflate(input, profile_in, profile_out);
	} else if (profile_in.bits > profile_out.bits) {
		output = vf_deflate(input, profile_in, profile_out, status);
	} else if (vf_profile_equal(profile_in, profile_out)) {
		output = vf_pass_thru(input, profile_in, profile_out, status);
	} else {
		output = vf_translate(input, profile_in, profile_out, status);
	}
	return output;
}

vf_base_t vf_base_from_float(float input, vf_profile_t profile, vf_status_t *status) {
	const uint32_t result = vf_cast(*(uint32_t*)&input, vf_profile_default(32), profile, status);
	return result;
}

vf_base_t vf_base_from_double(double input, vf_profile_t profile, vf_status_t *status) {
	const uint64_t result = vf_cast(*(uint64_t*)&input, vf_profile_default(64), profile, status);
	return result;
}

float vf_base_to_float(vf_base_t input, vf_profile_t profile, vf_status_t *status) {
	const uint32_t result = (uint32_t)vf_cast(input, profile, vf_profile_default(32), status);
	return *(float *)&result;
}

double vf_base_to_double(vf_base_t input, vf_profile_t profile, vf_status_t *status) {
	const uint64_t result = (uint64_t)vf_cast(input, profile, vf_profile_default(64), status);
	return *(double *)&result;
}

double vf_rounding_error(double input, vf_profile_t profile) {
	// Modified Veltkamp/Dekker algorithm
	static const int mnt_bits = 52;
	double x = input * (1 + pow(2.0, mnt_bits - profile.mnt_bits));
	double error = x - input - x + input;
	return error;
}

vf_profile_t vf_optimal_profile(double input, double epsilon, vf_profile_t *sorted_candidate_profiles, size_t count) {
	assert(sorted_candidate_profiles != NULL);
	assert(count > 0);
	
	vf_profile_t output_profile = vf_profile_default(64);
	size_t min = 0;
	size_t max = count - 1;
	while (max >= min) {
		size_t mid = min + ((max - min) / 2);
		vf_profile_t *target_profile = sorted_candidate_profiles + mid;
		if (output_profile.bits == target_profile->bits) {
			break;
		}
		double error = vf_rounding_error(input, *target_profile);
		if (fabs(error) <= fabs(epsilon)) {
			output_profile = *target_profile;
			if (mid == min) {
				break;
			}
			max = mid - 1;
		} else {
			if (mid == max) {
				break;
			}
			min = mid + 1;
		}
	}
	return output_profile;
}

vf_base_t vf_optimize(double input, double epsilon, vf_profile_t *sorted_candidate_profiles, size_t count, vf_profile_t *profile_out) {
	*profile_out = vf_optimal_profile(input, epsilon, sorted_candidate_profiles, count);
	vf_profile_t profile_in = vf_profile_default(64);
	vf_base_t base_input = vf_base_from_double(input, profile_in, NULL);
	vf_base_t output = vf_cast(base_input, profile_in, *profile_out, NULL);
	return output;
}

vf_base_t vf_min(vf_profile_t profile) {
	// pattern: 0 00001 0000000000
	return (vf_base_t)1 << profile.mnt_bits;
}

vf_base_t vf_max(vf_profile_t profile) {
	// pattern: 0 11110 1111111111
	return (((vf_base_t)1 << (profile.bits - 1)) - 1) & ~(vf_min(profile));
}

vf_base_t vf_lowest(vf_profile_t profile) {
	// pattern: 1 11110 1111111111
	return profile.sgn_mask | vf_max(profile);
}

//vf_base_t vf_epsilon(vf_profile_t profile) {
//	assert(0); // FIXME!!
//	return 0;
//}

//vf_base_t vf_round_error(vf_profile_t profile) {
//	assert(0); // FIXME!!
//	return 0;
//}

vf_base_t vf_qnan(vf_profile_t profile) {
	// pattern: 0 11111 1000000000
	return profile.exp_mask | (vf_base_t)1 << (profile.mnt_bits - 1);
}

vf_base_t vf_snan(vf_profile_t profile) {
	// pattern: 0 11111 0100000000
	return profile.exp_mask | (vf_base_t)1 << (profile.mnt_bits - 2);
}

vf_base_t vf_denorm_min(vf_profile_t profile) {
	// pattern: 0 00000 0000000001
	return (vf_base_t)1;
}

int vf_digits(vf_profile_t profile) {
	return profile.bits - profile.exp_bits;
}

int vf_digits10(vf_profile_t profile) {
	return (int)(0.301 * vf_digits(profile));
}

int vf_max_digits10(vf_profile_t profile) {
	return vf_digits10(profile) / 2;
}

int vf_min_exponent(vf_profile_t profile) {
	return -(int)profile.bias;
}

int vf_min_exponent10(vf_profile_t profile) {
	return (int)(0.301 * vf_min_exponent(profile));
}

int vf_max_exponent(vf_profile_t profile) {
	return (int)profile.bias;
}

int vf_max_exponent10(vf_profile_t profile) {
	return (int)(0.301 * vf_max_exponent(profile));
}

vf_base_t vf_nan(vf_profile_t profile) {
	// pattern: 0 11111 1000000000
	return vf_qnan(profile);
}

vf_base_t vf_inf(vf_profile_t profile) {
	// pattern: 0 11111 0000000000
	return profile.exp_mask;
}

vf_base_t vf_negate(vf_base_t value, vf_profile_t profile) {
	return value ^ profile.sgn_mask;
}

vf_base_t vf_abs(vf_base_t value, vf_profile_t profile) {
	return value & ~(profile.sgn_mask);
}

int vf_compare(vf_base_t lhs, vf_base_t rhs, vf_profile_t profile) {
	size_t shift = (sizeof(vf_base_t) * CHAR_BIT) - profile.bits;
	vf_signed_base_t sgned_lhs = lhs << shift;
	vf_signed_base_t sgned_rhs = rhs << shift;
	return (sgned_lhs < sgned_rhs) ? 1 : (sgned_lhs > sgned_rhs) ? -1 : 0;
}

int vf_isnormal(vf_base_t value, vf_profile_t profile) {
	return ((value & profile.exp_mask) != profile.exp_mask) && ((value & profile.exp_mask) != 0);
}

int vf_isfinite(vf_base_t value, vf_profile_t profile) {
	return ((value & profile.exp_mask) != profile.exp_mask);
}

int vf_isnan(vf_base_t value, vf_profile_t profile) {
	return ((value & profile.exp_mask) == profile.exp_mask) && ((value & profile.mnt_mask) != 0);
}

int vf_isinf(vf_base_t value, vf_profile_t profile) {
	return ((value & profile.exp_mask) == profile.exp_mask) && ((value & profile.mnt_mask) == 0);
}

int vf_iszero(vf_base_t value, vf_profile_t profile) {
	return vf_abs(value, profile) == 0;
}
