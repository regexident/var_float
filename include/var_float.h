#ifndef VAR_FLOAT_VAR_FLOAT_H
#define VAR_FLOAT_VAR_FLOAT_H

#include <stdint.h>
#include <float.h>

#ifdef __cplusplus
extern "C" {
#endif
	
	typedef uint64_t vf_base_t;
	typedef int64_t vf_signed_base_t;
	typedef double vf_float_base_t;
	
	typedef int vf_status_t;
	
	static const vf_status_t VF_STATUS_OK = 0;
	
	static const vf_status_t VF_STATUS_UNSUPPORTED = -1;
	static const vf_status_t VF_STATUS_INVALID_PROFILE = -2;
	
	static const vf_status_t VF_STATUS_OVERFLOW = -3;
	static const vf_status_t VF_STATUS_UNDERFLOW = -4;
	
	typedef struct {
		vf_base_t sgn_mask;
		vf_base_t exp_mask;
		vf_base_t mnt_mask;
		vf_base_t bias;
		vf_base_t mask;
		uint16_t bits;
		uint16_t sgn_bits;
		uint16_t exp_bits;
		uint16_t mnt_bits;
	} vf_profile_t;
	
	vf_profile_t vf_profile_default(uint16_t bits);
	vf_profile_t vf_profile_custom(uint16_t bits, uint16_t exp, vf_base_t bias);
	int vf_profile_equal(vf_profile_t lhs, vf_profile_t rhs);
	vf_status_t vf_validate(vf_profile_t profile);
	
	float vf_base_to_float(vf_base_t input, vf_profile_t profile, vf_status_t *status);
	double vf_base_to_double(vf_base_t input, vf_profile_t profile, vf_status_t *status);
	
	vf_base_t vf_base_from_float(float input, vf_profile_t profile, vf_status_t *status);
	vf_base_t vf_base_from_double(double input, vf_profile_t profile, vf_status_t *status);
	
	vf_base_t vf_inflate(vf_base_t input, vf_profile_t vf_in, vf_profile_t vf_out);
	vf_base_t vf_deflate(vf_base_t input, vf_profile_t vf_in, vf_profile_t vf_out, vf_status_t *status);
	vf_base_t vf_cast(vf_base_t value, vf_profile_t profile_in, vf_profile_t profile_out, vf_status_t *status);
	
	double vf_rounding_error(double input, vf_profile_t profile);
	vf_profile_t vf_optimal_profile(double input, double epsilon, vf_profile_t *sorted_candidate_profiles, size_t count);
	vf_base_t vf_optimize(double input, double epsilon, vf_profile_t *sorted_candidate_profiles, size_t count, vf_profile_t *profile_out);
	
	vf_base_t vf_min(vf_profile_t profile);
	vf_base_t vf_max(vf_profile_t profile);
	vf_base_t vf_lowest(vf_profile_t profile);
//	vf_base_t vf_epsilon(vf_profile_t profile);
//	vf_base_t vf_round_error(vf_profile_t profile);
	vf_base_t vf_qnan(vf_profile_t profile);
	vf_base_t vf_snan(vf_profile_t profile);
	vf_base_t vf_denorm_min(vf_profile_t profile);
	
	int vf_digits(vf_profile_t profile);
	int vf_digits10(vf_profile_t profile);
	int vf_max_digits10(vf_profile_t profile);
	
	int vf_min_exponent(vf_profile_t profile);
	int vf_min_exponent10(vf_profile_t profile);
	int vf_max_exponent(vf_profile_t profile);
	int vf_max_exponent10(vf_profile_t profile);
	
	vf_base_t vf_nan(vf_profile_t profile);
	vf_base_t vf_inf(vf_profile_t profile);
	
	vf_base_t vf_negate(vf_base_t value, vf_profile_t profile);
	vf_base_t vf_abs(vf_base_t value, vf_profile_t profile);
	
	int vf_compare(vf_base_t lhs, vf_base_t rhs, vf_profile_t profile);
	
	int vf_isnormal(vf_base_t value, vf_profile_t profile);
	int vf_isfinite(vf_base_t value, vf_profile_t profile);
	int vf_isnan(vf_base_t value, vf_profile_t profile);
	int vf_isinf(vf_base_t value, vf_profile_t profile);
	int vf_iszero(vf_base_t value, vf_profile_t profile);
	
#ifdef __cplusplus
}
#endif

#endif
