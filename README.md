# var\_float

**var\_float** implements a universal representation of IEEE 754 floating-point numbers convertible between arbitrary precisions.

## Features

* **Out-of-the-box** representations for floating point formats of sizes: **`8bit`**, **`16bit`** (half-float), **`24bit`**, **`32bit`** (float), **`40bit`**, **`48bit`**, **`56bit`** as well as **`64bit`** (double).
* **Custom formats** (arbitrary type sizes, exponent sizes & biases).
* **IEEE 754 compliant conversions** between formats of **arbitrary precision**.
* Smart detection of **most compact lossless representation** for any given real number.
* Smart detection of **most compact lossy representation** for any given real number within a specified **rounding error threshold**.

## Installation

Just copy the files in `"include/..."` and `"src/..."` into your project.

## Example usage

Packing `64bit` floating-point value into `16bit` floating-point representation:

	const double input = 42.0;
	const vf_profile_t profile = vf_profile_default(16);
	vf_status_t status; // used for checking for over/underflows
	const uint16_t output = (uint16_t) vf_base_from_double(input, profile, &status);
	// output == 20800 == "0 10100 0101000000" == 42.0 in 16bit

Unpacking `16bit` floating-point value to `64bit` floating-point representation:

	// 20800 == "0 10100 0101000000" == 42.0 in 16bit
	const uint16_t input = 20800;
	const vf_profile_t profile = vf_profile_default(16);
	vf_status_t status; // used for checking for over/underflows
	const uint16_t output = (uint16_t) vf_base_from_double(input, profile, &status);

Packing `64bit` floating-point value into most compact floating-point representation:

	const vf_profile_t profiles[] {
	    profile_8,
	    profile_16,
	    profile_32,
	};
	const size_t profiles_count = 3;
	const double input = 42.0;
	const double epsilon = 0.001;
	vf_profile_t profile; // optimal compact profile
	const uint16_t output = vf_optimize(input, epsilon, (vf_profile_t *)profiles, profiles_count, &profile);

There actually **is a lot more** to **var\_float** though. (There are a total of about ~40 functions for dealing with variable sized floats actually.) So feel encouraged to check them out in `var_float.h` yourself! :)

What **var\_float** does **not provide** however are implementations of floating-point arithmetic operations. Patches (if properly tested) are very much welcome though. ;)

## Dependencies

[Catch][1] for running unit tests. (optional)

## License

**var\_float** is available under the **BSD-3 clause license**.

[1]:	https://github.com/philsquared/Catch