#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <map>


struct point_t {
	point_t() = default;
	point_t(double x, double y) : x(x), y(y) {}

	double x, y;
};

struct line_t {
	line_t() = default;
	line_t(point_t p, double k) : p{p}, k{k} {}
	line_t(double x, double y, double k) : p{x,y}, k{k} {}

	point_t p;
	double k;
};


using ray_list = std::vector<line_t>;


/*!
 * \brief The base_curve represents piece-wise linear curve f,
 *  which is compliant to the requirements of network calculus.
 *
 * Network calculus implies the curve f has the following properties:
 *  - left-continuity [ lim {x -> -x1} f(x) = f(x1) ];
 *  - non-decreasing [ x' >= x" => f(x') >= f(x") ];
 *  - causality [ f(0) == 0 ].
 *
 * Curve constists of a non-empty ray sequence.
 * Each ray defines a half-line, which:
 *  - starts at the reference point (x,y), and
 *  - moves away from the origin with a slope k.
 *
 * A valid ray shall comply to the following restriction:
 *	0 <= x < inf && 0 <= y < inf && 0 <= k <= inf.
 *  (=> curves are non-decreasing)
 *
 * First ray starts at the origin (=> curves are causal).
 * For each subsequent ray its reference point belongs to half-line,
 *  defined by the previos ray (=> curves are continuous).
 *
 * However, the ray sequence must never include any redundancy:
 *  subsequent rays always have different reference points and slopes.
 */
class base_curve {
friend class ray_handle;
protected:
	ray_list rays;

public:

	/// using SFINAE to choose a proper handler for initializer list
	template<typename InputIt, typename ValueType> using iterates_over =
	std::enable_if<std::is_same<typename std::iterator_traits<InputIt>::value_type, ValueType>::value>;

	template<typename InputIt, typename iterates_over<InputIt, line_t>::type* = nullptr>
	static void validate_ray_list(InputIt first, InputIt last);

	template<typename InputIt, typename iterates_over<InputIt, line_t>::type* = nullptr>
	base_curve(InputIt first, InputIt last) : rays(first, last) {
		validate_ray_list(rays.cbegin(), rays.cend());
	}

	base_curve(const std::initializer_list<line_t> &ray_list) :
		base_curve(ray_list.begin(), ray_list.end()) {}

	template<typename OutputIt>
	OutputIt copy_ray_list(OutputIt it) const;

	bool is_convex() const {
		auto cmp = [](auto &a, auto &b) {return a.k < b.k;};
		return std::is_sorted(rays.cbegin(), rays.cend(), cmp);
	}

	bool is_concave() const {
		auto cmp = [](auto &a, auto &b) {return a.k > b.k;};
		return std::is_sorted(rays.cbegin(), rays.cend(), cmp);
	}

	virtual const base_curve&
	operator = (const base_curve &curve) {
		rays = curve.rays;
		return *this;
	}
};


struct rate_latency {
	rate_latency() = default;
	rate_latency(double rate, double latency) :
		rate(rate), latency(latency) {}

	double rate, latency;
};


/*!
 * \brief The convex_curve is a base_curve formed by
 *  a sequence of rays ordered by increasing slopes.
 *
 * Convex curves are build from either:
 *  - list of rays (which is similar to base curve), or
 *  - list of rate-latency pairs (which is a unique feature of convex_curve).
 *
 * Elements (R,L) of a valid rate-latency pair
 *  shall belong to the following domain:
 *  0 < R <= inf && 0 <= L < inf
 *
 * Pair (R,L) corresponds to a bi-segment curve:
 *  - the first one starts at origin and coinsides with 0X axis;
 *  - the second one starts at (L,0) and has a slope equal to R.
 *
 * Convex curve specified by a set of rate-latency pairs
 *  corresponds to a point-wise maximum among their curves.
 */
class convex_curve : public base_curve {
public:

	template<typename InputIt>
	static void validate_ray_list(InputIt first, InputIt last);

	convex_curve() : base_curve{{0.0, 0.0, 0.0}} {}
	convex_curve(const base_curve& curve) : base_curve(curve) {
		if (!curve.is_convex()) throw std::runtime_error("curve is not convex");
	}

	template<typename InputIt, typename iterates_over<InputIt, rate_latency>::type* = nullptr>
	convex_curve(InputIt first, InputIt last);

	convex_curve(const std::initializer_list<rate_latency> &rate_latency_list) :
		convex_curve(rate_latency_list.begin(), rate_latency_list.end()) {}

	template<typename InputIt, typename iterates_over<InputIt, line_t>::type* = nullptr>
	convex_curve(InputIt first, InputIt last) : base_curve(first, last) {
		if (!is_convex()) throw std::runtime_error("curve is not convex");
	}

	convex_curve(const std::initializer_list<line_t> &ray_list) :
		convex_curve(ray_list.begin(), ray_list.end()) {}

	template<typename OutputIt>
	OutputIt copy_rate_latency_list(OutputIt it) const;

	const convex_curve& operator = (const convex_curve& curve) {
		base_curve::operator=(curve);
		return *this;
	}

	virtual const base_curve& operator = (const base_curve& curve) override {
		if (!curve.is_convex()) throw std::runtime_error("curve is not convex");
		return base_curve::operator=(curve);
	}
};


struct rate_burst {
	rate_burst() = default;
	rate_burst(double rate, double burst) :
		rate(rate), burst(burst) {}

	double rate, burst;
};


/*!
 * \brief The concave_curve is a base_curve formed by
 *  a sequence of rays ordered by decreasing slopes.
 *
 * Concave curves are build from either:
 *  - list of rays (which is similar to base curve), or
 *  - list of rate-burst pairs (which is a unique feature of concave_curve).
 *
 * Elements (R,B) of a valid rate-burst pair
 *  shall belong to the following domain:
 *  0 <= R < inf && 0 <= B < inf
 *
 * Pair (R,B) corresponds to a bi-segment curve:
 *  - the first one starts at origin and coinsides with 0Y axis;
 *  - the second one starts at (0,B) and has a slope equal to R.
 *
 * Concave curve specified by a set of rate-burst pairs
 *  corresponds to a point-wise minimum among their curves.
 */
class concave_curve : public base_curve {
public:

	template<typename InputIt>
	static void validate_ray_list(InputIt first, InputIt last);

	concave_curve() : base_curve{{0.0, 0.0, 1.0/0.0}} {}
	concave_curve(const base_curve& curve) : base_curve(curve) {
		if (!curve.is_concave()) throw std::runtime_error("curve is not concave");
	}

	template<typename InputIt, typename iterates_over<InputIt, rate_burst>::type* = nullptr>
	concave_curve(InputIt first, InputIt last);

	concave_curve(const std::initializer_list<rate_burst> &rate_burst_list) :
		concave_curve(rate_burst_list.begin(), rate_burst_list.end()) {}

	template<typename InputIt, typename iterates_over<InputIt, line_t>::type* = nullptr>
	concave_curve(InputIt first, InputIt last) : base_curve(first, last) {
		if (!is_concave()) throw std::runtime_error("curve is not concave");
	}

	concave_curve(const std::initializer_list<line_t> &ray_list) :
		concave_curve(ray_list.begin(), ray_list.end()) {}

	template<typename OutputIt>
	OutputIt copy_rate_burst_list(OutputIt it) const;

	const concave_curve& operator = (const concave_curve& curve) {
		base_curve::operator=(curve);
		return *this;
	}

	virtual const base_curve& operator = (const base_curve& curve) override {
		if (!curve.is_concave()) throw std::runtime_error("curve is not concave");
		return base_curve::operator=(curve);
	}
};


//! Pointwise maximum for a pair of convex curves.
convex_curve max(const convex_curve &conv1, const convex_curve &conv2);

//! Pointwise minimum for a pair of concave curves.
concave_curve min(const concave_curve &conc1, const concave_curve &conc2);


//! Pointwise sum for a pair of convex curves.
convex_curve sum(const convex_curve &conv1, const convex_curve &conv2);

//! Pointwise sum for a pair of concave curves.
concave_curve sum(const concave_curve &conc1, const concave_curve &conc2);


//! Pointwise difference between a convex curve and a concave curve,
//!  with all the negative values are replaced by zero.
convex_curve leftover(const convex_curve &conv, const concave_curve &conc);

//! Min-Plus convolution for a pair of convex curves.
convex_curve convolve(const convex_curve &conv1, const convex_curve &conv2);

//! Min-Plus deconvolution for a concave curve and a convex curve.
concave_curve deconvolve(const concave_curve &conc, const convex_curve &conv);


//! Maximal horizontal difference between a concave and a convex curve.
double max_x_dist(const concave_curve &conc, const convex_curve &conv);

//! Maximal vertical difference between a concave and a convex curve.
double max_y_dist(const concave_curve &conc, const convex_curve &conv);


bool operator == (const base_curve &curve1, const base_curve &curve2);

std::ostream& operator << (std::ostream& stream, const base_curve &curve);
std::ostream& operator << (std::ostream& stream, const convex_curve &curve);
std::ostream& operator << (std::ostream& stream, const concave_curve &curve);


class ray_handle {
	static const ray_list& rays(const base_curve& curve) {return curve.rays;}
	static ray_list& rays(base_curve& curve) {return curve.rays;}

	friend bool operator == (const base_curve &lhs, const base_curve &rhs);

	friend base_curve max(const base_curve &curve1, const base_curve &curve2);
	friend base_curve min(const base_curve &curve1, const base_curve &curve2);
	friend base_curve sum(const base_curve &curve1, const base_curve &curve2);

	friend convex_curve max(const convex_curve &conv1, const convex_curve &conv2);
	friend convex_curve sum(const convex_curve &conv1, const convex_curve &conv2);

	friend concave_curve min(const concave_curve &conc1, const concave_curve &conc2);
	friend concave_curve sum(const concave_curve &conc1, const concave_curve &conc2);

	friend convex_curve leftover(const convex_curve &conv, const concave_curve &conc);
	friend convex_curve convolve(const convex_curve &conv1, const convex_curve &conv2);
	friend concave_curve deconvolve(const concave_curve &conc, const convex_curve &conv);

	friend double max_x_dist(const concave_curve &conc, const convex_curve &conv);
	friend double max_y_dist(const concave_curve &conc, const convex_curve &conv);
};


namespace detail {
	void validate_ray(const line_t &ray);
	void validate_rate_latency(const rate_latency &rl);
	void validate_rate_burst(const rate_burst &rb);
}


template<typename InputIt, typename base_curve::iterates_over<InputIt, line_t>::type*>
void base_curve::validate_ray_list(InputIt first, InputIt last) {
	// check validity of rays in the list
	for (auto iter = first; iter != last; ++iter) {
		detail::validate_ray(*iter);
	}

	// check restrictions on size and front side
	if (first == last) throw std::runtime_error("ray list is empty");
	if (first->p.x != 0.0 || first->p.y != 0.0) throw std::runtime_error("curve is not causal");

	// check restrictions to pairs of adjustent rays
	for(auto prev = first++; first != last; prev = first++) {
		double x_diff = first->p.x - prev->p.x, y_diff = first->p.y - prev->p.y;
		if ((x_diff == 0.0 && y_diff == 0.0) || prev->k == first->k)
			throw std::runtime_error("curve contains redundant rays");
		if (std::isinf(prev->k) && x_diff == 0.0 && y_diff > 0.0) continue;
		if (std::isfinite(prev->k) && y_diff == x_diff * prev->k) continue;
		throw std::runtime_error("curve is not continious");
	}
}

template<typename OutputIt>
OutputIt base_curve::copy_ray_list(OutputIt it) const {
	for(auto &ray : rays) *it++ = ray;
	return it;
}


template<typename InputIt>
void convex_curve::validate_ray_list(InputIt first, InputIt last) {
	// check validity of rays in the list
	for (auto iter = first; iter != last; ++iter) {
		detail::validate_ray(*iter);
	}

	// check restrictions on size and front side
	if (first == last) throw std::runtime_error("ray list is empty");
	if (first->p.x != 0.0 || first->p.y != 0.0) throw std::runtime_error("curve is not causal");

	// check restrictions to pairs of adjustent rays
	for(auto prev = first++; first != last; prev = first++) {
		if (prev->k >= first->k) throw std::runtime_error("curve is not convex");
		double x_diff = first->p.x - prev->p.x, y_diff = first->p.y - prev->p.y;
		if (x_diff <= 0.0) throw std::runtime_error("curve contains redundant rays");
		// x_diff * prev->k produces a valid reuslt, since 0 < x_diff < inf && 0 <= prev->k < inf
		if (y_diff != x_diff * prev->k) throw std::runtime_error("curve is not continious");
	}
}

template<typename InputIt, typename base_curve::iterates_over<InputIt, rate_latency>::type*>
convex_curve::convex_curve(InputIt first, InputIt last) : convex_curve() {
	if (first == last) return;

	// order given rays by increasing bursts
	std::map<double, double> latency_to_rate;
	for(; first != last; ++first) {
		detail::validate_rate_latency(*first);
		auto res = latency_to_rate.emplace(first->latency, first->rate);
		if (!res.second && res.first->second < first->rate)
			res.first->second = first->rate;
	}

	// handle the front side of the curve
	auto it = latency_to_rate.begin(), it_end = latency_to_rate.end();
	if (it->first > 0.0) rays.emplace_back(it->first, 0.0, it->second);
	else rays.back() = {it->first, 0.0, it->second};

	// filter out non-meaningful rays
	// calculate the Xs for the inflexion points
	while (++it != it_end) {
		if (rays.back().k >= it->second) continue;

		// remove previosly added rays,
		// if they are covered by the current one
		for (;; rays.pop_back()) {
			double x_diff = it->first - rays.back().p.x;
			double k_diff = 1.0/rays.back().k - 1.0/it->second;
			double y_next = x_diff / k_diff;

			if (y_next > rays.back().p.y) {
				rays.emplace_back(it->first, y_next, it->second);
				break;
			}
		}
	}

	// calculate proper X for each inflexion point
	auto rit = rays.begin(), rend = rays.end();
	while (++rit != rend && rit->p.y == 0.0);
	for(; rit != rend; ++rit) rit->p.x += rit->p.y / rit->k;
}

template<typename OutputIt>
OutputIt convex_curve::copy_rate_latency_list(OutputIt it) const {
	auto ray_it = rays.cbegin(), ray_end = rays.cend();
	if (ray_it->k == 1.0/0.0) ++ray_it;
	for(; ray_it != ray_end; ++ray_it) {
		double latency = ray_it->p.x - ray_it->p.y / ray_it->k;
		*it++ = {ray_it->k, latency};
	}
	return it;
}


template<typename InputIt>
void concave_curve::validate_ray_list(InputIt first, InputIt last) {
	// check validity of rays in the list
	for (auto iter = first; iter != last; ++iter) {
		detail::validate_ray(*iter);
	}

	// check restrictions on size and front side
	if (first == last) throw std::runtime_error("ray list is empty");
	if (first->p.x != 0.0 || first->p.y != 0.0) throw std::runtime_error("curve is not causal");

	// check restrictions to pairs of adjustent rays
	for(auto prev = first++; first != last; prev = first++) {
		if (prev->k <= first->k) throw std::runtime_error("curve is not convex");
		double x_diff = first->p.x - prev->p.x, y_diff = first->p.y - prev->p.y;
		if (y_diff <= 0.0) throw std::runtime_error("curve contains redundant rays");
		// y_diff / prev->k produces a valid result, since 0 < y_diff < inf && 0 < prev->k <= inf
		if (y_diff / prev->k != x_diff) throw std::runtime_error("curve is not continious");
	}
}

template<typename InputIt, typename base_curve::iterates_over<InputIt, rate_burst>::type*>
concave_curve::concave_curve(InputIt first, InputIt last) : concave_curve() {
	if (first == last) return;

	// order given rays by increasing bursts
	std::map<double, double> burst_to_rate;
	for(; first != last; ++first) {
		detail::validate_rate_burst(*first);
		auto res = burst_to_rate.emplace(first->burst, first->rate);
		if (!res.second && res.first->second > first->rate)
			res.first->second = first->rate;
	}

	// handle the front side of the curve
	auto it = burst_to_rate.begin(), it_end = burst_to_rate.end();
	if (it->first > 0.0) rays.emplace_back(0.0, it->first, it->second);
	else rays.back() = {0.0, it->first, it->second};

	// filter out non-meaningful rays
	// calculate the Xs for the inflexion points
	while (++it != it_end) {
		if (rays.back().k <= it->second) continue;

		// remove previosly added rays,
		// if they are covered by the current one
		for (;; rays.pop_back()) {
			double y_diff = it->first - rays.back().p.y;
			double k_diff = rays.back().k - it->second;
			double x_next = y_diff / k_diff;

			if (x_next > rays.back().p.x) {
				rays.emplace_back(x_next, it->first, it->second);
				break;
			}
		}
	}

	// calculate proper Y for each inflexion point
	auto rit = rays.begin(), rend = rays.end();
	while (++rit != rend && rit->p.x == 0.0);
	for(; rit != rend; ++rit) rit->p.y += rit->p.x * rit->k;
}

template<typename OutputIt>
OutputIt concave_curve::copy_rate_burst_list(OutputIt it) const {
	auto ray_it = rays.cbegin(), ray_end = rays.cend();
	if (ray_it->k == 0.0) ++ray_it;
	for(; ray_it != ray_end; ++ray_it) {
		double burst = ray_it->p.y - ray_it->p.x * ray_it->k;
		*it++ = {ray_it->k, burst};
	}
	return it;
}
