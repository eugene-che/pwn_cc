#include "pwn_cc.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <iterator>
#include <map>
#include <stdexcept>


namespace detail {
	void validate_ray(const line_t &ray) {
		if (!std::isfinite(ray.p.x) || ray.p.x < 0.0) throw std::runtime_error("invalid x");
		if (!std::isfinite(ray.p.y) || ray.p.y < 0.0) throw std::runtime_error("invalid x");
		if (!std::isgreaterequal(ray.k, 0.0)) throw std::runtime_error("invalid k");
	}

	void validate_rate_latency(const rate_latency &rl) {
		if (std::isnan(rl.rate) || rl.rate <= 0.0) throw std::runtime_error("invalid rate");
		if (!std::isfinite(rl.latency) || rl.latency < 0.0) throw std::runtime_error("invalid latency");
	}

	void validate_rate_burst(const rate_burst &rb) {
		if (!std::isfinite(rb.rate) || rb.rate < 0.0) throw std::runtime_error("invalid rate");
		if (!std::isfinite(rb.burst) || rb.burst < 0.0) throw std::runtime_error("invalid burst");
	}
}


bool operator == (const line_t &lhs, const line_t &rhs) {
	return lhs.p.x == rhs.p.x && lhs.p.y == rhs.p.y && lhs.k == rhs.k;
}

bool operator == (const base_curve &lhs, const base_curve &rhs) {
	return ray_handle::rays(lhs) == ray_handle::rays(rhs);
}


std::ostream& operator << (std::ostream& stream, const line_t &l) {
	return stream << "[ X = " << l.p.x << ", Y = " << l.p.y << ", K = " << l.k << " ]";
}

std::ostream& operator << (std::ostream& stream, const base_curve &curve) {
	std::ostream_iterator<line_t> out_it(stream, ", ");
	curve.copy_ray_list(out_it);
	return stream;
}


std::ostream& operator << (std::ostream& stream, const rate_latency &rl) {
	return stream << "[ R = " << rl.rate << ", L = " << rl.latency << " ]";
}

std::ostream& operator << (std::ostream& stream, const convex_curve &curve) {
	std::ostream_iterator<rate_latency> out_it(stream, ", ");
	curve.copy_rate_latency_list(out_it);
	return stream;
}


std::ostream& operator << (std::ostream& stream, const rate_burst &rb) {
	return stream << "[ R = " << rb.rate << ", B = " << rb.burst << " ]";
}

std::ostream& operator << (std::ostream& stream, const concave_curve &curve) {
	std::ostream_iterator<rate_burst> out_it(stream, ", ");
	curve.copy_rate_burst_list(out_it);
	return stream;
}


static
point_t crosspoint(const line_t &l1, const line_t &l2) {
	assert(l1.k != l2.k);
	if (std::isinf(l1.k)) {
		return {l1.p.x, l2.p.y + (l1.p.x - l2.p.x)*l2.k};
	} else
	if (std::isinf(l2.k)) {
		return {l2.p.x, l1.p.y + (l2.p.x - l1.p.x)*l1.k};
	} else {
		double b1 = l1.p.y - l1.p.x * l1.k;
		double b2 = l2.p.y - l2.p.x * l2.k;;
		double x = (b2 - b1) / (l1.k - l2.k);
		double y = b1 + x * l1.k;
		return {x, y};
	}
}

static
void push_ray(std::vector<line_t> &rays, line_t line) {
	if (rays.empty()) return rays.push_back(line);
	double x_diff = rays.back().p.x - line.p.x;
	double y_diff = rays.back().p.y - line.p.y;
	if (x_diff == 0.0 && y_diff == 0.0) return;
	double k_diff = rays.back().k - line.k;
	if (k_diff != 0.0 || y_diff != x_diff * line.k)
		return rays.push_back(line);
	rays.back().p.x = line.p.x;
	rays.back().p.y = line.p.y;
}

static
bool ends_above(const line_t &l1, const line_t &l2) {
	if (l1.k != l2.k) return l1.k > l2.k;
	if (std::isinf(l1.k)) return l1.p.x < l2.p.x;
	return l2.p.y < l1.p.y + (l2.p.x - l1.p.x) * l1.k;
}

void max(const ray_list &arg1, const ray_list &arg2, ray_list &res) {
	auto upper = std::make_pair(arg1.crbegin(), arg1.crend());
	auto lower = std::make_pair(arg2.crbegin(), arg2.crend());
	if (ends_above(*lower.first, *upper.first)) std::swap(upper, lower);

	res.clear();
	while (true) {
		auto &up = *upper.first, &lo = *lower.first;
		auto x_diff = up.p.x - lo.p.x;

		if (x_diff > 0.0) {
			if (up.p.y < lo.p.y + x_diff * lo.k) {
				push_ray(res, {crosspoint(lo, up), up.k});
				++upper.first;
				std::swap(upper, lower);
			} else push_ray(res, *upper.first++);
		} else
		if (x_diff < 0.0) {
			if (up.p.y < lo.p.y + x_diff * up.k) {
				push_ray(res, {crosspoint(lo, up), up.k});
				push_ray(res, *lower.first++);
				std::swap(upper, lower);
			} else ++lower.first;
		} else /* if (x_diff == 0.0) */ {
			if (up.p.y < lo.p.y) {
				push_ray(res, {crosspoint(lo, up), up.k});
				std::swap(upper, lower);
			}
			push_ray(res, *upper.first);
			// both curves start at the same point
			if (++lower.first == lower.second) break;
			if (++upper.first == upper.second) break;
		}
	}

	std::reverse(res.begin(), res.end());
}

void min(const ray_list &arg1, const ray_list &arg2, ray_list &res) {
	auto upper = std::make_pair(arg1.crbegin(), arg1.crend());
	auto lower = std::make_pair(arg2.crbegin(), arg2.crend());
	if (ends_above(*lower.first, *upper.first)) std::swap(upper, lower);

	res.clear();
	while (true) {
		auto &up = *upper.first, &lo = *lower.first;
		auto x_diff = up.p.x - lo.p.x;

		if (x_diff > 0.0) {
			if (up.p.y < lo.p.y + x_diff * lo.k) {
				push_ray(res, {crosspoint(lo, up), lo.k});
				push_ray(res, *upper.first++);
				std::swap(upper, lower);
			} else ++upper.first;
		} else
		if (x_diff < 0.0) {
			if (up.p.y < lo.p.y + x_diff * up.k) {
				push_ray(res, {crosspoint(lo, up), lo.k});
				++lower.first;
				std::swap(upper, lower);
			} else push_ray(res, *lower.first++);
		} else /* if (x_diff == 0.0) */ {
			if (up.p.y < lo.p.y) {
				push_ray(res, {crosspoint(lo, up), lo.k});
				std::swap(upper, lower);
			}
			push_ray(res, *lower.first);
			// both curves start at the same point
			if (++lower.first == lower.second) break;
			if (++upper.first == upper.second) break;
		}
	}

	std::reverse(res.begin(), res.end());
}

void sum(const ray_list &arg1, const ray_list &arg2, ray_list &res) {
	auto it1 = arg1.crbegin(), end1 = arg1.crend();
	auto it2 = arg2.crbegin(), end2 = arg2.crend();

	if (std::isinf(it1->k)) while (it2->p.x > it1->p.x) ++it2;
	if (std::isinf(it2->k)) while (it1->p.x > it2->p.x) ++it1;

	res.clear();
	while (true) {
		double x_diff = it1->p.x - it2->p.x;
		double y_sum = it1->p.y + it2->p.y;
		double k_sum = it1->k + it2->k;

		if (x_diff > 0.0) {
			y_sum += x_diff * it2->k;
			push_ray(res, {it1->p.x, y_sum, k_sum});
			++it1;
		} else
		if (x_diff < 0.0) {
			y_sum -= x_diff * it1->k;
			push_ray(res, {it2->p.x, y_sum, k_sum});
			++it2;
		} else
		if (x_diff == 0.0) {
			push_ray(res, {it1->p.x, y_sum, k_sum});
			// both curves start at the same point
			if (++it1 == end1 || ++it2 == end2) break;
		}
	}

	std::reverse(res.begin(), res.end());
}


// maximum of two convex curve is always a convex curve
convex_curve max(const convex_curve &conv1, const convex_curve &conv2) {
	convex_curve curve;
	max(ray_handle::rays(conv1), ray_handle::rays(conv2), ray_handle::rays(curve));
	return curve;
}

// minimum of two concave curves is always a concave curve
concave_curve min(const concave_curve &conc1, const concave_curve &conc2) {
	concave_curve curve;
	min(ray_handle::rays(conc1), ray_handle::rays(conc2), ray_handle::rays(curve));
	return curve;
}

// sum of two convex curve is always a convex curve
convex_curve sum(const convex_curve &conv1, const convex_curve &conv2) {
	convex_curve curve;
	sum(ray_handle::rays(conv1), ray_handle::rays(conv2), ray_handle::rays(curve));
	return curve;
}

// minimum of two concave curves is always a concave curve
concave_curve sum(const concave_curve &conc1, const concave_curve &conc2) {
	concave_curve curve;
	sum(ray_handle::rays(conc1), ray_handle::rays(conc2), ray_handle::rays(curve));
	return curve;
}



/* Leftover of a covex curve against a concave curve.
 *
 * Calculates segments of the overtop curve right-to-left
 *  whereas their reference points are above OX axis.
 * Crops the curve by OX as soon as it descends below.
 */
convex_curve leftover(const convex_curve &convex, const concave_curve &concave) {
	if (ray_handle::rays(concave).back().k >= ray_handle::rays(convex).back().k) return {}; // min curve

	auto convex_it = ray_handle::rays(convex).crbegin(), convex_end = ray_handle::rays(convex).crend();
	auto concave_it = ray_handle::rays(concave).crbegin(), concave_end = ray_handle::rays(concave).crend();
	if (std::isinf(convex_it->k)) while (concave_it->p.x > convex_it->p.x) ++concave_it;

	convex_curve curve;
	auto& res = ray_handle::rays(curve);
	res.clear();

	while (true) {
		double x_diff = convex_it->p.x - concave_it->p.x;
		double y_diff = convex_it->p.y - concave_it->p.y;
		double k_diff = convex_it->k - concave_it->k;

		res.emplace_back();
		auto &back = res.back();

		if (x_diff > 0.0) {
			back = {convex_it->p.x, y_diff - x_diff * concave_it->k, k_diff};
			++convex_it;
		} else
		if (x_diff < 0.0) {
			back = {concave_it->p.x, y_diff - x_diff * convex_it->k, k_diff};
			++concave_it;
		} else {
			// both curves start at the same point
			back = {convex_it->p.x, y_diff, k_diff};
			if (++convex_it == convex_end) break;
			if (++concave_it == concave_end) break;
		}

		if (back.p.y <= 0.0) {
			// crop the negative part of the curve
			back.p = {back.p.x - back.p.y / back.k, 0.0};
			if (back.p.x > 0.0) res.emplace_back(0.0, 0.0, 0.0);
			break;
		}
	}

	std::reverse(res.begin(), res.end());
	return curve;
}


/* Convolution for a pair of convex curves.
 *
 * The implementation is straightforward due to the following
 * THEOREM: Given two convex curves f(t) and g(t), curve of their
 *  convolution w(t) = (f*g)(t) coinsides with the one constructed by
 *  chaining the segments of f(t) and g(t) in the order of increasing slopes
 *  while starting from the origin of the coordinate plane.
 */
convex_curve convolve(const convex_curve &conv1, const convex_curve &conv2) {
	auto convex_it1 = ray_handle::rays(conv1).cbegin(), convex_end1 = ray_handle::rays(conv1).cend();
	auto convex_it2 = ray_handle::rays(conv2).cbegin(), convex_end2 = ray_handle::rays(conv2).cend();

	convex_curve curve;
	auto &res = ray_handle::rays(curve);
	res.clear();

	while (true) {
		double k_min = std::min(convex_it1->k, convex_it2->k);
		double x_sum = convex_it1->p.x + convex_it2->p.x;
		double y_sum = convex_it1->p.y + convex_it2->p.y;
		res.emplace_back(x_sum, y_sum, k_min);

		if (convex_it1->k == k_min && ++convex_it1 == convex_end1) break;
		if (convex_it2->k == k_min && ++convex_it2 == convex_end2) break;
	}

	return curve;
}


/* Deconvolution for a concave and a convex curve.
 *
 * This function is based on the following
 * THEOREM: Consider a convex concave f(t) and a convex curve g(t),
 *  such as f(t) == g(t) for exactly two arguments: 0 and t0.
 * For any t > 0 deconvolution w(t) = (f/g)(t) of these curves has
 *  the same value as a curve w0(t) constructed by chaining the segments of
 *  f(t) and g0(t) = {g(t) | 0 <= t <= t0} in the order of increasing slopes
 *  while starting from a point (-t0, -g(t0)).
 *
 * The above theorem allows to build w(t) as follows:
 * 1. find t0 -- the second crosspoint of f(t) and g(t);
 * 2. enumerate the segments of f(t) and g0(t) sorted by slopes;
 * 3. cut off the segments in the left half of the plane.
 *
 * However, there is a better idea:
 * One may also construct w(t) by shifting segments of f(t) and g(t):
 * - each segment of f(t) shifts by segments of g(t) with lower slopes.
 * - the segments of g(t) take the place among the segments of f(t).
 *
 * The implemented algorithm uses these 'shifts' as follows:
 * If forms segments of w(t) from right to left by mixing segments
 *  from the back of f(t) and from the front of g(t) as far as
 *  the shifted segments stay at the right half of the plane.
 * When staring point of the shifted segment becomes negative,
 *  the algorithms stops iteration and forms the front of w(t).
 */
concave_curve deconvolve(const concave_curve &conc, const convex_curve &conv) {
	if (ray_handle::rays(conc).back().k >= ray_handle::rays(conv).back().k) return {}; // max curve
	if (ray_handle::rays(conc).front().k <= ray_handle::rays(conv).front().k) return conc;
	// from now on the curves have exactly two crosspoints

	auto concave_it = ray_handle::rays(conc).crbegin();
	auto convex_it = ray_handle::rays(conv).cbegin(), convex_end = ray_handle::rays(conv).cend();
	while (convex_it->k <= concave_it->k) ++convex_it;

	concave_curve curve;
	auto &res = ray_handle::rays(curve);
	res.clear();

	double k_min, x_diff, y_diff;
	while (true) {
		k_min = std::min(convex_it->k, concave_it->k);
		x_diff = concave_it->p.x - convex_it->p.x;
		y_diff = concave_it->p.y - convex_it->p.y;

		// enforce shifted segments of to stay
		//  in the right half of the plane
		if (concave_it->p.x <= convex_it->p.x) break;
		res.emplace_back(x_diff, y_diff, k_min);

		if (concave_it->k == k_min) ++concave_it;
		if (convex_it->k == k_min && ++convex_it == convex_end) {
			// reach the last (endless) segment of the convex curve
			// => The following segments of concave curve should be
			//  shifted to only a fraction of this segment.
			// => The next shifted segment of concave curve
			//  will land into the left part of the plane.
			k_min = (--convex_it)->k;
			auto iter = ray_handle::rays(conc).crbegin();
			while (iter->p.y > iter->p.x * k_min) ++iter;
			point_t shift {crosspoint(*convex_it, *iter)};
			x_diff = concave_it->p.x - shift.x;
			y_diff = concave_it->p.y - shift.y;
			break;
		}
	}

	// deconvolution curve always cross OY above the origin
	res.emplace_back(0.0, y_diff - x_diff * k_min, k_min);
	res.emplace_back(0.0, 0.0, 1.0/0.0);
	std::reverse(res.begin(), res.end());
	return curve;
}


/* Maximal horizontal difference between a concave and a convex curve.
 *
 * For a given concave curve f(t) and a convex curve g(t)
 * their maximal difference max{t' - t"|f(t') == g(t")} equals to:
 *  - infinity, if f(t) exceeds g(t) for each t > 0;
 *  - zero, if t(t) <= g(t) for each t > 0;
 *  - length of a longest horizontal cross-section for a shape P,
 *   upper-bounded by f(t), and lower-bounded by g(t).
 *
 * Due to the properties of f(t) and g(t), shape P is convex polygon.
 * Thereby, its longest cross-section goes to one of its vertices,
 *  and we can find the maximal difference by checking each of them.
 * However, it is possible to reduce the search space:
 *  the algoritm calculates cross-sections for the vertices of P
 *  by increasing Ys, until the width of P starts to fall off
 *  (it works the same way as linear programming).
 */
double max_x_dist(const concave_curve &conc, const convex_curve &conv) {
	if (ray_handle::rays(conc).front().k <= ray_handle::rays(conv).front().k) return 0.0;
	if (ray_handle::rays(conc).back().k >= ray_handle::rays(conv).back().k) return 1.0/0.0;

	auto concave_it = ray_handle::rays(conc).cbegin(), concave_end = ray_handle::rays(conc).cend();
	auto convex_it = ray_handle::rays(conv).cbegin(), convex_end = ray_handle::rays(conv).cend();

	do {
		double y_diff = convex_it->p.y - concave_it->p.y;
		if (y_diff >= 0.0 && ++concave_it == concave_end) {
			for(--concave_it; concave_it->k > convex_it->k; ++convex_it);
			break;
		}
		if (y_diff <= 0.0 && ++convex_it == convex_end) {
			for(--convex_it; concave_it->k > convex_it->k; ++concave_it);
			break;
		}
	} while (concave_it->k > convex_it->k);

	double x_diff = convex_it->p.x - concave_it->p.x;
	double y_diff = convex_it->p.y - concave_it->p.y;

	// slopes of the segments can never be zero
	if (y_diff > 0.0) return x_diff - y_diff / concave_it->k;
	if (y_diff < 0.0) return x_diff - y_diff / convex_it->k;
	return x_diff;
}


/* Maximal vertical difference between a concave and a convex curve.
 *
 * For a given concave curve f(t) and a convex curve g(t)
 * their maximal difference max{f(t') - g(t")|t' == t"} equals to:
 *  - infinity, if f(t) exceeds g(t) for each t > 0;
 *  - zero, if t(t) <= g(t) for each t > 0;
 *  - length of a longest vertical cross-section for a shape P,
 *   upper-bounded by f(t), and lower-bounded by g(t).
 *
 * Due to the properties of f(t) and g(t), shape P is convex polygon.
 * Thereby, its longest cross-section goes to one of its vertices,
 *  and we can find the maximal difference by checking each of them.
 * However, it is possible to reduce the search space:
 *  the algoritm calculates cross-sections for the vertices of P
 *  by increasing Xs, until the height of P starts to fall off
 *  (it works the same way as linear programming).
 */
double max_y_dist(const concave_curve &conc, const convex_curve &conv) {
	if (ray_handle::rays(conc).front().k <= ray_handle::rays(conv).front().k) return 0.0;
	if (ray_handle::rays(conc).back().k >= ray_handle::rays(conv).back().k) return 1.0/0.0;

	auto concave_it = ray_handle::rays(conc).cbegin(), concave_end = ray_handle::rays(conc).cend();
	auto convex_it = ray_handle::rays(conv).cbegin(), convex_end = ray_handle::rays(conv).cend();

	do {
		double x_diff = convex_it->p.x - concave_it->p.x;
		if (x_diff >= 0.0 && ++concave_it == concave_end) {
			for(--concave_it; concave_it->k > convex_it->k; ++convex_it);
			break;
		}
		if (x_diff <= 0.0 && ++convex_it == convex_end) {
			for(--convex_it; concave_it->k > convex_it->k; ++concave_it);
			break;
		}
	} while (concave_it->k > convex_it->k);

	double x_diff = concave_it->p.x - convex_it->p.x;
	double y_diff = concave_it->p.y - convex_it->p.y;

	// slopes of the segments can never be infinite
	if (x_diff > 0.0) return y_diff - x_diff * convex_it->k;
	if (x_diff < 0.0) return y_diff - x_diff * concave_it->k;
	return y_diff;
}
