#define BOOST_TEST_MODULE pwncc
#include <boost/test/unit_test.hpp>

#include "pwn_cc.hpp"


BOOST_AUTO_TEST_SUITE(pwncc_convex)
BOOST_AUTO_TEST_CASE(construct_from_rate_latency_list) {
	BOOST_CHECK_THROW((convex_curve{{0.0, 1.0}}), std::runtime_error);
	BOOST_CHECK_THROW((convex_curve{{1.0, 1.0/0.0}}), std::runtime_error);
	BOOST_CHECK_THROW((convex_curve{{1.0, 0.0/0.0}}), std::runtime_error);
}
BOOST_AUTO_TEST_CASE(construct_from_ray_list) {
	BOOST_CHECK_THROW((convex_curve{{1.0, 0.0, 0.0}}), std::runtime_error);
	BOOST_CHECK_THROW((convex_curve{{0.0, 0.0, 0.0/0.0}}), std::runtime_error);
	BOOST_CHECK_THROW((convex_curve{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}), std::runtime_error);
	BOOST_CHECK_THROW((convex_curve{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}}), std::runtime_error);
	BOOST_CHECK_THROW((convex_curve{{0.0, 0.0, 2.0}, {1.0, 2.0, 1.0}}), std::runtime_error);
}
BOOST_AUTO_TEST_CASE(construction_equality) {
	convex_curve rl_constructed {{1.0/0.0, 9.0}, {1.0, 2.0}, {2.0, 2.0}, {20.0, 20.0}};
	convex_curve ray_constructed {{0.0, 0.0, 0.0}, {2.0, 0.0, 2.0}, {9.0, 14.0, 1.0/0.0}};
	BOOST_CHECK_EQUAL(rl_constructed, ray_constructed);
}
BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(pwncc_concave)
BOOST_AUTO_TEST_CASE(construct_from_rate_burst_list) {
	BOOST_CHECK_THROW((concave_curve{{1.0/0.0, 0.0}}), std::runtime_error);
	BOOST_CHECK_THROW((concave_curve{{0.0, 1.0/0.0}}), std::runtime_error);
	BOOST_CHECK_THROW((concave_curve{{1.0, 0.0/0.0}}), std::runtime_error);
}
BOOST_AUTO_TEST_CASE(construct_from_ray_list) {
	BOOST_CHECK_THROW((concave_curve{{1.0, 0.0, 0.0}}), std::runtime_error);
	BOOST_CHECK_THROW((concave_curve{{0.0, 0.0, 0.0/0.0}}), std::runtime_error);
	BOOST_CHECK_THROW((concave_curve{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}), std::runtime_error);
	BOOST_CHECK_THROW((concave_curve{{0.0, 0.0, 2.0}, {1.0, 1.0, 1.0}}), std::runtime_error);
	BOOST_CHECK_THROW((concave_curve{{0.0, 0.0, 2.0}, {1.0, 2.0, 3.0}}), std::runtime_error);
}
BOOST_AUTO_TEST_CASE(construction_equality) {
	concave_curve rb_constructed {{2.0, 0.0}, {1.0, 5.0}, {3.0, 0.0}, {0.0, 12.0}};
	concave_curve ray_constructed {{0.0, 0.0, 2.0}, {5.0, 10.0, 1.0}, {7.0, 12.0, 0.0}};
	BOOST_CHECK_EQUAL(rb_constructed, ray_constructed);
}
BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(pwncc_max)
BOOST_AUTO_TEST_CASE(case1) {
	convex_curve c1 {{0.0, 0.0, 0.0}, {2.0, 0.0, 2.0}, {9.0, 14.0, 1.0/0.0}};
	convex_curve c2 {{0.0, 0.0, 0.0}, {1.0, 0.0, 1.0}, {4.0, 3.0, 3.0}};
	convex_curve expected {{0.0, 0.0, 0.0}, {1.0, 0.0, 1.0}, {3.0, 2.0, 2.0}, {5.0, 6.0, 3.0}, {9.0, 18.0, 1.0/0.0}};
	BOOST_CHECK_EQUAL(expected, max(c1, c2));
}
BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(pwncc_min)
BOOST_AUTO_TEST_CASE(case1) {
	concave_curve c1 {{0.0, 0.0, 2.0}, {5.0, 10.0, 1.0}, {7.0, 12.0, 0.0}};
	concave_curve c2 {{0.0, 0.0, 1.0/0.0}, {0.0, 4.0, 1.0}, {5.0, 9.0, 0.5}};
	concave_curve expected {{0.0, 0.0, 2.0}, {4.0, 8.0, 1.0}, {5.0, 9.0, 0.5}, {11.0, 12.0, 0.0}};
	BOOST_CHECK_EQUAL(expected, min(c1, c2));
}
BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(pwncc_sum)
BOOST_AUTO_TEST_CASE(case1) {
	convex_curve c1 {{0.0, 0.0, 0.0}, {2.0, 0.0, 2.0}, {9.0, 14.0, 1.0/0.0}};
	convex_curve c2 {{0.0, 0.0, 0.0}, {1.0, 0.0, 1.0}, {4.0, 3.0, 3.0}};
	convex_curve expected {{0.0, 0.0, 0.0}, {1.0, 0.0, 1.0}, {2.0, 1.0, 3.0}, {4.0, 7.0, 5.0}, {9.0, 32.0, 1.0/0.0}};
	BOOST_CHECK_EQUAL(expected, sum(c1, c2));
}
BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(pwncc_leftover)
BOOST_AUTO_TEST_CASE(case1) {
	convex_curve convex {{0.0, 0.0, 0.0}, {1.0, 0.0, 1.0}, {5.0, 4.0, 5.0}, {7.0, 14.0, 1.0/0.0}};
	concave_curve concave {{0.0, 0.0, 1.0/0.0}, {0.0, 1.0, 2.0}, {2.0, 5.0, 1.0}};
	convex_curve expected {{0.0, 0.0, 0.0}, {6.0, 0.0, 4.0}, {7.0, 4.0, 1.0/0.0}};
	BOOST_CHECK_EQUAL(expected, leftover(convex, concave));
}
BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(pwncc_deconvolve)
BOOST_AUTO_TEST_CASE(case1) {
	concave_curve concave {{0.0, 0.0, 1.0/0.0}, {0.0, 1.0, 2.0}, {2.0, 5.0, 1.0}};
	convex_curve convex {{0.0, 0.0, 0.0}, {1.0, 0.0, 1.0}, {5.0, 4.0, 5.0}, {7.0, 14.0, 1.0/0.0}};
	concave_curve expected {{0.0, 0.0, 1.0/0.0}, {0.0, 4.0, 1.0}};
	BOOST_CHECK_EQUAL(expected, deconvolve(concave, convex));
}
BOOST_AUTO_TEST_CASE(trivial_convex_curve) {
	concave_curve concave {{0.0, 0.0, 1.0/0.0}, {0.0, 4.0, 2.0}, {2.0, 8.0, 1.0}};
	convex_curve convex {{0.0, 0.0, 2.5}};
	concave_curve expected {{0.0, 0.0, 1.0/0.0}, {0.0, 4.0, 2.0}, {2.0, 8.0, 1.0}};
	BOOST_CHECK_EQUAL(expected, deconvolve(concave, convex));
}
BOOST_AUTO_TEST_CASE(trivial_concave_curve) {
	concave_curve concave {{0.0, 0.0, 2.0}};
	convex_curve convex {{0.0, 0.0, 0.0}, {1.0, 0.0, 1.0}, {5.0, 4.0, 5.0}, {7.0, 14.0, 1.0/0.0}};
	concave_curve expected {{0.0, 0.0, 1.0/0.0}, {0.0, 6.0, 2.0}};
	BOOST_CHECK_EQUAL(expected, deconvolve(concave, convex));
}
BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(pwncc_max_x_dist)
BOOST_AUTO_TEST_CASE(case1) {
	concave_curve concave {{0.0, 0.0, 1.0/0.0}, {0.0, 1.0, 2.0}, {2.0, 5.0, 1.0}};
	convex_curve convex {{0.0, 0.0, 0.0}, {1.0, 0.0, 1.0}, {5.0, 4.0, 5.0}, {7.0, 14.0, 1.0/0.0}};
	BOOST_CHECK_EQUAL(3.5, max_x_dist(concave, convex));
}
BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(pwncc_max_y_dist)
BOOST_AUTO_TEST_CASE(case1) {
	concave_curve concave {{0.0, 0.0, 1.0/0.0}, {0.0, 1.0, 2.0}, {2.0, 5.0, 1.0}};
	convex_curve convex {{0.0, 0.0, 0.0}, {1.0, 0.0, 1.0}, {5.0, 4.0, 5.0}, {7.0, 14.0, 1.0/0.0}};
	BOOST_CHECK_EQUAL(4.0, max_y_dist(concave, convex));
}
BOOST_AUTO_TEST_SUITE_END()
