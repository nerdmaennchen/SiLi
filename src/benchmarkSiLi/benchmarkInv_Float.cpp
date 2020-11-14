#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include <catch2/catch.hpp>

#include <SiLi/SiLi.h>
#include <armadillo>
#include <eigen3/Eigen/Dense>
#include <random>
#include <string>

TEST_CASE("1x1 Matrix inverse (float)", "[benchmark][inverse][float]") {
	SECTION("SiLi") {
		auto m1 = SiLi::Matrix{{{float{2}}}};

		BENCHMARK("inverting matrices") {
			auto z  = inv(m1);
			return z;
		};
	}
	SECTION("armadillo") {
		auto m1 = arma::Mat<float>(1, 1);
		m1(0, 0) = float{2};

		BENCHMARK("inverting matrices") {
			auto z  = arma::Mat<float>{inv(m1)};
			return z;
		};
	}
	SECTION("Eigen3") {
		using Matrix = Eigen::Matrix<float, 1, 1, 0, 1, 1>;
		auto m1 = Matrix{};
		m1(0, 0) = float{2};

		BENCHMARK("inverting matrices") {
			auto z  = Matrix{m1.inverse()};
			return z;
		};
	}
}

TEST_CASE("2x2 Matrix inverse (float)", "[benchmark][inverse][float]") {
	SECTION("SiLi") {
		auto m1 = SiLi::Matrix{{{float{2}, float{3}},
								 {float{4}, float{5}}}};

		BENCHMARK("inverting matrices") {
			auto z = inv(m1);
			return z;
		};
	}
	SECTION("armadillo") {
		auto m1 = arma::Mat<float>(2, 2);
		m1 = {{float{2}, float{3}},
		      {float{4}, float{5}}};

		BENCHMARK("inverting matrices") {
			auto z  = arma::Mat<float>{inv(m1)};
			return z;
		};
	}
	SECTION("Eigen3") {
		using Matrix = Eigen::Matrix<float, 2, 2, 0, 2, 2>;
		auto m1 = Matrix{};
		m1(0, 0) = float{2};
		m1(0, 1) = float{3};
		m1(1, 0) = float{4};
		m1(1, 1) = float{5};

		BENCHMARK("inverting matrices") {
			auto z  = Matrix{m1.inverse()};
			return z;
		};
	}
}

TEST_CASE("3x3 Matrix inverse (float)", "[benchmark][inverse][float]") {
	SECTION("SiLi") {
		auto m1 = SiLi::Matrix{{{float{  2}, float{  3}, float{  4}},
								 {float{  4}, float{  5}, float{  5}},
								 {float{100}, float{200}, float{300}}}};


		BENCHMARK("inverting matrices") {
			auto z  = inv(m1);
			return z;
		};
	}
	SECTION("armadillo") {
		auto m1 = arma::Mat<float>(3, 3);
		m1 = {{float{  2}, float{  3}, float{  4}},
		      {float{  4}, float{  5}, float{  5}},
		      {float{100}, float{200}, float{300}}};

		BENCHMARK("inverting matrices") {
			auto z  = arma::Mat<float>{inv(m1)};
			return z;
		};

	}
	SECTION("Eigen3") {
		using Matrix = Eigen::Matrix<float, 3, 3, 0, 3, 3>;
		auto m1 = Matrix{};
		m1(0, 0) = float{2};
		m1(0, 1) = float{3};
		m1(0, 2) = float{4};
		m1(1, 0) = float{4};
		m1(1, 1) = float{5};
		m1(1, 2) = float{5};
		m1(2, 0) = float{100};
		m1(2, 1) = float{200};
		m1(2, 2) = float{300};

		BENCHMARK("inverting matrices") {
			auto z  = Matrix{m1.inverse()};
			return z;
		};
	}
}

TEST_CASE("4x4 Matrix inverse (float)", "[benchmark][inverse][float]") {
	constexpr int N = 4;
	SECTION("SiLi") {
		auto m1 = SiLi::Matrix<N, N, float>{};

		auto gen = std::mt19937{N};
		auto dist = std::uniform_real_distribution<float>{-1000., 1000.};
		for (int row{0}; row < N; ++row) {
			for (int col{0}; col < N; ++col) {
				m1(row, col) = dist(gen);
			}
		}

		BENCHMARK("inverting matrices") {
			auto z  = inv(m1);
			return z;
		};
	}
	SECTION("armadillo") {
		auto m1 = arma::Mat<float>(N, N);

		auto gen = std::mt19937{N};
		auto dist = std::uniform_real_distribution<float>{-1000., 1000.};
		for (int row{0}; row < N; ++row) {
			for (int col{0}; col < N; ++col) {
				m1(row, col) = dist(gen);
			}
		}

		BENCHMARK("inverting matrices") {
			auto z  = arma::Mat<float>{inv(m1)};
			return z;
		};
	}
	SECTION("Eigen3") {
		using Matrix = Eigen::Matrix<float, N, N, 0, N, N>;
		auto m1 = Matrix{};

		auto gen = std::mt19937{N};
		auto dist = std::uniform_real_distribution<float>{-1000., 1000.};
		for (int row{0}; row < N; ++row) {
			for (int col{0}; col < N; ++col) {
				m1(row, col) = dist(gen);
			}
		}

		BENCHMARK("inverting matrices") {
			auto z  = Matrix{m1.inverse()};
			return z;
		};
	}
}

TEST_CASE("5x5 Matrix inverse (float)", "[benchmark][inverse][float]") {
	constexpr int N = 5;
	SECTION("SiLi") {
		auto m1 = SiLi::Matrix<N, N, float>{};

		auto gen = std::mt19937{N};
		auto dist = std::uniform_real_distribution<float>{-1000., 1000.};
		for (int row{0}; row < N; ++row) {
			for (int col{0}; col < N; ++col) {
				m1(row, col) = dist(gen);
			}
		}

		BENCHMARK("inverting matrices") {
			auto z  = inv(m1);
			return z;
		};
	}
	SECTION("armadillo") {
		auto m1 = arma::Mat<float>(N, N);

		auto gen = std::mt19937{N};
		auto dist = std::uniform_real_distribution<float>{-1000., 1000.};
		for (int row{0}; row < N; ++row) {
			for (int col{0}; col < N; ++col) {
				m1(row, col) = dist(gen);
			}
		}

		BENCHMARK("inverting matrices") {
			auto z  = arma::Mat<float>{inv(m1)};
			return z;
		};
	}
	SECTION("Eigen3") {
		using Matrix = Eigen::Matrix<float, N, N, 0, N, N>;
		auto m1 = Matrix{};

		auto gen = std::mt19937{N};
		auto dist = std::uniform_real_distribution<float>{-1000., 1000.};
		for (int row{0}; row < N; ++row) {
			for (int col{0}; col < N; ++col) {
				m1(row, col) = dist(gen);
			}
		}

		BENCHMARK("inverting matrices") {
			auto z  = Matrix{m1.inverse()};
			return z;
		};
	}
}

TEST_CASE("10x10 Matrix inverse (float)", "[benchmark][inverse][float]") {
	constexpr int N = 10;
	SECTION("SiLi") {
		auto m1 = SiLi::Matrix<N, N, float>{};

		auto gen = std::mt19937{N};
		auto dist = std::uniform_real_distribution<float>{-1000., 1000.};
		for (int row{0}; row < N; ++row) {
			for (int col{0}; col < N; ++col) {
				m1(row, col) = dist(gen);
			}
		}

		BENCHMARK("inverting matrices") {
			auto z  = inv(m1);
			return z;
		};
	}
	SECTION("armadillo") {
		auto m1 = arma::Mat<float>(N, N);

		auto gen = std::mt19937{N};
		auto dist = std::uniform_real_distribution<float>{-1000., 1000.};
		for (int row{0}; row < N; ++row) {
			for (int col{0}; col < N; ++col) {
				m1(row, col) = dist(gen);
			}
		}

		BENCHMARK("inverting matrices") {
			auto z  = arma::Mat<float>{inv(m1)};
			return z;
		};
	}
	SECTION("Eigen3") {
		using Matrix = Eigen::Matrix<float, N, N, 0, N, N>;
		auto m1 = Matrix{};

		auto gen = std::mt19937{N};
		auto dist = std::uniform_real_distribution<float>{-1000., 1000.};
		for (int row{0}; row < N; ++row) {
			for (int col{0}; col < N; ++col) {
				m1(row, col) = dist(gen);
			}
		}

		BENCHMARK("inverting matrices") {
			auto z  = Matrix{m1.inverse()};
			return z;
		};
	}
}

TEST_CASE("20x20 Matrix inverse (float)", "[benchmark][inverse][float]") {
	constexpr int N = 20;
	SECTION("SiLi") {
		auto m1 = SiLi::Matrix<N, N, float>{};

		auto gen = std::mt19937{N};
		auto dist = std::uniform_real_distribution<float>{-1000., 1000.};
		for (int row{0}; row < N; ++row) {
			for (int col{0}; col < N; ++col) {
				m1(row, col) = dist(gen);
			}
		}

		BENCHMARK("inverting matrices") {
			auto z  = inv(m1);
			return z;
		};
	}
	SECTION("armadillo") {
		auto m1 = arma::Mat<float>(N, N);

		auto gen = std::mt19937{N};
		auto dist = std::uniform_real_distribution<float>{-1000., 1000.};
		for (int row{0}; row < N; ++row) {
			for (int col{0}; col < N; ++col) {
				m1(row, col) = dist(gen);
			}
		}

		BENCHMARK("inverting matrices") {
			auto z  = arma::Mat<float>{inv(m1)};
			return z;
		};
	}
	SECTION("Eigen3") {
		using Matrix = Eigen::Matrix<float, N, N, 0, N, N>;
		auto m1 = Matrix{};

		auto gen = std::mt19937{N};
		auto dist = std::uniform_real_distribution<float>{-1000., 1000.};
		for (int row{0}; row < N; ++row) {
			for (int col{0}; col < N; ++col) {
				m1(row, col) = dist(gen);
			}
		}

		BENCHMARK("inverting matrices") {
			auto z  = Matrix{m1.inverse()};
			return z;
		};
	}
}
