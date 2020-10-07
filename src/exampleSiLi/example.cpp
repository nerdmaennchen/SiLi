#include <SiLi/SiLi.h>
#include <SiLi/ostream.h>

using namespace SiLi;
int main() {
	{
		auto m1 = Matrix{{{1., 0.},
		                  {0., 1.},
		                  {3., 4.}}};

		auto m2 = Matrix{{{1., 0.},
		                  {0., 1.},
		                  {3., 4.}}};

		auto m3 = m1 * view_trans(m2);
		auto m4 = view_trans(m2) * m1;

		for (auto& v : m1) {
			std::cout << "m1: " << v << "\n";
		}
		std::cout << "\n\n";
		for (auto& v : view<1, 0, 2, 2>(m1)) {
			std::cout << "m1: " << v << "\n";
		}
		std::cout << "\n\n";
		std::cout << m1 << "\n";
		std::cout << "\n\n";

		std::cout << view<1, 0, 3, 2>(m1) << "\n";


	}
}
