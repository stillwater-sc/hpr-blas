#include <iostream>
#include <iomanip>

#define MTL_WITH_INITLIST 1
#include <boost/numeric/mtl/mtl.hpp>
#include <universal/number/posit/posit>

template<typename Real>
int countValue(const std::vector<Real>& V, const Real& ref) {
    int c = 0;
    for (auto v: V) {
        if (v == ref) ++c;
    }
    return c;
}

int main() {
    using namespace std;

    using Vector = mtl::dense_vector<float>;

    Vector v(10);

    v = 1;
    cout << v << endl;

    return EXIT_SUCCESS;
}
