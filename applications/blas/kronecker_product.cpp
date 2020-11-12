
#include <boost/numeric/mtl/mtl.hpp>

int main()
{
    using namespace mtl;
    
    dense2D<int> A(2, 2), B(2, 2), C(4, 4);
    
    for (size_t r= 0; r < 2; ++r)
        for (size_t c= 0; c < 2; ++c) {
            A[r][c]= (r+1) * 10 + c+1;
            B[r][c]= (r+1) * 1000 + (c+1) * 100;
        }
        
    //C= kron(A, B);
    //std::cout << "kron(A, B) is\n" << C;
    
    MTL_THROW_IF(C[0][0] != 12100, mtl::runtime_error("Wrong value in C[0][0]"));
    MTL_THROW_IF(C[3][3] != 48400, mtl::runtime_error("Wrong value in C[3][3]"));

    return 0;
}