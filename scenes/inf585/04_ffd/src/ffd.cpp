#include "ffd.hpp"



using namespace vcl;

// Helper to compute permutation(n, k) - k among n
//   avoids the recursive formula (too costly)
int binomial_coeff(int n, int k)
{
    int res = 1;
    if(k>n-k)
        k = n - k;
    for(int i=0; i<k; ++i) {
        res *= (n-i);
        res /= (i+1);
    }
    return res;
}

float b(int i, int N ,float u) //ith bernstein polynom applied to u
{
    return binomial_coeff(N,i)*std::pow(u,i)*std::pow(1-u,N-i);
}

buffer<grid_3D<float>> fill_weights(buffer<vec3> const position, grid_3D<vec3> const grid){

    buffer<grid_3D<float>> w_mijk;

    size_t const Nx = grid.dimension.x;
    size_t const Ny = grid.dimension.y;
    size_t const Nz = grid.dimension.z;
    size_t const N_vertex = position.size();

    float xm;
    float ym;
    float zm;
    float wm;
    for (size_t m = 0; m < N_vertex; ++m)
    {
        xm = position[m][0];
        ym = position[m][1];
        zm = position[m][2];

        w_mijk.push_back(grid_3D<float>(Nx,Ny,Nz));

            for (size_t i = 0; i < Nx; ++i) {
                    for (size_t j = 0; j < Ny; ++j) {
                            for (size_t k = 0; k < Nz; ++k) {
                                wm = b(i,Nx-1,xm)*b(j,Ny-1,ym)*b(k,Nz-1,zm);
                                w_mijk[m](i,j,k) = wm;
                            }
                    }
            }
    }

    return w_mijk;

};

// Computation of the FFD deformation on the position with respect to the grid
void ffd_deform(buffer<vec3>& position, grid_3D<vec3> const& grid, buffer<grid_3D<float>> w_mijk)
{

    // Get dimension of the grid
    size_t const Nx = grid.dimension.x;
    size_t const Ny = grid.dimension.y;
    size_t const Nz = grid.dimension.z;

    // Number of position to deform
    size_t const N_vertex = position.size();
    vec3 qm;
    for (size_t m = 0; m < N_vertex; ++m)
    {
        qm = {0,0,0};
        // Loop over all grid points
        for (size_t i = 0; i < Nx; ++i) {
                for (size_t j = 0; j < Ny; ++j) {
                        for (size_t k = 0; k < Nz; ++k) {

                                // TO DO: Should do some computations depending on the relative coordinates of the vertex k (you may need extra parameters to handle this), and the grid position.
                                // Note: A grid position can be accessed as grid(kx,ky,kz)
                                qm += w_mijk[m](i,j,k)*grid(i,j,k);
                        }
                }
        }
        position[m] = qm;

    }
}
