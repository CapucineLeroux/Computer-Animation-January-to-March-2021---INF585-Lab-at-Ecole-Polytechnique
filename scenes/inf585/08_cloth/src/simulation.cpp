#include "simulation.hpp"

using namespace vcl;
using namespace std;


// Fill value of force applied on each particle
// - Gravity
// - Drag
// - Spring force
// - Wind force
void compute_forces(grid_2D<vec3>& force, grid_2D<vec3> const& position, grid_2D<vec3> const& velocity, grid_2D<vec3>& normals, float side_length, simulation_parameters const& parameters, float wind_magnitude)
{
    size_t const N = force.size();        // Total number of particles of the cloth Nu x Nv
    int const N_dim = int(force.dimension[0]); // Number of particles along one dimension (square dimension)

    float const K  = parameters.K;
    float const m  = parameters.mass_total/N;
    float const mu = parameters.mu;
    float const	L0 = side_length/(N_dim-1.0f);

    // Gravity
    const vec3 g = {0,0,-9.81f};
    for(size_t k=0; k<N; ++k)
        force[k] = m*g;

    // Drag
    for(size_t k=0; k<N; ++k)
        force[k] += -mu*m*velocity[k];

    //wind
    vec3 w_direction = {wind_magnitude,0,0};
    float area = (norm(position(0,0)-position(0,N_dim-1))*norm(position(0,0)-position(N_dim-1,0)))/float(N);
    for(size_t k=0; k<N; ++k){
        if (k%(N) != N-1){
            //area = pow(norm(position[k+1]-position[k]),2);
            force[k] += area*dot(w_direction,normals[k])*normals[k];
        }
        else {
            //area = pow(norm(position[k-1]-position[k]),2);
            force[k] += area*dot(w_direction,normals[k])*normals[k];
        }
    }



    // TO DO: Add spring forces ...

    vec3 p, p_right, p_left, p_down, p_up, p1, p2, p3, p4, F_right, F_left, F_down, F_up, F1, F2, F3, F4;
    //1 = up right diagonal
    //2 = down right diagonal
    //3 = down left diagonal
    //4 = up left diagonal

    //bigger interior
    for (int j=1; j<N_dim-1; ++j){
        for (int i=1; i<N_dim-1; ++i){
            p = position(i,j);
            p_right = position(i,j+1);
            p_left = position(i,j-1);
            p_down = position(i+1,j);
            p_up = position(i-1,j);
            p1 = position(i-1,j+1);
            p2 = position(i+1,j+1);
            p3 = position(i+1,j-1);
            p4 = position(i-1,j-1);
            F_right = K*(norm(p_right-p)-L0)*(p_right-p)/norm(p_right-p);
            F_left = K*(norm(p_left-p)-L0)*(p_left-p)/norm(p_left-p);
            F_down = K*(norm(p_down-p)-L0)*(p_down-p)/norm(p_down-p);
            F_up = K*(norm(p_up-p)-L0)*(p_up-p)/norm(p_up-p);
            F1 = K*(norm(p1-p)-pow(2,0.5)*L0)*(p1-p)/norm(p1-p);
            F2 = K*(norm(p2-p)-pow(2,0.5)*L0)*(p2-p)/norm(p2-p);
            F3 = K*(norm(p3-p)-pow(2,0.5)*L0)*(p3-p)/norm(p3-p);
            F4 = K*(norm(p4-p)-pow(2,0.5)*L0)*(p4-p)/norm(p4-p);
            force(i,j) += F_right + F_left + F_down + F_up + F1 + F2 + F3 + F4;
        }
    }

    //corners
    p = position(0,0);
    p_right = position(0,1);
    p_down = position(1,0);
    p2 = position(1,1);
    F_right = K*(norm(p_right-p)-L0)*(p_right-p)/norm(p_right-p);
    F_down = K*(norm(p_down-p)-L0)*(p_down-p)/norm(p_down-p);
    F2 = K*(norm(p2-p)-pow(2,0.5)*L0)*(p2-p)/norm(p2-p);
    force(0,0) += F_right + F_down + F2;

    p = position(0,N_dim-1);
    p_left = position(0,N_dim-2);
    p_down = position(1,N_dim-1);
    p3 = position(1,N_dim-2);
    F_left = K*(norm(p_left-p)-L0)*(p_left-p)/norm(p_left-p);
    F_down = K*(norm(p_down-p)-L0)*(p_down-p)/norm(p_down-p);
    F3 = K*(norm(p3-p)-pow(2,0.5)*L0)*(p3-p)/norm(p3-p);
    force(0,N_dim-1) += F_left + F_down + F3;

    p = position(N_dim-1,0);
    p_right = position(N_dim-1,1);
    p_up = position(N_dim-2,0);
    p1 = position(N_dim-2,1);
    F_right = K*(norm(p_right-p)-L0)*(p_right-p)/norm(p_right-p);
    F_up = K*(norm(p_up-p)-L0)*(p_up-p)/norm(p_up-p);
    F1 = K*(norm(p1-p)-pow(2,0.5)*L0)*(p1-p)/norm(p1-p);
    force(N_dim-1,0) += F_right + F_up + F1;

    p = position(N_dim-1,N_dim-1);
    p_left = position(N_dim-1,N_dim-2);
    p_up = position(N_dim-2,N_dim-1);
    p4 = position(N_dim-2,N_dim-2);
    F_left = K*(norm(p_left-p)-L0)*(p_left-p)/norm(p_left-p);
    F_up = K*(norm(p_up-p)-L0)*(p_up-p)/norm(p_up-p);
    F4 = K*(norm(p4-p)-pow(2,0.5)*L0)*(p4-p)/norm(p4-p);
    force(N_dim-1,N_dim-1) += F_left + F_up + F4;

    //up borders
    for (int j=1; j<N_dim-1; ++j){
        p = position(0,j);
        p_right = position(0,j+1);
        p_left = position(0,j-1);
        p_down = position(1,j);
        p2 = position(1,j+1);
        p3 = position(1,j-1);
        F_right = K*(norm(p_right-p)-L0)*(p_right-p)/norm(p_right-p);
        F_left = K*(norm(p_left-p)-L0)*(p_left-p)/norm(p_left-p);
        F_down = K*(norm(p_down-p)-L0)*(p_down-p)/norm(p_down-p);
        F2 = K*(norm(p2-p)-pow(2,0.5)*L0)*(p2-p)/norm(p2-p);
        F3 = K*(norm(p3-p)-pow(2,0.5)*L0)*(p3-p)/norm(p3-p);
        force(0,j) += F_right + F_left + F_down + F2 + F3;
    }

    //down borders
    for (int j=1; j<N_dim-1; ++j){
        p = position(N_dim-1,j);
        p_right = position(N_dim-1,j+1);
        p_left = position(N_dim-1,j-1);
        p_up = position(N_dim-2,j);
        p1 = position(N_dim-2,j+1);
        p4 = position(N_dim-2,j-1);
        F_right = K*(norm(p_right-p)-L0)*(p_right-p)/norm(p_right-p);
        F_left = K*(norm(p_left-p)-L0)*(p_left-p)/norm(p_left-p);
        F_up = K*(norm(p_up-p)-L0)*(p_up-p)/norm(p_up-p);
        F1 = K*(norm(p1-p)-pow(2,0.5)*L0)*(p1-p)/norm(p1-p);
        F4 = K*(norm(p4-p)-pow(2,0.5)*L0)*(p4-p)/norm(p4-p);
        force(N_dim-1,j) += F_right + F_left + F_up + F1 + F4;
    }

    //left borders
    for (int i=1; i<N_dim-1; ++i){
        p = position(i,0);
        p_right = position(i,1);
        p_down = position(i+1,0);
        p_up = position(i-1,0);
        p1 = position(i-1,1);
        p2 = position(i+1,1);
        F_right = K*(norm(p_right-p)-L0)*(p_right-p)/norm(p_right-p);
        F_down = K*(norm(p_down-p)-L0)*(p_down-p)/norm(p_down-p);
        F_up = K*(norm(p_up-p)-L0)*(p_up-p)/norm(p_up-p);
        F1 = K*(norm(p1-p)-pow(2,0.5)*L0)*(p1-p)/norm(p1-p);
        F2 = K*(norm(p2-p)-pow(2,0.5)*L0)*(p2-p)/norm(p2-p);
        force(i,0) += F_right + F_down + F_up + F1 + F2;
    }

    //right borders
    for (int i=1; i<N_dim-1; ++i){
        p = position(i,N_dim-1);
        p_left = position(i,N_dim-2);
        p_down = position(i+1,N_dim-1);
        p_up = position(i-1,N_dim-1);
        p3 = position(i+1,N_dim-2);
        p4 = position(i-1,N_dim-2);
        F_left = K*(norm(p_left-p)-L0)*(p_left-p)/norm(p_left-p);
        F_down = K*(norm(p_down-p)-L0)*(p_down-p)/norm(p_down-p);
        F_up = K*(norm(p_up-p)-L0)*(p_up-p)/norm(p_up-p);
        F3 = K*(norm(p3-p)-pow(2,0.5)*L0)*(p3-p)/norm(p3-p);
        F4 = K*(norm(p4-p)-pow(2,0.5)*L0)*(p4-p)/norm(p4-p);
        force(i,N_dim-1) += F_left + F_down + F_up + F3 + F4;
    }

    //bending springs

    //smaller interior
    //2-neighbourhood
    for (int j=2; j<N_dim-2; ++j){
        for (int i=2; i<N_dim-2; ++i){
            p = position(i,j);
            p_right = position(i,j+2);
            p_left = position(i,j-2);
            p_down = position(i+2,j);
            p_up = position(i-2,j);
            F_right = K*(norm(p_right-p)-2.f*L0)*(p_right-p)/norm(p_right-p);
            F_left = K*(norm(p_left-p)-2.f*L0)*(p_left-p)/norm(p_left-p);
            F_down = K*(norm(p_down-p)-2.f*L0)*(p_down-p)/norm(p_down-p);
            F_up = K*(norm(p_up-p)-2.f*L0)*(p_up-p)/norm(p_up-p);
            force(i,j) += F_right + F_left + F_down + F_up;
        }
    }

    //corners
    for (int i=0 ; i<=1 ; i++){
        for (int j=0 ; j<=1 ; j++){
            p = position(i,j);
            p_right = position(i,j+2);
            p_down = position(i+2,j);
            F_right = K*(norm(p_right-p)-2.f*L0)*(p_right-p)/norm(p_right-p);
            F_down = K*(norm(p_down-p)-2.f*L0)*(p_down-p)/norm(p_down-p);
            force(i,j) += F_right + F_down;
        }
    }

    for (int i=0 ; i<=1 ; i++){
        for (int j=N_dim-2 ; j<=N_dim-1 ; j++){
            p = position(i,j);
            p_left = position(i,j-2);
            p_down = position(i+2,j);
            F_left = K*(norm(p_left-p)-2.f*L0)*(p_left-p)/norm(p_left-p);
            F_down = K*(norm(p_down-p)-2.f*L0)*(p_down-p)/norm(p_down-p);
            force(i,j) += F_left + F_down;
        }
    }

    for (int i=N_dim-2 ; i<=N_dim-1 ; i++){
        for (int j=0 ; j<=1 ; j++){
            p = position(i,j);
            p_right = position(i,j+2);
            p_up = position(i-2,j);
            F_right = K*(norm(p_right-p)-2.f*L0)*(p_right-p)/norm(p_right-p);
            F_up = K*(norm(p_up-p)-2.f*L0)*(p_up-p)/norm(p_up-p);
            force(i,j) += F_right + F_up;
        }
    }

    for (int i=N_dim-2 ; i<=N_dim-1 ; i++){
        for (int j=N_dim-2 ; j<=N_dim-1 ; j++){
            p = position(i,j);
            p_left = position(i,j-2);
            p_up = position(i-2,j);
            F_left = K*(norm(p_left-p)-2.f*L0)*(p_left-p)/norm(p_left-p);
            F_up = K*(norm(p_up-p)-2.f*L0)*(p_up-p)/norm(p_up-p);
            force(i,j) += F_left + F_up;
        }
    }

    //up double borders
    for (int j=2; j<N_dim-2; ++j){
        for (int i=0 ; i<=1 ; i++){
            p = position(i,j);
            p_right = position(i,j+2);
            p_left = position(i,j-2);
            p_down = position(i+2,j);
            F_right = K*(norm(p_right-p)-2.f*L0)*(p_right-p)/norm(p_right-p);
            F_left = K*(norm(p_left-p)-2.f*L0)*(p_left-p)/norm(p_left-p);
            F_down = K*(norm(p_down-p)-2.f*L0)*(p_down-p)/norm(p_down-p);
            force(i,j) += F_right + F_left + F_down;
        }
    }

    //down double borders
    for (int j=2; j<N_dim-2; ++j){
        for (int i=N_dim-2 ; i<=N_dim-1 ; i++){
            p = position(i,j);
            p_right = position(i,j+2);
            p_left = position(i,j-2);
            p_up = position(i-2,j);
            F_right = K*(norm(p_right-p)-2.f*L0)*(p_right-p)/norm(p_right-p);
            F_left = K*(norm(p_left-p)-2.f*L0)*(p_left-p)/norm(p_left-p);
            F_up = K*(norm(p_up-p)-2.f*L0)*(p_up-p)/norm(p_up-p);
            force(i,j) += F_right + F_left + F_up;
        }
    }

    //left double borders
    for (int i=2; i<N_dim-2; ++i){
        for (int j=0 ; j<=1 ; j++){
            p = position(i,j);
            p_right = position(i,j+2);
            p_down = position(i+2,j);
            p_up = position(i-2,j);
            F_right = K*(norm(p_right-p)-2.f*L0)*(p_right-p)/norm(p_right-p);
            F_down = K*(norm(p_down-p)-2.f*L0)*(p_down-p)/norm(p_down-p);
            F_up = K*(norm(p_up-p)-2.f*L0)*(p_up-p)/norm(p_up-p);
            force(i,j) += F_right + F_down + F_up;
        }
    }

    //right double borders
    for (int i=2; i<N_dim-2; ++i){
        for (int j=N_dim-2 ; j<=N_dim-1 ; j++){
            p = position(i,j);
            p_left = position(i,j-2);
            p_down = position(i+2,j);
            p_up = position(i-2,j);
            F_left = K*(norm(p_left-p)-2.f*L0)*(p_left-p)/norm(p_left-p);
            F_down = K*(norm(p_down-p)-2.f*L0)*(p_down-p)/norm(p_down-p);
            F_up = K*(norm(p_up-p)-2.f*L0)*(p_up-p)/norm(p_up-p);
            force(i,j) += F_left + F_down + F_up;
        }
    }


}

void numerical_integration(grid_2D<vec3>& position, grid_2D<vec3>& velocity, grid_2D<vec3> const& force, float mass, float dt)
{
    size_t const N = position.size();

    for(size_t k=0; k<N; ++k)
    {
        velocity[k] = velocity[k] + dt*force[k]/mass;
        position[k] = position[k] + dt*velocity[k];
    }

}


bool sphere_intersects_edge(vec3 p1, vec3 p2, float r, vec3 c){

    float beta = dot(c-p1,p2-p1)/dot(p2-p1,p2-p1);
    float alpha = 1.f - beta;
    vec3 orth_proj = alpha*p1 + beta*p2;

    if (0<=alpha && alpha<=1){
        return (norm(orth_proj-c)<r);
    }

    else if (norm(p1-c) <= norm(p2-c)){
        return (norm(p1-c)<r);
    }

    else if (norm(p2-c) <= norm(p1-c)){
        return (norm(p2-c)<r);
    }

    //return false;

}

/**vec3 minimal_point_face(vec3 p0, vec3 p1, vec3 p2, float r, vec3 c){

    vec3 p_min; //point of the triangle with minimal distance to c

    if ( fabs(dot(p1-p0,c-p0)) < pow(10,-2) ) {
        p_min = p0;
    }

    else {

        vec3 u1 = p1-p0;
        vec3 u2 = p2-p0;
        vec3 u = p0-c;

        float alpha = ( pow(dot(u1,u2),2) - pow(norm(u1),2)*pow(norm(u2),2) )/dot(u, pow(norm(u2),2)*u1 - dot(u1,u2)*u2 );
        float beta = -(dot(u,u2) + alpha*dot(u1,u2))/pow(norm(u2),2);

        if (0<=alpha && alpha<=1){
            p_min = p0 + alpha*u1 + beta*u2;
        }

        else {

            float l0 = norm(p0-c);
            float l1 = norm(p1-c);
            float l2 = norm(p2-c);

            if (l0<=l1 && l0<=l2){
                p_min = p0;
            }

            else if (l1<=l0 && l1<=l2){
                p_min = p1;
            }

            else if (l2<=l0 && l2<=l1){
                p_min = p2;
            }
        }

    }

    return p_min;

}*/


void detect_self_collision(grid_2D<vec3>& position, grid_2D<vec3>& velocity, size_t k){

    float epsilon = 5e-4f;
    vec3 p = position[k];
    vec3 v = velocity[k];
    int N_cloth = (int)sqrt((float)position.size());
    int ik = (int)k%N_cloth;
    int jk = int((int(k)-ik)/N_cloth);

    vec3 p1,p2,p3,n,proj,a;
    float a1,a2,a3;
    for (int i=0 ; i<N_cloth-4 ; i+=4){
        for (int j=0 ; j<N_cloth-4 ; j+=4){

            if ((i==ik && j==jk) or (i==ik-1 && j==jk-1)) {
                //ignore
            }
            else if (i==ik-1 && j==jk){
                //ignore lower triangle
                //upper triangle
                p1 = position(i,j);
                p2 = position(i+4,j+4);
                p3 = position(i,j+4);

                n = cross(p2-p1,p3-p1)/norm(cross(p2-p1,p3-p1));
                proj = p + dot(p1-p,n)*n;

                a = cross(p2-p1,p3-p1);
                a = a/pow(norm(a),2);
                a1 = dot(cross(p3-p2,proj-p2),a);
                a2 = dot(cross(p1-p3,proj-p3),a);
                a3 = dot(cross(p2-p1,proj-p1),a);

                if (norm(p-proj)<=epsilon && a1>=0 && a2>=0 && a3>=0 && fabs(a1+a2+a3-1)<=epsilon){
                    position[k] = proj + epsilon*n;
                    velocity[k] = v - dot(v,n)*n;
                }
            }
            else if (i==ik && j==jk-1){
                //ignore upper triangle
                //lower triangle
                p1 = position(i,j);
                p2 = position(i+4,j);
                p3 = position(i+4,j+4);

                n = cross(p2-p1,p3-p1)/norm(cross(p2-p1,p3-p1));
                proj = p + dot(p1-p,n)*n;

                a = cross(p2-p1,p3-p1);
                a = a/pow(norm(a),2);
                a1 = dot(cross(p3-p2,proj-p2),a);
                a2 = dot(cross(p1-p3,proj-p3),a);
                a3 = dot(cross(p2-p1,proj-p1),a);

                if (norm(p-proj)<=epsilon && a1>=0 && a2>=0 && a3>=0 && fabs(a1+a2+a3-1)<=epsilon){
                    position[k] = proj + epsilon*n;
                    velocity[k] = v - dot(v,n)*n;
                }
            }
            else {

                //upper triangle
                p1 = position(i,j);
                p2 = position(i+4,j+4);
                p3 = position(i,j+4);

                n = cross(p2-p1,p3-p1)/norm(cross(p2-p1,p3-p1));
                proj = p + dot(p1-p,n)*n;

                a = cross(p2-p1,p3-p1);
                a = a/pow(norm(a),2);
                a1 = dot(cross(p3-p2,proj-p2),a);
                a2 = dot(cross(p1-p3,proj-p3),a);
                a3 = dot(cross(p2-p1,proj-p1),a);

                if (norm(p-proj)<=epsilon && a1>=0 && a2>=0 && a3>=0 && fabs(a1+a2+a3-1)<=epsilon){
                    position[k] = proj + epsilon*n;
                    velocity[k] = v - dot(v,n)*n;
                }

                //lower triangle
                p1 = position(i,j);
                p2 = position(i+4,j);
                p3 = position(i+4,j+4);

                n = cross(p2-p1,p3-p1)/norm(cross(p2-p1,p3-p1));
                proj = p + dot(p1-p,n)*n;

                a = cross(p2-p1,p3-p1);
                a = a/pow(norm(a),2);
                a1 = dot(cross(p3-p2,proj-p2),a);
                a2 = dot(cross(p1-p3,proj-p3),a);
                a3 = dot(cross(p2-p1,proj-p1),a);

                if (norm(p-proj)<=epsilon && a1>=0 && a2>=0 && a3>=0 && fabs(a1+a2+a3-1)<=epsilon){
                    position[k] = proj + epsilon*n;
                    velocity[k] = v - dot(v,n)*n;
                }

            }
        }
    }

}


void apply_constraints(grid_2D<vec3>& position, grid_2D<vec3>& normal, grid_2D<vec3>& velocity, std::map<size_t, vec3> const& positional_constraints, obstacles_parameters const& obstacles)
{
    // Fixed positions of the cloth
    for(const auto& constraints : positional_constraints)
        position[constraints.first] = constraints.second;

    // To do: apply external constraints

    size_t const N = position.size();
    float r = obstacles.sphere_radius;
    vec3 p_sphere = obstacles.sphere_center;
    float z_min = obstacles.z_ground;
    vec3 p1 = obstacles.p1;
    vec3 p2 = obstacles.p2;
    float r_cylinder = obstacles.cylinder_radius;

    vec3 p, v,n,u,proj_position;
    float epsilon_ground = 5e-3f;
    float epsilon_sphere = 2*epsilon_ground;
    float epsilon_cylinder = 2*epsilon_ground;
    float t;
    for (size_t k=0 ; k<N ; k++){

            p = position[k];
            v = velocity[k];

            //collision with ground
            if (p.z<z_min+epsilon_ground){
                position[k].z = z_min+epsilon_ground;
                velocity[k].z = 0.f;
            }

            //collision with sphere
            if (norm(p-p_sphere)<r+epsilon_sphere) {

                n = (p-p_sphere)/norm(p-p_sphere);
                position[k] = p_sphere + (r+epsilon_sphere)*n;
                velocity[k] = v-dot(v,n)*n;

            }


            //collision with cylinder

            t = dot(p-p1,p2-p1)/dot(p2-p1,p2-p1);
            proj_position = (1.f-t)*p1 + t*p2;

            //check if vertex inside the bone cylinder
            if (-(r_cylinder+epsilon_cylinder)/norm(p1-p2)<=t && t<=1+(r_cylinder+epsilon_cylinder)/norm(p1-p2)) {

                if (0.f<=t && t<=1.f && norm(p-proj_position)<r_cylinder+epsilon_cylinder){
                    n = (p-proj_position)/norm(p-proj_position);
                    position[k] = proj_position + (r_cylinder+epsilon_cylinder)*n;
                    velocity[k] = v - dot(v,n)*n;
                }


                if (t<0.f && norm(p-p1)<r_cylinder+epsilon_cylinder){
                    u = (p-p1)/norm(p-p1);
                    position[k] = p1 + (r_cylinder+epsilon_cylinder)*u;
                    velocity[k] = v - dot(v,u)*u;
                }
                if (t>1.f && norm(p-p2)<r_cylinder+epsilon_cylinder){
                    u = (p-p2)/norm(p-p2);
                    position[k] = p2 + (r_cylinder+epsilon_cylinder)*u;
                    velocity[k] = v - dot(v,u)*u;
                }

            }

            //self_collision
            //detect_self_collision(position,velocity,k);

    }

}


void initialize_simulation_parameters(simulation_parameters& parameters, float L_cloth, size_t N_cloth)
{
	parameters.mass_total = 0.8f;
	parameters.K  = 5.0f;
	parameters.mu = 10.0f;
}

bool detect_simulation_divergence(grid_2D<vec3> const& force, grid_2D<vec3> const& position)
{
    bool simulation_diverged = false;
    const size_t N = position.size();
    for(size_t k=0; simulation_diverged==false && k<N; ++k)
    {
        const float f = norm(force[k]);
        const vec3& p = position[k];

        if( std::isnan(f) ) // detect NaN in force
        {
            std::cout<<"NaN detected in forces"<<std::endl;
            simulation_diverged = true;
        }

        if( f>600.0f ) // detect strong force magnitude
        {
            std::cout<<" **** Warning : Strong force magnitude detected "<<f<<" at vertex "<<k<<" ****"<<std::endl;
            simulation_diverged = true;
        }

        if( std::isnan(p.x) || std::isnan(p.y) || std::isnan(p.z) ) // detect NaN in position
        {
            std::cout<<"NaN detected in positions"<<std::endl;
            simulation_diverged = true;
        }
    }

    return simulation_diverged;

}
