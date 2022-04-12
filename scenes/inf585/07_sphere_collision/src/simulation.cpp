#include "simulation.hpp"

using namespace vcl;
using namespace std;


void simulate(std::vector<particle_structure>& particles, float dt)
{
	vec3 const g = {0,0,-9.81f};
	size_t const N = particles.size();
        float const alpha = 0.5f;
        float const beta = 0.5f;


	for (size_t k = 0; k < N; ++k)
	{
		particle_structure& particle = particles[k];

		vec3 const f = particle.m * g;

		particle.v = (1-0.9f*dt)*particle.v + dt*f;
		particle.p = particle.p + dt*particle.v;

                float r = particle.r;
                vec3 p = particle.p;

                vec3 a_surface;
                vec3 n_surface;
                float d;
                vec3 new_p;
                vec3 v_o;
                vec3 v_p;
                vec3 new_v;

                if (p.x-r<-1.f){
                    //cout<<"x < -1"<<endl;

                    a_surface = {-1,0,0};
                    n_surface = {1,0,0};

                    //cout<<particle.p<<endl;
                    d = r - dot(particle.p - a_surface,n_surface);
                    new_p = particle.p + d*n_surface;
                    particle.p = new_p;
                    //cout<<particle.p<<"\n"<<endl;

                    v_o = (dot(particle.v,n_surface))*n_surface;
                    v_p = particle.v - v_o;
                    new_v = alpha*v_p - beta*v_o;
                    particle.v = new_v;
                }

                if (p.x+r>1.f){
                    //cout<<"x > 1"<<endl;
                    a_surface = {1,0,0};
                    n_surface = {-1,0,0};

                    //cout<<particle.p<<endl;
                    d = r - dot(particle.p - a_surface,n_surface);
                    new_p = particle.p + d*n_surface;
                    particle.p = new_p;
                    //cout<<particle.p<<"\n"<<endl;

                    v_o = (dot(particle.v,n_surface))*n_surface;
                    v_p = particle.v - v_o;
                    new_v = alpha*v_p - beta*v_o;
                    particle.v = new_v;
                }

                if (p.y-r<-1.f){
                    //cout<<"y < -1"<<endl;

                    a_surface = {0,-1,0};
                    n_surface = {0,1,0};

                    //cout<<particle.p<<endl;
                    d = r - dot(particle.p - a_surface,n_surface);
                    new_p = particle.p + d*n_surface;
                    particle.p = new_p;
                    //cout<<particle.p<<"\n"<<endl;

                    v_o = (dot(particle.v,n_surface))*n_surface;
                    v_p = particle.v - v_o;
                    new_v = alpha*v_p - beta*v_o;
                    particle.v = new_v;
                }

                if (p.y+r>1.f){
                    //cout<<"y > 1"<<endl;

                    a_surface = {0,1,0};
                    n_surface = {0,-1,0};

                    //cout<<particle.p<<endl;
                    d = r - dot(particle.p - a_surface,n_surface);
                    new_p = particle.p + d*n_surface;
                    particle.p = new_p;
                    //cout<<particle.p<<"\n"<<endl;

                    v_o = (dot(particle.v,n_surface))*n_surface;
                    v_p = particle.v - v_o;
                    new_v = alpha*v_p - beta*v_o;
                    particle.v = new_v;
                }

                if (p.z-r<-1.f){
                    //cout<<"z < -1"<<endl;
                    a_surface = {0,0,-1};
                    n_surface = {0,0,1};

                    //cout<<particle.p<<endl;
                    d = r - dot(particle.p - a_surface,n_surface);
                    new_p = particle.p + d*n_surface;
                    particle.p = new_p;
                    //cout<<particle.p<<"\n"<<endl;

                    v_o = (dot(particle.v,n_surface))*n_surface;
                    v_p = particle.v - v_o;
                    new_v = alpha*v_p - beta*v_o;
                    particle.v = new_v;
                }

                if (p.z+r>1.f){
                    //cout<<"z > 1"<<endl;
                    a_surface = {0,0,1};
                    n_surface = {0,0,-1};

                    //cout<<particle.p<<endl;
                    d = r - dot(particle.p - a_surface,n_surface);
                    new_p = particle.p + d*n_surface;
                    particle.p = new_p;
                    //cout<<particle.p<<"\n"<<endl;

                    v_o = (dot(particle.v,n_surface))*n_surface;
                    v_p = particle.v - v_o;
                    new_v = alpha*v_p - beta*v_o;
                    particle.v = new_v;
                }

	}

	// To do :
	//  Handle collision ...

        const float mu = 0.5f;
        const float epsilon = 0.1f;
        float r1, r2;
        vec3 p1, p2;
        vec3 u;
        float m1, m2;
        float j;
        vec3 J;
        vec3 v1,v2;
        float d;
        for (size_t i1 = 0; i1 < N; ++i1)
        {
            for (size_t i2 = i1+1; i2 < N; ++i2){

                particle_structure& particle1 = particles[i1];
                particle_structure& particle2 = particles[i2];
                r1 = particle1.r;
                r2 = particle2.r;
                p1 = particle1.p;
                p2 = particle2.p;

                if (norm(p1-p2)<=r1+r2){

                    u = (p1-p2)/norm(p2-p1);
                    m1 = particle1.m;
                    m2 = particle2.m;
                    v1 = particle1.v;
                    v2 = particle2.v;
                    j = (2.f*(m1*m2)/(m1+m2))*dot(v2-v1,u);
                    J = j*u;

                    if (norm(v2-v1) > epsilon){
                        particle1.v = alpha*v1 + beta*J/m1;
                        particle2.v = alpha*v2 - beta*J/m2;
                    }

                    else {
                        particle1.v = mu*v1;
                        particle2.v = mu*v2;
                    }

                    d = r1 + r2 - norm(p1-p2);
                    particle1.p = p1 + (d/2.f)*u;
                    particle2.p = p2 - (d/2.f)*u;
                }

            }
        }
}
