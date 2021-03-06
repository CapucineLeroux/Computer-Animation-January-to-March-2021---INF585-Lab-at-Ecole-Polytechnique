#pragma once

#include "vcl/vcl.hpp"

struct simulation_parameters
{
    float mass_total; // total mass of the cloth
    float K; // stiffness
    float mu; // damping    
};

struct obstacles_parameters
{
    float z_ground = 0.0f;
    vcl::vec3 sphere_center = {0.15f,0.2f,0.5f};
    float sphere_radius = 0.15f;
    vcl::vec3 p1 = {0.3f,0.5f,0.2f};
    vcl::vec3 p2 = {0.2f,0.5f,0.5f};
    float cylinder_radius = 0.08f;
};

void initialize_simulation_parameters(simulation_parameters& parameters, float L_cloth, size_t N_cloth);
void compute_forces(vcl::grid_2D<vcl::vec3>& force, vcl::grid_2D<vcl::vec3> const& position, vcl::grid_2D<vcl::vec3> const& velocity, vcl::grid_2D<vcl::vec3>& normals, float side_length, simulation_parameters const& parameters, float wind_magnitude);

void numerical_integration(vcl::grid_2D<vcl::vec3>& position, vcl::grid_2D<vcl::vec3>& velocity, vcl::grid_2D<vcl::vec3> const& forces, float mass, float dt);

void detect_self_collision(vcl::grid_2D<vcl::vec3>& position, vcl::grid_2D<vcl::vec3>& velocity, size_t k);
void apply_constraints(vcl::grid_2D<vcl::vec3>& position, vcl::grid_2D<vcl::vec3>& normal, vcl::grid_2D<vcl::vec3>& velocity, std::map<size_t, vcl::vec3> const& positional_constraints, obstacles_parameters const& obstacles);
bool detect_simulation_divergence(vcl::grid_2D<vcl::vec3> const& force, vcl::grid_2D<vcl::vec3> const& position);
