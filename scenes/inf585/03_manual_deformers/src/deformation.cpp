#include "deformation.hpp"



using namespace vcl;

void apply_deformation(mesh& shape, // The position of shape are the one to be deformed
	vec2 const& tr,                 // Input gesture of the user in the 2D-screen coordinates - tr must be converted into a transformation applied to the positions of shape
	buffer<vec3> const& position_before_deformation,  // Initial reference position before the deformation
	buffer<vec3> const& normal_before_deformation,    // Initial reference normals before the deformation
	gui_widget const& widget,                         // Current values of the GUI widget
	picking_parameters const& picking,                // Information on the picking point
	rotation const& camera_orientation)               // Current camera orientation - allows to convert the 2D-screen coordinates into 3D coordinates
{

	float const r = widget.falloff; // radius of influence of the deformation
	size_t const N = shape.position.size();
	for (size_t k = 0; k < N; ++k)
	{
		vec3& p_shape = shape.position[k];                             // position to deform
		vec3 const& p_shape_original = position_before_deformation[k]; // reference position before deformation
		vec3 const& p_clicked = picking.p_clicked;                     // 3D position of picked position
                vec3 const& n_clicked = picking.n_clicked;                     // normal of the surface (before deformation) at the picked position

		float const dist = norm( p_clicked - p_shape_original );       // distance between the picked position and the vertex before deformation

		// TO DO: Implement the deformation models
		// **************************************************************** //
		// ...

		if (widget.deformer_type == deform_translate) // Case of translation
		{

                    if (widget.deformer_direction == direction_view_space) // Case of translation in the view space direction
                    {

                        // Hint: You can convert the 2D translation in screen space into a 3D translation in the view plane in multiplying
			//       camera_orientation * (tr.x, tr.y, 0)
			vec3 const translation = camera_orientation*vec3(tr,0.0f);

			// Fake deformation (linear translation in the screen space) 
			//   the following lines should be modified to get the expected smooth deformation
                        //if (dist < r)
                            //p_shape = p_shape_original + (1-dist/r)*translation;
                            p_shape = p_shape_original + (std::exp(-std::pow(dist/r,2)))*translation;
                    }


                    else if (widget.deformer_direction == direction_surface_normal) // Case of translation in the surface normal direction
                    {
                        vec3 const t = camera_orientation*vec3(tr,0.0f);
                        float signe = (t[0]*n_clicked[0] + t[1]*n_clicked[1] + t[2]*n_clicked[2])/std::fabs(t[0]*n_clicked[0] + t[1]*n_clicked[1] + t[2]*n_clicked[2]);
                        vec3 const translation = signe*norm(t)*n_clicked;

                        //if (dist < r)
                            p_shape = p_shape_original + (std::exp(-std::pow(dist/r,2)))*translation;


                    }


		}


		if (widget.deformer_type == deform_twist)
		{
			// Deformation to implement
                    if (widget.deformer_direction == direction_view_space) // Case of twist in the view space direction
                    {

                        vec3 R_axis = camera_orientation.matrix_row_z();

                        float w = std::exp(-std::pow(dist/r,2));
                        float R_angle;

                        if (tr.x >= 0){
                            R_angle = w*std::atan(tr.y/tr.x);
                        }
                        else {
                            if (tr.y >= 0){
                                R_angle = w*(std::atan(tr.y/tr.x)+3.14);
                            }
                            else {
                                R_angle = w*(std::atan(tr.y/tr.x)-3.14);
                            }
                        }

                        rotation const R(R_axis,R_angle);
                        p_shape = R*(p_shape_original-p_clicked) + p_clicked;

                    }


                    else if (widget.deformer_direction == direction_surface_normal) // Case of twist in the surface normal direction
                    {
                        vec3 R_axis = n_clicked;

                        float w = std::exp(-std::pow(dist/r,2));
                        float R_angle;

                        if (tr.x >= 0){
                            R_angle = w*std::atan(tr.y/tr.x);
                        }
                        else {
                            if (tr.y >= 0){
                                R_angle = w*(std::atan(tr.y/tr.x)+3.14);
                            }
                            else {
                                R_angle = w*(std::atan(tr.y/tr.x)-3.14);
                            }
                        }

                        rotation const R(R_axis,R_angle);
                        p_shape = R*(p_shape_original-p_clicked) + p_clicked;


                    }

		}


		if (widget.deformer_type == deform_scale)
		{
                    // Deformation to implement
                    float const scale = 1 + std::exp(-std::pow(dist/r,2))*tr.x;
                    p_shape = scale*(p_shape_original-p_clicked) + p_clicked;

		}

	}


}
