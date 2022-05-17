#include <nanogui/nanogui.h>

#include "../clothMesh.h"
#include "../misc/sphere_drawing.h"
#include "sphere.h"

using namespace nanogui;
using namespace CGL;

void Sphere::collide(PointMass &pm) {
  // TODO (Part 3): Handle collisions with spheres.
  if ((pm.position - origin).norm() > radius) return;
  Vector3D tangent_point = (pm.position - origin).unit() * radius + origin;
  Vector3D correction_vector = tangent_point - pm.last_position;
  pm.position = pm.last_position + (1-friction)*correction_vector;
}

bool Sphere::collide_bool(PointMass &pm) {
    if ((pm.position - origin).norm() > radius) return false;
    return true;
}

void Sphere::drop(double gravity_y) {
    force_simu = 0.1 * -pow(abs(gravity_y), log(step));
    origin.y += 0.00001 * force_simu;
}

void Sphere::render(GLShader &shader) {
  // We decrease the radius here so flat triangles don't behave strangely
  // and intersect with the sphere when renderer
  if (!collided) {
      m_sphere_mesh.draw_sphere(shader, origin, radius * 0.92);
  }
}
