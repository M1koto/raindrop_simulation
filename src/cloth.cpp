#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "cloth.h"
#include "collision/plane.h"
#include "collision/sphere.h"

using namespace std;

Cloth::Cloth(double width, double height, int num_width_points,
             int num_height_points, float thickness) {
  this->width = width;
  this->height = height;
  this->num_width_points = num_width_points;
  this->num_height_points = num_height_points;
  this->thickness = thickness;

  buildGrid();
  buildClothMesh();
}

Cloth::~Cloth() {
  point_masses.clear();
  springs.clear();

  if (clothMesh) {
    delete clothMesh;
  }
}

void Cloth::buildGrid() {
  // TODO (Part 1): Build a grid of masses and springs.

  // build masses
  double fixed = orientation == HORIZONTAL ? 1 : (2.0/1000) * (double)rand()/RAND_MAX + -1.0/1000;
  for (int j = 0; j < num_height_points; j++) {
    for (int i = 0; i < num_width_points; i++) {
      vector<int> v = { i, j};
      bool is_pinned = find(pinned.begin(), pinned.end(), v) != pinned.end();
      if (orientation == HORIZONTAL) {
        point_masses.emplace_back(Vector3D(width*i/num_width_points, 0, height*j/(num_height_points - 1)), is_pinned);
      } else {
        point_masses.emplace_back(Vector3D(width*i/num_width_points, height*j/(num_height_points - 1), (2.0/1000) * (double)rand()/RAND_MAX + -1.0/1000), is_pinned);
      }
    }
  }
  for (PointMass& p : point_masses) {
    p.last_position = p.position;
  }

  // build springs
  for (int i = 0; i < num_width_points; i++) {
    for (int j = 0; j < num_height_points; j++) {
      // structural
      // left
      if (i > 0) {
        springs.emplace_back(&point_masses[j*num_height_points + i], &point_masses[j*num_height_points + i - 1], STRUCTURAL);
      }
      // above
      if (j > 0) {
        springs.emplace_back(&point_masses[j*num_height_points + i], &point_masses[(j - 1)*num_height_points + i], STRUCTURAL);
      }

      // shearing
      // upper left
      if (i > 0 && j > 0) {
        springs.emplace_back(&point_masses[j*num_height_points + i], &point_masses[(j - 1)*num_height_points + i - 1], SHEARING);
      }
      // upper right
      if (i < num_width_points - 1 && j > 0) {
        springs.emplace_back(&point_masses[j*num_height_points + i], &point_masses[(j - 1)*num_height_points + i + 1], SHEARING);
      }

      // bending
      // two left
      if (i > 1) {
        springs.emplace_back(&point_masses[j*num_height_points + i], &point_masses[j*num_height_points + i - 2], BENDING);
      }
      // two above
      if (j > 1) {
        springs.emplace_back(&point_masses[j*num_height_points + i], &point_masses[(j - 2)*num_height_points + i], BENDING);
      }
    }
  }
}

int Cloth::idx(int y, int x) {
    return y * this->num_width_points + x;
}

void Cloth::waterDrop(int y, int x, double displacement) {
    // FIXME: using this function as rain drop initializer
    // FIXME: should probably use collision later on
    PointMass *pm = &this->point_masses[idx(y, x)];
    pm->start_position.y = displacement;
}

void Cloth::reset_raindrop() {
    // FIXME: using this function as rain drop initializer
    // FIXME: should probably use collision later on
    for (int i = 0; i < num_width_points; i++) {
      for (int j = 0; j < num_height_points; j++) {
        PointMass *pm = &this->point_masses[idx(i, j)];
        pm->start_position.y = 0;
      }
    }
}

bool Cloth::OutOfBound(int y, int x) {
    if (y < 0 || y >= this->num_height_points) return true;
    if (x < 0 || x >= this->num_width_points) return true;
    return false;
}

void Cloth::variedWaterDrop(int y, int x, int layers, double displacement, double ratio) {
    if (!OutOfBound(y, x)) {
        PointMass *pm = &this->point_masses[idx(y, x)];
        pm->start_position.y = displacement;
    }
    int move[4][2] = {{0, 1}, {1, 0}, {0, -1}, {-1, 0}};
    for (int layer = 1; layer < layers; layer++) {
        int row = y - layer;
        int column = x - layer;
        for (int m = 0; m < 4; m++) {
            for (int i = 0; i < 2 * layer; i++) {
                row += move[m][0];
                column += move[m][1];
                if (OutOfBound(row, column)) continue;
                PointMass *pm = &this->point_masses[idx(row, column)];
                pm->start_position.y = displacement * ratio;
            }
        }
        ratio *= ratio;
    }
}

void Cloth::rippleSpread(double frames_per_sec, double simulation_steps, ClothParameters *cp) {

    // TODO WaterRipple
    // fixme: only change the position on for 1 - boundary - 1 for now
    // fixme: might need change to if error occur
    for (int j = 0; j < this->num_height_points; j++) {
        for (int i = 0; i < this->num_width_points; i++) {
            int index = j * this->num_width_points + i;
            PointMass *p = &this->point_masses[index];
            p->last_position = p->position;
        }
    }
    double a = 0.5;
    double b = 1 - 4*a;
    // amplitude
    for (int j = 1; j < this->num_height_points - 1; j++) {
        for (int i = 1; i < this->num_width_points - 1; i++) {
            PointMass *p = &this->point_masses[idx(j, i)];
            PointMass *p_left = &this->point_masses[idx(j, i-1)];
            PointMass *p_right = &this->point_masses[idx(j, i+1)];
            PointMass *p_up = &this->point_masses[idx(j-1, i)];
            PointMass *p_down = &this->point_masses[idx(j+1, i)];
            // FIXME right shift ??? red underline on IDE
            p->position.y = ((p_left->last_position.y + p_right->last_position.y +
                              p_up->last_position.y + p_down->last_position.y) * a) + b*p->last_position.y;
            // attenuation
            // FIXME: rightshift
            p->position.y = p->position.y / 2;
        }
    }
}

void Cloth::ripple(PointMass& pm, double displacement) {
    pm.position.y -= displacement; // ?? pm.position or pm.start_position
}

void Cloth::stimulate(double frames_per_sec, double simulation_steps, ClothParameters *cp,
                     vector<Vector3D> external_accelerations,
                     vector<CollisionObject *> *collision_objects, double gravity) {
  double mass = width * height * cp->density / num_width_points / num_height_points;

  double delta_t = 1.0f / frames_per_sec / simulation_steps;

  double k = 1;
  double a = 10;
  double b = a/2;
  for (PointMass& p : point_masses) {
    p.new_position = p.position;
  }
  // calculate force for each point mass
  for (int j = 2; j < this->num_height_points - 2; j++) {
    for (int i = 2; i < this->num_width_points - 2; i++) {
      PointMass *p = &this->point_masses[idx(j, i)];
      PointMass *b1 = &this->point_masses[idx(j-1, i-1)];
      PointMass *a1 = &this->point_masses[idx(j, i-1)];
      PointMass *b3 = &this->point_masses[idx(j+1, i-1)];
      PointMass *a3 = &this->point_masses[idx(j-1, i)];
      PointMass *a4 = &this->point_masses[idx(j+1, i)];
      PointMass *b2 = &this->point_masses[idx(j-1, i+1)];
      PointMass *a2 = &this->point_masses[idx(j, i+1)];
      PointMass *b4 = &this->point_masses[idx(j+1, i+1)];

      double accel = b*((b1->position.y + b2->position.y + b3->position.y + b4->position.y) - 4*p->position.y)
          + a*((a1->position.y + a2->position.y + a3->position.y + a4->position.y) - 4*p->position.y);
      double vel = p->velocity(delta_t).y + accel*delta_t;
      double new_y = p->position.y + (1 - (cp->damping/100))*vel*delta_t;
      p->new_position = p->position;
      p->new_position.y = new_y;

//      double total_weight = 0;
//      double accel = 0;
//      for (int l = i - 2; l < i + 3; l++) {
//        for (int m = j - 2; m < j + 3; m++) {
//          if (l == i && m == j) continue;
//          PointMass *penis = &this->point_masses[idx(m, l)];
//          double dist2 = (penis->position - p->position).norm2();
//          accel += 1/dist2 * (penis->position.y - p->position.y);
//          total_weight += 1/dist2;
//        }
//      }
//      accel /= total_weight;
//      double vel = p->velocity(delta_t).y + accel*delta_t;
//      double new_y = p->position.y + (1 - (cp->damping/100))*vel*delta_t;
//      p->new_position = p->position;
//      p->new_position.y = new_y;
    }
  }
  for (CollisionObject* c : *collision_objects) {
      Sphere *s = (Sphere *) c;
      if (s->collided) continue;
      s->drop(gravity);
      s->step += 1;
      for (PointMass& pm: point_masses) {
          if (s->collide_bool(pm)) {
              double raindrop_displacement = -(pow(s->radius, 3) * s->force_simu);
              this->ripple(pm, raindrop_displacement);
              // hardcoded displacement first
              // TODO: change function parameter;
              s->collided = true;
              s->step = 0;
              break;
          }
      }
  }

  for (PointMass& p : point_masses) {
    p.last_position = p.position;
    p.position = p.new_position;
  }
}

void Cloth::stimulate_sine(double frames_per_sec, double simulation_steps, ClothParameters *cp,
                     vector<Vector3D> external_accelerations,
                     vector<CollisionObject *> *collision_objects, double gravity) {
  double delta_t = 1.0f / frames_per_sec / simulation_steps;

  double damping = 0.997;
  double wave_speed = 0.1;
  double freq = 10.0;
  for (PointMass& p : point_masses) {
    double dist = Vector3D(p.position.x - width / 2.0, 0.0, p.position.z - height / 2.0).norm();
    double displ = 0.0;
    if (wave_speed * p.time > dist) {
      displ = -p.amplitude * sin(freq * (p.time - dist / wave_speed));
    }
    p.position.y = displ;
    p.amplitude *= damping;
    p.time += delta_t;
  }
}

void Cloth::simulate(double frames_per_sec, double simulation_steps, ClothParameters *cp,
                     vector<Vector3D> external_accelerations,
                     vector<CollisionObject *> *collision_objects) {
  double mass = width * height * cp->density / num_width_points / num_height_points;
  double delta_t = 1.0f / frames_per_sec / simulation_steps;

  // TODO (Part 2): Compute total force acting on each point mass.
  Vector3D total_external = Vector3D();
  for (const Vector3D& v : external_accelerations) {
    total_external += v;
  }

  // total external force
  for (PointMass& p : point_masses) {
    p.forces = mass * total_external;
  }

  // spring correction forces
  for (const Spring& s : springs) {
    // spring type disabled
    if ((s.spring_type == CGL::STRUCTURAL && !cp->enable_structural_constraints)
      || (s.spring_type == CGL::SHEARING && !cp->enable_shearing_constraints)
      || (s.spring_type == CGL::BENDING && !cp->enable_bending_constraints)) {
      continue;
    }

    double Fs_magnitude = cp->ks * ((s.pm_a->position - s.pm_b->position).norm() - s.rest_length);
    if (s.spring_type == BENDING) Fs_magnitude *= 0.2;

    s.pm_a->forces += (s.pm_b->position - s.pm_a->position).unit() * Fs_magnitude;
    s.pm_b->forces += (s.pm_a->position - s.pm_b->position).unit() * Fs_magnitude;
  }


  // TODO (Part 2): Use Verlet integration to compute new point mass positions
  for (PointMass& p : point_masses) {
    if (p.pinned) continue;
    Vector3D new_position = p.position + (1 - (cp->damping/100))*(p.position - p.last_position)
        + pow(delta_t, 2)*p.forces/mass;
    p.last_position = p.position;
    p.position = new_position;
  }


  // TODO (Part 4): Handle self-collisions.
  build_spatial_map();
  for (PointMass& p : point_masses) {
    self_collide(p, simulation_steps);
  }


  // TODO (Part 3): Handle collisions with other primitives.
  for (PointMass& p : point_masses) {
    for (CollisionObject* c : *collision_objects) {
      c->collide(p);
    }
  }


  // TODO (Part 2): Constrain the changes to be such that the spring does not change
  // in length more than 10% per timestep [Provot 1995].
  for (Spring& s : springs) {
    if (s.pm_a->pinned && s.pm_b->pinned) continue;
    double correction_distance = (s.pm_a->position - s.pm_b->position).norm() - s.rest_length*1.1;
    if (correction_distance <= 0) continue;

    if (s.pm_a->pinned) {
      s.pm_b->position += (s.pm_a->position - s.pm_b->position).unit() * correction_distance;
    } else if (s.pm_b->pinned) {
      s.pm_a->position += (s.pm_b->position - s.pm_a->position).unit() * correction_distance;
    } else {
      s.pm_a->position += (s.pm_b->position - s.pm_a->position).unit() * correction_distance/2;
      s.pm_b->position += (s.pm_a->position - s.pm_b->position).unit() * correction_distance/2;
    }
  }

}

void Cloth::build_spatial_map() {
  for (const auto &entry : map) {
    delete(entry.second);
  }
  map.clear();

  // TODO (Part 4): Build a spatial map out of all of the point masses.
  for (PointMass& p: point_masses) {
    if (map.find(hash_position(p.position)) == map.end()) {
      vector<PointMass *> *v = new vector<PointMass *>;
      v->push_back(&p);
      map.insert({hash_position(p.position), v});
    }
    else {
      map.at(hash_position(p.position))->push_back(&p);
    }
  }

}

void Cloth::self_collide(PointMass &pm, double simulation_steps) {
  // TODO (Part 4): Handle self-collision for a given point mass.
  float key = hash_position(pm.position);
  // if (map.find(key) == map.end()) {
  //   return;
  // }
  vector<PointMass *> v = *map.at(key);
  Vector3D correction_vector = Vector3D();
  int num_vectors = 0;
  for (PointMass* &p : v) {
    if (p == &pm) {
      continue;
    }
    double dist = (p->position - pm.position).norm();
    if (dist <= 2 * thickness) {
      num_vectors++;
      correction_vector += (2 * thickness - dist) * (pm.position - p->position) / dist;
    }
  }
  if (num_vectors) correction_vector /= num_vectors;
  pm.position = pm.position + correction_vector / simulation_steps;

}

float Cloth::hash_position(Vector3D pos) {
  // TODO (Part 4): Hash a 3D position into a unique float identifier that represents membership in some 3D box volume.
  float w = width * 3.0 / num_width_points;
  float h = height * 3.0 / num_height_points;
  float t = w > h ? w : h;
  Vector3D truncated_pos = Vector3D(pos.x - fmod(pos.x, w), pos.y - fmod(pos.y, h), pos.z - fmod(pos.z, t));

  return (float) (truncated_pos.x + truncated_pos.y * width * height + truncated_pos.z * width * height * width * height);
}

///////////////////////////////////////////////////////
/// YOU DO NOT NEED TO REFER TO ANY CODE BELOW THIS ///
///////////////////////////////////////////////////////

void Cloth::reset() {
  PointMass *pm = &point_masses[0];
  double init_amplitude = 0.01;
  double init_radius = 0.02;
  for (int i = 0; i < point_masses.size(); i++) {
    pm->position = pm->start_position;
    pm->last_position = pm->start_position;
    pm->time = 0.0;
    double dist = Vector3D(pm->position.x - width / 2.0, 0.0, pm->position.z - height / 2.0).norm();
    pm->amplitude = dist < init_radius ? init_amplitude : init_amplitude / sqrt(dist / init_radius);
    pm++;
  }
}

void Cloth::buildClothMesh() {
  if (point_masses.size() == 0) return;

  ClothMesh *clothMesh = new ClothMesh();
  vector<Triangle *> triangles;

  // Create vector of triangles
  for (int y = 0; y < num_height_points - 1; y++) {
    for (int x = 0; x < num_width_points - 1; x++) {
      PointMass *pm = &point_masses[y * num_width_points + x];
      // Get neighboring point masses:
      /*                      *
       * pm_A -------- pm_B   *
       *             /        *
       *  |         /   |     *
       *  |        /    |     *
       *  |       /     |     *
       *  |      /      |     *
       *  |     /       |     *
       *  |    /        |     *
       *      /               *
       * pm_C -------- pm_D   *
       *                      *
       */

      float u_min = x;
      u_min /= num_width_points - 1;
      float u_max = x + 1;
      u_max /= num_width_points - 1;
      float v_min = y;
      v_min /= num_height_points - 1;
      float v_max = y + 1;
      v_max /= num_height_points - 1;

      PointMass *pm_A = pm                       ;
      PointMass *pm_B = pm                    + 1;
      PointMass *pm_C = pm + num_width_points    ;
      PointMass *pm_D = pm + num_width_points + 1;

      Vector3D uv_A = Vector3D(u_min, v_min, 0);
      Vector3D uv_B = Vector3D(u_max, v_min, 0);
      Vector3D uv_C = Vector3D(u_min, v_max, 0);
      Vector3D uv_D = Vector3D(u_max, v_max, 0);


      // Both triangles defined by vertices in counter-clockwise orientation
      triangles.push_back(new Triangle(pm_A, pm_C, pm_B,
                                       uv_A, uv_C, uv_B));
      triangles.push_back(new Triangle(pm_B, pm_C, pm_D,
                                       uv_B, uv_C, uv_D));
    }
  }

  // For each triangle in row-order, create 3 edges and 3 internal halfedges
  for (int i = 0; i < triangles.size(); i++) {
    Triangle *t = triangles[i];

    // Allocate new halfedges on heap
    Halfedge *h1 = new Halfedge();
    Halfedge *h2 = new Halfedge();
    Halfedge *h3 = new Halfedge();

    // Allocate new edges on heap
    Edge *e1 = new Edge();
    Edge *e2 = new Edge();
    Edge *e3 = new Edge();

    // Assign a halfedge pointer to the triangle
    t->halfedge = h1;

    // Assign halfedge pointers to point masses
    t->pm1->halfedge = h1;
    t->pm2->halfedge = h2;
    t->pm3->halfedge = h3;

    // Update all halfedge pointers
    h1->edge = e1;
    h1->next = h2;
    h1->pm = t->pm1;
    h1->triangle = t;

    h2->edge = e2;
    h2->next = h3;
    h2->pm = t->pm2;
    h2->triangle = t;

    h3->edge = e3;
    h3->next = h1;
    h3->pm = t->pm3;
    h3->triangle = t;
  }

  // Go back through the cloth mesh and link triangles together using halfedge
  // twin pointers

  // Convenient variables for math
  int num_height_tris = (num_height_points - 1) * 2;
  int num_width_tris = (num_width_points - 1) * 2;

  bool topLeft = true;
  for (int i = 0; i < triangles.size(); i++) {
    Triangle *t = triangles[i];

    if (topLeft) {
      // Get left triangle, if it exists
      if (i % num_width_tris != 0) { // Not a left-most triangle
        Triangle *temp = triangles[i - 1];
        t->pm1->halfedge->twin = temp->pm3->halfedge;
      } else {
        t->pm1->halfedge->twin = nullptr;
      }

      // Get triangle above, if it exists
      if (i >= num_width_tris) { // Not a top-most triangle
        Triangle *temp = triangles[i - num_width_tris + 1];
        t->pm3->halfedge->twin = temp->pm2->halfedge;
      } else {
        t->pm3->halfedge->twin = nullptr;
      }

      // Get triangle to bottom right; guaranteed to exist
      Triangle *temp = triangles[i + 1];
      t->pm2->halfedge->twin = temp->pm1->halfedge;
    } else {
      // Get right triangle, if it exists
      if (i % num_width_tris != num_width_tris - 1) { // Not a right-most triangle
        Triangle *temp = triangles[i + 1];
        t->pm3->halfedge->twin = temp->pm1->halfedge;
      } else {
        t->pm3->halfedge->twin = nullptr;
      }

      // Get triangle below, if it exists
      if (i + num_width_tris - 1 < 1.0f * num_width_tris * num_height_tris / 2.0f) { // Not a bottom-most triangle
        Triangle *temp = triangles[i + num_width_tris - 1];
        t->pm2->halfedge->twin = temp->pm3->halfedge;
      } else {
        t->pm2->halfedge->twin = nullptr;
      }

      // Get triangle to top left; guaranteed to exist
      Triangle *temp = triangles[i - 1];
      t->pm1->halfedge->twin = temp->pm2->halfedge;
    }

    topLeft = !topLeft;
  }

  clothMesh->triangles = triangles;
  this->clothMesh = clothMesh;
}
