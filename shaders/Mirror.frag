#version 330


uniform vec3 u_cam_pos;
uniform samplerCube u_texture_cubemap;

in vec4 v_position;
in vec4 v_normal;
in vec4 v_tangent;
out vec4 out_color;

void main() {
  // YOUR CODE HERE
  vec3 v_pos = vec3(v_position.x, v_position.y, v_position.z);
  vec3 wo = v_pos - u_cam_pos;

  vec3 n = normalize(vec3(v_normal.x, v_normal.y, v_normal.z));
  vec3 wi = reflect(wo, n);
  vec3 wr = refract(wo, n, 0.75);
  out_color = texture(u_texture_cubemap, wi) * 0.3 + texture(u_texture_cubemap, wr) * 0.7;
  out_color.a = 1;
}

