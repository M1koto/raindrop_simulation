#version 330

uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

uniform vec4 u_color;

uniform sampler2D u_texture_2;
uniform vec2 u_texture_2_size;

uniform float u_normal_scaling;
uniform float u_height_scaling;

in vec4 v_position;
in vec4 v_normal;
in vec4 v_tangent;
in vec2 v_uv;

out vec4 out_color;

float h(vec2 uv) {
  // You may want to use this helper function...
  return texture(u_texture_2, uv).x;
}

void main() {
  // YOUR CODE HERE
  vec3 v_pos = vec3(v_position.x, v_position.y, v_position.z);
  vec3 n = normalize(vec3(v_normal.x, v_normal.y, v_normal.z));
  vec3 t = vec3(v_tangent.x, v_tangent.y, v_tangent.z);
  vec3 b = cross(n, t);
  mat3 TBN;
  TBN[0] = t;
  TBN[1] = b;
  TBN[2] = n;

  float dU = (h(vec2(v_uv.x + 1 / u_texture_2_size.x, v_uv.y)) - h(v_uv)) * u_normal_scaling * u_height_scaling;
  float dV = (h(vec2(v_uv.x, v_uv.y + 1 / u_texture_2_size.y)) - h(v_uv)) * u_normal_scaling * u_height_scaling;

  vec3 no = vec3(-dU, -dV, 1);
  vec3 nd = TBN * no;


  float r_sq = dot(u_light_pos - v_pos, u_light_pos - v_pos);
  vec3 v = normalize(u_cam_pos - v_pos);
  vec3 l = normalize(u_light_pos - v_pos);
  vec3 h = normalize(l + v);

  vec3 ld = (u_light_intensity / r_sq) * max(0.0, dot(nd, l));
  vec3 ls = (u_light_intensity / r_sq) * pow(max(0.0, dot(nd, h)), 64);
  vec3 ltotal = vec3(0.01, 0.01, 0.01) + ld + ls;
  // (Placeholder code. You will want to replace it.)
  out_color = vec4(ltotal, 1.0) * u_color;
  out_color.a = 1.0;
}
