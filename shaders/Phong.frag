#version 330

uniform vec4 u_color;
uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

in vec4 v_position;
in vec4 v_normal;
in vec2 v_uv;

out vec4 out_color;

void main() {
  // YOUR CODE HERE
  vec3 v_pos = vec3(v_position.x, v_position.y, v_position.z);
  vec3 n = normalize(vec3(v_normal.x, v_normal.y, v_normal.z));
  float r_sq = dot(u_light_pos - v_pos, u_light_pos - v_pos);
  vec3 v = normalize(u_cam_pos - v_pos);
  vec3 l = normalize(u_light_pos - v_pos);
  vec3 h = normalize(l + v);

  vec3 ld = (u_light_intensity / r_sq) * max(0.0, dot(n, l));
  vec3 ls = (u_light_intensity / r_sq) * pow(max(0.0, dot(n, h)), 64);
  vec3 ltotal = vec3(0.01, 0.01, 0.01) + ld + ls;
  // (Placeholder code. You will want to replace it.)
  out_color = vec4(ltotal, 1.0) * u_color;
  out_color.a = 1.0;
}
