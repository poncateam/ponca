#version 150

in highp vec3 vertex;
in highp vec3 normal;

out vec3 _position;
out vec3 _normal;

uniform mat4 transform;

void main(void)
{
    mat4 MVI    = transpose(inverse(transform));
    _normal     = normalize(mat3(MVI) * normal);

    gl_Position = transform * vec4(vertex,1.0);
    _position   = gl_Position.xyz;
}
