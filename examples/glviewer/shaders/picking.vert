#version 330 core

layout(location = 0) in highp vec3 vertex;
layout(location = 1) in highp vec3 normal;
layout(location = 2) in highp uint ids;

uniform mat4 transform, projection;

flat out highp uint _id;

void main(void)
{
    _id = ids;
    gl_Position = projection * transform * vec4(vertex,1.0);
}
