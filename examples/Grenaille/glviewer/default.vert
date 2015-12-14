attribute highp vec3 vertex;
attribute highp vec3 normal;

varying vec3 _normal;

uniform mat4 transform;

void main(void)
{

   gl_Position = transform * vec4(vertex,1.0);
   _normal     = normalize(normal);
}
