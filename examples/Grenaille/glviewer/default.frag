#version 150

in vec3 _position;
in vec3 _normal;

out vec4 FragColor;

vec3 lightPos = vec3(0.,0., 1.);

void main(void)
{
    float ddot = max(0., dot(_normal, normalize(lightPos - _position)));
    FragColor = vec4(ddot, ddot, ddot, 1.0);
}
