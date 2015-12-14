
varying vec3 _normal;


void main(void)
{
    gl_FragColor = vec4(_normal.xyz, 1.0);
}
