#version 330 core

flat in highp uint _id;

out uint elementID;

void main(void)
{
    elementID  = _id;

    // get distance from barycentric coordinates
    //vec3 dists = mod(_ids, 1);


    //elementID = uint(dists.x < dists.y ?
    //   (dists.x < dists.z ? int(_ids.x) : int(_ids.z)) :
    //   (dists.y < dists.z ? int(_ids.y) : int(_ids.z)));
}
