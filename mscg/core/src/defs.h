#ifndef DEFS_H
#define DEFS_H

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>

typedef float vec3f[3];
typedef int   vec2i[2];
typedef int   vec3i[3];
typedef int   vec4i[4];

#define vector_sub(c, a, b) \
    float c##x = a[0] - b[0]; \
    float c##y = a[1] - b[1]; \
    float c##z = a[2] - b[2]; \
    if(c##x>box[0]*0.5) c##x-=box[0]; else if(c##x<-box[0]*0.5) c##x+=box[0]; \
    if(c##y>box[1]*0.5) c##y-=box[1]; else if(c##y<-box[1]*0.5) c##y+=box[1]; \
    if(c##z>box[2]*0.5) c##z-=box[2]; else if(c##z<-box[2]*0.5) c##z+=box[2]

#endif