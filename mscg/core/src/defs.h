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

#define RAD2DEG (180.0/M_PI)

#define vector_sub(c, a, b) \
    c[0] = a[0] - b[0]; \
    c[1] = a[1] - b[1]; \
    c[2] = a[2] - b[2]; \
    if(c[0]>box[0]*0.5) c[0]-=box[0]; else if(c[0]<-box[0]*0.5) c[0]+=box[0]; \
    if(c[1]>box[1]*0.5) c[1]-=box[1]; else if(c[1]<-box[1]*0.5) c[1]+=box[1]; \
    if(c[2]>box[2]*0.5) c[2]-=box[2]; else if(c[2]<-box[2]*0.5) c[2]+=box[2]

#define vector_add(c, a, b) \
    c[0] = a[0] + b[0]; \
    c[1] = a[1] + b[1]; \
    c[2] = a[2] + b[2];

#define vector_dot(a, b) (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])

#define vector_prod(c, a, b) \
    c[0] = a[1] * b[2] - a[2] * b[1]; \
    c[1] = a[2] * b[0] - a[0] * b[2]; \
    c[2] = a[0] * b[1] - a[1] * b[0]; 

#define vector_norm(a) { \
    float __s = 1.0 / sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]); \
    a[0] *= __s; a[1] *= __s; a[2] *= __s; \
}

#define vector_scale(b, a, s) \
    b[0] = a[0] * (s); \
    b[1] = a[1] * (s); \
    b[2] = a[2] * (s);

#endif
