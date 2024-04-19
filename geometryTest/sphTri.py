from pylab import *

def mag(x):
    "The maginutde of a vector x"
    return sqrt(dot(x,x))

def sphTriAngle(a,b,c):
    "The solid angle of a spherical triangle defined by points a,b and c"
    a = a/mag(a)
    b = b/mag(b)
    c = c/mag(c)
    A = 2*arcsin(0.5*mag(b-c))
    B = 2*arcsin(0.5*mag(a-c))
    C = 2*arcsin(0.5*mag(a-b))
    s = 0.5*(A + B + C)
    t = tan(0.5*s)*tan(0.5*(s-A))*tan(0.5*(s-B))*tan(0.5*(s-C));
    return 4*arctan(sqrt(abs(t)))

def greatCircDist(a,b):
    "Great circle distance between points a and b on a sphere"
    a = a/mag(a)
    b = b/mag(b)
    return 2*arcsin(0.5*mag(a-b))

def greatTriangleArea(a,b,c):
    "area of a triangle on the surface of a sphere defined wrongly"
    ab = greatCircDist(b,a)*(b-a)/mag(b-a)
    ac = greatCircDist(c,a)*(c-a)/mag(c-a)
    return 0.5*mag(cross(ab,ac))

# define 3 points on the surface of a sphere
a = array([1,0,0])
b = array([1,1,1])
c = array([0.8,0.2,0])
a = a/mag(a)
b = b/mag(b)
c = c/mag(c)

print   'solid angle       = ', sphTriAngle(a,b,c), \
      '\nor                  ', sphTriAngle(c,a,b), \
      '\nor                  ', sphTriAngle(b,c,a), \
      '\ngreatTriangleArea = ', greatTriangleArea(a,b,c), \
      '\nor                  ', greatTriangleArea(c,a,b), \
      '\nor                  ', greatTriangleArea(b,c,a), \
