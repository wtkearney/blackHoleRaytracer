//
//  blackhole.h
//  
//
//  Created by Will Kearney
//  Spring 2015
//

#ifndef ____blackhole__
#define ____blackhole__

typedef struct {
    double location[3];         // camera location (r, theta, phi)
    double beta;                // camera speed
    double motion[3];           // camera components of motion relative to FIDO
    double orientation[3][3];   // in camera reference frame, right-handed set of three orthonormal vectors {ex, ey, ez}
    
    double localSky[2];         // spherical polar coord system for the camera local sky (for directions of incoming rays)
    double lensAngle;
} Camera;

typedef struct {
    double r;
    double a;
} BlackHole;

typedef struct {
    // ray incoming direction on camera's local sky
    double theta;
    double phi;
    
    // Cartesian components of the unit vector N that points in the direction of the incoming ray
    double Nx;
    double Ny;
    double Nz;
    
    // Direction of motion of incoming ray (using equations for relativistic aberration, as measured by the FIDO in Cartesian coordinates aligned with those of the camera)
    double Fx;
    double Fy;
    double Fz;
    double Fr;
    double Ftheta;
    double Fphi;
    
    // Canonical momenta (covariante coordinate components of its four-momentum)
    double Pt;
    double Pr;
    double Ptheta;
    double Pphi;
    
    // Axial angular momentum and Carter constant
    double b;
    double q;
} Ray;

typedef struct {
    double rho;
    double delta;
    double sigma;
    double alpha;
    double omega;
    double nu;
    
    double a;
} KerrMetric;

void raytrace(Image *src, Camera *camera, BlackHole *blackhole, KerrMetric *kerr);

/* create camera */
Camera *create_camera(void);

/* create black hole */
BlackHole *create_black_hole(void);

/* create Kerr metric */
KerrMetric *create_kerr_metric(void);

/* Set camera location on the celestial sphere; direction of motion relative to the FIDO there (a unit vector B in the camera reference frame (if equatorial = 1, cameria is in circular, equatorial geodesic orbit); camera speed beta relative to the FIDO; right-handed set of three orthonormal basis vectors {ex, ey, ez} with ey=B along the direction of the camera's motion; and initialize spherical polar system for the camera's local sky */
void set_camera_parameters(Camera *camera, KerrMetric *kerr, double r, double theta, double phi, double B_r, double B_theta, double B_phi, double angle, int equatorial);

/* normalize the three-component vector */
void vector_normalize(double *vector);

/* sets ray location on camera's local sky point */
void set_ray_location(Ray *ray, double theta, double phi);

/* computes, in camera's proper reference frame, the Cartesian components of the unit vector N that points in the direction of the incoming ray */
void set_ray_normal(Ray *ray);

/* using equations for relativistic aberration, compute the direction of motion of the incoming ray as measured by the FIDO in Cartesian coordinates aligned with those of the camera. Then, compute components on the FIDO's spherical orthonormal basis */
void set_ray_direction(Ray *ray, Camera *camera);

/* compute the ray's canonical momenta (covariant coordinate components of its four-momentum) with its conserved energy -Pt set to unity as a convenient convention. Then, compute axial angular momentum b and the Carter constant q (these are also conserved quantities) */
void set_ray_momenta(Ray *ray, KerrMetric *kerr);

void set_black_hole(BlackHole *BH, double r, double a);

/* Express Kerr metric with supplied Boyer-Lindquist coordinates */
void set_kerr_metric(KerrMetric *kerr, double a, double r, double theta, double phi);



#endif /* defined(____blackhole__) */
