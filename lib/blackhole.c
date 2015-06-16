//
//  blackhole.c
//
//
//  Created by Will Kearney Spring 2015
//
//
//
// Contains functions to implement ray tracing in curved spacetime.
//
//
//

#include <stdio.h>
#include <math.h>
#include "graphics.h"


/* Do the actual raytracing */
void raytrace(Image *src, Camera *camera, BlackHole *blackhole, KerrMetric *kerr) {
    
    // TODO: better to walk across local sky with spherical coordinates? I think this might warp the view because have effectively made the viewing plane curved...
    double arcLengthHoriz = camera->lensAngle;
    double arcLengthVert = ((double)src->rows/(double)src->cols)*arcLengthHoriz;
    
    double horizStep = arcLengthHoriz/src->cols;
    double vertStep = arcLengthVert/src->rows;
    
    double rayPhi, rayTheta;
    
    printf("arcLengthHoriz: %f\n", arcLengthHoriz);
    printf("arcLengthVert: %f\n", arcLengthVert);

    printf("Raytracing...\n\n");

    // walk along the local sky, shooting a ray for each pixel
    for (int row=0;row<src->rows;row++) {
        rayTheta = (M_PI/2) - arcLengthVert/2 + row*vertStep;
        for (int col=0;col<src->cols;col++) {
            
            rayPhi = (M_PI) + arcLengthHoriz/2 - col*horizStep;
            
            Ray ray;
            set_ray_location(&ray, rayTheta, rayPhi);   // direction of incoming ray in camera spherical coord
            set_ray_normal(&ray);                       // direction of incoming ray in camera Cartesian coord
            set_ray_direction(&ray, camera);            // direction of motion of incoming ray in camera Cartesian coord
            set_ray_momenta(&ray, kerr);
            
//            printf("Ray info:\n");
//            printf("    Theta, Phi: (%f, %f)\n", ray.theta, ray.phi);
//            printf("    Nx, Ny, Nz: (%f, %f, %f)\n", ray.Nx, ray.Ny, ray.Nz);
//            printf("    Fx, Fy, Fz: (%f, %f, %f)\n", ray.Fx, ray.Fy, ray.Fz);
//            printf("    Fr, Ftheta, Fphi: (%f, %f, %f)\n", ray.Fr, ray.Ftheta, ray.Fphi);
//            printf("    Pt, Pr, Ptheta, Pphi: (%f, %f, %f, %f)\n", ray.Pt, ray.Pr, ray.Ptheta, ray.Pphi);
//            printf("    Axial angular momentum (b): %f\n", ray.b);
//            printf("    Carter constant (q): %f\n\n", ray.q);

            // check, using the ray constants b and q, whether it came from the horizon or the celestial sphere
            // first though, need to calculate some condition requirements
            double b1, b2 = 0.0;
            double r1, r2 = 0.0;
            double q0 = 0.0;
            double r0 = 0.0;
            // localize some things, also precalculate to keep things pretty & readable
            double b = ray.b;
            double a = kerr->a;
            double q = ray.q;
            double cubeRoot2 = pow(2, 1/3.);    // 1.259921049894873...
            
            // calculate the radius bounds of unstably trapped photons
            r1 = 2.0*(1.0 + cos((2.0/3.0)*acos(-a)));
            r2 = 2.0*(1.0 + cos((2.0/3.0)*acos(a)));
            // printf("Radius interval of unstably trapped photons [r1, r2]: %f, %f\n", r1, r2);
            
            // calculate a constant of motion for photons on the radius bounds
            b2 = -((r1*r1*r1) - 3*r1*r1 + a*a*r1 + a*a) / (a*(r1 - 1));
            b1 = -((r2*r2*r2) - 3*r2*r2 + a*a*r2 + a*a) / (a*(r2 - 1));
            printf("b1, b2: %f, %f\n", b1, b2);
            
            
            // get the value of r0 from parametric equations b(r0) in order to calculate q(r0)
            double tmpQ = (a*a + a*b - 3);
            double tmp = pow(sqrt(108.0*tmpQ*tmpQ*tmpQ + (54 - 54.0*a*a)*(54 - 54.0*a*a)) - 54.0*a*a + 54.0, 1/3.);
            r0 = -(cubeRoot2*tmpQ / tmp) + (tmp / 3.0*cubeRoot2) + 1.0;
            
            // printf("r0: %f\n", r0);
            double r03 = r0*r0*r0;
            
            // calculate another constant of motion
            q0 = -(r03*(r03 - 6.0*r0*r0 + 9.0*r0 - 4.0*a*a)) / (a*a*(b - 1)*(b - 1));
            
            // printf("q0(r0): %f\n", q0);
            
            if ((b1 < b && b < b2) && q < q0) {
                // there are no radial turning points for this {b, q}
                if (ray.Pr > 0) {
                    // ray came from horizon
                    src->data[row][col].rgb[0] = 1.0;
                    src->data[row][col].rgb[1] = 0.0;
                    src->data[row][col].rgb[2] = 0.0;
                } else {
                    // ray came from celestial sphere
                    src->data[row][col].rgb[0] = 0.0;
                    src->data[row][col].rgb[1] = 1.0;
                    src->data[row][col].rgb[2] = 0.0;
//                    src->data[row][col].rgb[0] = 176.0/255;
//                    src->data[row][col].rgb[1] = 196.0/255;
//                    src->data[row][col].rgb[2] = 222.0/255;
                }
                
            } else {
                // there are two radial turning points...
                // need to calculate some roots of equations for a null geodesic
                double rUp = 0.0;
                double delta = kerr->delta;
            
                double P = sqrt(delta*((b-a)*(b-a) + q));
                printf("P: %f\n", P);
                printf("a: %f\n", a);
                // the two positive real roots of R(r) = 0; we want the bigger one
                double rUp1 = -a*a - P + a*b;
                double rUp2 = -a*a + P + a*b;
                
                printf("rUp1: sqrt(%f)\n", rUp1);
                printf("rUp2: sqrt(%f)\n\n", rUp2);
                
                if (rUp1 > rUp2) {
                    rUp = rUp1;
                } else {
                    rUp = rUp2;
                }
                           
                // printf("rUp: %f\n", rUp);
                

                if (camera->location[0] >= rUp) {
                    // ray came from celestial sphere
                    src->data[row][col].rgb[0] = 0.0;
                    src->data[row][col].rgb[1] = 0.0;
                    src->data[row][col].rgb[2] = 1.0;
//                    src->data[row][col].rgb[0] = 176.0/255;
//                    src->data[row][col].rgb[1] = 196.0/255;
//                    src->data[row][col].rgb[2] = 222.0/255;
                } else {
                    // ray came from horizon
                    src->data[row][col].rgb[0] = 0.5;
                    src->data[row][col].rgb[1] = 0.0;
                    src->data[row][col].rgb[2] = 0.5;
                }
            }
            

            
        }
        
    }
  
}


/* create camera */
Camera *create_camera(void) {
    Camera *camera = NULL;
    camera = (Camera *)malloc(sizeof(Camera));
    
    return(camera);
}

/* create black hole */
BlackHole *create_black_hole(void) {
    BlackHole *blackhole = NULL;
    blackhole = (BlackHole *)malloc(sizeof(BlackHole));
    
    return(blackhole);
}

/* create Kerr metric */
KerrMetric *create_kerr_metric(void) {
    KerrMetric *kerr = NULL;
    kerr = (KerrMetric *)malloc(sizeof(KerrMetric));
    
    return(kerr);
}

/* Set camera location on the celestial sphere; direction of motion relative to the FIDO there (a unit vector B in the camera reference frame (if equatorial = 1, cameria is in circular, equatorial geodesic orbit); camera speed beta relative to the FIDO; right-handed set of three orthonormal basis vectors {ex, ey, ez} with ey=B along the direction of the camera's motion; and initialize spherical polar system for the camera's local sky */
void set_camera_parameters(Camera *camera, KerrMetric *kerr, double r, double theta, double phi, double B_r, double B_theta, double B_phi, double angle, int equatorial) {
    
    // set location
    camera->location[0] = r;
    camera->location[1] = theta;
    camera->location[2] = phi;
    
    // set speed (dictated by Kerr metric)
    double bigOmega = 1 / (kerr->a + pow(camera->location[0], 1.5));     // this is handy to compute
    camera->beta = (kerr->nu / kerr->alpha)*(bigOmega - kerr->omega);
    
    // set direction of motion, checking if in circular, equatorial geodesic orbit
    if (equatorial == 1) {
        camera->motion[0] = 0;
        camera->motion[1] = 0;
        camera->motion[2] = 1;
    } else {
        camera->motion[0] = B_r;
        camera->motion[1] = B_theta;
        camera->motion[2] = B_phi;
    }

    // set orthonormal basis vectors {ex, ey, ez} with ey=B
    //memcpy(camera->orientation[1], camera->motion, sizeof(camera->motion));
    //memcpy(camera->orientation[0], camera->location, sizeof(camera->location));  // this means camera will look at black hole
    
    // use cross product to calculate third orthogonal vector
//    camera->orientation[2][0] = camera->orientation[0][1]*camera->orientation[1][2] - camera->orientation[0][2]*camera->orientation[1][1];
//    camera->orientation[2][1] = camera->orientation[0][2]*camera->orientation[1][0] - camera->orientation[0][0]*camera->orientation[1][2];
//    camera->orientation[2][2] = camera->orientation[0][0]*camera->orientation[1][1] - camera->orientation[0][1]*camera->orientation[1][0];
    
    // initialize spherical polar system for local sky
    camera->localSky[0] = 0.0;
    camera->localSky[1] = 0.0;
    
    camera->lensAngle = angle;

}

/* normalize the three-component vector */
void vector_normalize(double *vector) {
    double length = sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
    
    vector[0] = vector[0] / length;
    vector[1] = vector[1] / length;
    vector[2] = vector[2] / length;
}

/* sets ray location on camera's local sky point */
void set_ray_location(Ray *ray, double theta, double phi) {
    ray->theta = theta;
    ray->phi = phi;
}

/* computes, in camera's proper reference frame, the Cartesian components of the unit vector N that points in the direction of the incoming ray */
void set_ray_normal(Ray *ray) {

    ray->Nx = sin(ray->theta)*cos(ray->phi);
    ray->Ny = sin(ray->theta)*sin(ray->phi);
    ray->Nz = cos(ray->theta);
    
    // normalize
    double length = sqrt(ray->Nx*ray->Nx + ray->Ny*ray->Ny + ray->Nz*ray->Nz);
    ray->Nx = ray->Nx / length;
    ray->Ny = ray->Ny / length;
    ray->Nz = ray->Nz / length;
}

/* using equations for relativistic aberration, compute the direction of motion of the incoming ray as measured by the FIDO in Cartesian coordinates aligned with those of the camera. Then, compute components on the FIDO's spherical orthonormal basis */
void set_ray_direction(Ray *ray, Camera *camera) {

    // localize some terms because cleanliness is next to godliness
    double beta = camera->beta;
    double Br = camera->motion[0];
    double Btheta = camera->motion[1];
    double Bphi = camera->motion[2];
    double Nx = ray->Nx;
    double Ny = ray->Ny;
    double Nz = ray->Nz;
    
    // compute direction of motion in Cartesian coord aligned with camera
    double Fx = (-sqrt(1 - beta*beta)*Nx)/(1 - beta*Ny);
    double Fy = (-Ny + beta)/(1 - beta*Ny);
    double Fz = (-sqrt(1 - beta*beta)*Nz)/(1 - beta*Ny);
    
    // normalize
    // TODO: not sure if I should normalize these. But I feel like I should because they comprise a direction of motion vector
    double length = sqrt(Fx*Fx + Fy*Fy + Fz*Fz);
    Fx = Fx / length;
    Fy = Fy / length;
    Fz = Fz / length;
    
    // compute components of ray on the FIDO's spherical orthonormal basis
    double k = sqrt(1 - Btheta*Btheta);
    
    double Fr = (Bphi/k)*Fx + Br*Fy + (Br*Btheta/k)*Fz;
    double Ftheta = Btheta*Fy - k*Fz;
    double Fphi = -(Br/k)*Fx + Bphi*Fy + (Btheta*Bphi/k)*Fz;
    
    // commit to ray
    ray->Fx = Fx;
    ray->Fy = Fy;
    ray->Fz = Fz;

    ray->Fr = Fr;
    ray->Ftheta = Ftheta;
    ray->Fphi = Fphi;
}

/* compute the ray's canonical momenta (covariant coordinate components of its four-momentum) with its conserved energy -Pt set to unity as a convenient convention. Then, compute axial angular momentum b and the Carter constant q (these are also conserved quantities) */
void set_ray_momenta(Ray *ray, KerrMetric *kerr) {
    
    // more localization
    double Fr = ray->Fr;
    double Ftheta = ray->Ftheta;
    double Fphi = ray->Fphi;
    double theta = ray->theta;
    
    double rho = kerr->rho;
    double delta = kerr->delta;
    double nu = kerr->nu;
    double alpha = kerr->alpha;
    double omega = kerr->omega;
    double a = kerr->a;
    
    // this is the energy measured by the FIDO, needed for measurements
    double E = 1 / (alpha + omega*nu*Fphi);
    
    double Pt = -1;
    double Pr = E*(rho/sqrt(delta))*Fr;
    double Ptheta = E*rho*Ftheta;
    double Pphi = E*nu*Fphi;
    
    // axial angular momentum b and Carter constant q
    // Note: if we had not set Pt = -1, then b = -Ptheta/Pt and q = (Carter constant)/(Pt*Pt)
    double b = Pphi;
    double q = Ptheta*Ptheta + cos(theta)*cos(theta)*((b*b)/(sin(theta)*sin(theta)) - a*a);
    
    // commit to ray
    ray->Pt = Pt;
    ray->Pr = Pr;
    ray->Ptheta = Ptheta;
    ray->Pphi = Pphi;
    
    ray->b = b;
    ray->q = q;
}




/* Set black hole parameters. */
void set_black_hole(BlackHole *BH, double r, double a) {
    BH->r = r;
    BH->a = a;
}


/* Express Kerr metric with supplied Boyer-Lindquist coordinates (don't actually use phi, and t coord excluded) */
void set_kerr_metric(KerrMetric *kerr, double a, double r, double theta, double phi) {
    
    // lets make things at least a little easier...
    double r2 = r*r;
    double a2 = a*a;
    double costheta = cos(theta);
    double sintheta = sin(theta);
    double cos2 = costheta*costheta;
    double sin2 = sintheta*sintheta;
    
    // calculate a bunch of lowercase greek letters
    kerr->rho = sqrt( r2 + a2*cos2 );
    kerr->delta = r2 - 2*r + a2;
    kerr->sigma = sqrt((r2 + a2)*(r2 + a2) - a2*kerr->delta*sin2);
    kerr->alpha = (kerr->rho*sqrt(kerr->delta))/kerr->sigma;
    kerr->omega = (2*a*r)/(kerr->sigma*kerr->sigma);
    kerr->nu = (kerr->sigma*sintheta)/kerr->rho;
    
    kerr->a = a;
    
}





