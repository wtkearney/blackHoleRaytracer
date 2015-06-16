//
//  proj4_scene1.c
//
//
//  Created by Will Kearney
//  Spring 2015
//

#include <stdio.h>
#include <time.h>
#include "graphics.h"



static const int HEIGHT = 600;
static const int WIDTH = 800;

static const double BLACK_HOLE_SPIN = 0.9;
static const double BLACK_HOLE_RADIUS = 5;


// camera specifications
static const double CAM_LOC_R = 20;
static const double CAM_LOC_THETA = M_PI/2;
static const double CAM_LOC_PHI = M_PI/3;

static const double CAM_DIR_R = 0;
static const double CAM_DIR_THETA = 0;
static const double CAM_DIR_PHI = 1;

static const double LENS_ANGLE = M_PI/2;


int main(void) {
    
    clock_t start, end;
    double cpu_time_used;
    
    start = clock();
    
    Image *src = NULL;

    Camera *camera = NULL;
    BlackHole *blackHole = NULL;
    KerrMetric *kerr = NULL;

    
    src = image_create(HEIGHT, WIDTH);
    camera = create_camera();
    blackHole = create_black_hole();
    kerr = create_kerr_metric();
    
    // set parameters of black hole
    set_black_hole(blackHole, BLACK_HOLE_RADIUS, BLACK_HOLE_SPIN);
    
    // express Kerr Metric with Boyer-Lindquist coordinates
    set_kerr_metric(kerr, blackHole->a, CAM_LOC_R, CAM_LOC_THETA, CAM_LOC_PHI);
    
    // Specify camera parameters for a nice relaxing circular, equatorial geodesic orbit looking at black hole
    set_camera_parameters(camera, kerr, CAM_LOC_R, CAM_LOC_THETA, CAM_LOC_PHI, CAM_DIR_R, CAM_DIR_THETA, CAM_DIR_PHI, LENS_ANGLE, 1);
    
    // Perform the actual raytracing
    raytrace(src, camera, blackHole, kerr);
    
    
    // Diagnostics
    printf("Black Hole info:\n");
    printf("    r: %f\n", blackHole->r);
    printf("    a: %f\n", blackHole->a);
           
    printf("Camera info:\n");
    printf("    Location: (%f, %f, %f)\n", camera->location[0], camera->location[1], camera->location[2]);
    printf("    Speed: %f\n", camera->beta);
    printf("    Motion components: (%f, %f, %f)\n", camera->motion[0], camera->motion[1], camera->motion[2]);
    
    printf("Kerr Metric info:\n");
    printf("    rho: %f\n", kerr->rho);
    printf("    delta: %f\n", kerr->delta);
    printf("    sigma: %f\n", kerr->sigma);
    printf("    alpha: %f\n", kerr->alpha);
    printf("    omega: %f\n", kerr->omega);
    printf("    nu: %f\n\n", kerr->nu);
    
    // write out image and convert
    image_write(src, "../images/blackhole.ppm");
    
    // it good to take care of memory
    image_free(src);
    
    free(camera);
    free(blackHole);
    free(kerr);

    
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Time: %f\n", cpu_time_used);
    
    
    printf("Converting and opening image...\n");
    system("convert ../images/blackhole.ppm ../images/blackhole.png");
    system("open ../images/blackhole.png");
    
    
    
}





