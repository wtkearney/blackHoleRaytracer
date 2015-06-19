# [Black Hole Ray Tracing in C](https://github.com/wtkearney/blackHoleRaytracer)

## Introduction

This began as the final project for a Advanced Computer Graphics class (Prof. Bruce Maxwell) in the Spring of 2015, but has since become a personal project as I continue to develop it. The following is a write-up describing the process and implementation; obviously, this is a work in progress, and thus there is no doubt that errors could still abound. The information herein should also not be taken as fact, as my physics knowledge is limited at best. Feel free to fork at your own discretion...!


## Gravitational Lensing
Gravitational lensing, on a simple level, occurs when light appears to bend as it travels through curved spacetime. Obvious examples from pop-media include the warping effect theoretically seen around black holes and high-mass neutron stars. Gravitational lensing, however, happens even as locally as around our own sun (though by a very small amount).

The geometry of spacetime around an object of extreme mass is defined mathematically by a metric tensor. A black hole can exist in four states, defined by two parameters: charge and rotation. The Schwarzschild metric describes the geometry around a Schwarzschild black hole, which carries zero charge, is spherically symmetrical, and is non-rotating. Similarly, the Kerr metric describes the geometry around a Kerr black hole, which is uncharged, has axial symmetry, and does have rotation. Obviously, these metrics are generalizations of each other; substituting zero as the angular momentum parameter in the Kerr metric will yield the Schwarzschild equations. These metrics are all highly non-linear, which makes exact solutions difficult. Numerical solving techniques must be used instead. The Kerr metric, written in Boyer-Lindquist coordinates, is:

These equations are used to define a photon’s geodesic, or world line, which is what we are interesting in. Particularly, we are concerned with a photon’s null geodesic, which is generally the path a massless particle (such as a photon) will take through curved spacetime. These geodesics are not time-like, as would be the case for objects of mass. When we talk about light “bending,” then, it is important to note that this isn’t particularly accurate language. Technically, the photon is traveling in a straight light as generalized to curved spacetime.

## Ray Tracing
Ray tracing is actually a useful context to understand how light behaves in a Kerr metric. Ordinarily, we would back-propagate light rays/photons from the camera in flat spacetime, and thus straight lines. Now, however, we must back-propagate light rays/photons from the camera in curved spacetime, using equations of motion derived from the Kerr metric to inform the photon’s world line.

There are a few implementations of ray tracing a black hole available online (see Antonelli), mostly people’s personal projects. Many of them use some clever algebra to reduce the equations of motion of a null geodesic around a Schwarzschild black hole into more manageable form. I got a bit lost when trying to follow these mathematical reductions and accompanying code, and thus decided to implement a ray tracer that wasn’t trying to be so clever. I thought it better to try and tackle something I better understood, and although the equations might be more complicated, I understood the mechanics much better and thus presumably would be able to debug code more effectively. After a bit of research, I settled on a paper written by the team at Double Negative, the team who did the visual effects for Christopher Nolan’s Interstellar. Their paper certainly went into great detail, but was also most accessible in terms of laying out the general concepts and accompanying equations.

## Gravitational Renderer
The process of rendering what one might see in the vicinity of a black hole is as follows:

1) Define parameters for the black hole (t, r, θ, φ). Create a locally non-rotating observer, or fiducial observer. The FIDO is at-rest in space, and has orthonormal basis vectors pointing along the spatial coordinate lines. Conceptually, the FIDO is what someone far away from the gravitational field would interpret as at rest, and does not move according to the spatial axis around the black hole.

2) Define parameters for the camera in orbit around the black hole: location (rc, θc, φc), direction of motion relative to the FIDO (a unit vector B), and speed (ß).

3) Create basis vectors in the camera’s reference frame (ex, ey, ez). ey is specified as the camera’s direction of motion B. A spherical polar system can then be set up around the camera to define the camera’s local sky. This is how we define the directions of incoming light rays.

4) Numerically integrate the null geodesic equation to propagate a specific light ray from the camera’s local sky to the celestial sphere. There are two possibilities: the ray hits the celestial sphere, as which point we can pull information from a texture map corresponding to the celestial sphere, or the ray falls below the proton sphere of the black hole, at which point we simply color that pixel black. (A third option exists if one is implementing an accretion disk, at which point a pixel can be pulled from a texture map of the accretion disk; we ignore this for now, along with phenomenon such as the Doppler effect and gravitational redshift. So did Christopher Nolan, so that’s totally legitimate.)

## Implementation
This renderer was built from scratch in C. The image data structure and functions as well as the color structure and functions were developed in previous projects and fine-tuned throughout the academic year; the ppm read/write library was provided by Bruce Maxwell in the fall of 2014 during the first semester of a computer graphics course.

At the core, the renderer is based on four data structures, corresponding to the camera, the black hole, a light ray, and the Kerr metric. The code logic is mostly easily followed from the source file. The black hole is initialized, and various Kerr metric parameters are solved for. Next, the camera parameters are specified and solved for based on the Kerr metric. It would be good to note that if we put our camera in a circular, equatorial geodesic orbit around the black hole, camera parameters tend to simplify nicely. The direction of motion (relative to the FIDO at that point) become Br = 0 Bθ = 0, Bφ = 1. The velocity ß also becomes easier to define. Next, the ray tracer method is called, taking the camera, the black hole, and the Kerr metric as arguments.

From here, a ray’s location, normal, direction, and momenta are calculated. A series of tests check the number of radial turning points the ray had, and then correspondingly if the ray originated from the celestial sphere or the horizon. If the ray originated from the celestial sphere, the point of origin is determined by numerically integrating the ray equations backwards in time.

## Results

First off, here is an example output for a single ray I am testing; the ray being tested is being back propagated along the normal of the camera plane, so it should land directly into the black hole.

    arcLengthHoriz: 0.523599
    arcLengthVert: 0.392699
    Raytracing...

    Ray info:
    Theta, Phi: (1.570796, 3.141593)
    Nx, Ny, Nz: (-1.000000, 0.000000, 0.000000)
    Fx, Fy, Fz: (0.942466, 0.334301, -0.000000)
    Fr, Ftheta, Fphi: (0.942466, 0.000000, 0.334301)
    Pt, Pr, Ptheta, Pphi: (-1.000000, 1.170099, 0.000000, 3.728645)
    Axial angular momentum (b): 3.728645
    Carter constant (q): 0.000000

    Radius interval of unstably trapped photons [r1, r2]: 2.188914, 3.629850
    b1, b2: -6.315649, 3.838494
    r0: 2.846423
    q0(r0): 11.812165
    Black Hole info:
    r: 2.000000
    a: 0.600000
    Camera info:
    Location: (10.000000, 1.570796, 0.000000)
    Speed: 0.334301
    Motion components: (0.000000, 0.000000, 1.000000)
    Kerr Metric info:
    rho: 10.000000
    delta: 80.360000
    sigma: 100.215767
    alpha: 0.894507
    omega: 0.001195
    nu: 10.021577

- arcLengthHoriz and arcLengthVert correspond to the viewing window size, or angle of the camera lens; 0.52 radians is ~30°. Because the viewing window is specified as 800x600, this means the vertical length is correspondingly 0.39 radians, or ~22.5°.

Ray Info

- theta and phi refer to the ray's location on the camera's local sky; (p/2, p) ≈ (1.5, 3.2), which means the ray is pointing directly into the black hole, or along the negative x-axis of the camera's proper reference frame (remember that our camera is in a nice relaxing circular, equatorial geodesic orbit).

- Nx, Ny, and Nz represent the normal vector in the camera's reference plane that points in the direction of the incoming ray. We can see that the vector points along the negative x-axis, as suspected. This is good.

- Fx, Fy, and Fz represent the ray's direction of motion, as measured by the FIDO in the camera's reference frame. Large positive x makes sense; we are back-propagating the ray, so it is traveling away from the camera. The positive y makes sense because we are in orbit moving along the positive y-direction = B.

- Fr, Ftheta and Fphi are the direction of motion on the FIDO's spherical orthonormal basis.

- Pt, Pr, Ptheta, and Pphi represent the ray's canonical momenta. Note that -Pt is set to unity, a convenient convention made possible by the fact that we are working with massless particles in relativistic contexts.

- b and q represent the ray's axial momentum and its Carter constant. These are constants of motion, terms that finally allow us to begin to determine if the ray originated from the celestial sphere or the horizon.

- r1 and r2 represent the radius bounds for a photon unstably trapped in a constant orbit around the black hole. We notice that the lower bound is just outside the radius of the actual black hole (r = 2, in this case).

- b1 is a function of R2 and b2 is a function of R1. b1 < b < b2 is one of the tests to determine the ray's point of origin.

- r0 is calculated from the angular momentum of the photon, and is then used to calculate the theoretical Carter constant at that radius (b and q can be defined parametrically with r).

- q0(b) can then be calculated once we have r0. q < q0(b) is another test to determine the ray's point of origin.

Black Hole Info

- r is the radius (or radius of the event horizon) of the black hole, and a is the angular momentum. If a = 0, we would of course have a Schwarzschild (non-rotating) black hole. The radius can tell us interesting things, such as where the proton sphere exists (1.5 times the Schwarzschild radius).

Camera Info

- The camera is located at (rc, θc, φc) = (10.0, 1.570796, 0.0) ≈ (10, p/2, 0.0), putting us nicely on the equator.
- As expected, our motion components (0.0, 0.0, 1.0) and speed indicate we are traveling in a circular, equatorial orbit. Speed is determined from location along with parameters from the Kerr metric.

Kerr Metric Info

- These are a bunch of parameters calculated from the Kerr Metric; they're used in many of the calculations, and thus useful to keep nearby.

Testing for Ray Origin

- The ray's constant's (b, q) are used to determine the origin of the ray with the following process...

(a) If b1 < b < b2 and q < q0(b), then there are no radial turning points. Therefore, if Pr > 0 (the ray's canonical momenta) at the camera's location then the ray comes from the horizon; otherwise, the ray comes from the celestial sphere.
(b) If condition (a) does not hold, then there are two radial turning points. If the camera radius is greater than the radius of the upper turning point, then the ray comes from the celestial sphere; otherwise, the ray comes from the horizon.

The radius of the upper turning point is the largest real root of R(r) = 0.

Compute Point of Origin

- This step is performed if the ray is confirmed as originating from the celestial sphere. The point of origin is computed by numerically integrating the ray equations, starting with the computed constants of motion, the ray momenta, and the camera location. I am planning on using a Runge-Kutta solver to do the numeric integration, but this section is still somewhat pending...


## Citations

Oliver James, Eugénie von Tunzelmann, Paul Franklin, and Kip S Thorne. "Gravitational lensing by spinning black holes in astrophysics, and in the movie Interstellar" Classical and Quantum Gravity 32.6 (2015).http://iopscience.iop.org/0264-9381/32/6/065001

Antonelli, Riccardo. http://rantonels.github.io/starless/

