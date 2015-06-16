//
//  Image.h
//  
//
//  Created by Will Kearney on 9/14/14.
//
//

#ifndef _image_h
#define _image_h

#include "color.h"

typedef struct {
    float rgb[3];
    float a;
    float z;
} FPixel;

typedef struct {
    /* pointer to space for storing Pixels */
    FPixel **data;
    
    /* number of rows and columns in the image; image size */
    int rows, cols, imagesize;
    
    /* char array to hold the filename of the image */
    char filename;
} Image;


/* function declarations */

// Constructors & destructors
Image *image_create(int rows, int cols);
void image_free(Image *src);
void image_init(Image *src);
int image_alloc(Image *src, int rows, int cols);
void image_dealloc(Image *src);

// I/O functions
Image *image_read(char *filename);
int image_write(Image *src, char *filename);
void fpixel_copy(FPixel *to, FPixel *from);

// Access
FPixel image_getf(Image *src, int c, int r);
float image_getc(Image *src, int c, int r, int b);
float image_geta(Image *src, int c, int r);
float image_getz(Image *src, int c, int r);
void image_setf(Image *src, int c, int r, FPixel val);
void image_setc(Image *src, int c, int r, int b, float val);
void image_seta(Image *src, int c, int r, float val);
void image_setz(Image *src, int c, int r, float val);

// Utility
void image_reset(Image *src);
void image_fill(Image *src, FPixel val);
void image_fillrgb(Image *src, float r, float g, float b);
void image_filla(Image *src, float a);
void image_setColor( Image *src, int c, int r, Color val);
Color image_getColor( Image *src, int c, int r, int level );

#endif
