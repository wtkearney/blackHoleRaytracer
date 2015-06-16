//
//  Image.c
//  
//
//  Updated by Will Kearney on 3/10/15
//
//


#include "graphics.h"


// Constructors & destructors

/* Allocates an Image structure and initializes the top level fields to appropriate values. Allocates space for an image of the specified size, unless either rows or cols is 0. Returns a pointer to the allocated Image structure. Returns a NULL pointer if the operation fails. */
Image *image_create(int rows, int cols) {
    Image *image = NULL;
    image = (Image *)malloc( sizeof(Image) );
    
    if (rows > 0 && cols > 0) {
        
        image->data = (FPixel **)malloc( sizeof(FPixel *)*rows );
        
        /* create rows */
        for (int i=0; i<rows;i++) {
            image->data[i] = (FPixel *)malloc(sizeof(FPixel)*cols);
        }

        
    } else {
        printf("Rows and columns must both be greater than 0.\n");
        image->rows = 0;
        image->cols = 0;
        image->data = NULL;
        return(image);
    }
    
    image->rows = rows;
    image->cols = cols;
    image->imagesize = rows * cols;
    image_reset(image);
    
    return(image);
}

/* De-allocates image data and frees the Image structure. */
void image_free(Image *src){
    
    if (src != NULL) {
        if (src->data != NULL) {
            
            /* free cols in each row */
            for (int i=0; i<src->rows;i++) {
                free( src->data[i] );
            }
            
            /* free rows */
            free(src->data);
        }
        free(src);
    }
}

/* Given an uninitialized Image structure, sets the rows and cols fields to zero and the data field to NULL. */
void image_init(Image *src){
    src->rows = 0;
    src->cols = 0;
    src->imagesize = 0;
    src->data = NULL;
}

/* Allocates space for the image data given rows and columns and initializes the image data to appropriate values, such as 0.0 for RGBA and 1.0 for Z. Returns 0 if the operation is successful. Returns a non-zero value if the operation fails. This function should probably free existing memory if rows and cols are both non-zero. */
int image_alloc(Image *src, int rows, int cols) {
    
    if (src->data != NULL) {
        
        /* free cols in each row */
        for (int i=0; i<src->rows;i++) {
            free( src->data[i] );
        }
        
        /* free rows */
        free(src->data);
    }
    
    
    
    if (rows > 0 && cols > 0) {
        src->rows = rows;
        src->cols = cols;
        
        src->data = (FPixel **)malloc( sizeof(FPixel *)*rows );
        
        /* create rows */
        for (int i=0; i<rows;i++) {
            src->data[i] = (FPixel *)malloc(sizeof(FPixel)*cols);
        }
        
        src->imagesize = rows*cols;
        
        image_reset(src);
        
        return(0);
    } else {
        printf("Rows and columns must both be greater than 0.\n");
        return(1);
    }
}

/* De-allocates image data and resets the Image structure fields. The function does not free the Image structure. */
void image_dealloc(Image *src){
    
    if (src->data != NULL) {
        
        /* free cols in each row */
        for (int i=0; i<src->rows;i++) {
            free( src->data[i] );
        }
        
        /* free rows */
        free(src->data);
    }

    src->data = NULL;
    src->rows = 0;
    src->cols = 0;
    src->imagesize = 0;
}



// I/O functions

/* reads a PPM image from the given filename. Initializes the alpha channel to 1.0 and the z channel to 1.0. Returns a NULL pointer if the operation fails. */
Image *image_read(char *filename) {
    int rows = 0;
    int cols = 0;
    int colors = 0;
    Pixel *image_array = readPPM(&rows, &cols, &colors, filename); // an array of pixels
    
    if (rows < 0 || cols < 0 ) {
        printf("rows: %d, cols: %d", rows, cols);
        return(NULL);
    }
    
    Image *image = image_create(rows, cols);
    for (int row=0;row<rows;row++) {
        for (int col=0;col<cols;col++) {

            image->data[row][col].rgb[0] = 256 - image_array[row*cols + col].r;
            image->data[row][col].rgb[1] = 256 - image_array[row*cols + col].g;
            image->data[row][col].rgb[2] = 256 - image_array[row*cols + col].b;
            image->data[row][col].a = 1.0;
            image->data[row][col].z = 1.0;
        }
    }
    return(image);
}

/* writes a PPM image to the given filename. Returns 0 on success. Optionally, you can look at the filename extension and write different file types. */
int image_write(Image *src, char *filename) {
    printf("Size: %d\n", src->imagesize);
    Pixel tempImage[ src->imagesize ]; // temp array to hold Pixels
    for (int row=0;row<src->rows;row++) {
        for (int col=0;col<src->cols;col++) {
            
            //printf("row: %d, col: %d\n", row, col);
            
            tempImage[row*src->cols + col].r = src->data[row][col].rgb[0] * 255;
            tempImage[row*src->cols + col].g = src->data[row][col].rgb[1] * 255;
            tempImage[row*src->cols + col].b = src->data[row][col].rgb[2] * 255;
            
        }
        
    }

    writePPM(tempImage, src->rows, src->cols, 255, filename);
    //free(tempImage);

    return(0);
}


//copies the FPixel data structure
void fpixel_copy(FPixel *to, FPixel *from) {
    to->rgb[0]=from->rgb[0];
    to->rgb[1]=from->rgb[1];
    to->rgb[2]=from->rgb[2];
    to->a=from->a;
    to->z=from->z;
}


// Access

/* returns the FPixel at (r, c). */
FPixel image_getf(Image *src, int r, int c) {
    return( src->data[r][c] );
}

/* returns the value of band b at pixel (r, c). */
float image_getc(Image *src, int r, int c, int b) {
    return( src->data[r][c].rgb[b] );
}

/* returns the alpha value at pixel (r, c). */
float image_geta(Image *src, int r, int c) {
    return( src->data[r][c].a );
}

/* returns the depth value at pixel (r, c). */
float image_getz(Image *src, int r, int c) {
    return( src->data[r][c].z );
}

/* sets the values of pixel (r, c) to the FPixel val. */
void image_setf(Image *src, int r, int c, FPixel val) {
    if (r < 0 || c < 0 || r >= src->rows || c >= src->cols) {
        return;
    } else {
        src->data[r][c] = val;
    }
}

/* sets the value of pixel (r, c) band b to val. */
void image_setc(Image *src, int r, int c, int b, float val) {
//    printf("Row: %d\n", r);
//    printf("Col: %d\n", c);
    if (val > 1) {
        val = 1.0;
    }
    if (r < 0 || c < 0 || r >= src->rows || c >= src->cols) {
        return;
    } else {
        src->data[r][c].rgb[b] = val;
    }
//    printf("Complete\n\n");

}

/* sets the alpha value of pixel (r, c) to val. */
void image_seta(Image *src, int r, int c, float val) {
    if (r < 0 || c < 0 || r >= src->rows || c >= src->cols) {
        return;
    } else {
        src->data[r][c].a = val;
    }
}

/* sets the depth value of pixel (r, c) to val. */
void image_setz(Image *src, int r, int c, float val) {
    if (r < 0 || c < 0 || r >= src->rows || c >= src->cols) {
        return;
    } else {
        src->data[r][c].z = val;
    }
}



// Utility

/* resets every pixel to a default value (e.g. Black, alpha value of 1.0, z value of 1.0). */
void image_reset(Image *src) {
    if (src == NULL) {
        printf("Here's your problem.\n");
    }
    
    for (int r=0;r<src->rows;r++) {
        for (int c=0;c<src->cols;c++) {
            src->data[r][c].rgb[0] = 0;
            src->data[r][c].rgb[1] = 0;
            src->data[r][c].rgb[2] = 0;
            src->data[r][c].a = 1.0;
            src->data[r][c].z = 1.0;
        }
    }

}

/* sets every FPixel to the given value. */
void image_fill(Image *src, FPixel val) {
    for (int r=0;r<src->rows;r++) {
        for (int c=0;c<src->cols;c++) {
            src->data[r][c] = val;
        }
    }
}

/* sets the (r, g, b) values of each pixel to the given color. */
void image_fillrgb(Image *src, float r, float g, float b) {
    for (int r=0;r<src->rows;r++) {
        for (int c=0;c<src->cols;c++) {
            src->data[r][c].rgb[0] = r;
            src->data[r][c].rgb[1] = g;
            src->data[r][c].rgb[2] = b;
        }
    }
}

/* sets the alpha value of each pixel to the given value. */
void image_filla(Image *src, float a) {
    for (int r=0;r<src->rows;r++) {
        for (int c=0;c<src->cols;c++) {
            src->data[r][c].a = 1.0;
        }
    }
}


/* copies the Color data to the proper pixel. */
void image_setColor( Image *src, int r, int c, Color val) {
    if (r < 0 || c < 0 || r >= src->rows || c >= src->cols || val.c[0] < 0 || val.c[1] < 0 || val.c[2] < 0) {
        return;
    } else {

        src->data[r][c].rgb[0] = val.c[0];
        src->data[r][c].rgb[1] = val.c[1];
        src->data[r][c].rgb[2] = val.c[2];
    }
}



/* returns a Color structure built from the pixel values. */
Color image_getColor( Image *src, int r, int c, int level ) {
    Color color;
    color.c[0] = src->data[r][c].rgb[0];
    color.c[1] = src->data[r][c].rgb[1];
    color.c[2] = src->data[r][c].rgb[2];
    return color;
}



