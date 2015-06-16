//  color.c
//  
//
//  Created by Hunter Standen on 9/22/14.
//

#include "graphics.h"

void color_copy(Color *to, Color *from) {
    to->c[0]=from->c[0];
    to->c[1]=from->c[1];
    to->c[2]=from->c[2];
    to->c[3]=from->c[3];
}

//sets the Color data
void color_set(Color *to, float r, float g, float b, float a) {
    to->c[0]=r;
    to->c[1]=g;
    to->c[2]=b;
    to->c[3]=a;
}