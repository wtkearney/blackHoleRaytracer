//  color.h
//  
//
//  Created by Hunter Standen and Will Kearney on 9/22/14.
//

#ifndef _color_h
#define _color_h

typedef struct{
    float c[4];
} Color;

//copies the Color data
void color_copy(Color *to, Color *from);

//sets the Color data
void color_set(Color *to, float r, float g, float b, float a);

#endif
