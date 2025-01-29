#ifndef MB_VIEWPORT_H
#define MB_VIEWPORT_H

#include "mb_types.h"

typedef struct
{
  f64 rmin, rmax;
  f64 imin, imax;
  i32 w, h;
} mb_viewport;

void mb_viewport_init (mb_viewport *vp, f64 r, f64 i, i32 w, i32 h);
void mb_viewport_zoom (mb_viewport *vp, f64 r, f64 i, f64 s);
void mb_viewport_to_plane (mb_viewport *vp, i32 x, i32 y, f64 *r, f64 *i);

#endif // MB_VIEWPORT_H

