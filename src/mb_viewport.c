#include "mb_viewport.h"

void
mb_viewport_init (mb_viewport *vp, f64 r, f64 i, i32 w, i32 h)
{
  const f64 zoom = 2.0;
  const f64 ratio = (double)w / (double)h;

  const f64 hw = zoom;
  const f64 hh = zoom / ratio;

  vp->rmin = r - hw;
  vp->rmax = r + hw;
  vp->imin = i - hh;
  vp->imax = i + hh;

  vp->w = w;
  vp->h = h;
}

void
mb_viewport_zoom (mb_viewport *vp, f64 r, f64 i, f64 s)
{
  const f64 dr = vp->rmax - vp->rmin;
  const f64 di = vp->imax - vp->imin;

  const f64 sr = dr * s;
  const f64 si = di * s;

  vp->rmin = r - (r - vp->rmin) / dr * sr;
  vp->rmax = vp->rmin + sr;

  vp->imin = i - (i - vp->imin) / di * si;
  vp->imax = vp->imin + si;
}

void
mb_viewport_to_plane (mb_viewport *vp, i32 x, i32 y, f64 *r, f64 *i)
{
  const f64 dr = vp->rmax - vp->rmin;
  const f64 di = vp->imax - vp->imin;
  *r = vp->rmin + (x / (f64)vp->w) * dr;
  *i = vp->imin + (y / (f64)vp->h) * di;
}


