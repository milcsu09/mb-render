#include "mb_viewport.h"
#include "mb_renderer.h"

int
main (void)
{
  srand(time(NULL));

  const int W = 1920;
  const int H = 1080;

  mb_viewport vp;
  mb_viewport_init (&vp, -.75, 0, W, H);

  // printf ("%f, %f\n", vp.rmin, vp.rmax);
  // printf ("%f, %f\n", vp.imin, vp.imax);

  mb_renderer r;
  mb_renderer_init (&r, &vp);

  mb_renderer_main (&r);

  mb_renderer_quit (&r);

  return 0;
}

