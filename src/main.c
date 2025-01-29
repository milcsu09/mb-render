#include <stdio.h>
#include "mb_viewport.h"
#include "mb_renderer.h"

int
main (void)
{
  mb_viewport vp;
  mb_viewport_init (&vp, -.75, 0, 1920, 1080);

  mb_renderer r;
  mb_renderer_init (&r, &vp);

  mb_renderer_main (&r);

  mb_renderer_quit (&r);

  return 0;
}

