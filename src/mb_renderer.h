#ifndef MB_RENDERER_H
#define MB_RENDERER_H

#include <SDL2/SDL.h>
#include <stdbool.h>

#include <pthread.h>
#include "mb_viewport.h"

typedef struct
{
  SDL_Window *window;
  SDL_Renderer *renderer;
  SDL_Texture *texture;
  mb_viewport *vp;
  i32 max_iter;
  bool active;
  bool recompute;

  bool dragging;
  struct
  {
    i32 sx, sy;
    f64 dx, dy;
  } drag;

  bool zooming;
  struct
  {
    i32 mx, my;
    f64 delta;
    f64 value;
  } zoom;
} mb_renderer;

void mb_renderer_init (mb_renderer *r, mb_viewport *vp);
void mb_renderer_quit (mb_renderer *r);

void mb_renderer_main (mb_renderer *r);

#endif // MB_RENDERER_H

