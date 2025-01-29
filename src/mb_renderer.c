#include "mb_renderer.h"
#include "mb_types.h"
#include <omp.h>

static void mb_renderer_poll (mb_renderer *r);
static void mb_renderer_recompute (mb_renderer *r);

void
mb_renderer_init (mb_renderer *r, mb_viewport *vp)
{
  SDL_Init (SDL_INIT_VIDEO);

  r->window = SDL_CreateWindow ("Mandelbrot-Render", SDL_WINDOWPOS_CENTERED,
                                SDL_WINDOWPOS_CENTERED, vp->w, vp->h,
                                SDL_WINDOW_SHOWN);
  r->renderer = SDL_CreateRenderer (r->window, -1, SDL_RENDERER_ACCELERATED);
  r->texture = SDL_CreateTexture (r->renderer, SDL_PIXELFORMAT_RGB888,
                                  SDL_TEXTUREACCESS_STREAMING, vp->w, vp->h);

  r->vp = vp;

  r->max_iter = 5000;
  r->active = true;
  r->recompute = true;

  r->dragging = false;
  r->drag.sx = r->drag.sy = 0;
  r->drag.dx = r->drag.dy = 0.0;

  r->zooming = false;
  r->zoom.mx = r->zoom.my = 0;
  r->zoom.delta = 1.0;
  r->zoom.value = 1.0;
}

void
mb_renderer_quit (mb_renderer *r)
{
  SDL_DestroyTexture (r->texture);
  SDL_DestroyRenderer (r->renderer);
  SDL_DestroyWindow (r->window);
  SDL_Quit ();
}

void
mb_renderer_main (mb_renderer *r)
{
  while (r->active)
    {
      mb_renderer_poll (r);

      if (r->zooming)
        {
          f64 mr, mi;
          mb_viewport_to_plane (r->vp, r->zoom.mx, r->zoom.my, &mr, &mi);
          mb_viewport_zoom (r->vp, mr, mi, r->zoom.delta);
          r->zoom.value /= r->zoom.delta;
          SDL_Delay (16);
        }

      if (r->recompute)
        {
          // r->max_iter = fmax (75 * log2(1.0 / (r->vp->rmax - r->vp->rmin)), 50);
          const u32 s = SDL_GetTicks ();
          mb_renderer_recompute (r);
          const u32 e = SDL_GetTicks ();

          printf ("(Recompute) %dms\n", e - s);
        }

      SDL_SetRenderDrawColor (r->renderer, 0, 0, 0, 255);
      SDL_RenderClear (r->renderer);

      SDL_Rect vr = {
        .x = r->zoom.mx - r->zoom.mx * r->zoom.value + r->drag.dx,
        .y = r->zoom.my - r->zoom.my * r->zoom.value + r->drag.dy,
        .w = r->vp->w * r->zoom.value,
        .h = r->vp->h * r->zoom.value,
      };

      SDL_RenderCopy (r->renderer, r->texture, NULL, &vr);
      SDL_RenderPresent(r->renderer);
    }
}

static void
mb_renderer_poll (mb_renderer *r)
{
  SDL_Event event;
  while (SDL_PollEvent (&event))
    switch (event.type)
      {
      case SDL_QUIT:
        r->active = false;
        break;
      case SDL_MOUSEBUTTONDOWN:
        if (event.button.button == SDL_BUTTON_MIDDLE)
          {
            r->dragging = true;
            r->drag.sx = event.button.x;
            r->drag.sy = event.button.y;
          }
        else if (event.button.button == SDL_BUTTON_LEFT)
          {
            r->zooming = true;
            r->zoom.delta = 0.9;
            r->zoom.mx = event.button.x;
            r->zoom.my = event.button.y;
            // SDL_GetMouseState (&r->zoom.mx, &r->zoom.my);
          }
        else if (event.button.button == SDL_BUTTON_RIGHT)
          {
            r->zooming = true;
            r->zoom.delta = 1.1;
            r->zoom.mx = event.button.x;
            r->zoom.my = event.button.y;
            // SDL_GetMouseState (&r->zoom.mx, &r->zoom.my);
          }
        break;
      case SDL_MOUSEBUTTONUP:
        if (event.button.button == SDL_BUTTON_MIDDLE)
          {
            mb_viewport *vp = r->vp;

            r->dragging = false;
            const f64 dx = r->drag.dx * (vp->rmax - vp->rmin) / vp->w;
            const f64 dy = r->drag.dy * (vp->imax - vp->imin) / vp->h;

            vp->rmin -= dx;
            vp->rmax -= dx;
            vp->imin -= dy;
            vp->imax -= dy;

            r->drag.dx = r->drag.dy = 0.0;
            r->recompute = true;
          }
        else if (event.button.button == SDL_BUTTON_LEFT)
          {
            r->zooming = false;
            r->zoom.value = 1.0;
            r->recompute = true;
          }
        else if (event.button.button == SDL_BUTTON_RIGHT)
          {
            r->zooming = false;
            r->zoom.value = 1.0;
            r->recompute = true;
          }
        break;
      case SDL_MOUSEMOTION:
        if (r->dragging)
          {
            r->drag.dx = event.motion.x - r->drag.sx;
            r->drag.dy = event.motion.y - r->drag.sy;
          }
        break;
      default:
        break;
      }
}

static void
mb_renderer_recompute (mb_renderer *r)
{
  u32 *pixels;
  int pitch;

  SDL_LockTexture (r->texture, NULL, (void **)&pixels, &pitch);

  SDL_PixelFormat *fmt = SDL_AllocFormat (SDL_PIXELFORMAT_RGB888);

  const mb_viewport *vp = r->vp;

  const u32 w = vp->w;
  const u32 h = vp->h;

  const f64 dr = vp->rmax - vp->rmin;
  const f64 di = vp->imax - vp->imin;

#pragma omp parallel for collapse(2) schedule(dynamic)
  for (u32 y = 0; y < h; ++y)
    for (u32 x = 0; x < w; ++x)
      {
        const f64 cr = vp->rmin + (x / (f64)w) * dr;
        const f64 ci = vp->imin + (y / (f64)h) * di;

        i32 i;
        f64 zr = cr;
        f64 zi = ci;

        for (i = 0; i < r->max_iter; ++i)
          {
            const f64 zr_sq = zr * zr;
            const f64 zi_sq = zi * zi;

            if (zr_sq + zi_sq > 4.0)
              break;

            const f64 t = zr_sq - zi_sq + cr;

            zi = 2.0 * zr * zi + ci;
            zr = t;
          }

        u32 color;

        if (i == r->max_iter)
          color = 0;
        else
          {
            const u8 r = (i * 9) % 256;
            const u8 g = (i * 8) % 256;
            const u8 b = (i * 6) % 256;
            color = SDL_MapRGB (fmt, r, g, b);
          }

        pixels[y * (pitch / 4) + x] = color;
      }

  SDL_FreeFormat (fmt);
  SDL_UnlockTexture (r->texture);

  r->recompute = false;
}


