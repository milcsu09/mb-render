#include "mb_renderer.h"
#include "mb_types.h"
#include <omp.h>

#define MB_USE_SIMD 1

static void mb_renderer_poll (mb_renderer *r);
static void mb_renderer_recompute (mb_renderer *r);

void
mb_renderer_init (mb_renderer *r, mb_viewport *vp)
{
  SDL_Init (SDL_INIT_VIDEO);

  r->window = SDL_CreateWindow ("Mandelbrot-Render", SDL_WINDOWPOS_CENTERED,
                                SDL_WINDOWPOS_CENTERED, vp->w, vp->h,
                                SDL_WINDOW_SHOWN);
  r->renderer = SDL_CreateRenderer (
      r->window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
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
  u32 max_time = 0;
  while (r->active)
    {
      mb_renderer_poll (r);

      if (r->zooming)
        {
          f64 mr, mi;
          mb_viewport_to_plane (r->vp, r->zoom.mx, r->zoom.my, &mr, &mi);
          mb_viewport_zoom (r->vp, mr, mi, r->zoom.delta);
          r->zoom.value /= r->zoom.delta;
        }

      if (r->recompute)
        {
          const u32 s = SDL_GetTicks ();
          mb_renderer_recompute (r);
          const u32 e = SDL_GetTicks ();
          const u32 t = e - s;

          if (t > max_time)
            max_time = t;

          printf ("(Recompute) %dms\n", t);
        }

      SDL_SetRenderDrawColor (r->renderer, 32, 24, 24, 255);
      SDL_RenderClear (r->renderer);

      SDL_Rect vr = {
        .x = r->zoom.mx - r->zoom.mx * r->zoom.value + r->drag.dx,
        .y = r->zoom.my - r->zoom.my * r->zoom.value + r->drag.dy,
        .w = r->vp->w * r->zoom.value,
        .h = r->vp->h * r->zoom.value,
      };

      SDL_RenderCopy (r->renderer, r->texture, NULL, &vr);
      SDL_RenderPresent (r->renderer);
    }

  printf ("(Maximum) %dms\n", max_time);
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
            r->zoom.delta = 0.95;
            r->zoom.mx = event.button.x;
            r->zoom.my = event.button.y;
          }
        else if (event.button.button == SDL_BUTTON_RIGHT)
          {
            r->zooming = true;
            r->zoom.delta = 1.05;
            r->zoom.mx = event.button.x;
            r->zoom.my = event.button.y;
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

#if !MB_USE_SIMD
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

#pragma omp parallel for schedule(dynamic)
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

#else

void
mb_renderer_recompute_iteration (f64 *real, f64 *imag, u32 *results, u32 nmax)
{
  __m256d zr = _mm256_setzero_pd ();
  __m256d zi = _mm256_setzero_pd ();
  __m256d cr = _mm256_loadu_pd (real);
  __m256d ci = _mm256_loadu_pd (imag);
  __m256d iter = _mm256_setzero_pd ();
  __m256d four = _mm256_set1_pd (4.0);
  __m256d one = _mm256_set1_pd (1.0);

  for (int i = 0; i < nmax; ++i)
    {
      __m256d zr2 = _mm256_mul_pd (zr, zr);
      __m256d zi2 = _mm256_mul_pd (zi, zi);
      __m256d mag2 = _mm256_add_pd (zr2, zi2);

      __m256d mask = _mm256_cmp_pd (mag2, four, _CMP_LT_OQ);
      if (_mm256_movemask_pd (mask) == 0)
        break;

      __m256d zrzi = _mm256_mul_pd (zr, zi);
      zr = _mm256_add_pd (_mm256_sub_pd (zr2, zi2), cr);
      zi = _mm256_add_pd (_mm256_add_pd (zrzi, zrzi), ci);

      iter = _mm256_add_pd (iter, _mm256_and_pd (mask, one));
    }

  __m128i int_iter = _mm256_cvttpd_epi32 (iter);
  _mm_storeu_si128 ((__m128i *)results, int_iter);
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

#pragma omp parallel for schedule(dynamic)
  for (u32 y = 0; y < h; ++y)
    for (u32 x = 0; x < w; x += 4)
      {
        f64 real[4], imag[4];
        for (u32 i = 0; i < 4; ++i)
          {
            const f64 cr = vp->rmin + ((x + i) / (f64)w) * dr;
            const f64 ci = vp->imin + (y / (f64)h) * di;
            real[i] = cr;
            imag[i] = ci;
          }

        u32 results[4];
        mb_renderer_recompute_iteration (real, imag, results, r->max_iter);

        for (int j = 0; j < 4; ++j)
          {
            int i = results[j];

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
            pixels[y * (pitch / 4) + (x + j)] = color;
          }
      }

  SDL_FreeFormat (fmt);
  SDL_UnlockTexture (r->texture);

  r->recompute = false;
}
#endif // MB_USE_SIMD

