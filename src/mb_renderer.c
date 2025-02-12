#include "mb_renderer.h"
#include "mb_types.h"
#include <omp.h>

/*
  0.25893945454588274, 0.00138771724246275
*/

#define MB_ZOOM_DELTA .015

static void mb_renderer_poll (mb_renderer *r);
static void mb_renderer_recompute (mb_renderer *r);
static void mb_renderer_recompute_aa (mb_renderer *r);

static inline void
mb_util_format_number (f64 num, char *buffer, size_t buffer_size)
{
  static const char *const suffixes[] = { "", "K", "M", "B", "T" };
  int suffix_index = 0;

  while (num >= 1000.0 && suffix_index < 4)
    {
      num /= 1000.0;
      suffix_index++;
    }

  snprintf (buffer, buffer_size, "%.1f%s", num, suffixes[suffix_index]);
}

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

  r->max_iter = 15000;
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
  u32 all_time = 0;
  u64 samples = 0;

  i32 pmx, pmy;
  SDL_GetMouseState (&pmx, &pmy);

  while (r->active)
    {
      mb_renderer_poll (r);

      const u8 *keys = SDL_GetKeyboardState (NULL);

      i32 mx, my;
      SDL_GetMouseState (&mx, &my);

      if (keys[SDL_SCANCODE_A])
        {
          mb_viewport_zoom (r->vp, 0.25893945454588274, 0.00138771724246275, .1);
          r->recompute = true;
        }

      if (keys[SDL_SCANCODE_X])
        {
          u32 *pixels;
          int pitch;

          SDL_LockTexture (r->texture, NULL, (void **)&pixels, &pitch);

          SDL_PixelFormat *fmt = SDL_AllocFormat (SDL_PIXELFORMAT_RGB888);

          const i32 dx = pmx - mx;
          const i32 dy = pmy - my;
    
          const i32 steps = abs(dx) > abs(dy) ? abs(dx) : abs(dy);

          const f64 xInc = dx / (float) steps;
          const f64 yInc = dy / (float) steps;

          f64 x = mx;
          f64 y = my;

          for (int i = 0; i <= steps; i++)
            {
              for (i32 j = -2; j < 2; ++j)
                for (i32 i = -2; i < 2; ++i)
                  {
                    const u64 idx = ((int)x + i) + ((int)y + j) * (pitch / 4);
                    pixels[idx] = SDL_MapRGB (fmt, 255, 0, 0);
                  }
              x += xInc;
              y += yInc;
            }

          SDL_FreeFormat (fmt);
          SDL_UnlockTexture (r->texture);
        }

      pmx = mx, pmy = my;

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

          all_time += t;
          if (t > max_time)
            max_time = t;

          samples++;

          char buffer[32];
          const f64 dr = r->vp->rmax - r->vp->rmin;
          mb_util_format_number (1.0 / dr, buffer, 32);
          printf ("(Recompute) %dms 1/%s\n", t, buffer);
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
  printf ("(Avarage) %.1fms\n", all_time / (float)samples);
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
            r->zoom.delta = 1.0 - MB_ZOOM_DELTA;
            r->zoom.mx = event.button.x;
            r->zoom.my = event.button.y;
          }
        else if (event.button.button == SDL_BUTTON_RIGHT)
          {
            r->zooming = true;
            r->zoom.delta = 1.0 + MB_ZOOM_DELTA;
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
      case SDL_KEYDOWN:
        switch (event.key.keysym.sym)
          {
          case SDLK_F2:
            mb_viewport_init (r->vp, -.75, 0, r->vp->w, r->vp->h);
            r->recompute = true;
            break;
          case SDLK_F3:
            {
              const f64 dr = (r->vp->rmin + r->vp->rmax) / 2;
              const f64 di = (r->vp->imin + r->vp->imax) / 2;
              printf ("%.17lf, %.17lf\n", dr, di);
            }
            break;
          case SDLK_F5:
            r->recompute = true;
            break;
          case SDLK_F6:
            mb_renderer_recompute_aa (r);
            break;
          default:
            break;
          }
        break;
      default:
        break;
      }
}

inline void
mandelbrot_SIMD (f64 const _x0[4], f64 const _y0[4], i32 _xs[4], i32 max)
{
  __m256d x0 = _mm256_loadu_pd (_x0);
  __m256d y0 = _mm256_loadu_pd (_y0);

  __m256d x1 = _mm256_setzero_pd ();
  __m256d y1 = _mm256_setzero_pd ();

  __m256d one = _mm256_set1_pd (1.0);
  __m256d four = _mm256_set1_pd (4.0);

  __m256d xs = _mm256_setzero_pd ();

  for (i32 i = 0; i < max; ++i)
    {
      __m256d x2 = _mm256_mul_pd (x1, x1);
      __m256d y2 = _mm256_mul_pd (y1, y1);

      __m256d s = _mm256_add_pd (x2, y2);

      __m256d mask = _mm256_cmp_pd (s, four, _CMP_LT_OQ);
      if (_mm256_movemask_pd (mask) == 0)
        break;

      __m256d d = _mm256_mul_pd (x1, y1);
      x1 = _mm256_add_pd (_mm256_sub_pd (x2, y2), x0);
      y1 = _mm256_add_pd (_mm256_add_pd (d, d), y0);

      xs = _mm256_add_pd (xs, _mm256_and_pd (mask, one));
    }

  __m128i ixs = _mm256_cvttpd_epi32 (xs);
  _mm_storeu_si128 ((__m128i *)_xs, ixs);
}

inline bool
mb_renderer_solid_guess (i32 const xs[4])
{
  const i32 a = xs[0];
  const i32 b = xs[1];
  const i32 c = xs[2];
  const i32 d = xs[3];
  return (a == b && b == c && c == d);
}

inline void
mb_renderer_rcompute_color (mb_renderer *r, i32 i, u8 *cr, u8 *cg, u8 *cb)
{
  (void) r;
  // -- Red and Cyan --
  // *cr = i * 6;
  // *cg = i * 3;
  // *cb = i * 3;
  *cr = i * 8;
  *cg = i * 4;
  *cb = i * 4;
  // *cr = i * 8;
  // *cg = i * 8;
  // *cb = i * 8;
}

inline void
mb_renderer_recompute (mb_renderer *r)
{
  u32 *pixels;
  int pitch;

  SDL_LockTexture (r->texture, NULL, (void **)&pixels, &pitch);

  SDL_PixelFormat *fmt = SDL_AllocFormat (SDL_PIXELFORMAT_RGB888);

  const mb_viewport *vp = r->vp;

  const i32 w = vp->w;
  const i32 h = vp->h;

  const f64 rmax = vp->rmax;
  const f64 rmin = vp->rmin;
  const f64 imax = vp->imax;
  const f64 imin = vp->imin;

  const f64 step_x = (rmax - rmin) / w;
  const f64 step_y = (imax - imin) / h;

#pragma omp parallel for schedule(dynamic, 1)
  for (i32 y = 0; y < h; y += 4)
    {
      for (i32 x = 0; x < w; x += 4)
        {
          i32 xs[4];

          const f64 x0[4] = {
            rmin + (x + 0) * step_x,
            rmin + (x + 3) * step_x,
            rmin + (x + 0) * step_x,
            rmin + (x + 3) * step_x,
          };

          const f64 y0[4] = {
            imin + (y + 0) * step_y,
            imin + (y + 0) * step_y,
            imin + (y + 3) * step_y,
            imin + (y + 3) * step_y,
          };

          mandelbrot_SIMD (x0, y0, xs, r->max_iter);

          if (mb_renderer_solid_guess (xs))
            {
              const i32 color = (xs[0] == r->max_iter) ? 0 : xs[0];

              for (i32 i = 0; i < 4; ++i)
                for (i32 j = 0; j < 4; ++j)
                  {
                    const u64 idx = (x + j) + (y + i) * (pitch / 4);
                    u8 cr, cg, cb;
                    mb_renderer_rcompute_color (r, color, &cr, &cg, &cb);
                    pixels[idx] = SDL_MapRGB (fmt, cr, cg, cb);
                  }
              continue;
            }

          for (i32 i = 0; i < 4; ++i)
            {
              const f64 x0[4] = {
                rmin + (x + 0) * step_x,
                rmin + (x + 1) * step_x,
                rmin + (x + 2) * step_x,
                rmin + (x + 3) * step_x,
              };

              const f64 im = imin + (y + i) * step_y;
              const f64 y0[4] = { im, im, im, im };

              i32 xs[4];
              mandelbrot_SIMD (x0, y0, xs, r->max_iter);

              for (i32 j = 0; j < 4; ++j)
                {
                  const i32 color = (xs[j] == r->max_iter) ? 0 : xs[j];
                  const u64 idx = (x + j) + (y + i) * (pitch / 4);
                  u8 cr, cg, cb;
                  mb_renderer_rcompute_color (r, color, &cr, &cg, &cb);
                  pixels[idx] = SDL_MapRGB (fmt, cr, cg, cb);
                }
            }
        }
    }

  SDL_FreeFormat (fmt);
  SDL_UnlockTexture (r->texture);

  r->recompute = false;
}

inline void
mb_renderer_recompute_aa (mb_renderer *r)
{
  const i32 samples = 16;

  u32 *pixels;
  int pitch;

  SDL_LockTexture (r->texture, NULL, (void **)&pixels, &pitch);

  SDL_PixelFormat *fmt = SDL_AllocFormat (SDL_PIXELFORMAT_RGB888);

  const mb_viewport *vp = r->vp;

  const i32 w = vp->w;
  const i32 h = vp->h;

  const f64 rmax = vp->rmax;
  const f64 rmin = vp->rmin;
  const f64 imax = vp->imax;
  const f64 imin = vp->imin;

  const f64 step_x = (rmax - rmin) / w;
  const f64 step_y = (imax - imin) / h;

#pragma omp parallel for schedule(dynamic, 1)
  for (i32 y = 0; y < h; y += 4)
    {
      for (i32 x = 0; x < w; x += 4)
        {
          i32 xs[4];

          const f64 x0[4] = {
            rmin + (x + 0) * step_x,
            rmin + (x + 3) * step_x,
            rmin + (x + 0) * step_x,
            rmin + (x + 3) * step_x,
          };

          const f64 y0[4] = {
            imin + (y + 0) * step_y,
            imin + (y + 0) * step_y,
            imin + (y + 3) * step_y,
            imin + (y + 3) * step_y,
          };

          mandelbrot_SIMD (x0, y0, xs, r->max_iter);

          if (mb_renderer_solid_guess (xs))
            {
              const i32 color = (xs[0] == r->max_iter) ? 0 : xs[0];

              for (i32 i = 0; i < 4; ++i)
                for (i32 j = 0; j < 4; ++j)
                  {
                    const u64 idx = (x + j) + (y + i) * (pitch / 4);
                    u8 cr, cg, cb;
                    mb_renderer_rcompute_color (r, color, &cr, &cg, &cb);
                    pixels[idx] = SDL_MapRGB (fmt, cr, cg, cb);
                  }
              continue;
            }

          for (i32 i = 0; i < 4; ++i)
            {
              u32 fr[4] = {0}, fg[4] = {0}, fb[4] = {0};

              for (i32 s = 0; s < samples; ++s)
                {
                  const f64 x0[4] = {
                    rmin + (x + 0 + (2.0 * rand() / RAND_MAX) - 1.0) * step_x,
                    rmin + (x + 1 + (2.0 * rand() / RAND_MAX) - 1.0) * step_x,
                    rmin + (x + 2 + (2.0 * rand() / RAND_MAX) - 1.0) * step_x,
                    rmin + (x + 3 + (2.0 * rand() / RAND_MAX) - 1.0) * step_x,
                  };

                  const f64 y0[4] = {
                    imin + (y + i + (2.0 * rand() / RAND_MAX) - 1.0) * step_y,
                    imin + (y + i + (2.0 * rand() / RAND_MAX) - 1.0) * step_y,
                    imin + (y + i + (2.0 * rand() / RAND_MAX) - 1.0) * step_y,
                    imin + (y + i + (2.0 * rand() / RAND_MAX) - 1.0) * step_y,
                  };

                  // const f64 im = imin + (y + i) * step_y;
                  // const f64 y0[4] = { im, im, im, im };

                  i32 xs[4];
                  mandelbrot_SIMD (x0, y0, xs, r->max_iter);

                  for (i32 j = 0; j < 4; ++j)
                    {
                      const i32 color = (xs[j] == r->max_iter) ? 0 : xs[j];
                      // const u64 idx = (x + j) + (y + i) * (pitch / 4);
                      u8 cr, cg, cb;
                      mb_renderer_rcompute_color (r, color, &cr, &cg, &cb);
                      fr[j] += cr;
                      fg[j] += cg;
                      fb[j] += cb;
                      // pixels[idx] = SDL_MapRGB (fmt, cr, cg, cb);
                    }
                }

              for (i32 j = 0; j < 4; ++j)
                {
                  const u64 idx = (x + j) + (y + i) * (pitch / 4);
                  pixels[idx] = SDL_MapRGB (fmt, fr[j] / samples, fg[j] / samples, fb[j] / samples);
                }
            }
        }
    }

  SDL_FreeFormat (fmt);
  SDL_UnlockTexture (r->texture);

  r->recompute = false;
}

