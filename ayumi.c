/* Author: Peter Sovietov */

#include <string.h>
#include "ayumi.h"
#include "ayumi_tables.c"

static void reset_segment(struct ayumi* const ay);
static void slide_up(struct ayumi* const ay) {
  ay->envelope += 1;
  if (ay->envelope > 31) {
    ay->envelope_segment ^= 1;
    reset_segment(ay);
  }
}
static void slide_down(struct ayumi* const ay) {
  ay->envelope -= 1;
  if (ay->envelope < 0) {
    ay->envelope_segment ^= 1;
    reset_segment(ay);
  }
}
static void hold_top(struct ayumi* const ay) { (void) ay; }
static void hold_bottom(struct ayumi* const ay) { (void) ay; }

static void (* const Envelopes[][2])(struct ayumi* const) = {
  {slide_down, hold_bottom},
  {slide_down, hold_bottom},
  {slide_down, hold_bottom},
  {slide_down, hold_bottom},
  {slide_up, hold_bottom},
  {slide_up, hold_bottom},
  {slide_up, hold_bottom},
  {slide_up, hold_bottom},
  {slide_down, slide_down},
  {slide_down, hold_bottom},
  {slide_down, slide_up},
  {slide_down, hold_top},
  {slide_up, slide_up},
  {slide_up, hold_top},
  {slide_up, slide_down},
  {slide_up, hold_bottom}
};

static void reset_segment(struct ayumi* const ay) {
  if (Envelopes[ay->envelope_shape][ay->envelope_segment] == slide_down
    || Envelopes[ay->envelope_shape][ay->envelope_segment] == hold_top) {
    ay->envelope = 31;
    return;
  }
  ay->envelope = 0;
}

static int update_noise(struct ayumi* const ay) {
  int bit0x3;
  ay->noise_counter += 1;
  if (ay->noise_counter >= (ay->noise_period << 1)) {
    ay->noise_counter = 0;
    bit0x3 = ((ay->noise ^ (ay->noise >> 3)) & 1);
    ay->noise = (ay->noise >> 1) | (bit0x3 << 16);
  }
  return ay->noise & 1;
}

static int update_envelope(struct ayumi* const ay) {
  ay->envelope_counter += 1;
  if (ay->envelope_counter >= ay->envelope_period) {
    ay->envelope_counter = 0;
    Envelopes[ay->envelope_shape][ay->envelope_segment](ay);
  }
  return ay->envelope;
}

static int update_tone(struct ayumi* ay, const int index) {
  struct tone_channel* ch = &ay->channels[index];
  ch->tone_counter += 1;
  if (ch->tone_counter >= ch->tone_period) {
    ch->tone_counter = 0;
    ch->tone ^= 1;
  }
  return ch->tone;
}

#define update_mixer_loop(i) \
  out = (update_tone(ay, i) | ay->channels[i].t_off) & \
        (noise | ay->channels[i].n_off); \
  out *= ay->channels[i].e_on ? envelope : ay->channels[i].volume * 2 + 1; \
  cur += dac_table[out]; \

static float update_mixer(struct ayumi* const ay) {
  int out;
  int noise = update_noise(ay);
  int envelope = update_envelope(ay);
  float cur = 0;
  update_mixer_loop(0);
  update_mixer_loop(1);
  update_mixer_loop(2);
  return ay->cur = cur;
}

int ayumi_configure(struct ayumi* const ay, 
                    const float clock_rate, const int sr) {
  int i;
  memset(ay, 0, sizeof(struct ayumi));
  ay->step = clock_rate / (sr * 8 * DECIMATE_FACTOR);
  ay->noise = 1;
  ayumi_set_envelope(ay, 1);
  for (i = 0; i < TONE_CHANNELS; i += 1) {
    ayumi_set_tone(ay, i, 1);
  }
  return ay->step < 1;
}

void ayumi_set_tone(struct ayumi* const ay, const int index,
                    const int period) {
  int period_masked = period & 0xffff;
  ay->channels[index].tone_period = (period_masked == 0) | period_masked;
}

void ayumi_set_noise(struct ayumi* const ay, const int period) {
  ay->noise_period = period & 0x1f;
}

void ayumi_set_mixer(struct ayumi* const ay, const int index, 
                     const int t_off, const int n_off, const int e_on) {
  ay->channels[index].t_off = t_off & 1;
  ay->channels[index].n_off = n_off & 1;
  ay->channels[index].e_on = e_on;
}

void ayumi_set_volume(struct ayumi* const ay,
                      const int index, const int volume) {
  ay->channels[index].volume = volume & 0xf;
}

void ayumi_set_envelope(struct ayumi* const ay, const int period) {
  int period_masked = period &  0xffff;
  ay->envelope_period = (period_masked == 0) | period_masked;
}

void ayumi_set_envelope_shape(struct ayumi* const ay, const int shape) {
  ay->envelope_shape = shape & 0xf;
  ay->envelope_counter = 0;
  ay->envelope_segment = 0;
  reset_segment(ay);
}

#define decimate_mac16 \
  asm volatile( \
  /* Load a[0-15], b[1-15] */ \
  "VLDMIA.32 %[a]!, { s0,  s1,  s2,  s3,  s4,  s5,  s6,  s7, " \
  "                   s8,  s9,  s10, s11, s12, s13, s14, s15 } \n\t" \
  "VLDMDB.32 %[b]!, {      s17, s18, s19, s20, s21, s22, s23, " \
  "                   s24, s25, s26, s27, s28, s29, s30, s31 } \n\t" \
  /* Decrement b because only loaded 15 floats */ \
  "SUB %[b], %[b], #4 \n\t" \
  /* Add a[1-7,9-15] = a[1-7,9-15] + b[1-7,9-15] */ \
  "VADD.F32 s1,  s1,  s17 \n\t" \
  "VADD.F32 s2,  s2,  s18 \n\t" \
  "VADD.F32 s3,  s3,  s19 \n\t" \
  "VADD.F32 s4,  s4,  s20 \n\t" \
  "VADD.F32 s5,  s5,  s21 \n\t" \
  "VADD.F32 s6,  s6,  s22 \n\t" \
  "VADD.F32 s7,  s7,  s23 \n\t" \
  "VADD.F32 s9,  s9,  s25 \n\t" \
  "VADD.F32 s10, s10, s26 \n\t" \
  "VADD.F32 s11, s11, s27 \n\t" \
  "VADD.F32 s12, s12, s28 \n\t" \
  "VADD.F32 s13, s13, s29 \n\t" \
  "VADD.F32 s14, s14, s30 \n\t" \
  "VADD.F32 s15, s15, s31 \n\t" \
  /* Load w[1-15] now that S17-S31 are free */ \
  "VLDMIA.32 %[w]!, {      s17, s18, s19, s20, s21, s22, s23, " \
  "                   s24, s25, s26, s27, s28, s29, s30, s31 } \n\t" \
  /* Increment w because only loaded 15 floats */ \
  "ADD %[w], %[w], #4 \n\t" \
  /* Multiply to compute w[1-7,9-15] * (a[1-7,9-15] + b[1-7,9-15]) \
     Interleaved with accumulating 
     w[1-7,9-15] * (a[1-7,9-15] + b[1-7,9-15]) \
     into y register */ \
  "VMUL.F32 s1,  s1,  s17 \n\t" \
  "VADD.F32 %[y], %[y], s1  \n\t" \
  "VMUL.F32 s2,  s2,  s18 \n\t" \
  "VADD.F32 %[y], %[y], s2  \n\t" \
  "VMUL.F32 s3,  s3,  s19 \n\t" \
  "VADD.F32 %[y], %[y], s3  \n\t" \
  "VMUL.F32 s4,  s4,  s20 \n\t" \
  "VADD.F32 %[y], %[y], s4  \n\t" \
  "VMUL.F32 s5,  s5,  s21 \n\t" \
  "VADD.F32 %[y], %[y], s5  \n\t" \
  "VMUL.F32 s6,  s6,  s22 \n\t" \
  "VADD.F32 %[y], %[y], s6  \n\t" \
  "VMUL.F32 s7,  s7,  s23 \n\t" \
  "VADD.F32 %[y], %[y], s7  \n\t" \
  "VMUL.F32 s9,  s9,  s25 \n\t" \
  "VADD.F32 %[y], %[y], s9  \n\t" \
  "VMUL.F32 s10, s10, s26 \n\t" \
  "VADD.F32 %[y], %[y], s10 \n\t" \
  "VMUL.F32 s11, s11, s27 \n\t" \
  "VADD.F32 %[y], %[y], s11 \n\t" \
  "VMUL.F32 s12, s12, s28 \n\t" \
  "VADD.F32 %[y], %[y], s12 \n\t" \
  "VMUL.F32 s13, s13, s29 \n\t" \
  "VADD.F32 %[y], %[y], s13 \n\t" \
  "VMUL.F32 s14, s14, s30 \n\t" \
  "VADD.F32 %[y], %[y], s14 \n\t" \
  "VMUL.F32 s15, s15, s31 \n\t" \
  : [a] "+r" (a), [b] "+r" (b), \
    [w] "+r" (w), [y] "+t" (y) \
  : \
  : "s0",  "s1",  "s2",  "s3",  "s4",  "s5",  "s6",  "s7", \
    "s8",  "s9",  "s10", "s11", "s12", "s13", "s14", "s15", \
           "s17", "s18", "s19", "s20", "s21", "s22", "s23", \
    "s24", "s25", "s26", "s27", "s28", "s29", "s30", "s31" \
);

static float decimate(float *x) {
  register float *a = x;
  register float *b = x+192;
  register const float *w = w;
  register float y = 0.0f;

  decimate_mac16;
  decimate_mac16;
  decimate_mac16;
  decimate_mac16;
  decimate_mac16;
  decimate_mac16;
  y += sinc_table[96] * x[96];

  memcpy(&x[FIR_SIZE - DECIMATE_FACTOR], x, DECIMATE_FACTOR * sizeof(float));
  return y;
}

static float* ayumi_process_internal(struct ayumi* const ay) {
  // Caller-save VFP registers (not preserved across loop iterations)
  register float y2_m_y0  __asm__("s2");
  register float fir_i    __asm__("s3");
  register float c0       __asm__("s4");  // These must be consecutive
  register float c1       __asm__("s5");  // to use VLDMIA and VSTMIA
  register float c2       __asm__("s6");  // .
  register float y0old    __asm__("s7");  // .
  register float y0       __asm__("s8");  // .
  register float y1       __asm__("s9");  // .
  register float y2       __asm__("s10"); // .
  register float cur      __asm__("s11"); // .

  // Callee-save VFP registers (preserved across loop iterations)
  register float half     __asm__("s17") = 0.5f;
  register float quarter  __asm__("s18") = 0.25f;
  register float ay_x     __asm__("s19") = ay->x;
  register float ay_step  __asm__("s20") = ay->step;
  register int i;
  register float* c = ay->interp.c;
  register float* y = ay->interp.y;
  register float* fir = &ay->fir[FIR_SIZE - ay->fir_index * DECIMATE_FACTOR];

  ay->fir_index = (ay->fir_index + 1) % (FIR_SIZE / DECIMATE_FACTOR - 1);

  for (i = DECIMATE_FACTOR - 1; i >= 0; i -= 1) {
    ay_x += ay_step;
    if (ay_x >= 1) {
      cur = update_mixer(ay);
      asm volatile(
        /* Decrement ay_x -= 1 */
        /* Load old y[0-3] */
        "VLDMIA.32 %[y], { %[y0old], %[y0], %[y1], %[y2] } \n\t"
        /* Compute y2_m_y0 = y[2] - y[0] */
        "VSUB.F32 %[y2_m_y0], %[y2], %[y0] \n\t"
        /* Compute c0 = 0.5f * y[1] + 0.25f * (y[0] + y[2]) */
        "VADD.F32 %[c0], %[y0], %[y2] \n\t"
        "VMUL.F32 %[c0], %[c0], %[quarter] \n\t"
        "VMUL.F32 %[c1], %[y1], %[half] \n\t" /* Use c1 to hold 0.5f * y[1] */
        "VADD.F32 %[c0], %[c0], %[c1] \n\t"
        /* Compute c1 = 0.5f * y2_m_y0 */
        "VADD.F32 %[c1], %[y2_m_y0], %[half] \n\t"
        /* Compute c2 = 0.25f * (y[3] - y[1] - y2_m_y0) */
        "VSUB.F32 %[c2], %[cur], %[y1] \n\t"
        "VSUB.F32 %[c2], %[c2], %[y2_m_y0] \n\t"
        "VMUL.F32 %[c2], %[c2], %[quarter] \n\t"
        /* Store computed c[0-2]. */
        /* Also tore old y[1-3] into y[0-2] and cur in y[3]. */
        /* Notice we store y0old in c[3]. c[3] is a don't care 
         * so this store doesn't do anything. */
        "VSTMIA.32 %[y], { %[c0], %[c1], %[c2], %[y0old], 
                           %[y0], %[y1], %[y2], %[cur] } \n\t"
        : [cur] "=t" (cur), [y2_m_y0] "=t" (y2_m_y0), 
          [c0] "=t" (c0), [c1] "=t" (c1), [c2] "=t" (c2), 
          [y0old] "=t" (y0old), [y0] "=t" (y0), [y1] "=t" (y1), [y2] "=t" (y2)
        : [y] "r" (y), [half] "t" (half), [quarter] "t" (quarter)
        :
      );
    }
    asm volatile(
      /* Compute and store fir[i] = (c[2] * ay->x + c[1]) * ay->x + c[0] 
         Also store ay->x back from register. */
      "VMUL.F32 %[fir_i], %[c2], %[ay_x] \n\t"
      "VADD.F32 %[fir_i], %[fir_i], %[c1] \n\t"
      "VMUL.F32 %[fir_i], %[fir_i], %[ay_x] \n\t"
      "VADD.F32 %[fir_i], %[fir_i], %[c0] \n\t"
      "VSTR.32 %[fir_i], [%[fir_i_ptr]] \n\t"
      "VSTR.32 %[ay_x], [%[ay_x_ptr]] \n\t"
      : [fir_i] "=t" (fir_i)
      : [c0] "t" (c0), [c1] "t" (c1), [c2] "t" (c2), [ay_x] "t" (ay_x),
        [fir_i_ptr] "r" (&fir[i]), [ay_x_ptr] "r" (&ay->x)
      :
    );
  }
  return fir;
}

void ayumi_process(struct ayumi* const ay) {
  ay->cur = decimate(ayumi_process_internal(ay));
}

static float dc_filter(struct dc_filter* const dc, 
                       const int index, const float x) {
  dc->sum += -dc->delay[index] + x;
  dc->delay[index] = x; 
  return x - dc->sum / DC_FILTER_SIZE;
}

void ayumi_remove_dc(struct ayumi* const ay) {
  ay->cur = dc_filter(&ay->dc, ay->dc_index, ay->cur);
  ay->dc_index = (ay->dc_index + 1) & (DC_FILTER_SIZE - 1);
}
