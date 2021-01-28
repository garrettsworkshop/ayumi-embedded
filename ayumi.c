/* Author: Peter Sovietov */

#include <string.h>
#include "ayumi.h"

static const float const dac_table[] = {
  0.0f, 0.0f,
  0.00999465934234f, 0.00999465934234f,
  0.0144502937362f, 0.0144502937362f,
  0.0210574502174f, 0.0210574502174f,
  0.0307011520562f, 0.0307011520562f,
  0.0455481803616f, 0.0455481803616f,
  0.0644998855573f, 0.0644998855573f,
  0.107362478065f, 0.107362478065f,
  0.126588845655f, 0.126588845655f,
  0.20498970016f, 0.20498970016f,
  0.292210269322f, 0.292210269322f,
  0.372838941024f, 0.372838941024f,
  0.492530708782f, 0.492530708782f,
  0.635324635691f, 0.635324635691f,
  0.805584802014f, 0.805584802014f,
  1.0f, 1.0f
};

static const float const sinc_table[] = {
   0.0f,
  -0.0000046183113992051936f,
  -0.0000111776164088722500f,
  -0.0000186102645020054320f,
  -0.0000251345861356310120f,
  -0.0000284942816906661970f,
  -0.0000263968287932751590f,
  -0.0000170942125588021560f,
   0.0f,
   0.0000237981935769668660f,
   0.0000512811602422021830f,
   0.0000776219782624342700f,
   0.0000967594266641204160f,
   0.0001024022930039340200f,
   0.0000893446142180771060f,
   0.0000548757001189491830f,
   0.0f,
  -0.0000698390822106801650f,
  -0.0001447966132360757000f,
  -0.0002115845291770830800f,
  -0.0002553506910655054400f,
  -0.0002622871437432210400f,
  -0.0002225880592702779900f,
  -0.0001332323049569570400f,
   0.0f,
   0.0001618257876705520600f,
   0.0003284617538509658100f,
   0.0004704561157618486300f,
   0.0005571385145753094400f,
   0.0005621256512151872600f,
   0.0004690191855396247800f,
   0.0002762486683895298600f,
   0.0f,
  -0.0003256417948683862200f,
  -0.0006518231028671038800f,
  -0.0009212778730931929800f,
  -0.0010772534348943575000f,
  -0.0010737727700273478000f,
  -0.0008855664539039263400f,
  -0.0005158189609076553400f,
   0.0f,
   0.0005954876719379527700f,
   0.0011803558710661009000f,
   0.0016527320270369871000f,
   0.0019152679330965555000f,
   0.0018927324805381538000f,
   0.0015481870327877937000f,
   0.0008947069583494130600f,
   0.0f,
  -0.0010178225878206125000f,
  -0.0020037400552054292000f,
  -0.0027874356824117317000f,
  -0.0032103299880219430000f,
  -0.0031540624117984395000f,
  -0.0025657163651900345000f,
  -0.0014750752642111449000f,
   0.0f,
   0.0016624165446378462000f,
   0.0032591192839069179000f,
   0.0045165685815867747000f,
   0.0051838984346123896000f,
   0.0050774264697459933000f,
   0.0041192521414141585000f,
   0.0023628575417966491000f,
   0.0f,
  -0.0026543507866759182000f,
  -0.0051990251084333425000f,
  -0.0072020238234656924000f,
  -0.0082672928192007358000f,
  -0.0081033739572956287000f,
  -0.0065831115395702210000f,
  -0.0037839040415292386000f,
   0.0f,
   0.0042781252851152507000f,
   0.0084176358598320178000f,
   0.0117256605746305500000f,
   0.0135504766477886720000f,
   0.0133881893699974960000f,
   0.0109795012423412590000f,
   0.0063812749416854130000f,
   0.0f,
  -0.0074212296041538880000f,
  -0.0148645630434021300000f,
  -0.0211435846221781040000f,
  -0.0250427505875860900000f,
  -0.0254735309425472010000f,
  -0.0216273100178821960000f,
  -0.0131043233832255430000f,
   0.0f,
   0.0170651339899804760000f,
   0.0369789192644519520000f,
   0.0582331806209395800000f,
   0.0790720120814059490000f,
   0.0976759987169523170000f,
   0.1123604593695093200000f,
   0.1217634357728773100000f,
   0.125f
};

static void reset_segment(struct ayumi* const ay);

static int update_tone(struct ayumi* ay, const int index) {
  struct tone_channel* ch = &ay->channels[index];
  ch->tone_counter += 1;
  if (ch->tone_counter >= ch->tone_period) {
    ch->tone_counter = 0;
    ch->tone ^= 1;
  }
  return ch->tone;
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

static void hold_top(struct ayumi* const ay) {
  (void) ay;
}

static void hold_bottom(struct ayumi* const ay) {
  (void) ay;
}

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

int update_envelope(struct ayumi* const ay) {
  ay->envelope_counter += 1;
  if (ay->envelope_counter >= ay->envelope_period) {
    ay->envelope_counter = 0;
    Envelopes[ay->envelope_shape][ay->envelope_segment](ay);
  }
  return ay->envelope;
}

static void update_mixer(struct ayumi* const ay) {
  int i;
  int out;
  int noise = update_noise(ay);
  int envelope = update_envelope(ay);
  ay->cur = 0;
  for (i = 0; i < TONE_CHANNELS; i += 1) {
    out = (update_tone(ay, i) | ay->channels[i].t_off) & (noise | ay->channels[i].n_off);
    out *= ay->channels[i].e_on ? envelope : ay->channels[i].volume * 2 + 1;
    ay->cur += ay->dac_table[out];
  }
}

int ayumi_configure(struct ayumi* const ay, const float clock_rate, const int sr) {
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

void ayumi_set_tone(struct ayumi* const ay, const int index, const int period) {
  int period_masked = period & 0xffff;
  ay->channels[index].tone_period = (period_masked == 0) | period_masked;
}

void ayumi_set_noise(struct ayumi* const ay, const int period) {
  ay->noise_period = period & 0x1f;
}

void ayumi_set_mixer(struct ayumi* const ay, const int index, const int t_off, const int n_off, const int e_on) {
  ay->channels[index].t_off = t_off & 1;
  ay->channels[index].n_off = n_off & 1;
  ay->channels[index].e_on = e_on;
}

void ayumi_set_volume(struct ayumi* const ay, const int index, const int volume) {
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

static float decimate(float* const x) {
  register float y = 0.0f;
  register const float const *a = x;
  register const float const *b = x+192;
  register const float const *w = w;

  for (int i = 0; i < 96; i+=16) {
    asm volatile(
      // Load a[0-15], b[1-15]
      "VLDMIA.32 %[a]!, { s0,  s1,  s2,  s3,  s4,  s5,  s6,  s7,  s8,  s9,  s10, s11, s12, s13, s14, s15 } \n\t"
      "VLDMDB.32 %[b]!, {      s17, s18, s19, s20, s21, s22, s23, s24, s25, s26, s27, s28, s29, s30, s31 } \n\t"
      // Decrement b because only loaded 15 floats
      "SUB %[b], %[b], #4 \n\t"
      // Add a[1-7,9-15] = a[1-7,9-15] + b[1-7,9-15]
      "VADD.F32 s1,  s1,  s17 \n\t"
      "VADD.F32 s2,  s2,  s18 \n\t"
      "VADD.F32 s3,  s3,  s19 \n\t"
      "VADD.F32 s4,  s4,  s20 \n\t"
      "VADD.F32 s5,  s5,  s21 \n\t"
      "VADD.F32 s6,  s6,  s22 \n\t"
      "VADD.F32 s7,  s7,  s23 \n\t"
      "VADD.F32 s9,  s9,  s25 \n\t"
      "VADD.F32 s10, s10, s26 \n\t"
      "VADD.F32 s11, s11, s27 \n\t"
      "VADD.F32 s12, s12, s28 \n\t"
      "VADD.F32 s13, s13, s29 \n\t"
      "VADD.F32 s14, s14, s30 \n\t"
      "VADD.F32 s15, s15, s31 \n\t"
      // Load w[1-15] now that S17-S31 are free
      "VLDMIA.32 %[w]!, { s17, s18, s19, s20, s21, s22, s23, s24, s25, s26, s27, s28, s29, s30, s31 } \n\t"
      // Increment w because only loaded 15 floats
      "ADD %[w], %[w], #4 \n\t"
      // Multiply to compute w[1-7,9-15] * (a[1-7,9-15] + b[1-7,9-15])
      "VMUL.F32 s1,  s1,  s17 \n\t"
      "VMUL.F32 s2,  s2,  s18 \n\t"
      "VMUL.F32 s3,  s3,  s19 \n\t"
      "VMUL.F32 s4,  s4,  s20 \n\t"
      "VMUL.F32 s5,  s5,  s21 \n\t"
      "VMUL.F32 s6,  s6,  s22 \n\t"
      "VMUL.F32 s7,  s7,  s23 \n\t"
      "VMUL.F32 s9,  s9,  s25 \n\t"
      "VMUL.F32 s10, s10, s26 \n\t"
      "VMUL.F32 s11, s11, s27 \n\t"
      "VMUL.F32 s12, s12, s28 \n\t"
      "VMUL.F32 s13, s13, s29 \n\t"
      "VMUL.F32 s14, s14, s30 \n\t"
      "VMUL.F32 s15, s15, s31 \n\t"
      // Accumulate w[1-7,9-15] * (a[1-7,9-15] + b[1-7,9-15]) into y register
      "VADD.F32 %[y], %[y], s1  \n\t"
      "VADD.F32 %[y], %[y], s2  \n\t"
      "VADD.F32 %[y], %[y], s3  \n\t"
      "VADD.F32 %[y], %[y], s4  \n\t"
      "VADD.F32 %[y], %[y], s5  \n\t"
      "VADD.F32 %[y], %[y], s6  \n\t"
      "VADD.F32 %[y], %[y], s7  \n\t"
      "VADD.F32 %[y], %[y], s9  \n\t"
      "VADD.F32 %[y], %[y], s10 \n\t"
      "VADD.F32 %[y], %[y], s11 \n\t"
      "VADD.F32 %[y], %[y], s12 \n\t"
      "VADD.F32 %[y], %[y], s13 \n\t"
      "VADD.F32 %[y], %[y], s14 \n\t"
      "VADD.F32 %[y], %[y], s15 \n\t"
      : [a] "+r" (a), [b] "+r" (b),
        [w] "+r" (w), [y] "+t" (y)
      :
      : "s0",  "s1",  "s2",  "s3",  "s4",  "s5",  "s6",  "s7",
        "s8",  "s9",  "s10", "s11", "s12", "s13", "s14", "s15",
               "s17", "s18", "s19", "s20", "s21", "s22", "s23",
        "s24", "s25", "s26", "s27", "s28", "s29", "s30", "s31"
    );
  }
  y += sinc_table[96] * x[96];

  memcpy(&x[FIR_SIZE - DECIMATE_FACTOR], x, DECIMATE_FACTOR * sizeof(float));
  return y;
}

static float* ayumi_process_internal(struct ayumi* const ay) {
  int i;
  float y1;
  float* c = ay->interp.c;
  float* y = ay->interp.y;
  float* fir = &ay->fir[FIR_SIZE - ay->fir_index * DECIMATE_FACTOR];
  ay->fir_index = (ay->fir_index + 1) % (FIR_SIZE / DECIMATE_FACTOR - 1);
  for (i = DECIMATE_FACTOR - 1; i >= 0; i -= 1) {
    ay->x += ay->step;
    if (ay->x >= 1) {
      ay->x -= 1;
      y[0] = y[1];
      y[1] = y[2];
      y[2] = y[3];
      update_mixer(ay);
      y[3] = ay->cur;
      y1 = y[2] - y[0];
      c[0] = 0.5f * y[1] + 0.25f * (y[0] + y[2]);
      c[1] = 0.5f * y1;
      c[2] = 0.25f * (y[3] - y[1] - y1);
    }
    fir[i] = (c[2] * ay->x + c[1]) * ay->x + c[0];
  }
}

void ayumi_process(struct ayumi* const ay) {
  ay->cur = decimate(ayumi_process_internal(ay));
}

static float dc_filter(struct dc_filter* const dc, const int index, const float x) {
  dc->sum += -dc->delay[index] + x;
  dc->delay[index] = x; 
  return x - dc->sum / DC_FILTER_SIZE;
}

void ayumi_remove_dc(struct ayumi* const ay) {
  ay->cur = dc_filter(&ay->dc, ay->dc_index, ay->cur);
  ay->dc_index = (ay->dc_index + 1) & (DC_FILTER_SIZE - 1);
}
