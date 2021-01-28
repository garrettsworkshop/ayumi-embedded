/* Author: Peter Sovietov */

#ifndef AYUMI_H
#define AYUMI_H

enum {
  TONE_CHANNELS = 3,
  DECIMATE_FACTOR = 8,
  FIR_SIZE = 192,
  DC_FILTER_SIZE = 1024
};

struct tone_channel {
  int tone_period;
  int tone_counter;
  int tone;
  int t_off;
  int n_off;
  int e_on;
  int volume;
};

struct interpolator {
  float c[4];
  float y[4];
};

struct dc_filter {
  float sum;
  float delay[DC_FILTER_SIZE];
};

struct ayumi {
  struct tone_channel channels[TONE_CHANNELS];
  int noise_period;
  int noise_counter;
  int noise;
  int envelope_counter;
  int envelope_period;
  int envelope_shape;
  int envelope_segment;
  int envelope;
  const float* dac_table;
  float step;
  float x;
  struct interpolator interp;
  float fir[FIR_SIZE * 2];
  int fir_index;
  struct dc_filter dc;
  int dc_index;
  float cur;
};

int ayumi_configure(struct ayumi* const ay, const float clock_rate, const int sr);
void ayumi_set_tone(struct ayumi* const ay, const int index, const int period);
void ayumi_set_noise(struct ayumi* const ay, const int period);
void ayumi_set_mixer(struct ayumi* const ay, const int index, const int t_off, const int n_off, const int e_on);
void ayumi_set_volume(struct ayumi* const ay, const int index, const int volume);
void ayumi_set_envelope(struct ayumi* const ay, const int period);
void ayumi_set_envelope_shape(struct ayumi* const ay, const int shape);
void ayumi_process(struct ayumi* const ay);
//void ayumi_process2(struct ayumi* const ay1, struct ayumi* const ay2);
void ayumi_remove_dc(struct ayumi* const ay);

#endif
