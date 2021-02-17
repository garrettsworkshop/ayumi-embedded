What is Ayumi-Embeded?
======================

Ayumi is an emulator of the Ayumi AY-3-8910 and YM2149 sound chips. Ayumi-Embedded is a port of Ayumi to the ARMv7-M microcontroller architecture.

The principal changes from Ayumi include the adoption of single-precision floating point, the elimination of the YM2149-compatible mode, and the elimination of stereo output. These changes help to make Ayumi compact enough to work well on 100 MHz-class ARM microcontrollers.

Ayumi API reference
===================

``` c
void ayumi_configure(struct ayumi* ay, int is_ym, float clock_rate, int sr)
```

*Configures the ayumi structure.*

**ay**: The pointer to the ayumi structure.

**clock_rate**: The clock rate of the chip.

**sr**: The output sample rate.

``` c
void ayumi_set_tone(struct ayumi* ay, int index, int period)
```

*Sets the tone period value for the specified sound channel.*

**ay**: The pointer to the ayumi structure.

**index**: The index of the sound channel.

**period**: The tone period value [0-4095].

``` c
void ayumi_set_noise(struct ayumi* ay, int period)
```

*Sets the noise period value.*

**ay**: The pointer to the ayumi structure.

**period**: The noise period value [0-31].

``` c
void ayumi_set_mixer(struct ayumi* ay, int index, int t_off, int n_off, int e_on)
```

*Sets the mixer value for the specified sound channel.*

**ay**: The pointer to the ayumi structure.

**index**: The index of the sound channel.

**t_off**: 1 if the tone is off.

**n_off**: 1 if the noise is off.

**e_on**: 1 if the envelope is on.

``` c
void ayumi_set_volume(struct ayumi* ay, int index, int volume)
```

*Sets the volume for the specified sound channel.*

**ay**: The pointer to the ayumi structure.

**index**: The index of the sound channel.

**volume**: The volume [0-15].

``` c
void ayumi_set_envelope(struct ayumi* ay, int period)
```

*Sets the envelope period value.*

**ay**: The pointer to the ayumi structure.

**period**: The envelope period value [0-65535].

``` c
void ayumi_set_envelope_shape(struct ayumi* ay, int shape)
```

*Sets the envelope shape value.*

**ay**: The pointer to the ayumi structure.

**shape**: The envelope shape index [0-15].

``` c
void ayumi_process(struct ayumi* ay)
```

*Renders the next stereo sample in ay->cur.*

**ay**: The pointer to the ayumi structure.

``` c
void ayumi_remove_dc(struct ayumi* ay)
```

*Removes the DC offset from the current sample.*

**ay**: The pointer to the ayumi structure.
