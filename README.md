# _Four Ways To Write A Pitch-Shifter_ example code

This is the code used to produce the audio demos for the  ADC22 presentation: _Four Ways To Write A Pitch-Shifter_.  It aimes for simplicity rather than performance.

## Organisation

All the different stretching methods are packed into one file.  It definitely could be tidier, and I hope to get to that soon - in the meantime, do get in contact if you have any questions about it.

Each class uses a similar streaming-style API.  You should be able to instantiate them simply:

```cpp
HybridPhaseStretch stretch;
```

Configure them with channel-count and block sizes:

```
double blockMs = 120; 
int blockSamples = int(sampleRate*0.001*blockMs);
stretch.configure(channels, blockSamples, blockSamples/4);
```

Then process blocks:

```
stretch.setTimeFactor(timeFactor);
stretch.setFreqFactor(freqFactor);

// For `outputBufferSize` output samples, how many input samples should we have?
int inputBufferSize = processor.samplesForOutput(outputBufferSize);
processor.process(inputBufferSize, inputSamples, outputBuffers, outputBufferSize);

processBlocks(stretch, timeFactor);
```

If you're not performing a time-stretch, then `inputBufferSize` will always match `outputBufferSize`.

###  Dependencies

It depends on [my DSP library](https://signalsmith-audio.co.uk/code/dsp/) for a couple of things (to get delay-buffers, FFT and Kaiser windows).

## Building

The main logic is in `shift-stretch.h`.

If you use `make`, it compiles `main.cpp` which produces a command-line interface, with built-in help.  (I wouldn't recommend focusing on that code - the main stuff is in `shift-stretch.h`).

## FFTs

It's currently written to use a Modified Real FFT, which has a half-bin offset (so we have _N_ complex bins/bands, instead of _N-1_ complex and 2 real ones).  I'd like to change this in future because it's probably more confusing than necessary.

It should be easy to write adapt this to a more typical Real FFT - the main difference should be in `bandToFreq()`, and perhaps handling the fact that the 0Hz and Nyquist bands will be real. 
