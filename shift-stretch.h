#ifndef SIGNALSMITH_EXAMPLE_SHIFT_STRETCH_H
#define SIGNALSMITH_EXAMPLE_SHIFT_STRETCH_H

#include "dsp/delay.h"
#include "dsp/windows.h"
#include <vector>

#ifndef LOG_EXPR
#	include <iostream>
#	define LOG_EXPR(expr) std::cout << #expr << " = " << (expr) << std::endl;
#endif

class OverlapAddStretch {
public:
	using Sample = double;

	OverlapAddStretch(bool isSpectral=false) : isSpectral(isSpectral) {}
	
	void configure(int channels, int blockSamples, int intervalSamples, int maxExtraInput=0) {
		this->channels = channels;
		this->blockSamples = blockSamples;
		this->intervalSamples = intervalSamples;
		this->maxSurplusInputSamples = maxExtraInput;

		inputHistory.resize(channels, blockSamples + maxExtraInput);
		summedOutput.resize(channels, blockSamples);
		blockBuffers.resize(blockSamples*channels);
		window.resize(blockSamples);
		if (isSpectral) {
			// Kaiser's a good window for spectral stuff, but not so great for time-domain
			auto kaiser = signalsmith::windows::Kaiser::withBandwidth(blockSamples*1.0/intervalSamples, true);
			kaiser.fill(window, blockSamples);
		} else {
			for (int i = 0; i < blockSamples; ++i) {
				double r = (i + 0.5)/blockSamples;
				window[i] = std::sin(r*M_PI); // sine window, becomes Hann when applied twice
			}
		}
		// Makes it add up nicely to 1 when applied twice
		signalsmith::windows::forcePerfectReconstruction(window, blockSamples, intervalSamples);
		
		intervalCounter = 0;
	}
	
	void reset() {
		inputHistory.reset();
		summedOutput.reset();
		intervalCounter = 0;
	}
	
	void setRate(double rate) {
		invTimeFactor = rate;
	}
	void setTimeFactor(double timeFactor) {
		invTimeFactor = 1/timeFactor;
	}
	/// How many input samples do we need to get this much output?
	int samplesForOutput(int outputSamples) const {
		double inputSamples = outputSamples*invTimeFactor - surplusInputSamples;
		return int(std::ceil(inputSamples));
	}
	
	void process(const Sample * const *inputs, int inputSamples, Sample **outputs, int outputSamples) {
		int inputFilledTo = 0;
		for (int o = 0; o < outputSamples; ++o) {
			if (++intervalCounter >= intervalSamples) {
				intervalCounter = 0;
				// Fill the block from the input
				int inputStart = int(std::round(o*invTimeFactor - surplusInputSamples - blockSamples));
				// For safety: don't go past the end of the block, or too far in the past
				inputStart = std::max(std::min(inputStart, inputSamples - blockSamples), -maxSurplusInputSamples - blockSamples);
				for (int c = 0; c < channels; ++c) {
					// Make sure we have enough input history
					auto input = inputs[c];
					auto history = inputHistory[c];
					for (int i = inputFilledTo; i < inputStart + blockSamples; ++i) {
						history[i] = input[i];
					}
					// Fill the block from history
					Sample *blockBuffer = channelBlock(c);
					for (int i = 0; i < blockSamples; ++i) {
						blockBuffer[i] = history[inputStart + i]*window[i];
					}
				}
				
				processBlock(inputStart - prevInputIndex);
				prevInputIndex = inputStart;
				
				// Add the block to the summed output
				for (int c = 0; c < channels; ++c) {
					Sample *blockBuffer = channelBlock(c);
					auto output = summedOutput[c];
					for (int i = 0; i < blockSamples; ++i) {
						output[i] += blockBuffer[i]*window[i];
					}
				}
			}
			for (int c = 0; c < channels; ++c) {
				outputs[c][o] = summedOutput[c][0];
				summedOutput[c][0] = 0;
			}
			++summedOutput;
		}
		
		// Copy in remaining input
		for (int c = 0; c < channels; ++c) {
			auto input = inputs[c];
			auto history = inputHistory[c];
			for (int i = inputFilledTo; i < inputSamples; ++i) {
				history[i] = input[i];
			}
		}
		inputHistory += inputSamples;
		prevInputIndex -= inputSamples;
		surplusInputSamples += inputSamples - outputSamples*invTimeFactor;
	}
	
	int inputLatency() const {
		return blockSamples/2;
	}
	int outputLatency() const {
		return blockSamples - inputLatency();
	}

protected:
	int channels = 0, blockSamples = 0;
	int intervalSamples = 0, intervalCounter = 0;
	double invTimeFactor = 1;
	
	Sample * channelBlock(int channel) {
		return blockBuffers.data() + channel*blockSamples;
	}
	
	virtual void processBlock(int inputIntervalSamples) {
		// Alter the blocks (for each channel) if we want to
		(void)inputIntervalSamples;
	}

	void scheduleNextBlock(int interval) {
		intervalCounter = intervalSamples - interval;
	}
private:
	bool isSpectral;

	// Multi-channel circular buffers
	signalsmith::delay::MultiBuffer<Sample> inputHistory, summedOutput;
	std::vector<Sample> blockBuffers, window;

	// Unused input samples, which may be fractional
	int maxSurplusInputSamples = 0;
	double surplusInputSamples = 0;
	int prevInputIndex = 0;
};

class WsolaStretch : public OverlapAddStretch {
public:
	void configure(int channels, int blockSamples, int intervalSamples, int searchSamples, int maxExtraInput=0) {
		OverlapAddStretch::configure(channels, blockSamples, intervalSamples, maxExtraInput);

		// Don't let our search reduce the crossover time by more than half
		searchSamples = std::min(searchSamples, (blockSamples - intervalSamples)/2);
		searchSamples = std::min(searchSamples, intervalSamples);
		this->searchSamples = searchSamples;

		previousBlocks.resize(channels*blockSamples);
	}

	void reset() {
		previousBlocks.assign(previousBlocks.size(), 0);
	}
protected:
	void processBlock(int) override {
		int bestOffset = 0;
		Sample bestDifferenceScore = -1;
		// Brute-force search to find the offset with minimum waveform difference
		for (int offset = -searchSamples; offset <= searchSamples; ++offset) {
			Sample sumDiff2 = 0, sum2 = 0;
			int startIndex = std::max(0, -offset), endIndex = std::min(this->blockSamples, this->blockSamples - offset);
			for (int i = startIndex; i < endIndex; ++i) {
				for (int c = 0; c < this->channels; ++c) {
					Sample prevSample = previousBlock(c)[i];
					Sample currentSample = this->channelBlock(c)[i + offset];
					Sample diff = prevSample - currentSample;
					sumDiff2 += diff*diff;
					sum2 += prevSample*prevSample + currentSample*currentSample;
				}
			}
			
			if (sum2 > 0) {
				Sample score = sumDiff2/sum2;
				if (bestDifferenceScore < 0 || score < bestDifferenceScore) {
					bestOffset = offset;
					bestDifferenceScore = score;
				}
			}
		}
		
		// Apply the offset
		if (bestOffset < 0) {
			for (int c = 0; c < this->channels; ++c) {
				Sample *block = this->channelBlock(c);
				for (int i = this->blockSamples - 1 - bestOffset; i >= 0; --i) {
					block[i - bestOffset] = block[i];
				}
				for (int i = 0; i < -bestOffset; ++i) {
					block[i] = 0;
				}
			}
		} else if (bestOffset > 0) {
			for (int c = 0; c < this->channels; ++c) {
				Sample *block = this->channelBlock(c);
				for (int i = bestOffset; i < this->blockSamples; ++i) {
					block[i - bestOffset] = block[i];
				}
				for (int i = this->blockSamples - bestOffset; i < this->blockSamples; ++i) {
					block[i] = 0;
				}
			}
		}
		
		// If our offset moved things earlier, request the next block earlier, etc.
		int nextInterval = this->intervalSamples - bestOffset/2;
		this->scheduleNextBlock(nextInterval);
		
		// Store the block for next time
		for (int c = 0; c < this->channels; ++c) {
			Sample *block = this->channelBlock(c), *prevBlock = this->previousBlock(c);
			for (int i = 0; i < this->blockSamples - nextInterval; ++i) {
				prevBlock[i] = block[i + nextInterval];
			}
			for (int i = this->blockSamples - nextInterval; i < this->blockSamples; ++i) {
				prevBlock[i] = 0;
			}
		}
	}
private:
	int searchSamples = 0;
	
	Sample *previousBlock(int channel) {
		return previousBlocks.data() + channel*this->blockSamples;
	}
	std::vector<Sample> previousBlocks;
};

#include "dsp/fft.h"
#include <complex>

class SpectralStretch : public OverlapAddStretch {
public:
	using Complex = std::complex<Sample>;

	SpectralStretch(bool kaiserWindow=true) : OverlapAddStretch(kaiserWindow) {}
	
	void configure(int channels, int blockSamples, int intervalSamples, double zeroPadding=1, int maxExtraInput=0) {
		OverlapAddStretch::configure(channels, blockSamples, intervalSamples, maxExtraInput);
		
		mrfft.setFastSizeAbove(blockSamples*zeroPadding);
		fftBuffer.resize(mrfft.size());
		bandCount = mrfft.size()/2;
		scalingFactor = 1.0/mrfft.size(); // the FFT round-trip scales things up, so we scale down again
		channelSpectra.resize(bandCount*channels);
	}
protected:
	virtual void processSpectrum(int inputIntervalSamples) {
		// Edit the spectrums using `channelSpectrum()`, `bands()` and `bandToFreq()`/`freqToBand()`
		(void)inputIntervalSamples;
	}

	Complex * channelSpectrum(int channel) {
		return channelSpectra.data() + channel*bandCount;
	}

	int bands() const {
		return bandCount;
	}
	int fftSize() const {
		return int(mrfft.size());
	}
	Sample bandToFreq(Sample band) const {
		return (band + 0.5f)/mrfft.size();
	}
	Sample freqToBand(Sample freq) const {
		return freq*mrfft.size() - 0.5f;
	}

	void timeShiftPhases(Sample shiftSamples, Complex *output) const {
		for (int b = 0; b < bandCount; ++b) {
			Sample phase = bandToFreq(b)*shiftSamples*(-2*M_PI);
			output[b] = {std::cos(phase), std::sin(phase)};
		}
	}

	void processBlock(int inputIntervalSamples) override final {
		for (int c = 0; c < this->channels; ++c) {
			Sample *block = this->channelBlock(c);
			Complex *spectrum = channelSpectrum(c);
			for (int i = 0; i < this->blockSamples; ++i) {
				fftBuffer[i] = block[i];
			}
			// Zero-padding
			for (int i = this->blockSamples; i < int(fftBuffer.size()); ++i) {
				fftBuffer[i] = 0;
			}
			mrfft.fft(fftBuffer, spectrum);
		}

		processSpectrum(inputIntervalSamples);

		for (int c = 0; c < this->channels; ++c) {
			Sample *block = this->channelBlock(c);
			Complex *spectrum = channelSpectrum(c);
			mrfft.ifft(spectrum, fftBuffer);
			for (int i = 0; i < this->blockSamples; ++i) {
				block[i] = fftBuffer[i]*scalingFactor;
			}
		}
	}

	static Complex generateComplex(Sample energy, Complex complexPhase) {
		Sample complexPhaseNorm = std::norm(complexPhase);
		if (complexPhaseNorm > 0) {
			return complexPhase*std::sqrt(energy/complexPhaseNorm);
		} else {
			Sample phase = Sample(2*M_PI)*rand()/RAND_MAX;
			Complex complexPhase = {std::cos(phase), std::sin(phase)};
			return std::sqrt(energy)*complexPhase;
		}
	}

private:
	signalsmith::fft::ModifiedRealFFT<Sample> mrfft{1};
	int bandCount = 0;
	Sample scalingFactor = 1;
	std::vector<Sample> fftBuffer;
	std::vector<Complex> channelSpectra;
};

class SpectralCutStretch : public SpectralStretch {
public:
	SpectralCutStretch(bool fixedPhase) : fixedPhase(fixedPhase) {}

	void configure(int channels, int blockSamples, int intervalSamples, double zeroPadding=2, int maxExtraInput=0) {
		SpectralStretch::configure(channels, blockSamples, intervalSamples, zeroPadding, maxExtraInput);
		
		energy.resize(this->bands());
		smoothedEnergy.resize(this->bands());
		newSpectra.resize(this->bands()*channels);
		prevSpectra.resize(this->bands()*channels);
		
		prevOutputRotations.resize(bands());
		timeShiftPhases(-intervalSamples, prevOutputRotations.data());
	}
	
	void reset() {
		prevSpectra.assign(prevSpectra.size(), 0);
	}

	void setFreqFactor(double factor) {
		freqFactor = factor;
	}
protected:
	virtual void processSpectrum(int) {
		for (int b = 0; b < this->bands(); ++b) {
			Sample e = 0;
			for (int c = 0; c < this->channels; ++c) {
				Complex bin = this->channelSpectrum(c)[b];
				e += std::norm(bin); // magnitude squared
			}
			energy[b] = smoothedEnergy[b] = e;
		}
		
		Sample smoothingFactor = 0.25; // Really this should depend on your overlap-ratio and stuff, but this whole thing's a bit approximate
		Sample smooth = energy[0];
		for (int b = 1; b < this->bands(); ++b) { // smooth upwards
			smooth += (smoothedEnergy[b] - smooth)*smoothingFactor;
			smoothedEnergy[b] = smooth;
		}
		for (int b = this->bands() - 1; b >= 0; --b) { // smooth downwards
			smooth += (smoothedEnergy[b] - smooth)*smoothingFactor;
			smoothedEnergy[b] = smooth;
		}
		
		int binIndex = 0;
		int prevSegmentStart = 0, prevSegmentEnd = 0;
		while (binIndex < this->bands()) {
			if (energy[binIndex] > smoothedEnergy[binIndex]) {
				if (prevSegmentEnd > 0) { // if it's not the first segment
					// backtrack until it's 6dB below the smoothed energy
					int segmentStart = binIndex;
					while (segmentStart > 0 && energy[segmentStart] > smoothedEnergy[segmentStart]*0.25f) {
						--segmentStart;
					}
				
					// extend this segment back and the previous one forwards
					int midPoint = (segmentStart + prevSegmentEnd)/2;
					// copy the previous segment across
					copySegmentToNew(prevSegmentStart, midPoint);
					prevSegmentStart = midPoint;
				}
				
				// and extend forward until it's 6dB below the smoothed energy
				int segmentEnd = binIndex + 1;
				while (segmentEnd < this->bands() && energy[segmentEnd] > smoothedEnergy[segmentEnd]*0.25f) {
					++segmentEnd;
				}
				prevSegmentEnd = binIndex = segmentEnd;
			} else {
				++binIndex;
			}
		}
		// Extend final band to the end, and copy it in
		copySegmentToNew(prevSegmentStart, this->bands());
		
		// Copy the new spectrum across
		for (int c = 0; c < this->channels; ++c) {
			Complex *spectrum = this->channelSpectrum(c);
			Complex *newSpectrum = newChannelSpectrum(c);
			Complex *prevSpectrum = prevChannelSpectrum(c);
			for (int b = 0; b < this->bands(); ++b) {
				spectrum[b] = newSpectrum[b];
				newSpectrum[b] = 0;
				prevSpectrum[b] = spectrum[b]*prevOutputRotations[b];
			}
		}
	}
private:
	bool fixedPhase;
	double freqFactor = 1;
	std::vector<Sample> energy, smoothedEnergy;
	std::vector<Complex> newSpectra, prevSpectra;
	Complex * newChannelSpectrum(int channel) {
		return newSpectra.data() + channel*this->bands();
	}
	Complex * prevChannelSpectrum(int channel) {
		return prevSpectra.data() + channel*this->bands();
	}
	std::vector<Complex> prevOutputRotations;
	
	// Copy a segment of the spectrum to the output spectrum, shifted in frequency
	void copySegmentToNew(int segmentStart, int segmentEnd) {
		// find centre of the segment by energy-weighted average
		double binTotal = 0, energyTotal = 0;
		for (int b = segmentStart; b < segmentEnd; ++b) {
			binTotal += b*energy[b];
			energyTotal += energy[b];
		}
		double binAverage = binTotal/(energyTotal + 1e-100);
		Sample centreFreq = this->bandToFreq(binAverage);
		Sample newCentreFreq = centreFreq*freqFactor;
		int binOffset = std::round(this->freqToBand(newCentreFreq) - binAverage);

		Complex phaseShift = 1;
		if (!fixedPhase) {
			Complex phaseShiftSum = 0;
			for (int c = 0; c < this->channels; ++c) {
				Complex *spectrum = this->channelSpectrum(c);
				Complex *prevSpectrum = prevChannelSpectrum(c);
				for (int b = segmentStart; b < segmentEnd; ++b) {
					int newB = b + binOffset;
					if (newB > 0 && newB < this->bands()) {
						phaseShiftSum += prevSpectrum[newB]*std::conj(spectrum[b]);
					}
				}
			}
			Sample norm = std::norm(phaseShiftSum);
			if (norm > 0) {
				phaseShift = phaseShiftSum/std::sqrt(norm);
			}
		}

		for (int c = 0; c < this->channels; ++c) {
			Complex *spectrum = this->channelSpectrum(c);
			Complex *newSpectrum = newChannelSpectrum(c);
			for (int b = segmentStart; b < segmentEnd; ++b) {
				int newB = b + binOffset;
				if (newB > 0 && newB < this->bands()) {
					newSpectrum[newB] += spectrum[b]*phaseShift;
				}
			}
		}
	}
};

class PhaseVocoderStretch : public SpectralStretch {
public:
	PhaseVocoderStretch(bool purePhase) : purePhase(purePhase) {}

	void configure(int channels, int blockSamples, int intervalSamples, double zeroPadding=2, int maxExtraInput=0) {
		SpectralStretch::configure(channels, blockSamples, intervalSamples, zeroPadding, maxExtraInput);

		prevInputSpectra.resize(bands()*channels);
		prevOutputSpectra.resize(bands()*channels);
		outputRotations.resize(bands()*channels);
		
		prevInputRotations.resize(bands());
		prevOutputRotations.resize(bands());
		timeShiftPhases(-intervalSamples, prevOutputRotations.data());
	}
	
	void reset() {
		prevInputSpectra.assign(prevInputSpectra.size(), 0);
		prevOutputSpectra.assign(prevInputSpectra.size(), 0);
		outputRotations.assign(prevInputSpectra.size(), 0);
	}
	
	double gain = 1;
protected:
	virtual void processSpectrum(int inputIntervalSamples) {
		// Scale phases by the ratio between our input and output steps
		Sample timeFactor = inputIntervalSamples > 0 ? intervalSamples/Sample(inputIntervalSamples) : 0;

		// Shift previous input/output back with appropriate phase
		timeShiftPhases(-inputIntervalSamples, prevInputRotations.data());
		for (int c = 0; c < channels; ++c) {
			Complex *prevInputBands = prevInputSpectrum(c);
			Complex *prevOutputBands = prevOutputSpectrum(c);
			for (int b = 0; b < bands(); ++b) {
				prevInputBands[b] *= prevInputRotations[b];
				prevOutputBands[b] *= prevOutputRotations[b];
			}
		}

		for (int c = 0; c < channels; ++c) {
			Complex *currentBands = channelSpectrum(c);
			Complex *prevInputBands = prevInputSpectrum(c);
			Complex *prevOutputBands = prevOutputSpectrum(c);
			for (int b = 0; b < bands(); ++b) {
				if (inputIntervalSamples > 0) {
					Complex rotation = currentBands[b]*std::conj(prevInputBands[b]);
					Sample rotationAbs = std::abs(rotation);
					Sample phase = std::arg(rotation)*timeFactor;
					outputRotations[b] = {rotationAbs*std::cos(phase), rotationAbs*std::sin(phase)};
					prevInputBands[b] = currentBands[b];
				}

				Sample outputEnergy = std::norm(currentBands[b]);
				Complex complexPhase = prevOutputBands[b]*outputRotations[b];
				if (!purePhase) {
					Sample existingEnergy = std::min(std::norm(prevOutputBands[b]), outputEnergy);
					Sample newEnergy = outputEnergy - existingEnergy;
					complexPhase = existingEnergy*complexPhase + newEnergy*currentBands[b];
				}
				currentBands[b] = generateComplex(outputEnergy, complexPhase);
				currentBands[b] *= gain;
				prevOutputBands[b] = currentBands[b];
			}
		}
	}
private:
	bool purePhase = true;
	std::vector<Complex> prevInputSpectra, prevOutputSpectra, outputRotations;
	std::vector<Complex> prevInputRotations, prevOutputRotations;
	Complex * prevInputSpectrum(int channel) {
		return prevInputSpectra.data() + channel*this->bands();
	}
	Complex * prevOutputSpectrum(int channel) {
		return prevOutputSpectra.data() + channel*this->bands();
	}
};

class PaulStretch : public SpectralStretch {
protected:
	void processSpectrum(int) override {
		Sample gain = std::sqrt(Sample(fftSize()/intervalSamples));
		for (int c = 0; c < channels; ++c) {
			Complex *spectrum = channelSpectrum(c);
			for (int b = 0; b < bands(); ++b) {
				Sample phase = Sample(2*M_PI)*rand()/RAND_MAX;
				Complex complexPhase = {std::cos(phase), std::sin(phase)};
				spectrum[b] *= complexPhase*gain;
			}
		}
	}
};

class VasePhocoderStretch : public SpectralStretch {
public:
	VasePhocoderStretch(bool stretchStride=true) : stretchStride(stretchStride) {}

	void configure(int channels, int blockSamples, int intervalSamples, double zeroPadding=2, int maxExtraInput=0) {
		SpectralStretch::configure(channels, blockSamples, intervalSamples, zeroPadding, maxExtraInput);

		newSpectrum.resize(bands());
		centreTimeRotations.resize(bands());
		timeShiftPhases(-blockSamples*0.5, centreTimeRotations.data());
	}

	void setFreqFactor(double factor) {
		freqFactor = factor;
	}

	double gain = 1;
protected:
	virtual void processSpectrum(int inputIntervalSamples) {
		Sample timeFactor = inputIntervalSamples > 0 ? intervalSamples/Sample(inputIntervalSamples) : 0;

		for (int c = 0; c < channels; ++c) {
			Complex *spectrum = channelSpectrum(c);
			
			// Rotate so the block is centered on t=0
			// This makes interpolation more sensible, as well as the phase-changes centred
			for (int b = 0; b < bands(); ++b) {
				spectrum[b] *= centreTimeRotations[b];
			}

			for (int b = 0; b < bands(); ++b) {
				Sample inputBin = freqToBand(bandToFreq(b)/freqFactor);

				Sample energy = getEnergy(spectrum, inputBin);
				Complex rotation;
				
				// Stretch vertical phase to expand time, either by scaling the phase, or by using a longer stride
				if (stretchStride) {
					Complex bin = getBin(spectrum, inputBin);
					Complex prevBin = getBin(spectrum, inputBin - timeFactor);
					rotation = bin*std::conj(prevBin);
				} else {
					Complex bin = getBin(spectrum, inputBin);
					Complex prevBin = getBin(spectrum, inputBin - 1);
					rotation = bin*std::conj(prevBin);
					Sample phaseRotation = std::arg(rotation)*timeFactor;
					rotation = {std::cos(phaseRotation), std::sin(phaseRotation)};
				}
				
				Complex phase = (b > 0) ? rotation*newSpectrum[b - 1] : inputBin;
				newSpectrum[b] = generateComplex(energy, phase);
			}
			
			// Rotate back again for output
			for (int b = 0; b < bands(); ++b) {
				spectrum[b] = newSpectrum[b]*std::conj(centreTimeRotations[b]);
				spectrum[b] *= gain;
			}
		}
	}
private:
	bool stretchStride;
	std::vector<Complex> newSpectrum;
	double freqFactor;

	Sample getEnergy(Complex *spectrum, double index) {
		int indexFloor = int(std::floor(index));
		Sample fractional = index - indexFloor;

		Complex lowBin = 0, highBin = 0;
		if (indexFloor >= 0 && indexFloor < bands()) {
			lowBin = spectrum[indexFloor];
		} else if (indexFloor < 0 && indexFloor >= -bands()){
			lowBin = std::conj(spectrum[1 - indexFloor]);
		}
		if (indexFloor + 1 >= 0 && indexFloor + 1 < bands()) {
			highBin = spectrum[indexFloor + 1];
		} else if (indexFloor + 1< 0 && indexFloor + 1>= -bands()){
			highBin = std::conj(spectrum[-indexFloor]);
		}
		Sample lowEnergy = std::norm(lowBin), highEnergy = std::norm(highBin);
		return lowEnergy + (highEnergy - lowEnergy)*fractional;
	}
	
	Complex getBin(Complex *spectrum, double index) {
		int indexFloor = int(std::floor(index));
		Sample fractional = index - indexFloor;

		Complex lowBin = 0, highBin = 0;
		if (indexFloor >= 0 && indexFloor < bands()) {
			lowBin = spectrum[indexFloor];
		} else if (indexFloor < 0 && indexFloor >= -bands()){
			lowBin = std::conj(spectrum[1 - indexFloor]);
		}
		if (indexFloor + 1 >= 0 && indexFloor + 1 < bands()) {
			highBin = spectrum[indexFloor + 1];
		} else if (indexFloor + 1< 0 && indexFloor + 1>= -bands()){
			highBin = std::conj(spectrum[-indexFloor]);
		}
		return lowBin + (highBin - lowBin)*fractional;
	}
	
	std::vector<Complex> centreTimeRotations;
};

class HybridPhaseStretch : public SpectralStretch {
public:
	HybridPhaseStretch(bool multipleTimeObservations=true, double pitchWeight=1, double timeWeight=2, double channelWeight=1, double maxWeight=1) : multipleTimeObservations(multipleTimeObservations), pitchWeight(pitchWeight), timeWeight(timeWeight), channelWeight(channelWeight), maxWeight(maxWeight) {}

	void configure(int channels, int blockSamples, int intervalSamples, double zeroPadding=2, int maxExtraInput=0) {
		SpectralStretch::configure(channels, blockSamples, intervalSamples, zeroPadding, maxExtraInput);

		prevInputSpectra.resize(bands()*channels);
		prevOutputSpectra.resize(bands()*channels);
		newOutputSpectra.resize(bands()*channels);
		
		horizontalRotations.resize(bands()*channels);

		centreTimeRotations.resize(bands());
		timeShiftPhases(-blockSamples*0.5, centreTimeRotations.data());
		prevInputRotations.resize(bands());
		prevOutputRotations.resize(bands());
		timeShiftPhases(-intervalSamples, prevOutputRotations.data());
	}
	
	void reset() {
		prevInputSpectra.assign(prevInputSpectra.size(), 0);
		prevOutputSpectra.assign(prevOutputSpectra.size(), 0);
		horizontalRotations.assign(horizontalRotations.size(), 0);
	}

	void setFreqFactor(double factor) {
		freqFactor = factor;
	}
protected:
	virtual void processSpectrum(int inputIntervalSamples) {
		Sample timeFactor = inputIntervalSamples > 0 ? intervalSamples/Sample(inputIntervalSamples) : 0;

		// Shift input and previous input/output with appropriate phase
		timeShiftPhases(-inputIntervalSamples, prevInputRotations.data());
		for (int c = 0; c < channels; ++c) {
			Complex *currentBands = channelSpectrum(c);
			Complex *prevInputBands = prevInputSpectrum(c);
			Complex *prevOutputBands = prevOutputSpectrum(c);
			for (int b = 0; b < bands(); ++b) {
				currentBands[b] *= centreTimeRotations[b]; // Rotate so the block is centered on t=0
				prevInputBands[b] *= prevInputRotations[b];
				prevOutputBands[b] *= prevOutputRotations[b];
			}
		}

		Complex *newSpectrum0 = newOutputSpectrum(0);
		for (int c = 0; c < channels; ++c) {
			Complex *currentBands = channelSpectrum(c);
			Complex *prevInputBands = prevInputSpectrum(c);
			Complex *prevOutputBands = prevOutputSpectrum(c);
			Complex *newSpectrum = newOutputSpectrum(c);
			
			Complex *horizontalRotations = channelHorizontalRotations(c);

			for (int b = 0; b < bands(); ++b) {
				Sample inputBin = freqToBand(bandToFreq(b)/freqFactor);

				Sample energy = getEnergy(currentBands, inputBin);
				energy /= freqFactor; // Keep total energy constant
				if (timeFactor > 1) energy *= std::sqrt(timeFactor);
				Complex phase = 0;
				Complex bin = getBin(currentBands, inputBin);
				
				Complex maxPrediction = 0;
				Sample maxNorm = 0;

				// Vase phocoder prediction
				Complex verticalPrediction1;
				{
					Complex prevBin = getBin(currentBands, inputBin - timeFactor);
					Complex prevOutput = (b >= 1) ? newSpectrum[b - 1] : 0;
					Complex verticalRotation = bin*std::conj(prevBin);
					verticalPrediction1 = prevOutput*verticalRotation;
					if (std::norm(verticalPrediction1) > maxNorm) {
						maxPrediction = verticalPrediction1;
						maxNorm = std::norm(verticalPrediction1);
					}
				}
				Sample weight = multipleTimeObservations ? 0.2 : 1;
				phase += verticalPrediction1*timeWeight*weight;
				if (multipleTimeObservations) {
					auto addTimeObservation = [&](int steps) {
						Complex verticalPrediction;
						{
							Complex prevBin = getBin(currentBands, inputBin - timeFactor*steps);
							Complex prevOutput = (b >= steps) ? newSpectrum[b - steps] : 0;
							Complex verticalRotation = bin*std::conj(prevBin);
							verticalPrediction = prevOutput*verticalRotation;
							if (std::norm(verticalPrediction) > maxNorm) {
								maxPrediction = verticalPrediction;
								maxNorm = std::norm(verticalPrediction);
							}
						}
						phase += verticalPrediction*timeWeight*weight;
					};
					addTimeObservation(2);
					addTimeObservation(4);
					addTimeObservation(8);
					addTimeObservation(16);
				}
				
				// Phase-vocoder (horizontal) predictions
				if (inputIntervalSamples > 0) {
					Complex prevBin = getBin(prevInputBands, inputBin);
					
					Complex rotation = bin*std::conj(prevBin);
					// Scale phase-rotation from input time-diff to output time-diff, and also by frequency
					rotation = scaleAngle(rotation, timeFactor*freqFactor);
					horizontalRotations[b] = rotation;
				}
				Complex horizontalPrediction = prevOutputBands[b]*horizontalRotations[b];
				if (std::norm(horizontalPrediction) > maxNorm) {
					maxPrediction = horizontalPrediction;
					maxNorm = std::norm(horizontalPrediction);
				}
				phase += horizontalPrediction*pitchWeight;

				// Channel predictions
				if (c > 0) {
					Complex prevBin = getBin(channelSpectrum(0), inputBin);
					Complex channelRotation = bin*std::conj(prevBin);
					Complex channelPrediction = newSpectrum0[b]*channelRotation;
					if (std::norm(channelPrediction) > maxNorm) {
						maxPrediction = horizontalPrediction;
						maxNorm = std::norm(channelPrediction);
					}
					phase += channelPrediction*channelWeight;
				}
				phase += maxPrediction*maxWeight;

				newSpectrum[b] = generateComplex(energy, phase);
			}

			for (int b = 0; b < bands(); ++b) {
				prevInputBands[b] = currentBands[b];
			}
		}
		for (int c = 0; c < channels; ++c) {
			Complex *currentBands = channelSpectrum(c);
			Complex *prevOutputBands = prevOutputSpectrum(c);
			Complex *newSpectrum = newOutputSpectrum(c);
			for (int b = 0; b < bands(); ++b) {
				prevOutputBands[b] = newSpectrum[b];
				currentBands[b] = newSpectrum[b]*std::conj(centreTimeRotations[b]);
			}
		}
	}
private:
	bool multipleTimeObservations;
	double freqFactor = 1;
	double pitchWeight, timeWeight, channelWeight, maxWeight;
	std::vector<Complex> newOutputSpectra;
	std::vector<Complex> prevInputSpectra, prevOutputSpectra, outputRotations;
	std::vector<Complex> horizontalRotations; // stored in case we get handed the same input twice
	Complex * prevInputSpectrum(int channel) {
		return prevInputSpectra.data() + channel*this->bands();
	}
	Complex * prevOutputSpectrum(int channel) {
		return prevOutputSpectra.data() + channel*this->bands();
	}
	Complex * newOutputSpectrum(int channel) {
		return newOutputSpectra.data() + channel*this->bands();
	}
	Complex * channelHorizontalRotations(int channel) {
		return horizontalRotations.data() + channel*this->bands();
	}
	std::vector<Complex> centreTimeRotations, prevInputRotations, prevOutputRotations;

	Sample getEnergy(Complex *spectrum, double index) {
		int indexFloor = int(std::floor(index));
		Sample fractional = index - indexFloor;

		Complex lowBin = 0, highBin = 0;
		if (indexFloor >= 0 && indexFloor < bands()) {
			lowBin = spectrum[indexFloor];
		} else if (indexFloor < 0 && indexFloor >= -bands()){
			lowBin = std::conj(spectrum[1 - indexFloor]);
		}
		if (indexFloor + 1 >= 0 && indexFloor + 1 < bands()) {
			highBin = spectrum[indexFloor + 1];
		} else if (indexFloor + 1< 0 && indexFloor + 1>= -bands()){
			highBin = std::conj(spectrum[-indexFloor]);
		}
		Sample lowEnergy = std::norm(lowBin), highEnergy = std::norm(highBin);
		return lowEnergy + (highEnergy - lowEnergy)*fractional;
	}
	
	Complex getBin(Complex *spectrum, double index) {
		int indexFloor = int(std::floor(index));
		Sample fractional = index - indexFloor;

		Complex lowBin = 0, highBin = 0;
		if (indexFloor >= 0 && indexFloor < bands()) {
			lowBin = spectrum[indexFloor];
		} else if (indexFloor < 0 && indexFloor >= -bands()){
			lowBin = std::conj(spectrum[1 - indexFloor]);
		}
		if (indexFloor + 1 >= 0 && indexFloor + 1 < bands()) {
			highBin = spectrum[indexFloor + 1];
		} else if (indexFloor + 1< 0 && indexFloor + 1>= -bands()){
			highBin = std::conj(spectrum[-indexFloor]);
		}
		return lowBin + (highBin - lowBin)*fractional;
	}
	
	Complex scaleAngle(Complex rotation, Sample multiplier) {
		Sample rotationAbs = std::abs(rotation);
		Sample phase = std::arg(rotation)*multiplier;
		return {rotationAbs*std::cos(phase), rotationAbs*std::sin(phase)};
	}
};

#endif // include guard
