#include "shift-stretch.h"

#include "util/simple-args.h"
#include "util/wav.h"

class WavCmd {
	std::string inputFile, outputFile;
	Wav inputWav, outputWav;
	
	double timeFactor = 1, freqFactor = 1;
	double blockMs = 80, overlapFactor = 4, searchMs = 10;
	bool trimLatency = false;

	template<class Processor>
	void processBlocks(Processor &processor, double stretchFactor) {
		int channels = inputWav.channels;
		int blockSize = 256;
		
		std::vector<std::vector<double>> inputBuffers(channels), outputBuffers(channels);
		std::vector<double *> inputPointers(channels), outputPointers(channels);
		for (auto &b : outputBuffers) b.resize(blockSize);
		
		outputWav.channels = inputWav.channels;
		int inputOffset = 0, outputOffset = 0;
		int inputLength = int(inputWav.length());
		int totalLatency = std::round(processor.inputLatency()*stretchFactor + processor.outputLatency());
		int outputLength = inputWav.length()*stretchFactor;
		while (outputOffset < outputLength + totalLatency*2) {
			// For `blockSize` output samples, how many input samples should we have?
			int inputSamples = processor.samplesForOutput(blockSize);
			// Make sure our input buffers are large enough
			if (inputSamples > int(inputBuffers[0].size())) {
				for (auto &b : inputBuffers) b.resize(inputSamples);
			}
			// Fill them up
			for (int c = 0; c < channels; ++c) {
				for (int i = 0; i < inputSamples; ++i) {
					// Fill input either from WAV, or with 0
					if (inputOffset + i < inputLength) {
						inputBuffers[c][i] = inputWav[c][inputOffset + i];
					} else {
						inputBuffers[c][i] = 0;
					}
				}
				inputPointers[c] = inputBuffers[c].data();
				outputPointers[c] = outputBuffers[c].data();
			}
			
			processor.process(inputPointers.data(), inputSamples, outputPointers.data(), blockSize);
			
			outputWav.samples.resize((outputOffset + blockSize)*channels);
			for (int c = 0; c < channels; ++c) {
				for (int i = 0; i < blockSize; ++i) {
					outputWav[c][outputOffset + i] = outputBuffers[c][i];
				}
			}
			
			inputOffset += inputSamples;
			outputOffset += blockSize;
		}
		if (trimLatency) {
			for (unsigned c = 0; c < outputWav.channels; ++c) {
				for (int i = 0; i < totalLatency; ++i) {
					outputWav[c][i + totalLatency] -= outputWav[c][totalLatency - 1 - i];
				}
			}
			outputWav.samples.erase(outputWav.samples.begin(), outputWav.samples.begin() + outputWav.channels*totalLatency);

			int end = outputLength;
			for (unsigned c = 0; c < outputWav.channels; ++c) {
				for (int i = 0; i < totalLatency; ++i) {
					outputWav[c][end - 1 - i] -= outputWav[c][end + i];
				}
			}
			outputWav.samples.erase(outputWav.samples.begin() + end*outputWav.channels, outputWav.samples.end());
		}
	}
public:
	enum class Mode{resample, overlapAdd, wsola, spectralCut, phaseVocoder, paulStretch, vasePhocoder, hybridPhase};
	Mode mode;

	WavCmd(SimpleArgs &args, Mode mode) : mode(mode) {
		// Collect mode-appropriate command-line arguments
		inputFile = args.arg<std::string>("input.wav", "Input WAV (must be 16-bit)");
		outputFile = args.arg<std::string>("output.wav", "Output WAV");
		freqFactor = args.flag<double>("freq", "Freq-scaling factor (e.g. 2 is +1 octave)", freqFactor);
		if (mode != Mode::resample) {
			timeFactor = args.flag<double>("time", "Time-scaling factor (e.g. 2 is twice as slow)", timeFactor);
		}
		bool fixedPhase = false, purePV = false, stretchPhase = false, singleTimeObservation = false;
		double pitchWeight = 1, timeWeight = 2, channelWeight = 1, maxWeight = 1;
		double zeroPadding = 2;
		double gain = 1;
		if (mode != Mode::resample) {
			if (mode == Mode::paulStretch) blockMs = 120;
			if (mode == Mode::phaseVocoder) blockMs = 120;
			if (mode == Mode::hybridPhase) blockMs = 120;
			blockMs = args.flag<double>("block", "Block length (ms)", blockMs);
			if (mode == Mode::overlapAdd) overlapFactor = 1.5;
			if (mode == Mode::wsola) overlapFactor = 2;
			overlapFactor = args.flag<double>("overlap", "Overlap factor for blocks", overlapFactor);
			if (mode == Mode::wsola) searchMs = args.flag<double>("search", "WSOLA search duration (ms)", searchMs);
			if (mode == Mode::spectralCut) fixedPhase = args.hasFlag("fixed-phase", "Don't phase-shift segments to match the previous block");
			if (mode == Mode::phaseVocoder) {
				purePV = args.hasFlag("pure", "Pure phase-vocoder (no transient/energy stuff)");
			}
			if (mode == Mode::phaseVocoder || mode == Mode::vasePhocoder) {
				gain = args.flag<double>("gain", "Gain factor (amplitude), default 1", 1);
			}
			if (mode == Mode::vasePhocoder) stretchPhase = args.hasFlag("stretch-phase", "Stretch phase (trigonmetrically) instead of using longer strides");
			if (mode == Mode::spectralCut || mode == Mode::phaseVocoder || mode == Mode::paulStretch || mode == Mode::vasePhocoder || mode == Mode::hybridPhase) {
				zeroPadding = args.flag<double>("zero-padding", "Zero-padding factor for spectral processing", zeroPadding);
			}
			if (mode == Mode::hybridPhase) {
				pitchWeight = args.flag<double>("pitch-weight", "Weighting for pitch-based phase information", pitchWeight);
				timeWeight = args.flag<double>("time-weight", "Weighting for timing-based phase information", timeWeight);
				channelWeight = args.flag<double>("channel-weight", "Weighting for inter-channel phase information", channelWeight);
				maxWeight = args.flag<double>("max-weight", "Additional weighting for the strongest phase predictor", maxWeight);
				if (mode == Mode::hybridPhase) {
					singleTimeObservation = args.hasFlag("single-vertical", "Use a single vertical (timing) observation, instead of combining multiple");
				}
			}
			trimLatency = args.hasFlag("trim", "Trim edges to remove processing latency");
		}
		args.errorExit();

		if (!inputWav.read(inputFile)) args.errorExit(inputWav.result.reason);
		outputWav.channels = inputWav.channels;
		outputWav.sampleRate = inputWav.sampleRate;

		int blockSamples = int(blockMs*0.001*inputWav.sampleRate + 0.5);
		int intervalSamples = int(blockSamples/overlapFactor);
		int searchSamples = int(searchMs*0.001*inputWav.sampleRate + 0.5);

		if (mode == Mode::resample) {
			outputWav = inputWav; // Just copy the waveform
			outputWav.sampleRate *= freqFactor;
		} else if (mode == Mode::overlapAdd) {
			OverlapAddStretch stretch;
			stretch.configure(inputWav.channels, blockSamples, intervalSamples);
			stretch.setTimeFactor(timeFactor*freqFactor);
			processBlocks(stretch, timeFactor*freqFactor);
			outputWav.sampleRate *= freqFactor;
		} else if (mode == Mode::wsola) {
			WsolaStretch stretch;
			stretch.configure(inputWav.channels, blockSamples, intervalSamples, searchSamples);
			stretch.setTimeFactor(timeFactor*freqFactor);
			processBlocks(stretch, timeFactor*freqFactor);
			outputWav.sampleRate *= freqFactor;
		} else if (mode == Mode::spectralCut) {
			SpectralCutStretch stretch(fixedPhase);
			stretch.configure(inputWav.channels, blockSamples, intervalSamples);
			stretch.setTimeFactor(timeFactor);
			stretch.setFreqFactor(freqFactor);
			processBlocks(stretch, timeFactor);
		} else if (mode == Mode::phaseVocoder) {
			PhaseVocoderStretch stretch(purePV);
			stretch.configure(inputWav.channels, blockSamples, intervalSamples, zeroPadding);
			stretch.gain = gain;
			stretch.setTimeFactor(timeFactor*freqFactor);
			processBlocks(stretch, timeFactor*freqFactor);
			outputWav.sampleRate *= freqFactor;
		} else if (mode == Mode::paulStretch) {
			PaulStretch stretch;
			stretch.configure(inputWav.channels, blockSamples, intervalSamples, zeroPadding);
			stretch.setTimeFactor(timeFactor*freqFactor);
			processBlocks(stretch, timeFactor*freqFactor);
			outputWav.sampleRate *= freqFactor;
		} else if (mode == Mode::vasePhocoder) {
			VasePhocoderStretch stretch(stretchPhase);
			stretch.configure(inputWav.channels, blockSamples, intervalSamples, zeroPadding);
			stretch.gain = gain;
			stretch.setTimeFactor(timeFactor);
			stretch.setFreqFactor(freqFactor);
			processBlocks(stretch, timeFactor);
		} else if (mode == Mode::hybridPhase) {
			HybridPhaseStretch stretch(!singleTimeObservation, pitchWeight, timeWeight, channelWeight, maxWeight);
			stretch.configure(inputWav.channels, blockSamples, intervalSamples, zeroPadding);
			stretch.setTimeFactor(timeFactor);
			stretch.setFreqFactor(freqFactor);
			processBlocks(stretch, timeFactor);
		}

		if (!outputWav.write(outputFile)) args.errorExit(outputWav.result.reason);
	}
};

int main(int argc, char **argv) {
	SimpleArgs args(argc, argv);

	if (args.command("rate", "Changes the sampling-rate, altering pitch and time together")) {
		WavCmd cmd(args, WavCmd::Mode::resample);
	} else if (args.command("overlap-add", "Overlap-add time-stretching")) {
		WavCmd cmd(args, WavCmd::Mode::overlapAdd);
	} else if (args.command("wsola", "WSOLA-ish time-stretch (nudges blocks around to get a good cross-fade)")) {
		WavCmd cmd(args, WavCmd::Mode::wsola);
	} else if (args.command("spectral-cut", "Pitch-shift by moving spectral segments up or down")) {
		WavCmd cmd(args, WavCmd::Mode::spectralCut);
	} else if (args.command("phase-vocoder", "Phase-vocoder time-stretching")) {
		WavCmd cmd(args, WavCmd::Mode::phaseVocoder);
	} else if (args.command("paul-stretch", "PaulStretch (random phase)")) {
		WavCmd cmd(args, WavCmd::Mode::paulStretch);
	} else if (args.command("vase-phocoder", "Time-synchronised phase vocoder")) {
		WavCmd cmd(args, WavCmd::Mode::vasePhocoder);
	} else if (args.command("hybrid-phase", "Blends several phase sources (phase-vocoder, vase-phocoder, inter-channel)")) {
		WavCmd cmd(args, WavCmd::Mode::hybridPhase);
	}
	args.errorCommand();
}
