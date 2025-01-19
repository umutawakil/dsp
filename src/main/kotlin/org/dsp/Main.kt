package org.dsp

import org.dsp.analysis.DiscreteFourierTransform.Companion.synthesizePartialDft
import org.dsp.config.Constants
import org.dsp.filters.FilterUtils
import org.dsp.modulation.WaveformEffect.Companion.compressSignal
import org.dsp.modulation.WaveformEffect.Companion.normalize
import org.dsp.signals.SignalGenerator.Companion.sineByWaveLength
import org.dsp.signals.SignalGenerator.Companion.synthesizeSineByVibratoHalfWavelengths
import org.dsp.tools.add

import org.dsp.wave.WaveUtils
import kotlin.math.*

fun main() {
    val resourceDirectory  = "/Users/umutawakil/Documents/Git/dsp/src/main/resources"
    val waveData = WaveUtils.getFileData(
        baseHalfPeriod    = 27.0, //TODO: This is kludgy as hell
        resourceDirectory = resourceDirectory,
        fileName          = "normalized.wav"//"harmonic-2.wav"//"test-normalized.wav"//"harmonic-2.wav"//"harmonic-2.wav"//"normalized.wav"//"harmonic-2.wav"//"normalized.wav"//"filter-out-1.wav"
    )

    val baseFrequency         = 444.0
    val baseWavelength        = (Constants.SAMPLE_RATE/baseFrequency)
    val targetHarmonic        = 2
    val targetFrequency       = targetHarmonic * baseFrequency
    val targetWaveLength      = baseWavelength / targetHarmonic
    val length                = Constants.SAMPLE_RATE//waveData.rawFileData.size//Constants.SAMPLE_RATE + 1
    val impulseResponseLength = length//* 2
    val sBandwidth            = baseFrequency//44.0//baseFrequency//44.0//baseFrequency/4 //(baseFrequency/4) May be all we need or perhaps even (baseFrequency/k) or some function of this.
    val frequencyDistance     = 1.0
    val estimatedHalfPeriods  = (length / (targetWaveLength/2)).toInt()

    /*val fOut1 = FilterUtils.bandpassFIR(
        input                 = waveData.rawFileData,//.map { it * 1000.0},
        impulseResponseLength = impulseResponseLength + 1,
        minFreqHz             = targetFrequency - (sBandwidth/2),
        maxFreqHz             = targetFrequency + (sBandwidth/2)
    )

    WaveUtils.writeSamplesToFile(
        fileName = "$resourceDirectory/output.wav",
        input    =  normalize(scale = 5000.0, input = fOut1)//postFilterInnerFormant.add(outerFormant))//.add(b).add(resc))
    )
    return */



    /*val out = FilterUtils.bandpassFIR(
        input                 = st,
        impulseResponseLength = impulseResponseLength + 1,
        minFreqHz             = targetFrequency - (sBandwidth/2) - attackAndDecayFrequencyChange,
        maxFreqHz             = targetFrequency + (sBandwidth/2) + attackAndDecayFrequencyChange
    )*/

    /*WaveUtils.writeSamplesToFile(
        fileName = "$resourceDirectory/output.wav",
        input    =  normalize(scale = 5000.0, input = st)
    )
    return */

    /*val dft = partialDft(
        x                 = normalize(scale = 5000.0, input = waveData.rawFileData),
        baseFrequency     = targetFrequency,
        sampleRate        = Constants.SAMPLE_RATE,
        bandwidth         = sBandwidth.toInt(),
        frequencyDistance = frequencyDistance
    )

    val idft = inversePartialDftPolar(
        magnitude         = dft.magnitude,
        phaseInPercent    = dft.phaseInPercent,
        baseFrequency     = targetFrequency,
        bandwidth         = sBandwidth.toInt(),
        frequencyDistance = frequencyDistance,
        sampleRate        = Constants.SAMPLE_RATE,
        length            = length
    )
    val wls = getWaves(data = idft).map { it.size.toDouble() } */

    /*val formant: List<Double> = synthesizeFormant(
            targetFrequency       = targetFrequency,
            bandwidth             = sBandwidth.toInt(),
            impulseResponseLength = impulseResponseLength,
            length                = length
    )*/

    val ocarina: List<Double> = synthesizeOcarina(
        frequency             = targetFrequency,
        bandwidth             = sBandwidth.toInt(),
        impulseResponseLength = impulseResponseLength,
        length                = length
    )

    WaveUtils.writeSamplesToFile(
        fileName = "$resourceDirectory/output.wav",
        input    =  normalize(scale = 5000.0, input = ocarina)
    )
}

fun synthesizeOcarina(
    frequency: Double,
    bandwidth: Int,
    impulseResponseLength: Int,
    length: Int
): List<Double> {
    val ocarinaSpectrum = ocarinaDft(baseFrequency = frequency)
    val lengthBuffer    = length * 2
    val buffer: MutableList<Double> = MutableList(lengthBuffer) { 0.0 }
    var maxLengthFormant = 0

    val formant: List<Double> = synthesizeFormant(
        targetFrequency       = frequency,
        bandwidth             = bandwidth,
        impulseResponseLength = impulseResponseLength,
        length                = length
    )

    /*for(i in ocarinaSpectrum.indices) {
        val currentBin = ocarinaSpectrum[i]
    }*/

   // return buffer
    return formant
}

fun synthesizeFormant(targetFrequency: Double, bandwidth: Int, impulseResponseLength: Int, length: Int) : List<Double> {
    val innerOuterFormantRatio        = 2 // Innerformant(central)/Outer(secondary)
    val envelopeChangePercent         = 0.05/2.0
    val attackAndDecayFrequencyChange = envelopeChangePercent * targetFrequency
    println("Attack/Decay Change: $attackAndDecayFrequencyChange")

    val vibratoSignal = generateAtomicVibratoSignal(
        baseWavelength = Constants.SAMPLE_RATE / targetFrequency,
        envelopeSlope  = envelopeChangePercent,
        timeInSamples  = length
    )
    val innerFormant: List<Double> = synthesizeSineByVibratoHalfWavelengths(halfWavelengths = vibratoSignal)

    val baseScalingAmplitude = 5000.0
    val postFilterInnerFormant = normalize(
        scale = innerOuterFormantRatio * baseScalingAmplitude,
        input = FilterUtils.bandpassFIR(
            input                 = innerFormant,
            impulseResponseLength = impulseResponseLength + if(impulseResponseLength % 2 == 0) { 1 } else { 0},
            minFreqHz             = targetFrequency - (bandwidth/2),
            maxFreqHz             = targetFrequency + (bandwidth/2)
        )
    )

    val outerFormant = normalize(
        scale = baseScalingAmplitude,
        input = synthesizeOuterFormant(
            amplitude     = 1.0,
            coefficient   = -1.0/64, //(16*decayFactor),
            baseFrequency = targetFrequency,
            bandWidth     = bandwidth,
            length        = length
        )
    )
    return postFilterInnerFormant.add(outerFormant)
}

fun generateAtomicVibratoSignal(
    baseWavelength: Double,
    envelopeSlope: Double,
    timeInSamples: Int
) : List<Double> {
    val vibratoRange = 0.5
    val options: List<Double> = listOf(
        (baseWavelength/2),
        (baseWavelength/2) + vibratoRange
    )

    val o: MutableList<Double> = mutableListOf()
    var count                  = 0.0
    while(count < timeInSamples) {
        val chosenValue = options[options.indices.random()]
        o.add(chosenValue)
        count += chosenValue
    }

    /*** Add transients to envelope (attack and decay) ***/
    val attackDecayLength = (Constants.SAMPLE_RATE * 0.2)
    var sum = 0.0
    var p = 0
    while (sum < attackDecayLength) {
        sum += o[p]
        p++
    }
    /** 5% worked well on the fundamental **/
    val initialPitchChange = baseWavelength * envelopeSlope
    val slope              = initialPitchChange / p

    val a: MutableList<Double> = mutableListOf()
    /** Attack **/
    for(i in 0 until p) {
        val line = initialPitchChange - (slope*i)
        a.add(line)
        o[i] += line
    }

    val b: MutableList<Double> = mutableListOf()
    /** Decay **/
    for(i in 0 until p) {
        val line = (slope*i)
        b.add(line)
        o[o.size - p + i] += line
    }

    return o
}

/**
 * ModulationIndex = DesiredFrequencyDeviation / Desired Modulation Frequency
 *
 * Frequency Deviation: Every fm seconds the carrier frequency fc
 * changes by +/- (Modulation Index * Modulation frequency)
 *  c(t) = Asin(fc*i + Isin(fm*i))
 *  **/
fun fmSynthesis(fc: Double, modulationIndex: Double, modSignal: List<Double>, length: Int) : List<Double>{
    val o: MutableList<Double> = mutableListOf()
    for(i in 0 until length) {
        val modulator = modulationIndex * modSignal[i]
        val value     = sin(((2*Math.PI*i) / (Constants.SAMPLE_RATE / fc)) + modulator )
        o.add(value)
    }
    return o
}

fun ocarinaDft(baseFrequency: Double) : List<FreqPair> {
    val numHarmonics = (Constants.SAMPLE_RATE / baseFrequency).toInt()
    val keyHarmonics: List<FreqPair> = listOf(
        FreqPair(amp = 0.0, freq = baseFrequency*0),
        FreqPair(amp = 168.0, freq = baseFrequency*2),
        FreqPair(amp = 136.0, freq = baseFrequency*3),
        FreqPair(amp = 40.0,  freq = baseFrequency*4),
        FreqPair(amp = 46.0,  freq = baseFrequency*5),
        FreqPair(amp = 53.0,  freq = baseFrequency*6),
        FreqPair(amp = 20.0,  freq = baseFrequency*7),
        FreqPair(amp = 40.0,  freq = baseFrequency*8),
        FreqPair(amp = 7.0,   freq = baseFrequency*9),
        FreqPair(amp = 17.0,  freq = baseFrequency*10)
    )

    val upperHarmonics: MutableList<FreqPair> = mutableListOf()
    var nextUpperHarmonic = keyHarmonics.size + 1
    val initialAmp = 1.0
    for(i in nextUpperHarmonic  until numHarmonics) {
        val currentAmplitude = initialAmp - (initialAmp/(numHarmonics)*i)
        upperHarmonics.add(FreqPair(amp = currentAmplitude, freq = i*baseFrequency))
    }

    return (keyHarmonics + upperHarmonics)//.subList(0, 20)
}

fun ocarinaSynth(baseFrequency: Double, length: Int) : List<Double> {
    val baseWavelength = Constants.SAMPLE_RATE / baseFrequency
    val numHarmonics    = ((baseWavelength / 2) + 1).toInt()
    //TODO: Need a DFT list of the first 10 harmonics (set fundamental to zero)
    val keyHarmonics: List<FreqPair> = listOf(
        FreqPair(amp = 0.0, freq = baseFrequency*0),
        FreqPair(amp = 168.0, freq = baseFrequency*2),
        FreqPair(amp = 136.0, freq = baseFrequency*3),
        FreqPair(amp = 40.0,  freq = baseFrequency*4),
        FreqPair(amp = 46.0,  freq = baseFrequency*5),
        FreqPair(amp = 53.0,  freq = baseFrequency*6),
        FreqPair(amp = 20.0,  freq = baseFrequency*7),
        FreqPair(amp = 40.0,  freq = baseFrequency*8),
        FreqPair(amp = 7.0,   freq = baseFrequency*9),
        FreqPair(amp = 17.0,  freq = baseFrequency*10)
    )

    val upperHarmonics: MutableList<FreqPair> = mutableListOf()
    var nextUpperHarmonic = keyHarmonics.size + 1
    for(i in nextUpperHarmonic  until numHarmonics) {
        val currentAmplitude = 6.0 - (5.0/(numHarmonics)*i)
        upperHarmonics.add(FreqPair(amp = currentAmplitude, freq = i*baseFrequency))
    }

    println("\r\nkey harmonics:")
    for(i in keyHarmonics.indices) {
        println("A: ${keyHarmonics[i].amp}, F: ${keyHarmonics[i].freq}")
    }
    println("\r\nupper harmonics:")
    for(i in upperHarmonics.indices) {
        println("$i, A: ${upperHarmonics[i].amp}, F: ${upperHarmonics[i].freq}")
    }

    //TODO: Loop through the DFT and generate corresponding formant and add it to the main buffer (called underWaveform)
    val frequencyBinData = (keyHarmonics + upperHarmonics).subList(0, 20)
    //var underWaveform: List<Double> = List(length) {0.0}
    var k = 0

    /**TODO: Is it worth the performance gains to do this procedurally? **/
    val underWaveform: List<Double> = frequencyBinData.fold((0..length).map {0.0}) { acc: List<Double>, currentBin ->
        acc.zip(
            synthesizeOuterFormant(
                amplitude     = currentBin.amp,
                coefficient   = -(1.0/16),//-(1.0/32),
                bandWidth     = 444,
                baseFrequency = currentBin.freq,
                length        = length
            )
        ) { a: Double, b: Double ->  a + b }


       /* println("Formant: $k, of ${frequencyBinData.size}")
        k++
        val formant: List<Double> = createFormant(
            amplitude     = currentBin.amp,
            coefficient   = -(1.0/16),//-(1.0/32),
            formantSize   = 444,
            baseFrequency = currentBin.freq,
            length        = length
        )
        return acc.zip(formant) { a: Double, b: Double ->  a + b }*/
    }
    println("Under-waveform: ${underWaveform.size}")

    /*for(currentBin in frequencyBinData) {
        println("Formant: $k, of ${frequencyBinData.size}")
        val currentFormant = formantToSpectrum(
            amplitude     = currentBin.amp,
            coefficient   = -(1.0/16),//-(1.0/32),
            formantSize   = 444,
            baseFrequency = currentBin.freq,
            length        = length,
            iterations    = 1
        )
        underWaveform = underWaveform.zip(currentFormant) { a, b -> a + b }
        k++
    }*/

    val normalizedUnderWaveform = normalize(scale = 600.0, input = underWaveform)
    val uavg = normalizedUnderWaveform.map {abs(it)}.average()
    println("UAVG: $uavg")
    //val vib = Vibrato.genericVibratoSignal(baseHalfWavelength = baseWavelength/2.0, durationInSamples = length)
    //listData(data = vib)
    val fundamental = sineByWaveLength(amplitude = 5000.0, waveLength = baseWavelength, size = length)
   // val result =  fundamental.zip(normalizedUnderWaveform) { a, b -> a + b }

    val result = compressSignal(targetAmplitude = 5000.0, baseWavelength = baseWavelength, input = normalizedUnderWaveform)

    println("Result: ${result.size}")

    return result
}

fun synthesizeOuterFormant(
    amplitude: Double,
    coefficient: Double,
    bandWidth: Int,
    baseFrequency: Double,
    length: Int,
    distanceBetweenPoints: Double = 1.0
): List<Double> {
    val spectrumAmplitudeEnvelope = createOuterFormantEnvelope(
        amplitude   = amplitude,
        coefficient = coefficient,
        length      = bandWidth
    )
    return synthesizePartialDft(
        magnitude         = spectrumAmplitudeEnvelope,
        baseFrequency     = baseFrequency,
        bandwidth         = bandWidth,
        frequencyDistance = distanceBetweenPoints,
        length            = length,
        sampleRate        = Constants.SAMPLE_RATE
    )
}

//TODO: Why length + 1
private fun createOuterFormantEnvelope(amplitude: Double, coefficient: Double, length: Int) : List<Double> {
    val functionalLength = (length / 2) + 1
    val signal = expF(amplitude = amplitude, coefficient = coefficient, length = functionalLength )
    val left   = signal.reversed().subList(0, signal.size - 1).map { it }
    val middle = 0.00//amplitude //At the moment this is left empty because the central frequency is created and modulated elsewhere
    val right  = signal.subList(1, signal.size)

    return left + middle + right
}

private fun expF(amplitude: Double, coefficient: Double,  length: Int) : List<Double> {
    val o: MutableList<Double> = mutableListOf()
    for(i in 0 until length) {
        val value = amplitude* Math.E.pow(coefficient * i.toDouble())
        o.add(value)
    }
    return o
}

class FreqPair(val amp: Double, val freq: Double)
