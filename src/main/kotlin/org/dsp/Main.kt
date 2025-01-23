package org.dsp

import org.dsp.analysis.DiscreteFourierTransform
import org.dsp.analysis.DiscreteFourierTransform.Companion.synthesizeInharmonicDft
import org.dsp.analysis.WaveformAnalyzer.Companion.getHistogramStats
import org.dsp.analysis.WaveformAnalyzer.Companion.getWaves
import org.dsp.config.Constants
import org.dsp.filters.FilterUtils
import org.dsp.modulation.WaveformEffect.Companion.normalize
import org.dsp.signals.SignalGenerator.Companion.singleSine
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
    val targetHarmonic        = 1
    val targetFrequency       = targetHarmonic * baseFrequency
    val targetWaveLength      = baseWavelength / targetHarmonic
    val length                = Constants.SAMPLE_RATE//waveData.rawFileData.size//Constants.SAMPLE_RATE + 1
    val impulseResponseLength = length//* 2
    val bandwidth             = baseFrequency//44.0//baseFrequency//44.0//baseFrequency/4 //(baseFrequency/4) May be all we need or perhaps even (baseFrequency/k) or some function of this.
    val frequencyDistance     = 1.0

    println("Organic wavelength histogram:")
    val realWaveLengths = waveData.fullWaves.map { it.size }.filter { it in 107..109 }
    getHistogramStats(start = 0, length = realWaveLengths.size, signal = realWaveLengths.map { it.toDouble()})
    println()


    /*val signal = expF(
        amplitude = 100.0,
        coefficient = -1/4.0,
        length = 54
    )
    listData(data = signal)
    return*/

    /*val fout1 = FilterUtils.bandpassFIR(
        input                 = waveData.rawFileData,
        impulseResponseLength = impulseResponseLength + 1,
        minFreqHz             = (2*baseFrequency) - (bandwidth/2),
        maxFreqHz             = (3*baseFrequency) + (bandwidth/2)
    )

    WaveUtils.writeSamplesToFile(
        fileName = "$resourceDirectory/output.wav",
        input    =  normalize(scale = 5000.0, input = fout1)
    )
    return */



    /*val testFrequecy = 0.5*baseFrequency
    val s = sineByFrequency(amplitude = 1.0, frequency = testFrequecy, size = (Constants.SAMPLE_RATE/testFrequecy).toInt())
    val cout = FilterUtils.convolution(input = waveData.rawFileData, impulseResponse = s) */

    val ocarina: List<Double> = synthesizeOcarina(
        frequency             = baseFrequency,
        bandwidth             = bandwidth.toInt(),
        frequencyDistance     = frequencyDistance,
        length                = length
    )
    val hf = getWaves(data = ocarina).map { it.size }
    val fullWaves = mutableListOf<Double>()
    for(i in 1 until hf.size -1 step 2) {
        fullWaves.add((hf[i] + hf[i - 1]).toDouble())
    }
    getHistogramStats(start = 0, length = fullWaves.size, signal = fullWaves)

    /*val fout2 = FilterUtils.bandpassFIR(
        input                 = waveData.rawFileData,
        impulseResponseLength = impulseResponseLength + 1,
        minFreqHz             = (2*baseFrequency) - (bandwidth/2),
        maxFreqHz             = (3*baseFrequency) + (bandwidth/2)
    )*/

    val dfto = DiscreteFourierTransform.partialDft(input = normalize(scale = 5000.0, input = waveData.rawFileData), fundamentalFrequency = baseFrequency)
    val dfts = DiscreteFourierTransform.partialDft(input = normalize(scale = 5000.0, input = ocarina             ), fundamentalFrequency = baseFrequency)

    val organic   = normalize(scale = 5000.0, input = dfto.magnitude)
    val synthetic = normalize(scale = 5000.0, input = dfts.magnitude)

    for(i in organic.indices) {
        println("O: ${organic[i]}, S: ${synthetic[i]}")
    }

    WaveUtils.writeSamplesToFile(
        fileName = "$resourceDirectory/output.wav",
        input    =  normalize(scale = 5000.0, input = ocarina)
    )
}

fun trueFmSynthesis(fc: Double, modAmp: Double, fm: Double, length: Int) : List<Double> {
    val o: MutableList<Double> = mutableListOf()
    for (i in 0 until length) {
        val modulator            = modAmp * sin((2*Math.PI*i)/ (Constants.SAMPLE_RATE/fm))
        val functionalWavelength = (Constants.SAMPLE_RATE/fc) + modulator
        o.add(sin((2*Math.PI*i)/functionalWavelength))
    }
    return o
}

fun synthesizeOcarina(
    frequency: Double,
    bandwidth: Int,
    frequencyDistance: Double,
    length: Int
): List<Double> {
    val ocarinaSpectrum             = ocarinaDft(baseFrequency = frequency)
    val lengthBuffer                = length * 2
    var buffer: MutableList<Double> = MutableList(lengthBuffer) { 0.0 }
    var maxLengthFormant            = 0

    /** Create fhe fundamental formant with its jitter signal and added to the buffer**/
    val fundamentalWavelength:Int = Constants.SAMPLE_RATE / frequency.toInt()
    //val rrange = ((fundamentalWavelength * 0.03)/2).toInt()
    //val fundamentalVibratoSignal = (0 until(length / fundamentalWavelength)).map { fundamentalWavelength + (-rrange..rrange).random()}
    val fundamentalWavelengthSignal: List<Int> = generateFundamentalWavelengthSignal(
        baseWavelength = fundamentalWavelength.toDouble(), jitterRange = 1.0, timeInSamples = length
    ).map { it.toInt()}
    val fundamentalFormant: List<Double> = normalize(
        scale = ocarinaSpectrum[0], input = fundamentalWavelengthSignal.map { singleSine(amplitude = 1.0, size = it)}.flatten()
    )
    for(f in fundamentalFormant.indices) {
        buffer[f] += fundamentalFormant[f]
    }

    for(spectrumIndex in 1 until ocarinaSpectrum.size) {
        val harmonicNumber = spectrumIndex + 1
        println("Calculating formant ${harmonicNumber}....")
        val amplitude = ocarinaSpectrum[spectrumIndex]
        val formant   = normalize(
            scale = amplitude,
            input = synthesizeFormant(
                fundamentalWavelengthSignal = fundamentalWavelengthSignal,
                harmonicNumber           = harmonicNumber,
                centerFrequency          = frequency*harmonicNumber,
                bandwidth                = bandwidth,
                frequencyDistance        = frequencyDistance,
                length                   = length
            )
        )
        for(f in formant.indices) {
            buffer[f] += formant[f]
        }
        maxLengthFormant = Math.max(maxLengthFormant, formant.size)
    }

    println("Buffer size: ${buffer.size}, maxLengthFormant: $maxLengthFormant")
    return buffer.subList(0, maxLengthFormant)
}

fun synthesizeFormant(
    fundamentalWavelengthSignal: List<Int>,
    harmonicNumber: Int,
    centerFrequency: Double,
    bandwidth: Int,
    frequencyDistance: Double,
    length: Int
) : List<Double> {
    val outerInnerFormantRatio        = 0.5 // Innerformant(central)/Outer(secondary)
    val envelopeChangePercent         = 0.05/2.0
    val attackAndDecayFrequencyChange = envelopeChangePercent * centerFrequency
    println("Synthesizing formant at $centerFrequency...")
    println("Attack/Decay Change: $attackAndDecayFrequencyChange")

    val waveLengthSignal = periodBoundJitter(
        fundamentalWavelengthSignal = fundamentalWavelengthSignal,
        harmonicNumber           = harmonicNumber
    )
    val innerFormant: List<Double> = waveLengthSignal.map { singleSine(amplitude = 1.0, size = it)}.flatten()

    val baseScalingAmplitude = 1.0
    val postFilterInnerFormant = normalize(
        scale = baseScalingAmplitude,
        input = FilterUtils.inharmonicDftFilter(
            input                 = innerFormant,
            centerFrequency       = centerFrequency,
            bandwidth             = bandwidth,
            frequencyDistance     = frequencyDistance,
            length                = length
        )
    )

    return postFilterInnerFormant

    val outerFormant = normalize(
        scale = outerInnerFormantRatio * baseScalingAmplitude,
        input = synthesizeOuterFormant(
            amplitude       = 1.0,
            coefficient     = -1.0/16, //-1.0/64//(16*decayFactor),
            centerFrequency = centerFrequency,
            bandWidth       = bandwidth,
            length          = length
        )
    )
    return postFilterInnerFormant.add(outerFormant)
}

class Period(var modified: Boolean, var size: Int)

fun periodBoundJitter(
    fundamentalWavelengthSignal: List<Int>,
    harmonicNumber: Int
) : List<Int> {
    return fundamentalWavelengthSignal.map {
        periodBoundJitterHelper(
            fundamentalWavelength = it,
            harmonicNumber        = harmonicNumber,
            timeInSamples         = it
        )
    }.flatten()
}


/**TODO: This should probably be added to a fundamental signal as a baseline so (fundamental + jitter). This better
 * captures whats actually going on
 */
/** TODO: This really should not need a timeInSamples variable any longer now that it's been tailored to only execute
 * for a single period driven by a wavelength signal of the fundamental harmonic
  */
fun periodBoundJitterHelper(
    fundamentalWavelength: Int,
    harmonicNumber: Int,
    timeInSamples: Int
) : List<Int> {
    val numPeriods: Int  = (timeInSamples / fundamentalWavelength)
    val periods: MutableList<Period> = mutableListOf()
    /** Organize each fundamental period into sequences of integer wavelengths periods of the harmonic in question so
     * that the total samples utiled by the sequence of the harmonic is equal to one period of the fundamental.
     *
     * Note: The rounding leftover is distributed randomly. This so far produces no noise and makes sounds that are identical to if
     * the wavelength was not rounded. But if you do NOT distribute the rounding leftover randomly then the perception of a larger
     * overall inharmonic waveform takes shape and the sound is metallic.
     */
    for(p in 0 until numPeriods) {
        val roundedWaveLength = (fundamentalWavelength.toDouble() / harmonicNumber).toInt()
        val leftOver          = fundamentalWavelength - (roundedWaveLength * harmonicNumber)

        repeat(harmonicNumber) {
            periods.add(Period(modified = false, size = roundedWaveLength))
        }

        var m = 0
        while(m != leftOver) {
            val period = periods[periods.size - 1 - (0 until harmonicNumber).random()]
            if(!period.modified) {
                period.size += 1
                period.modified = true
                m++
            }
        }

        /** Reset the modified flag (TODO: a more functional approach is needed **/
        for(j in 0 until harmonicNumber) {
            periods[periods.size - 1 - j].modified = false
        }

        //TODO: How much wavelength variation should be distributed per cycle?

        /** TODO: increase the wavelength of x random periods **/
        m = 0
        while(m != ceil(harmonicNumber/4.0).toInt()) {
            val period = periods[periods.size - 1 - (0 until harmonicNumber).random()]
            if(!period.modified) {
                period.size += 4
                period.modified = true
                m++
            }
        }

        /** TODO: Decrease the wavelength of y random periods **/
        m = 0
        while(m != ceil(harmonicNumber/4.0).toInt()) {
            val period = periods[periods.size - 1 - (0 until harmonicNumber).random()]
            if(!period.modified) {
                period.size -= 4
                period.modified = true
                m++
            }
        }

        for(j in 0 until harmonicNumber) {
            val lastSum = periods.subList(periods.size - harmonicNumber, periods.size).map{it.size}.sum()
            if(fundamentalWavelength != lastSum) {
                throw RuntimeException("Overflow!!! lastSum: $lastSum, fundamentalWavelength: $fundamentalWavelength")
            }
        }
    }
    return periods.map { it.size }
}


fun generateFundamentalWavelengthSignal(
    baseWavelength: Double,
    jitterRange: Double,
    timeInSamples: Int
) : List<Double> {

    /*val options: List<Double> =  listOf(
        13.0,
        14.0
    )*/
   /* val options: List<Double> =  listOf(
        (baseWavelength),
        (baseWavelength) + jitterRange
    )*/
    val options: MutableList<Double> = mutableListOf()
    repeat(48) {
        options.add(108.0)
    }
    repeat(41) {
        options.add(109.0)
    }
    repeat(10) {
        options.add(107.0)
    }

    val o: MutableList<Double> = mutableListOf()
    var count                  = 0.0
    while(count < timeInSamples) {
        val chosenValue = options[options.indices.random()]
        o.add(chosenValue)
        count += chosenValue
    }
    return o
/*
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

    return o */
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

//TODO: This will need an equation or interpolation for larger wavelengths
fun ocarinaDft(baseFrequency: Double) : List<Double> {
    /*return listOf(
        0.0,
        184.94136177017273,
        140.46091762627321
    )

    return listOf(10000.0) + expF(
        amplitude = 100.0,
        coefficient = -1/8.0, //we know 16 is probably too slow and 8 is too high so trying 12
        length = 53
    )*/

    return listOf(
        5000.0,
        184.94136177017273,
        140.46091762627321,
        17.283687719655802,
        19.56940324244725,
        16.471444299042613,
        7.645258488700816,
        13.155503081204106,
        3.1738433883470223,
        11.090018522589068,
        3.377728871538592,
        2.4105179938835897,
        2.6326739544744324,
        1.3535911399238596,
        1.325863224012884,
        1.9896012499274944,
        1.245440038107859,
        1.8455239538957742,
        1.5639680504707987,
        0.39465191515382453,
        0.46833289648269744,
        0.5298771810521437,
        1.9716322980630892,
        0.9816853213814685,
        2.1737145189673246,
        1.3080943163112606,
        1.0971473138370247,
        0.4371143657228572,
        0.5579105666475643,
        0.6365539641994094,
        0.7777038170934596,
        0.6624235768601714,
        0.8050031102392452,
        0.22751403114466967,
        0.25274331757475743,
        0.719038328568466,
        0.5003373504325695,
        1.1950639641280625,
        0.9188029950171505,
        0.2833869410091478,
        1.2693203907689667,
        1.0651773286036699,
        0.8537571148560699,
        0.7924304142226949,
        1.0700755357613827,
        0.3595929619239763,
        0.32098046957076676,
        0.44539751244265574,
        0.3351952066898757,
        0.3437698351189163,
        0.6563719826362677,
        0.48182549855966167,
        0.8198157090588679,
        0.8094969820504015
    )
    /*val keyHarmonics: List<FreqPair> = listOf(
        FreqPair(amp = 5000.0, freq = baseFrequency*1),
        FreqPair(amp = 168.0, freq = baseFrequency*2),
        FreqPair(amp = 136.0, freq = baseFrequency*3),
        FreqPair(amp = 40.0,  freq = baseFrequency*4),
        FreqPair(amp = 46.0,  freq = baseFrequency*5),
        FreqPair(amp = 53.0,  freq = baseFrequency*6),
        FreqPair(amp = 20.0,  freq = baseFrequency*7),
        FreqPair(amp = 40.0,  freq = baseFrequency*8),
        FreqPair(amp = 7.0,   freq = baseFrequency*9),
        FreqPair(amp = 17.0,  freq = baseFrequency*10)
    )*/
    /*val keyHarmonics: List<FreqPair> = listOf(
        FreqPair(amp = 3232.0256977414124, freq  = baseFrequency*1),
        FreqPair(amp = 71.61919166328113, freq   = baseFrequency*2),
        FreqPair(amp = 95.18529931957727,  freq  = baseFrequency*3),
        FreqPair(amp = 28.001202204999867,  freq = baseFrequency*4),
        FreqPair(amp = 26.158342224967402,  freq = baseFrequency*5),
        FreqPair(amp = 21.218082603339695,  freq = baseFrequency*6),
        FreqPair(amp = 13.677313608548639,  freq = baseFrequency*7),
        FreqPair(amp = 13.322686613835051,  freq = baseFrequency*8),
        FreqPair(amp = 5.805614515338418,   freq = baseFrequency*9),
        FreqPair(amp = 7.975981155678555,  freq  = baseFrequency*10),
        FreqPair(amp = 4.298937191983834,  freq  = baseFrequency*11),
        FreqPair(amp = 2.830145339630758,  freq  = baseFrequency*12),
        FreqPair(amp = 2.085193674698621,  freq  = baseFrequency*13),
        FreqPair(amp = 1.8777036438052082,  freq = baseFrequency*14),
        FreqPair(amp = 1.7718408647317523,  freq = baseFrequency*15),
        FreqPair(amp = 1.7164821739120368,  freq = baseFrequency*16),
        FreqPair(amp = 1.4628793034302057,  freq = baseFrequency*17),
        FreqPair(amp = 1.3708961989390205,  freq = baseFrequency*18),
        FreqPair(amp = 1.2377887839757293,  freq = baseFrequency*19)
    )*/
}

fun synthesizeOuterFormant(
    amplitude: Double,
    coefficient: Double,
    bandWidth: Int,
    centerFrequency: Double,
    length: Int,
    distanceBetweenPoints: Double = 1.0
): List<Double> {
    val spectrumAmplitudeEnvelope = createOuterFormantEnvelope(
        amplitude   = amplitude,
        coefficient = coefficient,
        length      = bandWidth
    )
    return synthesizeInharmonicDft(
        magnitude         = spectrumAmplitudeEnvelope,
        baseFrequency     = centerFrequency,
        bandwidth         = bandWidth,
        frequencyDistance = distanceBetweenPoints,
        length            = length
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
