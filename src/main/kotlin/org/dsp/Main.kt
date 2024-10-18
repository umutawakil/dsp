package org.dsp


import org.dsp.config.Constants
import org.dsp.modulation.WaveformEffect.Companion.normalize
import org.dsp.wave.WaveUtils
import kotlin.math.cos
import kotlin.math.pow

fun main() {
    val resourceDirectory  = "/Users/umutawakil/Documents/Git/dsp/src/main/resources"
    val waveData = WaveUtils.getFileData(
        baseHalfPeriod    = 27.0,
        resourceDirectory = resourceDirectory,
        fileName          = "harmonic-2.wav"//"test-normalized.wav"//"harmonic-2.wav"//"harmonic-2.wav"//"normalized.wav"//"harmonic-2.wav"//"normalized.wav"//"filter-out-1.wav"
    )

    val fOut = formantToSpectrum(
        amplitude     = 100.0,
        coefficient   = -(1.0/16),//-(1.0/32),
        formantSize   = 444,
        baseFrequency = 888.0,
        length        = 48000,
        iterations    = 1
    )

    /** TODO: Even and odd decomposition signals **/

    WaveUtils.writeSamplesToFile(
        fileName = "$resourceDirectory/output.wav",
        input    =  normalize(scale = 5000.0, input = fOut)
    )
    return
}
fun formantToSpectrum(
    amplitude: Double,
    coefficient: Double,
    formantSize: Int,
    baseFrequency: Double,
    length: Int,
    distanceBetweenPoints: Double = 1.0,
    iterations: Int
) : List<Double> {
    var o: MutableList<Double> = MutableList(length) {0.0}

    repeat(iterations) {
        val fOut = formantToSpectrumHelper(
            amplitude     = amplitude,
            coefficient   = coefficient,
            formantSize   = formantSize,
            baseFrequency = baseFrequency,
            distanceBetweenPoints = distanceBetweenPoints,
            length                = length
        )
        o = o.zip(fOut) { a, b -> a + b}.toMutableList() //TODO: Very inefficient
    }
    return o
}
fun formantToSpectrumHelper(
    amplitude: Double,
    coefficient: Double,
    formantSize: Int,
    baseFrequency: Double,
    length: Int,
    distanceBetweenPoints: Double = 1.0
): List<Double> {
    //val distanceBetweenPoints = 5.0
    val f1 = formant1(amplitude = amplitude, coefficient = coefficient, length = formantSize)

    val spectra: MutableList<FreqPair> = mutableListOf()
    for(i in f1.indices) {
        spectra.add(
            FreqPair(
                amp = f1[i],
                freq = baseFrequency - (distanceBetweenPoints*formantSize) + (i*distanceBetweenPoints)
            )
        )
    }

    //listData(input = spectra.map {it.amp})

    return synthesizeSpectra(baseFrequency = baseFrequency, spectra = spectra, length = length)
}

fun formant1(amplitude: Double, coefficient: Double, length: Int) : List<Double> {
    val signal = expF(amplitude = amplitude, coefficient = coefficient, length = length + 1 )
    return signal.reversed().subList(0, signal.size - 1)+
            amplitude +
            signal.subList(1, signal.size)
}

/** Coefficient value of -0.25 seems **/
fun expF(amplitude: Double, coefficient: Double,  length: Int) : List<Double> {
    val o: MutableList<Double> = mutableListOf()
    for(i in 0 until length) {
        val value = amplitude* Math.E.pow(coefficient * i.toDouble())
        o.add(value)
    }
    return o
}

class FreqPair(val amp: Double, val freq: Double)
fun synthesizeSpectra(baseFrequency: Double, spectra: List<FreqPair>,length: Int) : List<Double> {
    val baseWaveLength = Constants.SAMPLE_RATE / baseFrequency

    val phaseMax = (baseWaveLength.toInt())/1
    val o: MutableList<Double> = MutableList(length) { 0.0 }

    for(k in spectra.indices) {
        val currentWaveLength = Constants.SAMPLE_RATE/ spectra[k].freq
        val amplitude         = spectra[k].amp
        val phaseShift = Math.PI
        val ps = (0..phaseMax).random()
        for(i in o.indices) {
            o[i] += amplitude*cos(((2*Math.PI*(i + ps))/currentWaveLength) + phaseShift)
        }
    }
    return o
}


