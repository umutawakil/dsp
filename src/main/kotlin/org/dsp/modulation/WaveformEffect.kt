package org.dsp.modulation

import org.dsp.analysis.WaveformAnalyzer
import kotlin.math.abs

@Suppress("unused")
class WaveformEffect {
    companion object {
        fun interpolate(waveform: List<Double>, delta: Int): List<Double> {
            if (delta == 0) return waveform

            val outputSize = waveform.size + delta
            if (outputSize <= 0) {
                println("Empty List WTF!!!: $outputSize")
                return emptyList()
            }

            val result = mutableListOf<Double>()
            val inputSize = waveform.size.toDouble()

            for (i in 0 until outputSize) {
                val inputIndex = i * (inputSize - 1) / (outputSize - 1)
                val lowerIndex = inputIndex.toInt()
                val upperIndex = minOf(lowerIndex + 1, waveform.lastIndex)

                if (lowerIndex == upperIndex) {
                    result.add(waveform[lowerIndex])
                } else {
                    val fraction = inputIndex - lowerIndex
                    val interpolatedValue = waveform[lowerIndex] * (1 - fraction) + waveform[upperIndex] * fraction
                    result.add(interpolatedValue)
                }
            }

            return result
        }

        fun normalize(scale: Double, input: List<Double>) : List<Double> {
            var max = 0.0
            for(i in input.indices) {
                if(abs(max) < abs(input[i])) {
                    max = abs(input[i])
                }
            }
            if(max == 0.0) { return  input }
            return input.map { scale * (it/max)}
        }

        fun normalizeWaves(input: List<Double>, scale: Double) : List<Double> {
            return WaveformAnalyzer.getWaves(data = input).map { normalize(scale = scale, input = it) }.flatten()
        }
    }
}