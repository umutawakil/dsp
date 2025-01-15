package org.dsp.modulation

import org.dsp.analysis.WaveformAnalyzer
import org.dsp.analysis.WaveformAnalyzer.Companion.findPeak
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
            /*var max = 0.0
            for(i in input.indices) {
                if(abs(max) < abs(input[i])) {
                    max = abs(input[i])
                }
            }*/
            val max = findPeak(input = input)
            val sign = if(max <0) { -1 } else { 1 }
            return input.map { sign * scale * (it/max)}
        }

        fun normalizeWaves(input: List<Double>, scale: Double) : List<Double> {
            return WaveformAnalyzer.getWaves(data = input).map { normalize(scale = scale, input = it) }.flatten()
        }

        fun toggleZero() : Int  {
            return if(((0..1000).random()) % 2 == 0) {
                1
            } else {
                0
            }
        }
        fun toggleSign() : Int {
            return if(((0..1000).random()) % 2 == 0) {
                1
            } else {
                -1
            }
        }
    }
}