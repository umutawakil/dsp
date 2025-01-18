package org.dsp.modulation

import org.dsp.analysis.WaveformAnalyzer
import org.dsp.analysis.WaveformAnalyzer.Companion.findPeak
import org.dsp.signals.SignalGenerator
import kotlin.math.abs

@Suppress("unused")
class WaveformEffect {
    companion object {

        /** Some of these compression functions may not actually be needed once a fundamental formant is created but the
         * idea is powerful and it has been proven to work well with a fake fundamental with no outerformant but driven up and
         * down by compression with an "underwaveform". Its not going to be used to synthesize an Ocarina but having two
         * different signals drive each other with one muffling the other is definitely a feature to revisit.
         */
        fun compressSignalWithWavelengths(
            input: List<Double>,
            waveLengths: List<Int>
        ) : List<Double> {
            val targetAmplitude = 5000.0
            val out             = mutableListOf<Double>()
            var direction       = 1
            var pos             = 0

            for (i in waveLengths.indices) {
                val compressedPeriod = compressHalfPeriodHelper(
                    targetAmplitude = targetAmplitude*direction,
                    input           = input,
                    start           = pos,
                    end             = pos + waveLengths[i]
                )
                out.addAll(compressedPeriod)
                direction *= -1
                pos += waveLengths[i]
            }
            return out
        }

        fun compressSignal(
            targetAmplitude: Double,
            baseWavelength: Double,
            input: List<Double>
        ) : List<Double> {
            val halfWavelength  = (baseWavelength / 2).toInt()

            val out           = mutableListOf<Double>()
            val numPeriods    = input.size / halfWavelength
            var direction     = 1

            for(p in 0 until numPeriods) {
                val compressedPeriod = compressHalfPeriodHelper(
                    targetAmplitude = direction*targetAmplitude,
                    input           = input,
                    start           = halfWavelength*p,
                    end             = (halfWavelength*p) + halfWavelength
                )
                out.addAll(compressedPeriod)
                direction *= -1
            }
            return out
        }

        private fun compressHalfPeriodHelper(
            input: List<Double>,
            targetAmplitude: Double,
            start: Int,
            end: Int
        ) : List<Double> {
            val currentUnderWaveform = input.subList(start, end)
            val iteration1           = compressHalfPeriod(amplitude = targetAmplitude, stepSize = 100.0, input = currentUnderWaveform)
            val iteration2           = compressHalfPeriod(amplitude = targetAmplitude, stepSize = 10.0, input  = iteration1)

            return compressHalfPeriod(amplitude = targetAmplitude, stepSize = 1.0, input   = iteration2)
        }
        private fun compressHalfPeriod(amplitude: Double, stepSize: Double, input: List<Double>) : List<Double> {
            var temp: List<Double>
            var currentAmplitude = 0.0//TODO: Start at zero for those moments its already in range  //direction * range
            var iteration        = 0
            var previousDistance = 0.0

            while(true) {
                //val inputMax = input.max()
                temp = input.zip(SignalGenerator.halfSine(amplitude = currentAmplitude, size = input.size)) { a, b -> a + b }
                val peak = findPeak(input = temp)
                val currentDistance = abs(peak - amplitude)

                //TODO: Need to check if the currentDistance is larger than the previous meaning an overshoot occurred

                if (currentDistance <= stepSize) {
                    break
                }
                if((currentDistance > previousDistance) && (previousDistance > 0.0)) {
                    println("Overshoot detected -> previous: $previousDistance, currentDistance: $currentDistance, target: $amplitude, current: $currentAmplitude, step: $stepSize")
                    break
                }
                previousDistance = currentDistance

                var dir = 1
                if(amplitude > peak) {
                    currentAmplitude += stepSize
                } else if(amplitude < peak ) {
                    currentAmplitude -= stepSize
                    dir = -1
                } else {
                    throw RuntimeException("This should never execute given the check above")
                }

                //if(iteration > 50) {
                //println("i $iteration, t: $amplitude, iMax: $inputMax, d: $currentDistance, cA: $currentAmplitude, p: $peak, range: $stepSize, d: $dir")
                //}
                iteration++

                /** So far this doesn't seem to get hit anymore **/
                if(iteration > 1000) {
                    println("iteration: $iteration, target: $amplitude, distance: $currentDistance, currentAmplitude: $currentAmplitude, peak: $peak, range: $stepSize, dir: $dir")
                    //listData(data = input)
                    //println()
                    //listData(data = temp)
                    throw RuntimeException("Computation limit reached")
                }
            }
            return temp
        }

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