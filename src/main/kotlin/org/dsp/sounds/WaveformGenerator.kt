package org.dsp.sounds

import org.dsp.modulation.WaveformEffect.Companion.normalize
import kotlin.math.abs
import kotlin.math.cos
import kotlin.math.sin
import kotlin.math.sqrt
import kotlin.random.Random

@Suppress("unused")
class WaveformGenerator {

    companion object {

        fun sineByWaveLength(amplitude: Double, waveLength: Double, size: Int) : List<Double> {
            return (0 until size).map { amplitude * sin((2*Math.PI*it)/waveLength)}
        }

        fun sineByCycles(amplitude: Double, cycles: Double, size: Int) : List<Double> {
            return (0 until size).map { amplitude * sin((2*Math.PI*cycles*it)/size)}
        }

        fun sineWave(amplitude: Double, k: Int, x: Array<Double>) {
            for(i in x.indices) {
                x[i] = amplitude * sin((2 * Math.PI * k*i) / x.size)
            }
        }

        fun cosineWave(amplitude: Double, k: Int, x: Array<Double>) {
            for(i in x.indices) {
                x[i] = amplitude * cos((2 * Math.PI * k*i) / x.size)
            }
        }

        fun halfSine(size: Int) : List<Double> {
            return (0 until size).map {
                sin(
                    (2.0*Math.PI*it)/(2.0*size)
                )
            }
        }

        fun pitchedNoiseWaves(amplitude: Double, waveLengths: List<Double>) : List<List<Double>> {
            if(waveLengths.size % 2 != 0) { throw RuntimeException("Number of wavelengths supplied must be even. Num waves: ${waveLengths.size}") }

            val o: MutableList<List<Double>> = mutableListOf()
            for(i in waveLengths.indices step 2) {
                if(i + 1 == waveLengths.size) { break }
                val wA = waveLengths[i]
                val wB = waveLengths[i + 1]
                val noiseA = whiteNoise(range = amplitude, size = wA.toInt())//.map { amplitude*(it + 1.0) }
                val noiseB = whiteNoise(range = amplitude, size = wB.toInt()).map { amplitude*(it + 1.0) }
                o.add(noiseA)
                o.add(noiseB.map {-1.0*it})
            }
            return o
        }
        fun pitchedNoiseWaves(waveLength: Int, periods: Int) : List<List<Double>> {
            val halfPeriod = waveLength/2
            val o: MutableList<List<Double>> = mutableListOf()
            while(o.size < periods) {
                val noise = whiteNoise(range = 1.0, size = halfPeriod)//.map { it + 1.0 }
                o.add(noise)
                o.add(noise.map {-1.0*it})
            }
            return o
        }
        fun pitchedNoise(waveLength: Int, size: Int) : List<Double> {
            val halfPeriod = waveLength/2
            val o: MutableList<Double> = mutableListOf()
            while(o.size < size) {
                o.addAll(whiteNoise(size = halfPeriod).map { it + 1.0 } )
                o.addAll(whiteNoise(size = halfPeriod).map {-1.0*(it + 1.0)})
            }
            return o
        }

        fun mirrorPhaser(waveLengths: List<Double>): List<Double> {
            /** Pretty cool phaser I must say **/
            val o: MutableList<Double> = mutableListOf()
            val sizer = 100*800
            var oi = 0
            var inc = 0
            while(oi < sizer) {
                if(inc + 1 >= waveLengths.size) { break  }
                val lenA = waveLengths[inc].toInt()
                val lenB = waveLengths[inc + 1].toInt()

                //TODO: Add the half sine multiplication after confirming all is well

                val pos:List<Double>  = normalize(scale = 5000.0, input = brownNoise(size = lenA).map{ abs(it) })
                val neg: List<Double> = normalize(scale = 5000.0, input = brownNoise(size = lenB).map { abs(it) *-1})
                o.addAll(pos)
                o.addAll(neg)
                oi += pos.size + neg.size
                inc += 2
            }
            return o
        }

        /**** Conventional waveforms ****/

        /** TODO: What is the real best way to represent and control white noise as well as brown? These various types I have below don't all sound the same **/
        @Suppress("MemberVisibilityCanBePrivate")
        fun whiteNoise(size: Int, mean: Double = 0.0, stdDev: Double = 1.0): List<Double> {
            val random = java.util.Random()
            return List(size) {
                random.nextGaussian() * stdDev + mean
            }
        }

        @Suppress("MemberVisibilityCanBePrivate")
        fun whiteNoise(range: Double, size: Int) : List<Double> {
            val o: MutableList<Double> = MutableList(size) {0.0}
            for(i in o.indices) {
                o[i] = (0..range.toInt()).random().toDouble()
            }
            val avg = o.average()

            return o.map { (it - avg) + range}
        }

        fun brownianNoise(size: Int): List<Int> {
            val noise = mutableListOf<Int>()
            var currentValue = 5 // Start in the middle of the range

            noise.add(currentValue)

            for (i in 1 until size) {
                // Generate a step of -1, 0, or 1
                val step = Random.nextInt(-1, 2)

                // Update the current value and constrain it to the range [0, 10]
                currentValue = (currentValue + step).coerceIn(0, 10)

                noise.add(currentValue)
            }

            return noise
        }

        fun brownianNoise(size: Int, stepSize: Double, minValue: Double = 0.0, maxValue: Double = 10.0): List<Double> {
            val noise = mutableListOf<Double>()
            var currentValue = (minValue + maxValue) / 2 // Start in the middle of the range

            noise.add(currentValue)

            for (i in 1 until size) {
                // Generate a step that's exactly -stepSize, 0, or +stepSize
                val step = when (Random.nextInt(3)) {
                    0 -> -stepSize
                    1 -> 0.0
                    else -> stepSize
                }

                // Calculate the new value
                var newValue = currentValue + step

                // If the new value is outside the range, clip it to the nearest boundary
                newValue = newValue.coerceIn(minValue, maxValue)

                currentValue = newValue
                noise.add(currentValue)
            }

            return noise
        }

        @Suppress("MemberVisibilityCanBePrivate")
        fun brownNoise(size: Int): List<Double> {
            val noise = mutableListOf<Double>()
            var lastValue = 0.0
            val scale = sqrt(1.0 / size)

            for (i in 0 until size) {
                // Generate a random value between -1 and 1
                val whiteNoise = Random.nextDouble(-1.0, 1.0)

                // Integrate the white noise
                lastValue += whiteNoise * scale

                // Ensure the values stay within a reasonable range
                lastValue = lastValue.coerceIn(-1.0, 1.0)

                noise.add(lastValue)
            }

            val avg = noise.average()

            return noise.map{it - avg}
        }


    }
}