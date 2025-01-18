package org.dsp.signals

import org.dsp.analysis.WaveformAnalyzer
import org.dsp.config.Constants
import org.dsp.modulation.WaveformEffect.Companion.normalize
import kotlin.math.abs
import kotlin.math.cos
import kotlin.math.sin
import kotlin.math.sqrt
import kotlin.random.Random

@Suppress("unused")
class SignalGenerator {

    companion object {

        /** This accounts for the fact that the actual wavelength and period are not the same thing when we talk of sampling
         * Many frequencies exist between 13 and 14 samples. This function allows for fractional periods which it will treat
         * the zero crossings in an integer fashion because we have no choice but will still compute them as close as possible.
         * We can then use the output to modulate a carrier frequency to get even better results perhaps
         *
         * **/
        fun synthesizeSineByVibratoHalfWavelengths(halfWavelengths: List<Double>) : List<Double> {
            var direction = 1.0
            val o: MutableList<Double> = mutableListOf()
            for(w in halfWavelengths) {
                var t = 0
                while(true) {
                    val data = direction * sin((Math.PI*(t + 0.5))/w)
                    if(!WaveformAnalyzer.hasOppositeSigns(a = data, b = direction)) {
                        o.add(data)
                        t++
                    } else {
                        break
                    }
                }
                direction *= -1.0
            }

            return o
        }

        fun synthesizeCosineByVibratoFullWavelengths(fullWavelengths: List<Double>) : List<Double> {
            val o: MutableList<Double> = mutableListOf()
            for (w in fullWavelengths) {
                var t    = 0
                var prev = 0.0

                while(true) {
                    val data = cos((2*Math.PI*(t + 0.5))/w)
                    if ((data > 0.0) && (prev > data) && (t > w/4)) {
                        break
                    }
                    o.add(data)
                    t++
                    prev = data
                }
            }

            return o
        }
        fun synthesizeSineByVibratoFullWavelengths(fullWavelengths: List<Double>) : List<Double> {
            val o: MutableList<Double> = mutableListOf()
            for (w in fullWavelengths) {
                var t               = 0
                var directionChange = 0
                var direction       = 1.0

                while(true) {
                    val data = sin((2*Math.PI*(t + 0.5))/w)
                    if (WaveformAnalyzer.hasOppositeSigns(a = data, b = direction)) {
                        direction *= -1.0
                        directionChange++
                        if(directionChange == 2) {
                            break
                        }
                    }
                    o.add(data)
                    t++
                }
            }
            return o
        }



        /*fun sineByVibrato(amplitude: Double, baseWavelength: Double, waveLengths: List<Double>) : List<Double> {
            val o: MutableList<Double> = mutableListOf()
            var direction = 1
            for(w in waveLengths) {
                o.addAll(halfSine(amplitude = direction*amplitude, size = (baseWavelength*w).toInt()/2))
                direction *= -1
            }
            return o
        }*/

        /*fun sineByFrequency(offset: Double, amplitude: Double, frequency: Double, size: Int): List<Double> {
            return (0 until size).map { i ->
                val t = i.toDouble() / Constants.SAMPLE_RATE
                amplitude * sin(2 * Math.PI * frequency * (t + offset / Constants.SAMPLE_RATE))
            }
        }*/

        fun sineByFrequency(offset: Double = 0.0, amplitude: Double, frequency: Double, size: Int): List<Double> {
            return (0 until size).map { amplitude * sin((2*Math.PI*(it + offset))/(Constants.SAMPLE_RATE/frequency)) }
        }
        fun cosineByFrequency(offset: Double = 0.0, amplitude: Double, frequency: Double, size: Int): List<Double> {
            return (0 until size).map { amplitude * cos((2*Math.PI*(it + offset))/(Constants.SAMPLE_RATE/frequency)) }
        }

        fun sineByWaveLength(amplitude: Double, waveLength: Double, size: Int) : List<Double> {
            //return (0 until size).map { amplitude * sin((2*Math.PI*(it + offset))/waveLength)}
            return (0 until size).map { amplitude * sin((2*Math.PI*it)/waveLength) }
        }
        fun cosineByWaveLength(amplitude: Double, waveLength: Double, size: Int) : List<Double> {
            return (0 until size).map { amplitude * cos((2*Math.PI*it)/waveLength)}
        }

        fun sineByCycles(amplitude: Double, cycles: Double, size: Int) : List<Double> {
            return (0 until size).map { amplitude * sin((2*Math.PI*cycles*it)/size)}
        }
        fun cosineByCycles(amplitude: Double, cycles: Double, size: Int) : List<Double> {
            return (0 until size).map { amplitude * cos((2*Math.PI*cycles*it)/size)}
        }

        /*fun tailoredHalfSine(amplitude: Double, size: Int): List<Double> {
            val line = 1.0 / size
            return (1 until size + 1).map { it ->
                amplitude * sin(
                    (Math.PI * (it - line)) / size
                )
            }
        }*/
        /** This approach tries to use values spread out more evenly without the useless 0 value a the begining
         * Wich in some contexts creates clipping. But in actuality I'm not sure anymore if this is solving anything. The
         * original fear I may have had was when both the start and end of my sine function had a zero value so the wave periods
         * would start with two adjacent zeros creating a clip. This may be irrelavant now.
         *
         * TODO: What are the benefits of this approach?
         * **/
        fun tailoredHalfSine(amplitude: Double, size: Int): List<Double> {
            return (0 until size).map { it ->
                val x = (it + 0.5) / size  // Distribute points evenly between 0 and 1
                amplitude * sin(x * Math.PI)
            }
        }

        /*fun sineWavesFromHalfPeriods(halfPeriods: List<Double>) : List<Double> {
            val o: MutableList<Double> = mutableListOf()
            var direction = 1
            var temp: MutableList<Double> = mutableListOf()

            for(i in halfPeriods.indices) {
                if(WaveformAnalyzer.hasOppositeSigns(a = data[i], b = direction)) {//|| (data[i] == -0.0 && direction == 1)) {
                    //println("Less than zero: ${coercedValue * direction}, data: ${data[i]}, direction: $direction")
                    //if(data[i] * direction <= 0 && (i > 0)) {
                    direction *= -1
                    o.add(temp)
                    temp = mutableListOf()
                }
                temp.add(data[i])
            }

            //TODO: This logic is suppose to pickup the last wave since theres no extra zero crossing to suggests it's there
            //not sure if this is the best way of dealing with this case. Could be.
            if(temp.size != 0) {
                //println("T: ${temp.size}, o: ${o.size}")
                o.add(temp)
            }
            return o
        }*/

        fun halfSine(amplitude: Double, size: Int) : List<Double> {
            return (0 until size).map {
                amplitude*sin(
                    (2.0*Math.PI*it)/(2.0*size)
                )
            }
        }

        fun hat(amplitude: Double, size: Int) : List<Double> {
            val amp = amplitude*0.5
            return (0 until size).map {
                -1*amp*cos(
                    (2.0*Math.PI*it)/size
                ) + amp
            }
        }

        fun halfHat(amplitude: Double, size: Int) : List<Double> {
            val amp = amplitude*0.5
            return (0 until size).map {
                -1*amp*cos(
                    (2.0*Math.PI*it)/(2*size)
                ) + amp
            }
        }
        fun halfCosine(amplitude: Double, size: Int) : List<Double> {
            return (0 until size).map {
                amplitude*cos(
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

        fun whiteNoise(length: Int): List<Double> {
            val noise = mutableListOf<Double>()
            var sum = 0.0

            // Generate random values between -1 and 1
            repeat(length) {
                noise.add(Random.nextDouble(-1.0, 1.0))
                sum += noise.last()
            }

            // Calculate the average
            val average = sum / length

            // Adjust values to ensure the average is zero
            return noise.map { it - average }
        }

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

            //return o.map { (it - avg) + range}
            val result = o.map { (it - avg)}
            println("White noise Average: ${result.average()}")
            return result
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