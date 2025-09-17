package org.dsp.analysis

import org.dsp.config.Constants
import org.dsp.signals.SignalGenerator
import kotlin.math.*

@Suppress("unused")
class DiscreteFourierTransform {
    companion object {
        class FrequencyResponse(
            val real: List<Double>,
            val imaginary: List<Double>,
            val magnitude: List<Double>,
            val phaseInPercent: List<Double>
        )

        fun inharmonicDft(
            x: List<Double>,
            baseFrequency: Double,
            bandwidth: Int,
            frequencyDistance: Double
        ) : FrequencyResponse {
            val r: MutableList<Double> = mutableListOf()
            val i: MutableList<Double> = mutableListOf()

            val start       = baseFrequency - (bandwidth/2.0)
            val stop        = baseFrequency + (bandwidth/2.0)
            var currentFreq = start
            while (currentFreq < stop) { //TODO: Should this be <= ?
                val sine   = SignalGenerator.sineByFrequency(amplitude = -1.0, frequency = currentFreq, size = x.size)
                val cosine = SignalGenerator.cosineByFrequency(amplitude = 1.0, frequency = currentFreq, size = x.size)
                i.add(x.zip(sine) { a, b -> a * b}.sum() * -1)
                r.add(x.zip(cosine) { a, b -> a * b}.sum())

                currentFreq += frequencyDistance
            }

            val real      = r.map { it / x.size}
            val imaginary = i.map {it / x.size}

            val magnitude = real.zip(imaginary) { a, b -> sqrt(a.pow(2.0) + b.pow(2.0))}
            val phase     = calculatePhaseInPercent(real = real, imaginary = imaginary)

            return FrequencyResponse(
                real           = real,
                imaginary      = imaginary,
                magnitude      = magnitude,
                phaseInPercent = phase
            )
        }

        //TODO: There exists duplication between this function and it's rectangular equivalent below
        fun inverseInharmonicDftPolar(
            phaseOffset: Int = 0,
            magnitude: List<Double>,
            phaseInPercent: List<Double>,
            baseFrequency: Double,
            bandwidth: Int,
            frequencyDistance:  Double,
            length:    Int
        ) : List<Double> {
            val o: MutableList<Double> = MutableList(length) { 0.0 }
            val start             = baseFrequency - (bandwidth/2.0)
            val stop              = baseFrequency + (bandwidth/2.0)
            var currentFreq       = start
            var frequencyBinIndex = 0

            println("Start frequency: $start, End frequency: $stop")

            while (currentFreq < stop) { //TODO: should this be <= ? If so it needs to change in the specialDft function too?
                val currentWaveLength:Double = Constants.SAMPLE_RATE / currentFreq
                val currentOffset = if(phaseOffset == 0) {
                    (phaseInPercent[frequencyBinIndex] / 100.0) * currentWaveLength
                } else {
                    phaseOffset.toDouble()
                }

                for(i in o.indices) {
                    o[i] += magnitude[frequencyBinIndex] * cos((2*Math.PI*(i - currentOffset)/currentWaveLength) )
                }
                currentFreq += frequencyDistance
                frequencyBinIndex++
            }
            return o
        }

        fun inverseInharmonicDftRectangular(
            imaginary: List<Double>,
            real:      List<Double>,
            baseFrequency: Double,
            bandwidth: Double,
            frequencyDistance: Double,
            length:    Int
        ) : List<Double> {
            val o: MutableList<Double> = MutableList(length) { 0.0 }
            val start          = baseFrequency - (bandwidth/2.0)
            val stop           = baseFrequency + (bandwidth/2.0)
            var currentFreq    = start
            var amplitudeIndex = 0

            println("Start frequency: $start, End frequency: $stop")

            while (currentFreq < stop) { //TODO: should this be <= ? If so it needs to change in the specialDft function too?
                val currentWaveLength:Double = Constants.SAMPLE_RATE / currentFreq

                for(i in o.indices) {
                    val currentReal      = real[amplitudeIndex]*cos((2*Math.PI*i)/currentWaveLength)
                    val currentImaginary = imaginary[amplitudeIndex]*sin((2*Math.PI*i)/currentWaveLength)
                    o[i]                += currentReal + currentImaginary

                }
                currentFreq += frequencyDistance
                amplitudeIndex++
            }
            return o
        }

        /**
         * TODO: Move to a synthesis class of some sort?
         *
         * This is more for audio synthesis purposes and not so much for analysis like the other inverse functions.
         * It will be tailored to suit musical purposes where the other IDFT functions are for verifying work based on
         * standard calculations which is why this version defaults to randomized phase which in any other
         * situation would be a riduclous baseline.**/
        fun synthesizeInharmonicDft(
            magnitude: List<Double>,
            phaseOffset: List<Int>,
            baseFrequency: Double,
            bandwidth: Double,
            frequencyDistance: Double,
            length:    Int
        ) : List<Double> {
            val o: MutableList<Double> = MutableList(length) { 0.0 }
            val start          = baseFrequency - (bandwidth/2.0)
            val stop           = baseFrequency + (bandwidth/2.0)
            var currentFreq    = start
            var amplitudeIndex = 0

            //println("Start frequency: $start, End frequency: $stop")
            var frequencyCount = 0

            while (currentFreq <= stop) { //TODO: should this be <= ? If so it needs to change in the specialDft function too?
                //println("current: $currentFreq, index: $amplitudeIndex / ${magnitude.size}")
                val currentWaveLength:Double = Constants.SAMPLE_RATE / currentFreq
                /*var phaseOffset: Int         = if(currentFreq != baseFrequency) {
                    ((-currentWaveLength/2).toInt()..(currentWaveLength/2).toInt()).random()
                } else {
                    0
                }*/
               // val phaseOffset: Int = (0..length).random()

                for(i in o.indices) {
                    o[i] += magnitude[amplitudeIndex]*sin((2*Math.PI*(i + phaseOffset[amplitudeIndex]))/currentWaveLength)
                }
                currentFreq = "%.4f".format(currentFreq + frequencyDistance).toDouble() //frequencyDistance
                amplitudeIndex++
                frequencyCount++
            }
            println("FrequencyCount: $frequencyCount, currentFreq: $currentFreq")
            return o
        }

        fun synthesizeFrequencyRange(
            magnitude: List<Double>,
            startFrequency: Double,
            frequencyRange: Double,
            spacing: Double,
            length:    Int
        ) : List<Double> {
            val o: MutableList<Double> = MutableList(length) { 0.0 }
            val stop           = startFrequency + frequencyRange
            var currentFreq    = startFrequency
            var amplitudeIndex = 0

            println("Start frequency: $startFrequency, End frequency: $stop, values: ${magnitude.size}")
            var frequencyCount = 0

            while (currentFreq < stop) { //TODO: should this be <= ? If so it needs to change in the specialDft function too?
                //println("Current frequency: $currentFreq")
                val currentWaveLength:Double = Constants.SAMPLE_RATE / currentFreq
                val phaseOffset: Int         = (
                        (-currentWaveLength/2).toInt()..(currentWaveLength/2).toInt()
                        ).random()


                for(i in o.indices) {
                    if(amplitudeIndex == magnitude.size) {
                        println("breaking early. amplitude index and magnitude mismatch")
                        break
                    }
                    o[i] += magnitude[amplitudeIndex]*sin((2*Math.PI*(i + phaseOffset))/currentWaveLength)
                }
                currentFreq += spacing
                amplitudeIndex++
                frequencyCount++
            }
            println("FrequencyCount: $frequencyCount")
            return o
        }

        //TODO: What are all the commented out lines I was using below
        fun dft(x: List<Double>) : FrequencyResponse {
            val harmonics                      = (x.size /2) + 1
            val real: MutableList<Double>      = MutableList(harmonics) { 0.0 }
            val imaginary: MutableList<Double> = MutableList(harmonics) { 0.0 }

            for (j in 0 until harmonics) {
                for (i in x.indices) {
                    imaginary[j] += x[i] * sin((2 * Math.PI * j * i) / x.size) * (-1)
                    real[j]      += x[i] * cos((2 * Math.PI * j * i) / x.size)
                }
            }
            for(i in real.indices) {
                if(i == 0 || i == real.size - 1) {
                    real[i] = real[i] / harmonics
                    //R[i] / R.size

                    //R[i] / (R.size / 2.0)
                } else {
                    real[i] = real[i] / harmonics
                //R[i] / (R.size / 2.0)
                    //R[i] = R[i] / R.size

                    //R[i] / (R.size / 2.0)
                }
            }

            for (i in imaginary.indices) {
                imaginary[i] = -1.0*imaginary[i] / imaginary.size
                //I[i] = (-1*I[i]) / (I.size / 1.0)
                //I[i] = (abs(-1*I[i])) / (I.size / 2.0)

                //I[i] = (abs(-1*I[i])) / (harmonics / 2.0)//(abs(-1*I[i])) / (I.size / 2.0)
            }
            /*for(i in I.indices) {
                if(i == 1) {
                    R[i] = 0.0
                    I[i] = 0.0
                }
            }*/
            imaginary[0] = 0.0

            val magnitude = real.zip(imaginary) { a, b -> sqrt(a.pow(2.0) + b.pow(2.0))}
            val phase     = calculatePhaseInPercent(real = real, imaginary = imaginary)

            return FrequencyResponse(
                real           = real,
                imaginary      = imaginary,
                magnitude      = magnitude,
                phaseInPercent = phase
            )
        }

        /**Useful for getting a general DFT idea over a block of data if you have some certainty about the fundamental frequency
         * such as multiple seconds of audio recording. Take with a grain of salt because the original scaling was meant for a
         * single cycle of the fundamental frequency.
         * **/
        fun partialDft(input: List<Double>, fundamental: Double) : FrequencyResponse {
            val fundamentalWavelength          = Constants.SAMPLE_RATE / fundamental
            val harmonics                      = ((fundamentalWavelength /2) + 1).toInt()
            val real: MutableList<Double>      = MutableList(harmonics) { 0.0 }
            val imaginary: MutableList<Double> = MutableList(harmonics) { 0.0 }

            for (h in 0 until harmonics) {
                for (i in input.indices) {
                    imaginary[h] += input[i] * sin((2 * Math.PI * h * i) / fundamentalWavelength) * (-1)
                    real[h]      += input[i] * cos((2 * Math.PI * h * i) / fundamentalWavelength)
                }
            }

            /**TODO: Some of this scaling may need to be changed since this is a partial DFT. Perhaps the other style
             * should be called the "IhharmonicDft" to differentiate it form this one?....
             * **/
            for(i in real.indices) {
                if(i == 0 || i == real.size - 1) {
                    real[i] = real[i] / harmonics
                    //R[i] / R.size

                    //R[i] / (R.size / 2.0)
                } else {
                    real[i] = real[i] / harmonics
                    //R[i] / (R.size / 2.0)
                    //R[i] = R[i] / R.size

                    //R[i] / (R.size / 2.0)
                }
            }

            for (i in imaginary.indices) {
                imaginary[i] = -1.0*imaginary[i] / imaginary.size
                //I[i] = (-1*I[i]) / (I.size / 1.0)
                //I[i] = (abs(-1*I[i])) / (I.size / 2.0)

                //I[i] = (abs(-1*I[i])) / (harmonics / 2.0)//(abs(-1*I[i])) / (I.size / 2.0)
            }
            /*for(i in I.indices) {
                if(i == 1) {
                    R[i] = 0.0
                    I[i] = 0.0
                }
            }*/
            imaginary[0] = 0.0

            val magnitude = real.zip(imaginary) { a, b -> sqrt(a.pow(2.0) + b.pow(2.0))}
            val phase     = calculatePhaseInPercent(real = real, imaginary = imaginary)

            return FrequencyResponse(
                real           = real,
                imaginary      = imaginary,
                magnitude      = magnitude,
                phaseInPercent = phase
            )
        }

        private fun calculatePhaseInPercent(real: List<Double>, imaginary: List<Double>): List<Double> {
            return real.zip(imaginary) { a, b ->
                atan2(b, a) * (180.0 / Math.PI)
            }.map{ (it/360)*100.0}//TODO// should this come back?.map { unwrapPhase(degrees = it)}
        }

        /** TODO: Phase unwrapping has been paused because I really wasn't encountering instances of say 500 degrees nor would that
         * have a negative impact on my downstream mathematics in terms of multiples of the wavelength. Should it be added back
         * or scrapped?
         */
        private fun unwrapPhase(degrees: Double): Double {
            // Normalize the angle to the range [0, 360)
            val normalized = degrees % 360

            // Adjust for negative values
            val adjusted = if (normalized < 0) normalized + 360 else normalized

            // Final adjustment to ensure it's within [-180, 180]
            return if (adjusted > 180) {
                adjusted - 360
            } else {
                adjusted
            }
        }

        /** Some minor interpolation occurs here because we force the length of the result
         * to be the same as the original and this at times will violate (N/2) + 1 when N is not even.
         * But the result is pretty damn close and inaudibly different.
         */
        //TODO: What are all the commented out lines I was using below
        fun inverseDftRectangular(real: List<Double>, imaginary: List<Double>,size: Int) : List<Double> {
            println("Size: $size")
            val actualLength = real.size * 2 - 2
            val x: MutableList<Double> = MutableList(size) { 0.0 }
            //TODO: check the actualLength usage
            for (i in x.indices) {
                for (j in real.indices) {
                    //x[i] = x[i] + R[j]* cos(((2 * Math.PI * j * i) / x.size)) + I[j]* sin(((2 * Math.PI * j * i) / x.size))
                    val realPart = real[j]* cos(((2 * Math.PI * j * i) / actualLength))
                    val imagPart = imaginary[j]* sin(((2 * Math.PI * j * i) / actualLength))
                    x[i]         = x[i] + realPart + imagPart
                }
                //x[i] = x[i] / 2
            }
            for(i in x.indices) {
                if(abs(x[i]) < 0.00000001) {
                    x[i] = 0.0
                }
            }
            return x
        }

        /** TODO: Should be noted that computing the IDFT without phase makes a world of difference. Should remove this
         * as soon as a solid algorithm for incorporating phase is determined.
         * **/
        @Deprecated(message = "Use an inverse rectangular function that will naturally incorporate phase.")
        fun inverseDftPolar(m: List<Double>) : List<Double> {
            val length  = m.size * 2 - 2
            val x: MutableList<Double> = MutableList(length) { 0.0 }
            for (i in 0 until length) {
                for (j in m.indices) {
                    x[i] = x[i] + m[j]* cos(((2 * Math.PI * j * i) / length))
                }
            }
            return x
        }

        fun synthesizeSpectrum(
            fundamentalFrequency: Double,
            magnitude           : List<Double>,
            phase               : List<Int>,
            length              : Int
        ) : List<Double> {
            val fundamentalWavelength  = Constants.SAMPLE_RATE / fundamentalFrequency
            val x: MutableList<Double> = MutableList(length) { 0.0 }
            val phaseIndex: List<Int>  = phase//if(phase != null) { phase } else { (magnitude.indices).map { (0 until length).random()}

                for(h in magnitude.indices) {
                val harmonicNumber = h + 1
                //val phaseIndex     = (0 until length).random()

                for (i in 0 until length) {
                    /*val basisFunction = if(h == 0) {
                        1.0
                    } else {
                        sin(((2 * Math.PI * h * (i + phaseIndex[h])) / fundamentalWavelength))
                    }*/
                    x[i] += magnitude[h] * sin(((2 * Math.PI * harmonicNumber * (i + phaseIndex[h])) / fundamentalWavelength))
                }
            }
            return x
        }
        fun synthesizeInverseDft(offset: Int = 0, m: List<Double>, length: Int) : List<Double> {
            val x: MutableList<Double> = MutableList(length) { 0.0 }
            for(j in m.indices) {
                if(j == 0) {
                    for (i in 0 until length) {
                        x[i] = x[i] + m[j]
                    }
                    continue
                }
                val phaseIndex = (0 until length).random() + offset
                for (i in 0 until length) {
                    x[i] = x[i] + m[j] * sin(((2 * Math.PI * j * (i + phaseIndex)) / length))
                }
            }
            return x
        }

        fun inversePolarWithRandomizedPhase(m: List<Double>, length: Int) : List<Double> {
            val x: MutableList<Double> = MutableList(length) { 0.0 }
            for(j in m.indices) {
                val phaseIndex = (0 until length).random()
                for (i in 0 until length) {
                    x[i] = x[i] + m[j] * cos(((2 * Math.PI * j * (i + phaseIndex)) / length))
                }
            }
            return x
        }
    }
}