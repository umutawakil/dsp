package org.dsp.signals

import org.dsp.analysis.WaveformAnalyzer
import kotlin.math.abs

@Suppress("unused")
class Vibrato {
    companion object {

        private const val INTERNAL_BASE: Double = 52.0 //TODO This may need some cleanup since fundamentalFormant doesn't use it

        //444.444 = 108 wavelength or 54 half wavelength
        fun fundamentalFormant(baseWavelength: Double, durationInSamples: Int): List<Double> {
            /*val histogramMap = mapOf(
                8.0 to 1,
                6.0 to 2,
                5.0 to 5,
                4.0 to 8,
                3.0 to 8,
                2.0 to 81,
                1.0 to 216,
                0.0 to 442,
                -1.0 to 29,
                -2.0 to 4,
                -3.0 to 1
            ) */
            val histogramMap = mapOf(
                12.0 to 2,
                10.0 to 5,
                8.0 to 8,
                6.0 to 8,
                4.0 to 81,
                2.0 to 216,
                0.0 to 442,
                -2.0 to 29,
                -4.0 to 4,
                -6.0 to 1
            )
            return generateVibratoSignal(
                internalBaseToUse  = 54.0,
                histogramMap       = histogramMap,
                baseWavelength     = baseWavelength,
                durationInSamples  = durationInSamples,
                stepSize           = 2.0
            )
        }

        fun fluidVibratoSignal(baseWavelength: Double, durationInSamples: Int): List<Double> {
            val histogramMap = mapOf(
                9.0 to 1,
                8.0 to 3,
                7.0 to 5,
                6.0 to 6,
                5.0 to 6,
                4.0 to 14,
                3.0 to 39,
                2.0 to 91,
                1.0 to 298,
                0.0 to 226,
                -1.0 to 90,
                -2.0 to 22,
                -3.0 to 9,
                -4.0 to 5,
                -5.0 to 1,
                -7.0 to 3
            )
            return generateVibratoSignal(
                histogramMap       = histogramMap,
                baseWavelength     = baseWavelength,
                durationInSamples  = durationInSamples
            )
        }

        fun genericVibratoSignal(baseWavelength: Double, durationInSamples: Int): List<Double> {
            val histogramMap = mapOf(
                16.0 to 1,
                15.0 to 2,
                14.0 to 1,
                13.0 to 1,
                12.0 to 1,
                11.0 to 3,
                10.0 to 1,
                9.0 to 5,
                8.0 to 20,
                7.0 to 35,
                6.0 to 32,
                5.0 to 123,
                4.0 to 144,
                3.0 to 23,
                2.0 to 4,
                1.0 to 1,
                0.0 to 1
            )

            return generateVibratoSignal(
                histogramMap       = histogramMap,
                baseWavelength     = baseWavelength,
                durationInSamples  = durationInSamples
            )
        }

        private fun generateVibratoSignal(
            internalBaseToUse: Double = INTERNAL_BASE,
            histogramMap: Map<Double, Int>,
            baseWavelength: Double,
            durationInSamples: Int,
            stepSize: Double = 1.0,
        ) : List<Double> {
            val baseHalfWavelength = baseWavelength / 2

            val paritySignal = generateParitySumSignal(
                internalBaseToUse  = internalBaseToUse,
                histogramMap       = histogramMap,
                baseHalfWavelength = baseHalfWavelength,
                durationInSamples  = durationInSamples/2,
                stepSize           = stepSize
            )

            println("Parity Sum:")
            WaveformAnalyzer.getHistogramStats(start = 0, length = paritySignal.size, signal = paritySignal)
            println()

            val diffSignal   = (paritySignal.indices).map { (0..1).random()}
            val output: MutableList<Double> = mutableListOf()

            for(i in paritySignal.indices) {
                if(diffSignal[i] == 0) {
                    val even = (paritySignal[i] / 2).toInt()
                    val odd  = (paritySignal[i] / 2).toInt()
                    output.add(even.toDouble())
                    output.add(odd.toDouble())
                } else {
                    var diff = 0
                    while(diff == 0) {
                        diff = (-1..1).random()
                    }
                    val even = (paritySignal[i] / 2).toInt() + diff
                    val odd  = paritySignal[i] - even
                    output.add(even.toDouble())
                    output.add(odd)
                }
            }

            for (i in output.indices) {
                //output[i] =  baseHalfWavelength*((output[i] + INTERNAL_BASE) / INTERNAL_BASE)
                //output[i] =  ((output[i] + internalBaseToUse) / internalBaseToUse)
            }

            return output
        }

        private fun generateParitySumSignal(
            internalBaseToUse: Double,
            histogramMap: Map<Double, Int>,
            baseHalfWavelength: Double,
            durationInSamples: Int,
            stepSize: Double
        ) : List<Double> {
            val histogramList = histogramMap.toList()

            val options: List<Double> = histogramList.flatMap { it1 ->
                (0..it1.second).map { it1.first }
            }

            val o: MutableList<Double> = mutableListOf()
            o.add(histogramList[0].first)
            var duration = 0.0

            while(duration < durationInSamples) {
                while(true) {
                    val index = (options.indices).random()
                    val result = options[index]

                    if (abs(result - o.last()) <= stepSize) {
                        o.add(result)

                        duration += ((result + internalBaseToUse)/internalBaseToUse)*baseHalfWavelength
                        break
                    }
                }
            }
            return o
        }


        /** Code below is from the legacy vibrato generation algorithm---------
         * It was replaced because it uses a series of fixed histograms and can only function for a set length across
         * each histogram
         * **/

        fun createPitchNote(windowSizeInSamples: Int, base: Double) : List<Double> {
            val atomicVibrato = calculateAtomicVibrato(
                scale               = 1.0,
                windowSizeInSamples = windowSizeInSamples,
                base                = base//52.0
            )
            val so: MutableList<Double> = mutableListOf()
            for(i in atomicVibrato.first.indices) {
                if(i >= atomicVibrato.second.size) {break } //TODO: Can they be equalized?
                so.add(atomicVibrato.first[i])
                so.add(atomicVibrato.second[i])
            }
            return so
        }

        private fun calculateAtomicVibrato(
            @Suppress("SameParameterValue") scale: Double,
            windowSizeInSamples: Int,
            base: Double
        ) : Pair<List<Double>, List<Double>> {
            val histogramsEven: List<Map<Double, Int>> = listOf(
                mapOf(5.0 to 32, 4.0 to 27, 3.0 to 29, 2.0 to 10, 1.0 to 3),
                mapOf(4.0 to 2, 3.0 to 23, 2.0 to 42, 1.0 to 23, 0.0 to 7, -1.0 to 3),
                mapOf(4.0 to 2, 3.0 to 18, 2.0 to 43, 1.0 to 28, 0.0 to 8, -1.0 to 1),
                mapOf(4.0 to 2, 3.0 to 16, 2.0 to 38, 1.0 to 27, 0.0 to 12, -1.0 to 4, -2.0 to 1)
            )

            val histogramsParitySum: List<Map<Double, Int>> = listOf(
                mapOf(10.0 to 10, 9.0 to 6, 8.0 to 19, 7.0 to 33, 6.0 to 26, 5.0 to 7),
                mapOf(6.0 to 8, 5.0 to 40, 4.0 to 48, 3.0 to 4 ),
                mapOf(6.0 to 4, 5.0 to 36, 4.0 to 48, 3.0 to 12),
                mapOf(6.0 to 3, 5.0 to 34, 4.0 to 41, 3.0 to 12, 2.0 to 3, 1.0 to 5, 0.0 to 1)
            )

            println("Calculating even signal pdf...")
            val evenSignal = calculateSignalFromPdf(
                initialValue        = 5.0,
                scale               = scale,
                windowSizeInSamples = windowSizeInSamples,
                base                = base,
                histograms          = histogramsEven
            ).toMutableList()

            var evenWindowSum = 0.0
            var p = 0
            while(evenWindowSum < windowSizeInSamples) {
                evenWindowSum += evenSignal[p] + base
                p++
                println("evenWindowSum: $evenWindowSum, windowSizeInSamples: $windowSizeInSamples")
            }
            println("Even signal number count: $p")

            println("\r\nCalculating parity Sum pdf...")
            val paritySumSignal = calculateSignalFromPdf(
                initialValue        = 10.0,
                scale               = scale,
                windowSizeInSamples = windowSizeInSamples,//(1.5 * windowSizeInSamples).toInt(),
                base                = base,
                histograms          = histogramsParitySum
            )
            /*var sumB = 0.0
            for(t in 0 until 100) {
                sumB += paritySumSignal[t]
            }
            println("SumB: $sumB")*/
            val oddSignal = paritySumSignal.zip(evenSignal) { a, b -> a - b}.toMutableList()

            /** TODO: Might make sense to move this logic into the signal generation code **/
            for (i in evenSignal.indices) {
                evenSignal[i] += base
            }
            for (i in oddSignal.indices) {
                oddSignal[i] += base
            }

            return Pair(first = evenSignal, second = oddSignal)
        }

        private fun calculateSignalFromPdf(
            initialValue: Double,
            scale: Double,
            windowSizeInSamples: Int,
            base: Double,
            histograms: List<Map<Double, Int>>
        ) : List<Double> {
            val maxDifference          = scale* 1.0 //2.0
            var prev                   = scale*initialValue//10.0//15.0
            val o: MutableList<Double> = mutableListOf()
            for (j in histograms.indices) {
                val histogram = histograms[j]
                val options: MutableList<Double> = mutableListOf()
                for(entry in histogram) {
                    repeat(entry.value) {
                        options.add(scale*entry.key)
                    }
                }
                options.sort() //TODO: Is this really necessary?

                var windowSum = 0.0
                var p = 0
                //while(if(usePeriods){ p < periods } else { windowSum < windowSizeInSamples}) {
                //while(windowSum < windowSizeInSamples) {
                repeat(windowSizeInSamples) {
                //for (i in 0 until windowSizeInSamples) {
                    var trials = 0
                    if(windowSum > 0) {
                        //if(i != 0) {
                        var index = (0 until options.size).random()
                        while (abs(options[index] - prev) > maxDifference) {
                            index = (0 until options.size).random()
                            trials++
                            if (trials > options.size) {
                                println("trial#:$trials, prev: $prev option: ${options[index]}, diff: ${abs(options[index] - prev)}")
                            }
                        }

                        //TODO: Cleanup this redundancy so the if statement can be fed to index as a single line
                        prev = options[index]
                        o.add(prev)
                        windowSum += prev + base

                    } else  { /** [Important to prevent audible waves]A new window is starting so chose the option closest to the last point in the last window!! **/
                    var smallestDiffIndex = 0
                        for (oi in options.indices) {
                            if(abs(options[smallestDiffIndex] - prev) > abs(options[oi] - prev)) {
                                smallestDiffIndex = oi
                            }
                        }
                        prev = options[smallestDiffIndex]
                        o.add(prev)
                        windowSum += prev + base
                    }
                    p++
                }
                //println("WindowSum: $windowSum, P: $p")
            }
            return o
        }
    }
}