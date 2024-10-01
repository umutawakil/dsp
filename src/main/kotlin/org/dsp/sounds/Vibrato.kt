package org.dsp.sounds

import kotlin.math.abs

@Suppress("unused")
class Vibrato {
    companion object {

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

        /** This is a more generalized form of the other because it allows for arbitrary values for size while the other
         * is hog tied to a collection of histograms of fixed length
         */
        fun generateVibratoSignal(noBase: Boolean, base: Double,size: Int) : List<Double> {
            val internalBase = 52.0
            val paritySignal = generateParitySumSignal(size = size/2)
            val diffSignal   = (0 until size).map { (0..1).random()}
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

            if(noBase) {
                for (i in output.indices) {
                    output[i] =  (output[i] / internalBase)
                }
            } else {
                for (i in output.indices) {
                    output[i] = (base * (output[i] / internalBase)) + base
                }
            }

            return output
        }

        //TODO: Should be able to choose what span of the histogram is desired if the initial envelope is not desired
//TODO: As well as a starting point
        private fun generateParitySumSignal(size: Int) : List<Double> {
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
            val histogramList = histogramMap.toList()

            val stepSize = 1.0
            val options: MutableList<Double> = mutableListOf()
            for(p in histogramList) {
                repeat(p.second) {
                    options.add(p.first)
                }
            }
            //var tempOptions = options.toMutableList()
            //var tempHistogramMap = histogramMap.toMutableMap()

            val o: MutableList<Double> = mutableListOf()
            o.add(histogramList[0].first)
            /*tempOptions.removeAt(index = 0)
            tempHistogramMap[histogramList[0].first] = histogramMap[histogramList[0].first]!! - 1*/
            //var check = false

            while(o.size != size) {
                while(true) {
                    val index = (0 until options.size).random()
                    val result = options[index]

                    if (abs(result - o.last()) <= stepSize) {
                        o.add(result)
                        /*if(check) {
                            tempHistogramMap[result] = tempHistogramMap[result]!! - 1
                            if (tempHistogramMap[result] == 0) {
                                tempOptions = tempOptions.subList(index + 1, tempOptions.size)
                            } else {
                                tempOptions.removeAt(index = index)
                            }
                            if (tempOptions.size == 0) {
                                //println("No more options")
                                tempOptions = options.toMutableList()
                                tempHistogramMap = histogramMap.toMutableMap()
                                check = false
                            }
                        }*/
                        break
                    }
                }
            }
            return o
        }


        /** Code below is from the legacy vibrato generation algorithm **/

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