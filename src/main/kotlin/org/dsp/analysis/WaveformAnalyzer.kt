package org.dsp.analysis

import org.dsp.modulation.WaveformEffect.Companion.normalize
import kotlin.math.abs
import kotlin.math.sign

@Suppress("unused")
class WaveformAnalyzer {
    companion object {

        fun standardDeviation(data: List<Double>): Double {
            // Ensure the list has at least one element to calculate standard deviation
            if (data.isEmpty()) return 0.0

            // Calculate the mean
            val mean = data.sum() / data.size

            // Calculate the sum of squared differences from the mean
            val sumOfSquaredDifferences = data.sumOf { (it - mean).let { diff -> diff * diff } }

            // Calculate and return the standard deviation
            return Math.sqrt(sumOfSquaredDifferences / data.size) // For population standard deviation
            // return Math.sqrt(sumOfSquaredDifferences / (numbers.size - 1)) // For sample standard deviation
        }
        fun sumOfSquaredDifferences(data: List<Double>): Double {
            // Ensure the list has at least two elements to calculate differences
            if (data.size < 2) return 0.0

            // Compute the sum of squared differences
            return data.zipWithNext { a, b -> (b - a).let { it * it } }.sum()
        }


        /** TODO: Move this and other related to WaveAnalyzer  and create a WaveGenerator class as well as WaveformEffects or Effects **/
//TODO: The organic waveLength extraction algorithm may be leaving off the last half period

        fun getWavesP(data: List<Double>): List<List<Double>> {
            val waves = mutableListOf<List<Double>>()
            if (data.isEmpty()) return waves

            var currentWave = mutableListOf<Double>()
            var lastSign = if (data[0] >= 0.0) 1 else -1

            for (value in data) {
                val currentSign = when {
                    value > 0.0 -> 1
                    value < 0.0 || value.toBits() == (-0.0).toBits() -> -1
                    else -> 0
                }

                if (currentSign != 0 && currentSign != lastSign) {
                    if (currentWave.isNotEmpty()) {
                        waves.add(currentWave)
                        currentWave = mutableListOf()
                    }
                    lastSign = currentSign
                }
                currentWave.add(value)
            }

            if (currentWave.isNotEmpty()) {
                waves.add(currentWave)
            }

            return waves
        }

        fun getFullWavelengths(data: List<Double>) : List<Double> {
            val halfWaves = getWaves(data = data)
            val o: MutableList<Double> = mutableListOf()
            for(i in 1 until halfWaves.size step 2) {
                o.add((halfWaves[i].size + halfWaves[i - 1].size).toDouble())
            }
            return o
        }
        fun getWaves(data: List<Double>) : List<List<Double>> {
            val o: MutableList<List<Double>> = mutableListOf()
            var direction = if (data[0] >= 0.0) { 1.0 } else { -1.0 }
            var temp: MutableList<Double> = mutableListOf()

            for(i in data.indices) {
                if(hasOppositeSigns(a = data[i], b = direction)) {//|| (data[i] == -0.0 && direction == 1)) {
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
        }

        /** -0.0 is a pain. You must not use == on it and was said somewhere to use equals on it's sign value **/
        fun hasOppositeSigns(a: Double, b: Double) : Boolean {
            val aVal = if((a < 0.0) || (a.sign.equals(-0.0))) { -1 } else { 1 }
            val bVal = if((b < 0.0) || (b.sign.equals(-0.0))) { -1 } else { 1 }
            return aVal * bVal < 0
        }

        fun findPeakPosition(input: List<Double>) : Int {
            var peakIndex = 0
            for(i in input.indices) {
                if(abs(input[peakIndex]) < abs(input[i])) {
                    peakIndex = i
                }
            }
            return peakIndex
        }
        fun findPeak(input: List<Double>) : Double {
            return input[findPeakPosition(input = input)]
        }
        fun findPeaks(input: List<Double>) : List<Double> {
            return getWaves(input).map {findPeak(it)}
        }
        fun listData(data: List<*>) {
            for(i in data.indices) {
                println("$i ${data[i]}")
            }
            println()
        }

        fun derivative(input: List<Double>) : List<Double> {
            val o: MutableList<Double> = mutableListOf()
            for(i in 1 until input.size) {
                o.add(input[i] - input[i - 1])
            }
            return o
        }

        fun integral(initialValue: Double, input: List<Double>) : List<Double> {
            val o: MutableList<Double> = mutableListOf()
            var buffer = initialValue
            o.add(buffer)
            for(i in 1 until input.size) {
                buffer += input[i]
                o.add(buffer)
            }
            return o
        }

        fun evenSignal(input: List<Double>) : List<Double> {
            val o: MutableList<Double> = mutableListOf()
            for(i in input.indices) {
                if(i == 0) {
                    o.add(input[i])
                    continue
                }
                o.add((input[i] + input[input.size - i])/2)
            }
            return o
        }
        fun oddSignal(input: List<Double>) : List<Double> {
            val o: MutableList<Double> = mutableListOf()
            for(i in input.indices) {
                if(i == 0) {
                    o.add(input[i])
                    continue
                }
                o.add((input[i] - input[input.size - i])/2)
            }
            return o
        }

        fun getHistogramStats(start: Int, length: Int, signal: List<Double>) : Map<Double, Int> {
            val histogram: MutableMap<Double, Int> = mutableMapOf()
            for (i in 0 until length) {
                if(i + start >= signal.size) {
                    println("Exiting early")
                    break
                }
                val currentSignalValue = signal[i + start]
                histogram[currentSignalValue] = histogram.getOrDefault(currentSignalValue, 0) + 1
            }
           // println("Histogram Entries: ${histogram.entries.size}")
           // println("Histogram Keys: ${histogram.keys.size}")
           // println("Histogram Values: ${histogram.values.size}")
            val total:Double = histogram.values.sum().toDouble()

            val sortedEntries = histogram.entries.sortedByDescending { it.value }
            for( entry in sortedEntries) {
                println("Count: ${entry.value}, Value: ${entry.key},  ${100.0*(entry.value/total)}%")
            }

            /** TODO: Not sure why I did it this way but this treats counts as unique **/
            /*val countMap: MutableMap<Int, Double> = mutableMapOf()
            for(histogramEntry in histogram.entries) {
                countMap[histogramEntry.value] = histogramEntry.key
            }

            for(c in countMap.keys.sortedDescending()) {
                println("Count: $c, Value: ${countMap[c]},  ${100.0*(c/total)}%")
            }*/

            return histogram
        }
        fun plotData(normalize: Boolean = true, scale: Double = 100.0, data: List<Double>) {
            val x = if(normalize) { normalize(scale = scale, input = data) } else { data }
            for(i in x.indices) {
                plot(index = i, value = x[i])
            }
        }
        private fun plot(index: Int, value: Double) {
            print("$index ")
            for(i in 0 until abs(value.toInt())) {
                if(value.toInt() < 0) {
                    print("-")
                }
                if(value.toInt() > 0) {
                    print("+")
                }
                if(value == 0.0) {
                    print("0")
                }
            }
            print(" $value")
            println()
        }

        fun plotSignals(signalA: List<Double>, signalB: List<Double>) {
            for(i in signalA.indices) {
                if(i == signalB.size) { break }
                plots(index = i, valueA = signalA[i], valueB = signalB[i])
            }
        }
        private fun plots(index: Int, valueA: Double, valueB: Double) {
            print("$index ")
            plotInline(value = valueA)
            print("  ")
            plotInline(value = valueB)
            println()
        }
        private fun plotInline(value: Double) {
            for(i in 0 until abs(value.toInt())) {
                if(value.toInt() < 0) {
                    print("-")
                }
                if(value.toInt() > 0) {
                    print("+")
                }
                if(value == 0.0) {
                    print(" ")
                }
            }
        }
    }
}