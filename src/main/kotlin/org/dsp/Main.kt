package org.dsp

import org.dsp.file.FileUtils
import org.dsp.wave.WaveUtils
import java.io.*
import kotlin.math.abs
import kotlin.math.sin

fun main() {
    val resourceDirectory     = "/Users/umutawakil/Documents/Git/dsp/src/main/resources"
    val readFile              = File("$resourceDirectory/test-sample.wav")
    val inputStream           = FileInputStream(readFile)
    val fileBuffer: ByteArray = inputStream.readBytes()
    val dataOffset            = 44
    val dataSizeLocation      = 40
    val dataSize              = fileBuffer.size - dataOffset
    println("FileSize: $dataSize")
    val dataSizeFromFile      = fourBytesToInt(bytes = fileBuffer, pos = dataSizeLocation)
    println("DataSizeFromFile: $dataSizeFromFile")

    var currentWaveLength             = 0
    val waveLengths: MutableList<Int> = mutableListOf()
    val waveLength: MutableList<Int>  = mutableListOf()
    val s:MutableList<Short>          = mutableListOf()
    var sign: Short                   = bytesToSample(bytes = fileBuffer, pos = dataOffset)
    var currentSample:Short

    for (i in 0 .. dataSize step 2) {
        if (i + dataOffset >= dataSize) break

        currentSample = bytesToSample(bytes = fileBuffer, pos = i + dataOffset)
        s.add(currentSample)

        if(currentSample * sign >= 0) {
            currentWaveLength++
        } else {
            waveLengths.add(currentWaveLength)
            currentWaveLength = 1
        }
        if(currentSample != 0.toShort()) {
            sign = currentSample
        }
    }

    for(i in 0 until waveLengths.size step 2) {
        if(i + 1 >= waveLengths.size) { break }
        waveLength.add(waveLengths[i] + waveLengths[i + 1])
    }

    val waves: MutableList<MutableList<Short>> = mutableListOf()
    var l = 0
    for(i in waveLengths.indices) {
        val temp: MutableList<Short> = mutableListOf()
        for(j in 0 until waveLengths[i]) {
            temp.add(s[l])
            l++
        }
        waves.add(temp)
    }

    val maxes: MutableList<Double> = mutableListOf()
    for(i in waves.indices) {
        var max = 0.0
        for(j in 0 until waves[i].size) {
            if(max < abs(waves[i][j].toDouble())) {
                max = abs(waves[i][j].toDouble())
            }
        }
        maxes.add(max)
    }

    val even: MutableList<Double> = mutableListOf()
    val odd: MutableList<Double>  = mutableListOf()
    val full: MutableList<Double> = mutableListOf()

    for(i in waveLengths.indices) {
        if(i % 2 == 0) {
            //println("e($i) v: ${waveLengths[i]}")
            even.add((waveLengths[i] - 52.0).toDouble())

        } else {
            odd.add((waveLengths[i] - 52.0).toDouble())
        }
        full.add(waveLengths[i] - 52.0)
    }

    val atomicVibrato = calculateAtomicVibrato(scale = 1.0, windowSize = 100, base = 52.0)
    val so: MutableList<Double> = mutableListOf()

    for(i in atomicVibrato.first.indices) {
        so.add(atomicVibrato.first[i])
        so.add(atomicVibrato.second[i])
    }

    val output: MutableList<Short> = mutableListOf()
    val input: List<Double> = waveLengthToOutput(peaks = maxes, vibratoSignal = so)
    input.forEach { output.add(it.toInt().toShort())}
    writeSamplesToFile(fileName = "$resourceDirectory/output.wav", samples = output)
}

fun calculateAtomicVibrato(
    scale: Double,
    windowSize: Int,
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

    val paritySumSignal = calculateSignalFromPdf(
        initialValue = 10.0,
        scale        = scale,
        windowSize   = windowSize,
        histograms   = histogramsParitySum
    )

    val evenSignal = calculateSignalFromPdf(
        initialValue = 5.0,
        scale        = scale,
        windowSize   = windowSize,
        histograms   = histogramsEven
    ).toMutableList()
    val oddSignal    = paritySumSignal.zip(evenSignal) { a, b -> a - b}.toMutableList()

    for (i in evenSignal.indices) {
        evenSignal[i] += base
    }
    for (i in oddSignal.indices) {
        oddSignal[i] += base
    }

    return Pair(first = evenSignal, second = oddSignal)
}

fun calculateSignalFromPdf(
    initialValue: Double,
    scale: Double,
    windowSize: Int,
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

        for (i in 0 until windowSize) {
            var trials = 0
            var index = (0 until options.size).random()
            while (abs(options[index] - prev) > maxDifference) {
                index = (0 until options.size).random()
                trials++
                if(trials > options.size) {
                    println("trial#:$trials, prev: $prev option: ${options[index]}, diff: ${abs(options[index] - prev)}")
                }
            }
            prev = options[index]
            o.add(prev)
        }

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
        val value = signal[i + start]
        histogram[value] = histogram.getOrDefault(value, 0) + 1
    }
    val total:Double = histogram.values.sum().toDouble()
    for(k in histogram.keys.sorted()) {
        val value:Int = histogram[k]!!
        println("value: $k, Count: $value,  ${100.0*(value/total)}%")
    }
    //println("Total: ${total}")
    return histogram
}
fun plotSignal(signal: List<Double>) {
    for(i in 0 until signal.size) {
         plot(index = i, value = signal[i])
    }
}
fun plot(index: Int, value: Double) {
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
    //print(" $value")
    println()
}

fun plotSignals(signalA: List<Double>, signalB: List<Double>) {
    for(i in signalA.indices) {
        if(i == signalB.size) { break }
        plots(index = i, valueA = signalA[i], valueB = signalB[i])
    }
}
fun plots(index: Int, valueA: Double, valueB: Double) {
    print("$index ")
    plotInline(value = valueA)
    print("  ")
    plotInline(value = valueB)
    println()
}
fun plotInline(value: Double) {
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

fun waveLengthToOutput(peaks: List<Double>, vibratoSignal: List<Double>) : List<Double> {
    val input: MutableList<Double> = mutableListOf()
    var direction = 1
    for (k in 0 until vibratoSignal.size) {
        val amp           = 5000.0//peaks[k]//5000//maxes[k]
        var previous      = 0.0
        var j             = 1

        while(true) {
            val halfPeriod = vibratoSignal[k]
            //val currentSampleValue = direction*amp * triangle(i = j , waveLength = 2*halfPeriod, amplitude = amp)//sin((2 * Math.PI * j) / (2*(halfPeriod)))
            val currentSampleValue = direction*amp * sin((2 * Math.PI * j) / (2*(halfPeriod)))

            //println("$previous $currentSampleValue")
            if (previous * currentSampleValue <= 0 && j != 1) {
                break
            }
            input.add(currentSampleValue)
            previous = currentSampleValue
            j++
        }
        direction *= -1
    }
    return input
}

fun writeSamplesToFile(fileName: String, samples: List<Short>) {
    val file              = File(fileName)
    file.delete()
    file.createNewFile()
    val os: OutputStream = FileOutputStream(file)
    val bos              = BufferedOutputStream(os)
    val outFile          = DataOutputStream(bos)

    WaveUtils.writeWaveHeader(
        outFile         = outFile,
        numberOfSamples = samples.size
    )
    for(i in samples.indices) {
        FileUtils.writeShortLE(
            out   = outFile,
            value = samples[i]
        )
    }
    outFile.flush()
    outFile.close()
}
fun fourBytesToInt(bytes: ByteArray, pos: Int): Int {
    return (bytes[pos].toInt() and 0xFF) or
            ((bytes[pos + 1].toInt() and 0xFF) shl 8) or
            ((bytes[pos + 2].toInt() and 0xFF) shl 16) or
            ((bytes[pos + 3].toInt() and 0xFF) shl 24)
}

fun bytesToSample(bytes: ByteArray, pos: Int) : Short {
    return ((bytes[pos].toInt() and 0xFF) or ((bytes[pos + 1].toInt() and 0xFF) shl 8)).toShort()
}
