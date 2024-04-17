package org.dsp

import org.dsp.file.FileUtils
import org.dsp.wave.WaveUtils
import java.io.*
import kotlin.math.abs
import kotlin.math.sin

fun main() {
    //data starts at 44
    //data size starts at 40

    val resourceDirectory     = "/Users/umutawakil/Documents/Git/dsp/src/main/resources"
    val readFile              = File("$resourceDirectory/test-sample.wav")
    val inputStream           = FileInputStream(readFile)
    val fileBuffer: ByteArray = inputStream.readBytes()
    val dataOffset            = 44
    val dataSizeLocation      = 40
    val dataSize              = fileBuffer.size - dataOffset
    println("FileSize: " + dataSize)
    val dataSizeFromFile      = fourBytesToInt(bytes = fileBuffer, pos = dataSizeLocation)
    println("DataSizeFromFile: $dataSizeFromFile")

    var currentWaveLength = 0
    var currentSample:Short = 0
    var max:Short = 0
    var previousSample: Short = bytesToSample(bytes = fileBuffer, pos = dataOffset)
    val waveLengths: MutableList<Int> = mutableListOf()
    val maxes: MutableList<Short> = mutableListOf()
    val s:MutableList<Short> = mutableListOf()

    for (i in 0 .. dataSize step 2) {
        if(i + dataOffset >= dataSize) break

        currentSample = bytesToSample(bytes = fileBuffer, pos = i + dataOffset)
        s.add(currentSample)

        if (abs(currentSample.toDouble()) > abs(max.toDouble())) {
            max = currentSample
        }
        if ((currentSample <= 0 && previousSample <= 0) || (currentSample >= 0 && previousSample >= 0)) {
            currentWaveLength++
        } else {
            waveLengths.add(currentWaveLength)
            maxes.add(max)
            currentWaveLength = 1
            max = 0
        }
        previousSample = currentSample
    }

     val maxPairs: MutableList<Pair<Short,Short>> = mutableListOf()
    for(i in 0 until maxes.size step 2) {
        maxPairs.add(Pair(maxes[i],maxes[i + 1]))
    }

    println("Max Pairs: ${maxPairs.size}")
    println("Wave halfs: ${waveLengths.size}")
    val averageWaveLength = waveLengths.average()
    println("Average waveHalf length: $averageWaveLength")

    val periodToPeriodEnvelops: MutableList<Double> = mutableListOf()
    for(i in 0..maxes.size step 2) {
        if((i + 2) >= maxes.size) break

        val delta = (1*(abs(maxes[i + 2].toDouble()) - abs(maxes[i].toDouble()))) / abs(maxes[i].toDouble())
        periodToPeriodEnvelops.add(1 + delta)
    }
    println("Maxes: ${maxes.size}")
    println("Period-to-period envelopes: ${periodToPeriodEnvelops.size}")

    val outputWaveLength      = 110 //Maybe 100 is better since the source sample-rate is 48k as opposed to 44.1k
    val cycles                = periodToPeriodEnvelops.size//198//200
    val x: MutableList<Short> = mutableListOf()

    for (k in 0 until cycles) {
        for (i in 0 until outputWaveLength) {
            //TODO: Add a wavelength signal to modulate by

            val currentValue = sin((2 * Math.PI * i) / outputWaveLength)
            var ampMod:Short
            if(currentValue < 0) {
                ampMod = abs(maxPairs[k].first.toDouble()).toInt().toShort()
            } else {
                ampMod = abs(maxPairs[k].second.toDouble()).toInt().toShort()
            }
            val finalValue = currentValue * ampMod
            x.add(finalValue.toInt().toShort())
        }
    }

    writeSamplesToFile(fileName = "$resourceDirectory/output.wav", samples = x)
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
