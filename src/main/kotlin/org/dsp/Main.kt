package org.dsp

import org.dsp.analysis.DiscreteFourierTransform
import org.dsp.file.FileUtils
import org.dsp.wave.WaveUtils
import java.io.*
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
    for(i in waveLength.indices) {
        val temp: MutableList<Short> = mutableListOf()
        for(j in 0 until waveLength[i]) {
            temp.add(s[l])
            l++
        }
        waves.add(temp)
    }

    println("Number of waves: ${waves.size}")

    val magnitude: MutableList<Array<Double>> = mutableListOf()
    val phase: MutableList<Array<Double>>  = mutableListOf()

    for(t in waves.indices) {
        val x: Array<Double>  = waves[t].map { it.toDouble()}.toTypedArray()
        val m: Array<Double> = Array((x.size/2) + 1) { 0.0}
        val p: Array<Double> = Array((x.size/2) + 1) { 0.0}
        DiscreteFourierTransform.dft(x = x, m = m, phase = p)
        magnitude.add(m)
        phase.add(p)
       //println("$t ${m[5]}")
    }

    val start    = 50
    val harmonic = 5
    val swap: MutableList<Double> = mutableListOf()
    for(l in start until magnitude.size) {
        swap.add(magnitude[l][harmonic])
    }
    //swap.forEach { println(it)}
    val m: Array<Double> = Array((swap.size/2) + 1) { 0.0}
    val p: Array<Double> = Array((swap.size/2) + 1) { 0.0}
    DiscreteFourierTransform.dft(x = swap.map { it.toDouble()}.toTypedArray(), m = m, phase = p)
    m.forEach { println(it) }


    val output: MutableList<Short> = mutableListOf()
    for (d in 0 until magnitude.size) {
        val periodLength = waveLength[d] //104
        val harmonics    = 11//(periodLength / 2) + 1

        for (t in 0 until periodLength) {
            var sampleValue = 0.0
            for(k in 0 until harmonics) { //52
                sampleValue += magnitude[d][k] * sin((2 * Math.PI * k * t) / periodLength)
            }
            output.add(sampleValue.toInt().toShort())
        }
    }
    writeSamplesToFile(fileName = "$resourceDirectory/output.wav", samples = output)
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
