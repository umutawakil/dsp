package org.dsp

import org.dsp.config.Constants
import org.dsp.file.FileUtils
import org.dsp.filters.FilterUtils
import org.dsp.sounds.SoundGenerator
import org.dsp.wave.WaveUtils
import java.io.*
import java.nio.ByteBuffer
import java.nio.ByteOrder
import kotlin.math.abs
import kotlin.math.sin

fun main() {
    //data starts at 44
    //data size starts at 40

    val resourceDirectory = "/Users/umutawakil/Documents/Git/dsp/src/main/resources"
    val readfile          = File("$resourceDirectory/test-sample.wav")
    val inputStream       = FileInputStream(readfile)

    val wavDataChunkSize = getWavDataChunkSize("$resourceDirectory/test-sample.wav")
    println("WAV data chunk size: $wavDataChunkSize bytes")

    val fileBuffer: ByteArray = inputStream.readBytes()
    val dataOffset            = 44
    val dataSizeLocation      = 40
    val dataSize              = fileBuffer.size - dataOffset
    println("FileSize: " + dataSize)

    val dataSizeFromFile      = byteArrayToInt(bytes = fileBuffer, pos = dataSizeLocation)
    println("DataSizeFromFile: $dataSizeFromFile")

    var currentWaveLength = 0
    var currentSample:Short = 0
    var max:Short = 0
    var previousSample: Short = bytesToSample(bytes = fileBuffer, pos = dataOffset)
    val waveLengths: MutableList<Int> = mutableListOf()
    val maxes: MutableList<Short> = mutableListOf()

    for(i in dataOffset .. (dataSize / 2) step 2) {
        currentSample = bytesToSample(bytes = fileBuffer, pos = i)
        if(abs(currentSample.toDouble()) > abs(max.toDouble())) {
            max = currentSample
        }
        //println(currentSample)
        if(currentSample * previousSample > 0) {
            currentWaveLength++
        } else {
            waveLengths.add(currentWaveLength)
            maxes.add(max)
            currentWaveLength = 1
            max = 0
        }
        previousSample = currentSample
    }
    println("Periods: ${waveLengths.size}")

    val outputWaveLength  = 50
    val cycles            = 397
    val outputNumSamples  = outputWaveLength * cycles
    val x: Array<Short>   = Array(outputNumSamples){0.toShort()}
    var currentPosition   = 0
    for (k in 0 until cycles) {
        for(i in 0 until outputWaveLength) {
            x[currentPosition] = (1000 * sin((2 * Math.PI * i) / outputWaveLength)).toInt().toShort()
            currentPosition++
        }
    }
    writeSamplesToFile(fileName = "$resourceDirectory/output.wav", samples = x)

   // waveLengths.forEach { println(it)}
    //maxes.forEach { println(it)}

    /*for(i in 0..(maxes.size/2) step 2) {
        val delta = (100*(abs(maxes[i + 1].toDouble()) - abs(maxes[i].toDouble()))) / abs(maxes[i].toDouble())
        println(delta)
        //println(abs(maxes[i + 1].toDouble()) - abs(maxes[i].toDouble()))
    }*/

    /*for(i in 0..(maxes.size/2) step 2) {
        val delta = (100*(abs(maxes[i + 2].toDouble()) - abs(maxes[i].toDouble()))) / abs(maxes[i].toDouble())
        println(delta)
    }*/

    /*val file              = File("$resourceDirectory/output.wav")
    file.delete()
    file.createNewFile()
    val os: OutputStream = FileOutputStream(file)
    val bos              = BufferedOutputStream(os)
    val outFile          = DataOutputStream(bos)
    val seconds          = 1 //
    val numberOfSamples  = 44100*seconds // 1436 seconds for 44100*5 (5 seconds takes 23 minutes worth of processing time)
    val pitch            = 441.0

    val y                = Array(numberOfSamples){0.0}
    val x                = Array(y.size) {0.0}

    WaveUtils.writeWaveHeader(
        outFile         = outFile,
        numberOfSamples = numberOfSamples
    )
    SoundGenerator.inharmonic(
        amplitude = 5000.0,
        fc        = pitch,
        output    = x
    )

    for(i in y.indices) {
        FileUtils.writeShortLE(
            out   = outFile,
            value = y[i].toInt().toShort()
        )
    }
    outFile.flush()
    outFile.close()*/
}

fun writeSamplesToFile(fileName: String, samples: Array<Short>) {
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
fun byteArrayToInt(bytes: ByteArray, pos: Int): Int {
    return (bytes[pos].toInt() and 0xFF) or
            ((bytes[pos + 1].toInt() and 0xFF) shl 8) or
            ((bytes[pos + 2].toInt() and 0xFF) shl 16) or
            ((bytes[pos + 3].toInt() and 0xFF) shl 24)
}

fun bytesToSample(bytes: ByteArray, pos: Int) : Short {
    return ((bytes[pos].toInt() and 0xFF) or ((bytes[pos + 1].toInt() and 0xFF) shl 8)).toShort()
}

fun getWavDataChunkSize(wavFilePath: String): Int {
    val file = File(wavFilePath)
    val bytes = file.readBytes()

    // Check if the file is a valid WAV file
    if (bytes.size < 12 || !isWavFile(bytes)) {
        throw IllegalArgumentException("$wavFilePath is not a valid WAV file")
    }

    // Find the start of the "data" chunk
    val dataChunkIndex = findDataChunk(bytes)
    if (dataChunkIndex == -1) {
        throw IllegalArgumentException("$wavFilePath does not contain a valid data chunk")
    }

    // Read the data chunk size
    val dataChunkSizeBytes = bytes.sliceArray(dataChunkIndex + 4..dataChunkIndex + 7)
    val dataChunkSize = ByteBuffer.wrap(dataChunkSizeBytes).order(ByteOrder.LITTLE_ENDIAN).int
    return dataChunkSize
}

fun isWavFile(bytes: ByteArray): Boolean {
    return bytes[0] == 'R'.code.toByte() &&
            bytes[1] == 'I'.code.toByte() &&
            bytes[2] == 'F'.code.toByte() &&
            bytes[3] == 'F'.code.toByte()
}

fun findDataChunk(bytes: ByteArray): Int {
    for (i in 0 until bytes.size - 8) {
        if (bytes[i] == 'd'.code.toByte() &&
            bytes[i + 1] == 'a'.code.toByte() &&
            bytes[i + 2] == 't'.code.toByte() &&
            bytes[i + 3] == 'a'.code.toByte()) {
            return i
        }
    }
    return -1
}