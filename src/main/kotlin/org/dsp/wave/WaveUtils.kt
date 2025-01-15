package org.dsp.wave

import org.dsp.analysis.WaveformAnalyzer
import org.dsp.config.Constants
import org.dsp.file.FileUtils
import java.io.*
import kotlin.math.abs

class WaveUtils {

    companion object {
        @Suppress("unused")
        class WaveData(
            val rawFileData: List<Double>,
            val waveLengths: List<Double>,
            val waves:       List<List<Double>>,
            val peaks:       List<Double>,
            val averagePeak: Double,
            val fullWaves:   List<List<Double>>,
            val maxes:       List<Double>,
            val even:        List<Double>,
            val odd:         List<Double>,
            val full:        List<Double>
        )

        /**
         * Important to recognize this is for analyzing waveforms with one zero crossing per period. If the signal you want to analyze has multiple you
         * probably need to do the following.
         *
         * 1) Zero average it and or low pass filter the signal to remove all DC components
         * 2) Take a "signed" magnitude DFT so that the sign on the harmonics is preserved which is necessary to transmit the tremolo normalization that occurs
         * with in-harmonics.
         *
         * TODO: Sample usage. addi it below
         *
         *
         *     val resourceDirectory  = "/Users/umutawakil/Documents/Git/dsp/src/main/resources"
         *     val waveData = WaveUtils.getFileData(
         *         baseHalfPeriod    = 27.0,
         *         resourceDirectory = resourceDirectory,
         *         fileName          = "harmonic-2.wav"//"test-normalized.wav"//"harmonic-2.wav"//"harmonic-2.wav"//"normalized.wav"//"harmonic-2.wav"//"normalized.wav"//"filter-out-1.wav"
         *     )
         *
         *
         */
        fun getFileData(baseHalfPeriod: Double, resourceDirectory: String, fileName: String) : WaveData {
            val readFile              = File("$resourceDirectory/$fileName")//normalized.wav,test-sample.wav, test.wave
            val inputStream           = FileInputStream(readFile)
            val fileBuffer: ByteArray = inputStream.readBytes()
            val dataOffset            = 44
            val dataSizeLocation      = 40
            val dataSize              = fileBuffer.size - dataOffset
            println("FileSize: $dataSize")
            val dataSizeFromFile      = fourBytesToInt(bytes = fileBuffer, pos = dataSizeLocation)
            println("DataSizeFromFile: $dataSizeFromFile")

            var currentWaveLength                = 0
            val waveLengths: MutableList<Double> = mutableListOf()
            val s:MutableList<Short>             = mutableListOf()
            var sign: Short                      = bytesToSample(bytes = fileBuffer, pos = dataOffset)
            var currentSample:Short

            for (i in 0 .. dataSize step 2) {
                if (i + dataOffset >= dataSize) break

                currentSample = bytesToSample(bytes = fileBuffer, pos = i + dataOffset)
                s.add(currentSample)

                if(currentSample * sign >= 0) {
                    currentWaveLength++
                } else {
                    waveLengths.add(currentWaveLength.toDouble())
                    currentWaveLength = 1
                }
                if(currentSample != 0.toShort()) {
                    sign = currentSample
                }
            }
            val rawFileData: List<Double> = s.map{it.toDouble()}

            val waves: MutableList<MutableList<Double>> = mutableListOf()
            var l = 0
            for(i in waveLengths.indices) {
                val temp: MutableList<Double> = mutableListOf()
                for(j in 0 until waveLengths[i].toInt()) {
                    temp.add(rawFileData[l])
                    l++
                }
                waves.add(temp)
            }

            val fullWaves: MutableList<List<Double>> = mutableListOf()
            for(i in 0 until waves.size step 2) {
                if(i + 1 >= waves.size) { break }
                fullWaves.add(waves[i] + waves[i + 1])
            }

            val maxes: MutableList<Double> = mutableListOf()
            for (i in waves.indices) {
                var max = 0.0
                for(j in 0 until waves[i].size) {
                    if(max < abs(waves[i][j])) {
                        max = abs(waves[i][j])
                    }
                }
                maxes.add(max)
            }

            val even: MutableList<Double> = mutableListOf()
            val odd: MutableList<Double>  = mutableListOf()
            //val full: MutableList<Double> = mutableListOf()

            //TODO: The 52 needs to be determined automagically
            for (i in waveLengths.indices) {
                if (i % 2 == 0) {
                    even.add((waveLengths[i] - baseHalfPeriod))

                } else {
                    odd.add((waveLengths[i] - baseHalfPeriod))//52.0))
                }
            }
            val full = even.zip(odd) { a, b -> a + b}

            val peaks         = waves.map { WaveformAnalyzer.findPeak(it)}
            val averagePeak   = peaks.average()

            return WaveData(
                rawFileData = rawFileData,
                waveLengths = waveLengths,
                waves       = waves,
                peaks       = peaks,
                averagePeak = averagePeak,
                fullWaves   = fullWaves,
                maxes       = maxes,
                odd         = odd,
                even        = even,
                full        = full
            )
        }

        private fun writeWaveHeader(outFile: OutputStream, numberOfSamples: Int) {
            val channels       = 1
            val sampleRate     = Constants.SAMPLE_RATE
            val numChannels    = 1
            val bitsPerSample  = 16
            val byteRate       =  sampleRate * numChannels * (bitsPerSample/8)
            val totalAudioLen  = (numberOfSamples * channels * bitsPerSample) / 8
            val totalDataLen   = totalAudioLen + 36

            val header = ByteArray(44)

            header[0] = 'R'.code.toByte() // RIFF/WAVE header

            header[1] = 'I'.code.toByte()
            header[2] = 'F'.code.toByte()
            header[3] = 'F'.code.toByte()
            header[4] = (totalDataLen and 0xff).toByte()
            header[5] = (totalDataLen shr 8 and 0xff).toByte()
            header[6] = (totalDataLen shr 16 and 0xff).toByte()
            header[7] = (totalDataLen shr 24 and 0xff).toByte()
            header[8] = 'W'.code.toByte()
            header[9] = 'A'.code.toByte()
            header[10] = 'V'.code.toByte()
            header[11] = 'E'.code.toByte()

            header[12] = 'f'.code.toByte() // 'fmt ' chunk
            header[13] = 'm'.code.toByte()
            header[14] = 't'.code.toByte()
            header[15] = ' '.code.toByte()
            header[16] = 16 // 4 bytes: size of 'fmt ' chunk

            header[17] = 0
            header[18] = 0
            header[19] = 0
            header[20] = 1 // format = 1

            header[21] = 0
            header[22] = channels.toByte()
            header[23] = 0
            header[24] = (sampleRate and 0xff).toByte()
            header[25] = (sampleRate shr 8 and 0xff).toByte()
            header[26] = (sampleRate shr 16 and 0xff).toByte()
            header[27] = (sampleRate shr 24 and 0xff).toByte()
            header[28] = (byteRate and 0xff).toByte()
            header[29] = (byteRate shr 8 and 0xff).toByte()
            header[30] = (byteRate shr 16 and 0xff).toByte()
            header[31] = (byteRate shr 24 and 0xff).toByte()
            header[32] = (16 / 8).toByte() //mono block align(2 * 16 / 8).toByte() // block align

            header[33] = 0
            header[34] = 16 // bits per sample

            header[35] = 0
            header[36] = 'd'.code.toByte()
            header[37] = 'a'.code.toByte()
            header[38] = 't'.code.toByte()
            header[39] = 'a'.code.toByte()
            header[40] = (totalAudioLen and 0xff).toByte()
            header[41] = (totalAudioLen shr 8 and 0xff).toByte()
            header[42] = (totalAudioLen shr 16 and 0xff).toByte()
            header[43] = (totalAudioLen shr 24 and 0xff).toByte()

            outFile.write(header, 0, 44)
        }

        @Suppress("unused")
        fun writeSamplesToFile(fileName: String, input: List<Double>) {
            val samples = input.map {it.toInt().toShort()}
            val file              = File(fileName)
            file.delete()
            file.createNewFile()
            val os: OutputStream = FileOutputStream(file)
            val bos              = BufferedOutputStream(os)
            val outFile          = DataOutputStream(bos)

            writeWaveHeader(
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

        @Suppress("SameParameterValue")
        private fun fourBytesToInt(bytes: ByteArray, pos: Int): Int {
            return (bytes[pos].toInt() and 0xFF) or
                    ((bytes[pos + 1].toInt() and 0xFF) shl 8) or
                    ((bytes[pos + 2].toInt() and 0xFF) shl 16) or
                    ((bytes[pos + 3].toInt() and 0xFF) shl 24)
        }

        private fun bytesToSample(bytes: ByteArray, pos: Int) : Short {
            return ((bytes[pos].toInt() and 0xFF) or ((bytes[pos + 1].toInt() and 0xFF) shl 8)).toShort()
        }

        fun readDataFile(fileName: String) : List<Double> {
            val file = File(fileName)
            return file.readLines().map { it.toDouble()}
        }
    }
}