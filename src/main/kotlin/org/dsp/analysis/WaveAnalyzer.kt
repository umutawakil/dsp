package org.dsp.analysis

import com.github.psambit9791.wavfile.WavFile
import java.io.*
import java.nio.ByteBuffer
import java.nio.ByteOrder
import java.util.*
import kotlin.collections.HashSet


class WaveAnalyzer {
    companion object {

        fun dataToSignal(fileName: String) : Array<Double>{
            val waveFile:WavFile = WavFile.openWavFile(File(fileName))
            val frames = waveFile.framesRemaining.toInt()
            val output = Array<Double>(frames) {0.0}

            println("Frames: $frames")
            var sampleBuffer: IntArray = IntArray(1)
            waveFile.readFrames(sampleBuffer, 1)

            var bytesRead = -1
            var i = 0
            while (bytesRead != 0) {
                //println("samples: $numSamples")
                bytesRead = waveFile.readFrames(sampleBuffer, 1)
                if (bytesRead == 0) {
                    println("BytesRead was 0")
                    break
                }
                output[i] = sampleBuffer[0].toDouble()
                //println("Out: ${output[i]}")
                i++
            }
            println("Samples read -> output: ${output.size}, i: $i")
            return output
        }
        fun getWaveEnvelope(fileName: String, envelope: LinkedList<Double>) {
            val waveFile = WavFile.openWavFile(File(fileName))
            var sampleBuffer: IntArray = IntArray(1)
            waveFile.readFrames(sampleBuffer, 1)

            val set = HashSet<Int>()
            var numSamples = 0
            var bytesRead = 1
            var sum = 0
            var collecting = false
            var wave = 0
            var averageDiff = 0.0
            var prev = 0.0
            var peak = 0.0
            while (bytesRead != 0) {
                //println("samples: $numSamples")
                bytesRead = waveFile.readFrames(sampleBuffer, 1)
                if(bytesRead == 0) {
                    println("BytesRead was 0")
                    break
                }
                if(sampleBuffer[0]> 0) {
                    collecting = true
                    sum += sampleBuffer[0]
                    wave++
                    if(sampleBuffer[0] > peak) {
                        peak = sampleBuffer[0].toDouble()
                    }
                } else {
                    /*if (numSamples == 0) {
                        prev = sum * 2.0
                        averageDiff = prev
                    } else {
                        averageDiff = Math.abs(((sum * 2) - averageDiff)) / 2
                    }*/
                   // if(wave != 0) { println("waveLength: ${wave*2}, sum: ${sum*2}, prev: $prev, diff: ${(sum*2) - prev }")}
                   // wave = 0

                }
                if(collecting && (sampleBuffer[0] < 0)) {
                    var diff = (peak - prev)
                    var perc = (diff / peak) * 100
                    set.add(wave*2)
                    //println("waveLength: ${wave*2}, peak: ${peak}, prev: $prev, diff: $diff, perc: $perc")
                    wave = 0
                    prev = peak

                    collecting = false
                    envelope.add(((sum * 2.0))*1000)
                    sum = 0
                    peak = 0.0
                }
                numSamples++
            }
            set.forEach{println("${(44100.0/it)}")}
            println("AverageDiff: $averageDiff")
            println("Num of samples: $numSamples")
            println("Periods: " + envelope.size)
        }
    }
}