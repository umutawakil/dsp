package org.dsp.analysis

import java.io.File
import kotlin.math.*

class DiscreteFourierTransform {
    companion object {
        fun synthesizeFromDftFile(fileName: String, start: Int, length: Int) : Array<Double> {
            val file = File(fileName)
            val lines = file.readLines().subList(start, start + length)
            val amplitude = Array<Double>(lines.size){0.0}
            val phase     = Array<Double>(lines.size){0.0}

            for(i in lines.indices) {
                amplitude[i] = (lines[i].split(":")[0]).toDouble()
                phase[i]     = (lines[i].split(":")[0]).toDouble()
            }

            val N = 39250 //What on earth was I doing here?
            val output = Array<Double>(N) {0.0}
            //results.forEach{println("$it")}
            for(k in amplitude.indices) {
                for(i in 0 until N) {
                    output[i] = output[i] + amplitude[k]*sin(((2*Math.PI*(start + k)*i)/N) + phase[k])
                }
            }
            return output
        }
        fun dft(x: Array<Double>, m: Array<Double>, phase: Array<Double>) {
            if(m.size != (x.size/2) + 1) {
                throw RuntimeException("Frequency domain input m must be ((N/2) + 1) in length. A rule of the DFT")
            }
            for (j in m.indices) {
                m[j]          = 0.0
                var imaginary = 0.0
                var real      = 0.0
                for (i in x.indices) {
                    imaginary += x[i] * sin((2 * Math.PI * j * i) / x.size) * (-1)
                    real      += x[i] * cos((2 * Math.PI * j * i) / x.size)
                }
                m[j] = sqrt(real*real + imaginary*imaginary)
                if(real == 0.0) {
                    real = 0.00000000001
                }
                phase[j] = atan((imaginary / real))

                if (j != 0 && (j != m.size - 1)) {
                    m[j] = m[j] / (x.size / 2)
                } else {
                    m[j] = m[j] / x.size
                }
                //println("k: $j, M[k]: ${m[j]}")
            }
        }

        fun dftStandard(x: Array<Double>, real: Array<Double>, imaginary: Array<Double>) {
            for (j in real.indices) {
                for (i in x.indices) {
                    imaginary[j] += x[i] * sin((2 * Math.PI * j * i) / x.size) * (-1)
                    real [j]     += x[i] * cos((2 * Math.PI * j * i) / x.size)
                }
                if(real[j] == 0.0) {
                    real[j] = 0.00000000001
                }

                if (j != 0 && (j != real.size - 1)) {
                    real[j] = real[j] / (x.size / 2)
                    imaginary[j] = imaginary[j] / (x.size / 2)
                } else {
                    real[j] = real[j] / x.size
                    imaginary[j] = imaginary[j] / x.size
                }
            }
        }

        fun inverseDFT(m: Array<Double>,phase: Array<Double>, x: Array<Double>) {
            for (i in x.indices) {
                for (j in m.indices) {
                    x[i] = x[i] + m[j]* cos(((2 * Math.PI * j * i) / x.size) + phase[j])
                }
            }
        }
        fun inverseDftStandard(real: Array<Double>,imaginary: Array<Double>, x: Array<Double>) {
            if((real.size != (x.size/2) + 1) || (imaginary.size != (x.size/2) + 1)) {
                throw RuntimeException("Frequency domain input m must be ((N/2) + 1) in length. A rule of the DFT")
            }
            for (i in x.indices) {
                for (j in real.indices) {
                    x[i] = x[i] + real[j]* cos((2 * Math.PI * j * i) / x.size) + imaginary[j]* cos((2 * Math.PI * j * i) / x.size)
                }
            }
        }
    }
}