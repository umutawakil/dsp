package org.dsp.analysis

import org.dsp.listData
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

        fun dft(x: List<Double>) : List<Double> {
            val m: MutableList<Double> = mutableListOf()
            repeat((x.size /2) + 1) {
                m.add(0.0)
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

                if (j != 0 && (j != m.size - 1)) {
                    m[j] = m[j] / (x.size / 2)
                } else {
                    m[j] = m[j] / x.size
                }
                //println("k: $j, M[k]: ${m[j]}")
            }
            return m
        }

        fun subHarmonicDft(x: List<Double>, base: Double, numHarmonics: Int) : List<Double> {
            val m: MutableList<Double>   = MutableList(numHarmonics) { 0.0 }
            val reX: MutableList<Double> = MutableList(numHarmonics) { 0.0 }
            val imX: MutableList<Double> = MutableList(numHarmonics) { 0.0 }

            for (k in 1 until numHarmonics) {
                m[k]          = 0.0
                var imaginary = 0.0
                var real      = 0.0
                for (i in x.indices) {
                    imaginary += x[i] * sin((2 * Math.PI * i) / (base*k)) * (-1)
                    real      += x[i] * cos((2 * Math.PI * i) / (base*k))
                }
                m[k] = sqrt(real*real + imaginary*imaginary)
                reX[k] = real
                imX[k] = imaginary
                /*if (k != 0 && (k != m.size - 1)) {
                    m[k] = m[k] / (x.size / 2)
                } else {
                    m[k] = m[k] / x.size
                }*/
            }
            println("Real")
            listData(input = reX)
            println("imaginary")
            listData(input = imX)
            println("magnitude")
            listData(input = m)
            return m
        }

        fun dftRectangular(x: List<Double>) : Pair<List<Double>, List<Double>> {
            val harmonics = (x.size /2) + 1

            val R: MutableList<Double> = mutableListOf()
            repeat(harmonics) {
                R.add(0.0)
            }
            val I: MutableList<Double> = mutableListOf()
            repeat(harmonics) {
                I.add(0.0)
            }

            for (j in 0 until harmonics) {
                for (i in x.indices) {
                    I[j] += x[i] * sin((2 * Math.PI * j * i) / x.size) * (-1)
                    R[j] += x[i] * cos((2 * Math.PI * j * i) / x.size)
                }
            }
            for(i in R.indices) {
                if(i == 0 || i == R.size - 1) {
                    R[i] = R[i] / harmonics
                    //R[i] / R.size

                    //R[i] / (R.size / 2.0)
                } else {
                    R[i] = R[i] / harmonics
                //R[i] / (R.size / 2.0)
                    //R[i] = R[i] / R.size

                    //R[i] / (R.size / 2.0)
                }
            }

            for (i in I.indices) {
                I[i] = -1.0*I[i] / I.size
                //I[i] = (-1*I[i]) / (I.size / 1.0)
                //I[i] = (abs(-1*I[i])) / (I.size / 2.0)

                //I[i] = (abs(-1*I[i])) / (harmonics / 2.0)//(abs(-1*I[i])) / (I.size / 2.0)
            }
            /*for(i in I.indices) {
                if(i == 1) {
                    R[i] = 0.0
                    I[i] = 0.0
                }
            }*/
            I[0] = 0.0

            return Pair(R, I)
        }

        /** Some minor interpolation occurs here because we force the length of the result
         * to be the same as the original and this at times will violate (N/2) + 1 when N is not even.
         * But the result is pretty damn close and inaudibly different.
         */
        fun inverseDftRectangular(size: Int, R: List<Double>, I: List<Double>) : List<Double> {
            val x: MutableList<Double> = mutableListOf()
            /*repeat(2*(R.size - 1)) {
                x.add(0.0)
            }*/

            repeat(size) {
                x.add(0.0)
            }

            for (i in x.indices) {
                for (j in 0 until R.size) {
                    //x[i] = x[i] + R[j]* cos(((2 * Math.PI * j * i) / x.size)) + I[j]* sin(((2 * Math.PI * j * i) / x.size))
                    val realPart = R[j]* cos(((2 * Math.PI * j * i) / x.size))
                    val imagPart = I[j]* sin(((2 * Math.PI * j * i) / x.size))
                    x[i]         = x[i] + realPart + imagPart
                }
                //x[i] = x[i] / 2
            }
            for(i in x.indices) {
                if(abs(x[i]) < 0.00000001) {
                    x[i] = 0.0
                }
            }
            return x
        }

        fun inverseDft(m: List<Double>, harmonics: Int) : List<Double> {
            val x: MutableList<Double> = mutableListOf()
            repeat(2*(m.size - 1)) {
                x.add(0.0)
            }

            for (i in x.indices) {
                for (j in 0 until harmonics) {
                    x[i] = x[i] + m[j]* cos(((2 * Math.PI * j * i) / x.size))
                }
            }
            return x
        }

        fun inverseDftWithLength(m: List<Double>, length: Int) : List<Double> {
            val x: MutableList<Double> = mutableListOf()
            repeat(length) {
                x.add(0.0)
            }

            for (i in 0 until length) {
                for (j in 0 until m.size) {
                    x[i] = x[i] + m[j]* sin(((2 * Math.PI * j * i) / length))
                }
            }
            return x
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

        fun inverseDFTRaw(m: List<Double>): List<Double> {
            val x: MutableList<Double> = MutableList(2*(m.size - 1)) {0.0}
            for (i in x.indices) {
                for (j in m.indices) {
                    x[i] = x[i] + m[j]* sin(((2 * Math.PI * j * i) / x.size))
                }
            }
            return x
        }

        fun inverseRaisedCosine(m: List<Double>, length: Int): List<Double> {
            val x: MutableList<Double> = MutableList(length) {0.0}
            for (i in x.indices) {
                //for (j in 0 until 3) {
                for (j in m.indices) {
                    if(j == 0) { continue }
                    x[i] = x[i] + (-1*m[j]* cos(((2 * Math.PI * j * (i + x.size)) / x.size))) + m[j]
                }
            }
            return x
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