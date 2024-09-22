package org.dsp.filters

import java.lang.Math.cos
import java.util.LinkedList
import java.util.Queue
import kotlin.math.abs
import kotlin.math.sin

class FilterUtils {
    companion object {
        fun feedbackCombFilter(delayCoefficient: Double,delay: Int, x: Array<Double>, y: Array<Double>) {
            for(i in x.indices) {
                if(i - delay < 0) {
                    y[i] = x[i]
                  // continue
                } else {
                    y[i] = x[i] + delayCoefficient * y[i - delay]
                    //println("i: $i, x: ${x[i]}, xdelayed: ${delayCoefficient * y[i - delay]}, y: ${y[i]}")
                }
            }
        }
        fun compress(amplitude: Double, output: Array<Double>) {
            var hillCount  = 0
            var peak       = 0.0
            for(i in output.indices) {
                if(i == 0){
                    peak = output[i]
                    hillCount++
                    continue
                }
                var currentSample = output[i]
                if((output[i-1] * output[i]) >= 0) {
                    println("CurrentSample: $currentSample")
                    if((currentSample*currentSample) > (peak*peak)) {
                        peak = currentSample
                    }
                    hillCount++
                }
                if(((i + 1 < output.size) && ((output[i + 1] * output[i]) < 0)) ||(i + 1 == output.size)) {
                    val coefficient = amplitude / abs(peak)
                    println("amplitude: $amplitude, Coeff: $coefficient, peak: $peak, hillCount: $hillCount")
                    for(j in 0 until hillCount) {
                        val index = i + j - (hillCount - 1)
                        output[index] = coefficient*output[index]
                    }
                    if((i + 1) < output.size) {
                        hillCount = 0
                        peak      = output[i + 1]
                    }
                }
            }
        }

        //TODO: M is not the length of h but h.size - 1
        //TODO: H must be odd so m can always be even
        fun windowedSinc(fc: Double, h: Array<Double>) {
            val m = h.size - 1
            for(i in h.indices) {
                if(i - (m/2) == 0) {
                    h[i] = 2*Math.PI*fc
                } else {
                    val window = 0.42 - (0.52*cos((2*Math.PI*i)/m)) + (0.08*cos((4*Math.PI*i)/m))
                    h[i] = (sin(2 * Math.PI * fc * (i - (m / 2))) / (i - (m / 2))) * window
                }
            }
            var sum = 0.0
            h.forEach { sum += it }
            for(i in h.indices) {
                h[i] = h[i]/sum
            }
        }
        fun convolve(h: Array<Double>, input: Array<Double>, output: Array<Double>) {
            for(i in output.indices) {
                for(j in h.indices) {
                    if((i - j) < 0) continue
                    if((i - j) > (input.size - 1)) continue
                    output[i] = output[i] + h[j]*input[i - j]
                }
            }
        }

        fun combFilterFeedForward(delay: Int, x: Array<Double>, y: Array<Double>) {
            for(i in x.indices) {
                if((i - delay) < 0) {
                    y[i] = 0.0
                } else {
                    y[i] = 0.5*x[i] + (0.5*x[i - delay])
                }
                //y[i] = x[i]
                //println("Y: ${y[i]}")
            }
        }

        /*fun bandPassFilterBank(frequencies: Array<Double>, amplitudes: Array<Double>, bw: Double, x: Array<Double>, y: Array<Double>) {
            //bandPassFilter(amplitude = amplitudes[0], f = frequencies[0], bw = bw, x = x, y = y)
            for(i in frequencies.indices) {
                val temp = Array(y.size){0.0}
                bandPassFilter(f = frequencies[i], bw = bw, x = x, y = temp)
                add(amplitude = amplitudes[i], temp, y)
                //bandPassFilterAdditive(amplitude = amplitudes[i], f = frequencies[i], bw = bw, x = x, y = y)
            }
        }*/

        fun add(amplitude: Double, x: Array<Double>, y: Array<Double>) {
            for(i in x.indices) {
                y[i] = y[i] + (amplitude*x[i])
                //println("Y: ${y[i]}")
            }
        }

        fun bandPassFilterAdditive(amplitude: Double, f: Double, bw: Double, x: List<Double>) : List<Double> {
            //y[n] ' a0x[n] % a1x[n&1] % a2x[n&2] % a3x[n&3] %  ̨ % b1y[n&1] % b2y[n&2] % b3y[n&3] %

            val r  = 1 - 3 * bw
            val k  = (1 - (2*r*cos(2*Math.PI*f)) + (r*r)) / (2 - (2*cos(2*Math.PI*f)))
            val a0 = 1 - k
            val a1 = 2*(k - r)*cos(2 * Math.PI * f)
            val a2 = r * r - k
            val b1 = 2*r*cos(2*Math.PI*f)
            val b2 = -1*r*r

            val y: MutableList<Double> = MutableList(x.size) {0.0}
            for (i in 2 until x.size) {
                y[i] = y[i] + amplitude*(a0*x[i] + a1*x[i - 1] + a2*x[i - 2] + b1*y[i - 1] + b2*y[i - 2])
             }
            return y
        }

        /**
         *     val wlength   = 100
         *     val frequency = 48000.0/wlength
         *     val fc        = frequency/48000.0
         *     val bw        = 75.0/48000.0
         */
        fun bandPassFilter(f: Double, bw: Double, x: List<Double>, iterations: Int): List<Double> {
            val r  = 1 - 3 * bw
            val k  = (1 - (2*r*cos(2*Math.PI*f)) + (r*r)) / (2 - (2*cos(2*Math.PI*f)))
            val a0 = 1 - k
            val a1 = 2*(k - r)*cos(2 * Math.PI * f)
            val a2 = r * r - k
            val b1 = 2*r*cos(2*Math.PI*f)
            val b2 = -1*r*r
            //println("Filter coefficients -> r: $r, k: $k, a0: $a0, a1: $a1, a2: $a2, b1: $b1, b2: $b2")

            val y: MutableList<Double> = MutableList(x.size) {0.0}
            repeat(iterations) {
                for (i in 2 until x.size) {
                    y[i] = a0 * x[i] + a1 * x[i - 1] + a2 * x[i - 2] + b1 * y[i - 1] + b2 * y[i - 2]
                }
            }
            return y
        }

        fun lowPassFilter(iterations: Int, fc: Double, x: List<Double>) : List<Double> {
            val y: MutableList<Double> = MutableList(x.size) {0.0}
            val xe = Math.pow(Math.E, -1*Math.PI*fc)
            val a0 = 1 - xe
            val b1 = xe

            repeat(iterations) {
                for (i in 1 until y.size) {
                    y[i] = a0 * x[i] + b1 * y[i - 1]
                }
            }

            return y
        }

        /** TODO: For the moment I see better results using negative gain and half the target wavelength then when
         * using positive gain and the exact target wavelength. Maybe a DSP forum can explain this.
         *
         * One thing I can indeed explain is the negative feedback has less DC offset issues. When using the positive
         * feedback the DC rises incredibly and that includes when the gain is a fraction.
         */
        fun delayLine(gain: Double, delay: Int, x: List<Double>, iterations: Int): List<Double> {
            //TODO: if(gain > 0) { throw RuntimeException("Gain must be negative for the effect you actually want!!!") }
            val y: MutableList<Double> = MutableList(x.size) {0.0}

            repeat(iterations) {
                val queue: Queue<Double> = LinkedList<Double>()
                repeat(delay) {
                    queue.add(0.0)
                }

                for (i in y.indices) {
                    y[i] = (gain * queue.poll()) + x[i]
                    queue.add(y[i])
                }
            }
            return y
        }

        /** Incredibly powerful high pass filter great for cleaning up brown noise output in particular **/
        fun dcRemoval(x: List<Double>) : List<Double> {
            val y: MutableList<Double> = MutableList(x.size) {0.0}
            val a = 0.9999
            for(i in 1 until x.size) {
                y[i] = x[i] - x[i-1] + a * y[i-1]
            }
            return y
        }

        fun deAmplify(factor: Double, x: Array<Double>) {
            for(i in x.indices) {
                x[i] = x[i] / (factor*10.0)
            }
        }
    }
}