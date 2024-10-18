package org.dsp.filters

import org.dsp.config.Constants
import kotlin.math.cos
import java.util.LinkedList
import java.util.Queue
import kotlin.math.pow
import kotlin.math.sin

@Suppress("unused")
class FilterUtils {
    companion object {

        fun bandpassFIR(input: List<Double>, impulseResponseLength: Int,minFreqHz: Double, maxFreqHz: Double) : List<Double> {
            val lowerBandKernel = windowedSinc(fc = minFreqHz/ Constants.SAMPLE_RATE, size = impulseResponseLength)
            val upperBandKernel = windowedSinc(fc = maxFreqHz/ Constants.SAMPLE_RATE, size = impulseResponseLength)
            val passBandKernel  = upperBandKernel.zip(lowerBandKernel) { a, b -> a - b }

            return convolution(input = input, impulseResponse = passBandKernel)
        }

         /**
          *
          * Powerful low pass filter. YOu can uncomment a particular line to use an ideal IDFT impulse response or build one from the windowed-sinc filter because
         * there is serious ripple and transients in the "ideal" signal. Make note that you are not to re-amplify the output by some large gain because the
         * cancelled signals are there just dropped in amplitude significantly. The output signal is meant to be used as is so if you are filtering just a sine wave
         * to block and test that the filter is blocking it will be there. Don't amplify the output unless you also have WANTED components present so the relative amplitude
         * between the pass band and stop band will be audible; otherwise, it will seem like your filter is not working regardless of how long you make the kernel.
         */
        fun lowPassFIR(input: List<Double>, impulseResponseLength:  Int, fc: Double ) : List<Double> {
            val harmonics = (impulseResponseLength/2 ) + 1
            val bwPerbin = (Constants.SAMPLE_RATE / 2) / harmonics

            var cutOffHarmonic = 1
            while(Constants.SAMPLE_RATE/(impulseResponseLength/cutOffHarmonic) < fc) { //TODO: I've gone back and fourth on whether it should be < fc or <= fc. At the moment I think if theres a harmonic equal to fc you want it in the passband and nothing else
                cutOffHarmonic++
            }
            val cutOffFrequency = Constants.SAMPLE_RATE/(impulseResponseLength/cutOffHarmonic)
            val transitionBandFrequency = (4.0* Constants.SAMPLE_RATE)/impulseResponseLength
            val numHarmonicsForTransition = transitionBandFrequency / bwPerbin
            println()
            println("Fundamental Filter frequency: ${Constants.SAMPLE_RATE/impulseResponseLength}HZ")
            println("Harmonic frequency is: $cutOffFrequency HZ, harmonic number: $cutOffHarmonic")
            println("Bandwidth Per Bin: $bwPerbin HZ")

            println("Number of Harmonics/Bins: $harmonics, cutOffFrequency: $cutOffFrequency HZ")
            println("Transition band estimate: $transitionBandFrequency HZ")
            println("Harmonics for transition: $numHarmonicsForTransition")

            val dft:MutableList<Double> = MutableList(harmonics) {0.0}
            for(h in 0 until cutOffHarmonic + 1) {
                dft[h] = 1.0
            }

            val impulseResponse = windowedSinc(fc = fc/ Constants.SAMPLE_RATE, size = impulseResponseLength)

            //val impulseResponse = DiscreteFourierTransform.inverseDFTRaw(m = dft)
            return convolution(input = input, impulseResponse = impulseResponse)
        }

        //TODO: M is not the length of h but h.size - 1
        //TODO: H must be odd so m can always be even
        private fun windowedSinc(fc: Double, size: Int) : List<Double> {
            val h: MutableList<Double> = MutableList(size) {0.0}
            val m = h.size - 1
            for(i in h.indices) {
                if(i - (m/2) == 0) {
                    h[i] = 2*Math.PI*fc
                } else {
                    val window = 0.42 - (0.52* cos((2 * Math.PI * i) / m)) + (0.08*cos((4*Math.PI*i)/m))
                    h[i] = (sin(2 * Math.PI * fc * (i - (m / 2))) / (i - (m / 2))) * window
                }
            }
            var sum = 0.0
            h.forEach { sum += it }
            for(i in h.indices) {
                h[i] = h[i]/sum
            }
            return h
        }

        fun convolve(input: List<Double>, impulseResponse: List<Double>) : List<Double> {
            val output: MutableList<Double> = MutableList(impulseResponse.size + input.size - 1) {0.0}
            for(i in output.indices) {
                for(j in impulseResponse.indices) {
                    if((i - j) < 0) continue
                    if((i - j) > (input.size - 1)) continue
                    output[i] = output[i] + impulseResponse[j]*input[i - j]
                }
            }
            return output
        }

        /** Cleaner convolution algorithm I created that just uses the idea that you turn the input signal into a series of echos
         * sliding the input signal over to each point in the impulse response and scale it by the corresponding value at that point.
         * The length of input + IR - 1 comes from the fact that at the end the IR sticks out completely minus one sample where h[0]
         * lines up with input[N-1]
         */
        fun convolution(input: List<Double>, impulseResponse: List<Double>) : List<Double> {
            val output: MutableList<Double> = MutableList(impulseResponse.size + input.size - 1) { 0.0 }
            for(h in impulseResponse.indices) {
                val currentAmplitude = impulseResponse[h]
                for(i in input.indices) {
                    output[h + i] += currentAmplitude*input[i]
                }
            }
            return output
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
         *     val length   = 100
         *     val frequency = 48000.0/length
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
            val xe = Math.E.pow(-1 * Math.PI * fc)
            val a0 = 1 - xe

            repeat(iterations) {
                for (i in 1 until y.size) {
                    y[i] = a0 * x[i] + xe * y[i - 1]
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
                val queue: Queue<Double> = LinkedList()
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