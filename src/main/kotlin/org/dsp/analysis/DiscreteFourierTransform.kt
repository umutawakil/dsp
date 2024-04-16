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

            val N = 39250
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

        /*1) How many of the harmonics outside of the center frequency can be killed to still maintain a decent sound?
        2) Whats the effect on sound if the higher harmonics are cut? They seem to only contribute noise.
        3) The lower harmonics might need to stay. Even if a fundemental is 440 hz it will still coincide with 220 hz and 110 hz ever so often. Not guranteed its even a factor but could be something
        interesting to add or experiment with if its not present. To my knowledge, subharmonics are more of a percussion sound thing?
        4) Consider the key harmonic values may be the problem. Perhaps you need to spread out the amplitude so that no one harmonic has such high values but they converge to a peak
        5) Consider what you take off needs to be redistibuted otherwise you make harnmonics too loud. I think there is a requirement to conserve total energy and re-distibute it.
        This can be done in a loop. Pick the harmonics and the quantity to chop off and a loop can go back and forth distributing the ampliitude across all harmonics and/or up to a cutoff.

        6) This wont yield results immediatley but i can establish rules about the effect of redistributing amplitude in total and over cuttoff frequence as well as can
        you effectively overcome the tremelo created by fractional inharmonics.

        *7)I suspect the secret is in distributing the amplitude in such a way that no one tone is audible
        by itself, that in order to be heard they must all be summed so nothing stands alone only the mass. ANd as such you can not hear any individual tone, only the collective
        union. SOme proof to this is how hundreds of harmonics between the fundementals multiples will just have .5 for an amplitude but over hundreds this is more magnitude than
        the few key harmonics that will have 80 or 90. You can not hear them individual but when summed their mass exceeds the handful of high points. The possiblity exists
        that removing the video game feel to these noise is simply a question of bringing more equilibrium to the Frequency domain. If true than what we perceive as a tone
        has a lot to do with how the frequency domain converges.*/

        fun inverseDFT(m: Array<Double>,phase: Array<Double>, x: Array<Double>) {
            for (i in x.indices) {
                for (j in m.indices) {
                    x[i] = x[i] + m[j]* cos(((2 * Math.PI * j * i) / x.size) + phase[j])
                }
            }
        }
    }
}