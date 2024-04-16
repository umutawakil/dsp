package org.dsp.sounds

import org.dsp.config.Constants
import org.dsp.file.FileUtils
import java.io.DataOutputStream
import kotlin.math.cos
import kotlin.math.pow
import kotlin.math.sin

class SoundGenerator {
    class SP(val a: Double, val f: Double)
    companion object {

        fun inharmonic(amplitude: Double, fc: Double, output: Array<Double>) {
            /*val numInharmonics = 2
            val harmonicDistance = 1
            //val temp = Array<Double>(output.size) {0.0}
            for(j in 0 until numInharmonics) {
                for(i in output.indices) {
                    val currentSinussoidValue = sin((2 * Math.PI * i * (fc + (j*harmonicDistance))) / Constants.SAMPLE_RATE)
                    output[i] = output[i] + currentSinussoidValue
                }
            }*/
            val d = Constants.SAMPLE_RATE / 3
            for(i in output.indices) {
                output[i] = output[i] - sin((2 * Math.PI * fc*(i + d)) / Constants.SAMPLE_RATE)
                output[i] = output[i] - sin((2 * Math.PI * (fc + 1)*(i + d)) / Constants.SAMPLE_RATE)
                output[i] = output[i] - sin((2 * Math.PI * (fc + 2)*(i + d)) / Constants.SAMPLE_RATE)
                //output[i] = output[i] + sin((2 * Math.PI * (fc + 2)*(i + 22050)) / Constants.SAMPLE_RATE)
            }

            for(i in output.indices) {
                output[i] *= amplitude
            }
        }

        fun dirtySignal(amplitude: Double, fc: Double, output: Array<Double>) {
            val waveLength        = (Constants.SAMPLE_RATE / fc).toInt()
            val periods           = output.size / waveLength
            val minAmpModPercent  = 0.1 // (minAmpMod/10) percent. Divide it again by 100 to use it as a coefficient or just 1000 from the start
            val maxAmpModPercent  = 1
            val ampModIncrement   = (maxAmpModPercent - minAmpModPercent) / 10000
            val minAmpModInterval = 2
            val maxAmpModInterval = 6
            var direction         = true
            var dirtyAmp          = 0.0
            var dirtyAmpInterval  = 0
            var ampModUnderway    = false

            val vibratoVariance   = 2//10
            var vibIntervalActive = false
            var vibratoCounter    = 0
            val vMin              = waveLength - (vibratoVariance/2)
            val vMax              = waveLength + (vibratoVariance/2)
            var workingWaveLength = 0
            var vibratoLength     = 0

            var i = 0
            var p = 0
            var ampIntervalCounter = 0
            while (i < output.size) {
                if(ampModUnderway && (ampIntervalCounter == dirtyAmpInterval)) {
                    direction          = !direction
                    ampModUnderway     = false
                    ampIntervalCounter = 0
                }

                //TODO: for increasing add 1 to the coefficient produced by multiplying ampModIncrement with
                //TODO: the random number output from 1 to 100 for the amp
                if(!ampModUnderway) {
                    dirtyAmpInterval = (minAmpModInterval .. maxAmpModInterval).random()
                    ampModUnderway = true
                }
                dirtyAmp = (1..100).random() * ampModIncrement
                if(direction)  {
                    dirtyAmp += 1
                } else {
                    dirtyAmp = 1 - dirtyAmp
                }

                if(vibIntervalActive && (vibratoCounter == vibratoLength)) {
                    vibIntervalActive = false
                    vibratoCounter = 0
                }
                //println("d: " + dirtyAmp
                if(!vibIntervalActive) {
                    vibIntervalActive = true
                    workingWaveLength = (vMin..vMax).random()
                    vibratoLength = 10//(2..10).random()
                }

                println("w: $workingWaveLength")
                for(j in 0 until workingWaveLength) {
                    output[i] = (dirtyAmp*amplitude)*sin((2*Math.PI*j)/workingWaveLength )
                    i++
                    if(i >= output.size) {
                        break
                    }
                }
                vibratoCounter++
                ampIntervalCounter++
                p++
            }
        }
        fun sineEnvelope(
            amplitude: Double,
            attack: Int,
            release: Int,
            fc: Double,
            output: Array<Double>
        ) {
            val waveLength     = (Constants.SAMPLE_RATE / fc).toInt()
            val attackPeriods  = (attack / waveLength)
            val releasePeriods = (release / waveLength)
            var samplePosition = 0

            samplePosition = loopByPeriod(
                increase   = true,
                periods    = attackPeriods,
                waveLength = waveLength,
                amplitude  = amplitude,
                position   = samplePosition,
                output     = output
            )

            val sustain = output.size - (attack + release)
            val sustainPeriods = sustain / waveLength
            val complexAmpEnvelope = dirtyEnvelope(
                amplitude = 25.0,
                size      = sustainPeriods
            )
            println("Sustain: $sustain, attack: $attack, release: $release")
            var p = 0
            var sustainCalculated = 0
            //TODO: calculatedEnvelope
            //TODO: length
            //TODO: length counter
            while(p < sustainPeriods) {
                var w = 0
                val envelope = complexAmpEnvelope[p]


                while(w < waveLength) {
                    val phaseData = (2 * Math.PI * w) / waveLength

                    output[samplePosition] = (amplitude + envelope) * sin(phaseData)
                    samplePosition++
                    w++
                    sustainCalculated++
                    if(sustainCalculated == sustain) {
                        break
                    }
                }
                p++
            }

            loopByPeriod(
                increase   = false,
                periods    = releasePeriods,
                waveLength = waveLength,
                amplitude  = amplitude,
                position   = samplePosition,
                output     = output
            )
        }

        private fun dirtyEnvelope(amplitude: Double,size: Int) : Array<Double>{
            val envelope = Array<Double>(size) {0.0}

            /*for(k in 1 until 22) {
                val localAmplitude = 105 - k//105 - (4.66*k)
                //println(localAmplitude)
                for(i in envelope.indices) {
                    var sinussoid = localAmplitude*sin((2*Math.PI*i*k)/envelope.size)
                    if(sinussoid < 0) {
                        //sinussoid *= -1
                    }
                    envelope[i] = envelope[i] + sinussoid
                }
            }*/
            val minLength = 20
            val maxLength = 100
            val minSlope  = 1
            val maxSlope  = 2
            var i = 0
            var direction = true
            while(i < size) {
                val length = (minLength..maxLength).random()
                var slope  = (minSlope .. maxSlope).random().toDouble()
                if(!direction) {
                    slope *= -1
                }
                direction = !direction

                i = createLine(
                    slope      = slope,
                    start      = i,
                    length     = length,
                    input      = envelope,
                    intercept  = if(i > 0) {envelope[i - 1]} else { 0.0}
                )
            }


            for(i in envelope.indices) {
                envelope[i] = amplitude * envelope[i]
            }
            return envelope
        }

        private fun createLine(
            slope: Double,
            intercept: Double,
            start: Int,
            length: Int,
            input: Array<Double>
        ): Int {
            for(i in 0 until length) {
                if(i + start >= input.size) {
                    break
                }
                //println("Y: $intercept, m: $slope, i: $i, result: ${intercept + slope*i}")
                input[i + start] = intercept + slope*i
            }
            //println("LastIndex: $lastIndex")
            return start + length
        }

        //TODO: 13 period attack
        private fun loopByPeriod(
            increase: Boolean,
            periods: Int,
            waveLength: Int,
            amplitude: Double,
            position: Int,
            output: Array<Double>
        ) : Int {
            var samplePosition = position
            val m =  amplitude / periods

            for(p in 0 until periods) {
                //y = mx = amplitude

                var w = 0
                while(w < waveLength) {
                    //println("samplePosition: $samplePosition, w: $w, p: $p")
                    val amp = if(increase) { m*p} else {amplitude - (m*p)}

                    output[samplePosition] = amp*sin((2*Math.PI*w)/waveLength)
                    w++
                    samplePosition++
                }
            }
            return samplePosition
        }

        fun spikeAmplitudes(leftCoeff: Double, outputLeft: Array<Double>, rightCoeff: Double, outputRight: Array<Double>) {
            for(i in outputLeft.indices) {
                outputLeft[i] = Math.pow(leftCoeff*2.75,i.toDouble())
            }
            for(i in outputRight.indices) {
                outputRight[i] = Math.pow(rightCoeff*2.75, i.toDouble())
            }
        }

        //TODO: The last element should be fc but seems to be off by one
        fun applySpike(
            fc: Double,
            leftBW: Double,
            leftAmp: Array<Double>,
            rightBW: Double,
            rightAmp: Array<Double>,
            output: Array<Double>
        ) {
            val leftIncrement  = leftBW / (leftAmp.size - 1)
            var leftFc = fc - leftBW
            println("Left increment amount -> $leftIncrement")
            for(j in leftAmp.indices) {
                println("Left Amp: ${leftAmp[j]}, frequency: ${leftFc}")
                for(i in output.indices) {
                    output[i] = output[i] + leftAmp[j]*sin((2*Math.PI*leftFc*i)/Constants.SAMPLE_RATE)
                    //println("Ol: ${output[i]}")
                }
                leftFc += leftIncrement
            }

            /** Now process the right side of the spike **/
            val rightIncrement = rightBW / rightAmp.size
            var rightFc = fc + rightBW

            for(j in rightAmp.indices) {
                println("Right Amp: ${rightAmp[j]}, frequency: ${rightFc}")
                for(i in output.indices) {
                    output[i] = output[i] + 100*rightAmp[j]*sin((2*Math.PI*rightFc*i)/Constants.SAMPLE_RATE)
                    //println("Or: ${output[i]}")
                }
                rightFc -= rightIncrement
            }
        }

        fun writeWhiteNoise(outFile: DataOutputStream, numSamples: Int) {
            val max = 2.0.pow(16.0).toInt() - 1

            repeat(numSamples) {
                val result = (0..max).random().toShort()
                //println("R: $result")
                FileUtils.writeShortLE(
                    out   = outFile,
                    value = result
                )
            }
        }

        fun generateWhiteNoise(input: Array<Double>) {
            val max = 2.0.pow(16.0).toInt() - 1

            for(i in input.indices) {
                input[i] = (0..max).random().toDouble()
            }
        }

        fun sineWave(amplitude: Double, k: Int, x: Array<Double>) {
            for(i in x.indices) {
                x[i] = amplitude * sin((2 * Math.PI * k*i) / x.size)
            }
        }

        fun cosineWave(amplitude: Double, k: Int, x: Array<Double>) {
            for(i in x.indices) {
                x[i] = amplitude * cos((2 * Math.PI * k*i) / x.size)
            }
        }

        fun additiveSynthesis(pairs: Array<SP>, seconds: Int, output: Array<Double>) {
            for(p in pairs.indices) {
                for(i in output.indices) {
                    output[i] = output[i] + pairs[p].a*sin(2*Math.PI*(pairs[p].f*seconds)*i / output.size)
                }
            }
        }
    }
}