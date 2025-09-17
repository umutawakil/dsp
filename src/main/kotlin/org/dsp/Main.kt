package org.dsp

import org.dsp.analysis.DiscreteFourierTransform.Companion.dft
import org.dsp.analysis.DiscreteFourierTransform.Companion.inharmonicDft
import org.dsp.analysis.DiscreteFourierTransform.Companion.inverseDftRectangular
import org.dsp.analysis.DiscreteFourierTransform.Companion.inverseInharmonicDftPolar
import org.dsp.analysis.DiscreteFourierTransform.Companion.inversePolarWithRandomizedPhase
import org.dsp.analysis.DiscreteFourierTransform.Companion.partialDft
import org.dsp.analysis.DiscreteFourierTransform.Companion.synthesizeFrequencyRange
import org.dsp.analysis.DiscreteFourierTransform.Companion.synthesizeInharmonicDft
import org.dsp.analysis.DiscreteFourierTransform.Companion.synthesizeSpectrum
import org.dsp.analysis.WaveformAnalyzer.Companion.findPeak
import org.dsp.analysis.WaveformAnalyzer.Companion.findPeakPosition
import org.dsp.analysis.WaveformAnalyzer.Companion.findPeaks
import org.dsp.analysis.WaveformAnalyzer.Companion.getHistogramStats
import org.dsp.analysis.WaveformAnalyzer.Companion.getWaves
import org.dsp.analysis.WaveformAnalyzer.Companion.listData
import org.dsp.analysis.WaveformAnalyzer.Companion.plotData
import org.dsp.analysis.WaveformAnalyzer.Companion.standardDeviation
import org.dsp.config.Constants
import org.dsp.filters.FilterUtils
import org.dsp.filters.FilterUtils.Companion.convolution
import org.dsp.modulation.WaveformEffect.Companion.normalize
import org.dsp.modulation.WaveformEffect.Companion.normalizeWaves
import org.dsp.signals.Ocarina
import org.dsp.signals.SignalGenerator.Companion.halfPeriodFrequenciesToWaves
import org.dsp.signals.SignalGenerator.Companion.sineByCycles
import org.dsp.signals.SignalGenerator.Companion.sineByFrequency
import org.dsp.signals.SignalGenerator.Companion.singleHalfPeriod
import org.dsp.signals.SignalGenerator.Companion.whiteNoise
import org.dsp.tools.add
import org.dsp.wave.WaveUtils
import java.lang.Math.pow
import kotlin.math.*

fun main() {
    val resourceDirectory  = "/Users/umutawakil/Documents/Git/dsp/src/main/resources"
    val waveData = WaveUtils.getFileData(
        baseHalfPeriod    = 27.0, //TODO: This is kludgy as hell
        resourceDirectory = resourceDirectory,
        fileName          = "truncate.wav"
    )
    /** Phase
     *
     * TODO: Test the effect of placing a constant phase signal on all the inharmonics of a simple formant
     * and synthesize it
     * **/

    /** Confirm these two cases before running final test
     * TODO: confirm first that the sum of the envelope signals form the complement
     * TODO: confirm the complement + the original base signal envelope form a constant
     */


    /*val testSignal        = listOf(0.0, 2.0, 4.0, 10.0, 2.0)
    val magnitudePercentage = listOf(25.0, 25.0, 25.0, 25.0)
    val results = compress(
        baseSignal          = testSignal,
        magnitudePercentage = magnitudePercentage
    )
    results.forEach {
        listData(data = it)
        println(("sum: ${it.sum()}"))
        println("")
    }

    val complement: MutableList<Double> = mutableListOf()
    for(i in testSignal.indices) {
        var sum = 0.0
        for(j in magnitudePercentage.indices) {
            sum += results[j][i]
        }
        complement.add(sum)
    }
    val balancedSignal = testSignal.add(complement)
    println("Complement")
    listData(data = complement)
    println("Balanced signal")
    listData(data = balancedSignal)
    return*/

    val baseFrequency         = 445.0//444.0
    val baseWavelength        = (Constants.SAMPLE_RATE/baseFrequency)
    val targetHarmonic        = 2
    val targetFrequency       = targetHarmonic * baseFrequency
    val targetWaveLength      = Constants.SAMPLE_RATE / targetFrequency
    val length                = waveData.rawFileData.size//Constants.SAMPLE_RATE//waveData.rawFileData.size
    val impulseResponseLength = Constants.SAMPLE_RATE + 1//* 2
    val bandwidth             = baseFrequency//44.0//baseFrequency//44.0//baseFrequency/4 //(baseFrequency/4) May be all we need or perhaps even (baseFrequency/k) or some function of this.
    val frequencyDistance     = 1.0

    /*val fbw = 0.0125
    val fo = FilterUtils.bandpassFIR(
        input                 = waveData.rawFileData,
        impulseResponseLength = 2*Constants.SAMPLE_RATE + 1,
        minFreqHz             = baseFrequency - fbw ,
        maxFreqHz             = baseFrequency + fbw,
    )
    val absAvg = fo.map { abs(it)}.average()

    println("AbsAvg: $absAvg ")
    //1296.3496148319305 for 444
    // 661.2011116691183 for 445
    // 80.484 for 445 with fbw = 0.125*/
    val space = 2.0
    val s0 = sineByFrequency(amplitude = 1.0, frequency  = 300.0 + (space*1), size = length)
    val s1 = sineByFrequency(amplitude = 1.0, frequency  = 300.0 + (space*2), size = length)
    val s2 = sineByFrequency(amplitude = 1.0, frequency  = 300.0 + (space*3), size = length)
    val s3 = sineByFrequency(amplitude = 1.0, frequency  = 300.0 + (space*4), size = length)
    val sn = s0.add(s1).add((s2.add(s3)).reversed())

    //s0.add(s1) starts off then decays to zero
    //s2.add(s3) decays then increases
    //sn has its dip offset closer to the right

    WaveUtils.writeSamplesToFile(
        fileName = "$resourceDirectory/output.wav",
        input    =  withAbsAverage(
            average = 5000.0,
            input   = sn//n//filtered
        )
    )
    return

    /*val csynth = compressionSynth(
        frequency = 444.0,
        length    = length
    )*/

    /*val magnitude: List<Double> = listOf(0.0) + exponentialFunction(
        amplitude   = 5000.0,
        coefficient = -1.0/192.0,//-1.0 / 192,//512,//192.0,
        length      = (baseFrequency/2.0).toInt() + 1
    )
    val phase : List<Int> = (0 until magnitude.size).map { (0 until length).random()}

    val epiFundamentalMod = synthesizeSpectrum(
        fundamentalFrequency = 0.5,
        magnitude            = magnitude,
        phase                = phase,
        length               = length
    ).map { abs(it)}
    val epiFundamental = sineByFrequency(
        amplitude = 1.0,
        frequency = baseFrequency,
        size      = length
    ).zip(epiFundamentalMod) { a,b -> a * b}

    //TODO: generate a matrix of modulation signals using the epi-fundamental envelope
    val percentages = Ocarina.getPercentages()
    println("Percentages sum (should be 100.0): ${percentages.sum()}")
    val harmonicFormantModDfts:List<List<Double>> = compress(
        baseSignal          = epiFundamental,//magnitude,
        magnitudePercentage = percentages
    )

    validateMatrix(baseSignal = epiFundamental, results = harmonicFormantModDfts)

    return*/

    /*val phase : List<Int> = (0 until 445).map { (0 until length).random()}
    val test = simpleFormant(
        amplitude = 5000.0,
        frequency = 444.0,
        bandwidth = bandwidth,
        booster   = 0.0,
        phase     = phase,
        length    = length
    ) */

    val noise = whiteNoise(length = (length / (targetWaveLength * 0.5)).toInt())//.map{abs(it) }

    val fout = FilterUtils.lowPassFIR(
        input                 = noise,
        fc                    = 2500.0,
        impulseResponseLength = (noise.size * 2) + 1
    ).map { abs(it)}
    //val waves = getWaves(data = fout)
    val peaks = fout//whiteNoise(length = (length / (targetWaveLength * 0.5)).toInt()).map{abs(it) }//findPeaks(input = fout).map { abs(it)}
    val delta = derivative(input = peaks)
    val sg    = signGroups(input = delta)
    val rsg   = realSignGroups(input = delta)

    listData(data = sg)
    listData(data = delta)
    listData(data = rsg)

    //TODO: Need an integrator that can take the sums over a given sign

    //getHistogramStats(start = 0, length = delta.size, signal = delta)

    var direction = 1.0
    val o: MutableList<Double> = mutableListOf()
    for(i in peaks.indices) {
        o.addAll(singleHalfPeriod(amplitude = direction*peaks[i], frequency = targetFrequency))
        direction *= -1
    }

   /* val filtered = FilterUtils.bandpassFIR(
        input                 = o,
        minFreqHz             = targetFrequency - (bandwidth/2.0),
        maxFreqHz             = targetFrequency + (bandwidth/2.0),
        impulseResponseLength = Constants.SAMPLE_RATE + 1
    )*/

    WaveUtils.writeSamplesToFile(
        fileName = "$resourceDirectory/output.wav",
        input    =  withAbsAverage(
            average = 5000.0,
            input   = o//filtered
        )
    )
}

fun realSignGroups(input: List<Double>) : List<Double> {
    val o: MutableList<Double> = mutableListOf()
    var runningValue           = 1.00//(input[0]/100.0) + 1.0
    var sum                    = input[0]

    for(value in input) {
        if(sum * value < 0) {
            o.add((runningValue - 1) * 100.0)
            sum          = 0.0
            runningValue = 1.0
        }
        val quantity = (runningValue * (value / 100.0))
        runningValue += quantity
        sum          += value
        //runningValue *= ((value/100.0) + 1.0)
    }
    o.add((runningValue - 1.0) * 100.0)
    return o
}
fun signGroups(input: List<Double>) : List<Double> {
    val o: MutableList<Double> = mutableListOf()
    var sum = input[0]
    for(value in input) {
        if(sum * value < 0) {
            o.add(sum)
            sum = 0.0
        }
        sum += value
    }
    o.add(sum)
    return o
}

fun derivative(input: List<Double>) : List<Double> {
    val o: MutableList<Double> = mutableListOf()
    for(i in 1 until input.size) {
        o.add(((input[i] - input[i - 1])/input[i]) * 100.0)
    }
    return o
}
fun filterFormant(frequency: Double, peaks: List<Double>) : List<Double> {
    var direction = 1.0
    val o: MutableList<Double> = mutableListOf()
    for(i in peaks.indices) {
        o.addAll(singleHalfPeriod(amplitude = direction*peaks[i], frequency = frequency))
        direction *= -1
    }
    return o
}

fun validateMatrix(
    baseSignal: List<Double>,
    results: List<List<Double>>
) {
    val complement: MutableList<Double> = MutableList(baseSignal.size) { 0.0}
    println("Number of results: ${results.size}")
    for (i in results.indices) {
        for(j in results[i].indices) {
            complement[j] += results[i][j]
        }
        // println("s: $sum")
    }
    val balancedSignal = baseSignal.add(complement)
    //println("Complement")
    //listData(data = complement)

    println("Balanced signal")
    plotData(data = balancedSignal)
}
fun simpleSynth(
    frequency: Double,
    length: Int
) : List<Double> {
    val underwaveform: MutableList<Double> = MutableList(length) { 0.0}
    val bins                               = ((frequency/1.0) + 1).toInt()

    val phase : List<Int> = (0 until bins).map { (0 until length).random()}
    for(h in Ocarina.dft.indices) {
        val harmonicNumber = h + 1
        println("Harmonic: $harmonicNumber...")
        val formant = withAbsAverage(
            average = Ocarina.dft[h],
            input   = simpleFormant(
                amplitude = 1.0,
                frequency = frequency * harmonicNumber,
                bandwidth = frequency,
                booster   = 4.0,
                phase     = phase,
                length    = length
            )
        )

        for(j in formant.indices) {
            underwaveform[j] += formant[j]
        }
    }

    val underwaveformAverage = underwaveform.map { abs(it)}.average()
    val fundamental = sineByFrequency(
        amplitude = underwaveformAverage * 30,
        frequency = frequency,
        size      = underwaveform.size
    )
    return fundamental.add(underwaveform)

   // return underwaveform
}
fun compressionSynth(
    frequency       : Double,
    length          : Int
) : List<Double> {
    val magnitude: List<Double> = exponentialFunction(
        amplitude   = 1.0,
        coefficient = -1.0/64.0,//-1.0/192.0,//-1.0 / 192,//512,//192.0,
        length      = (frequency/2.0).toInt() //+ 1
    )
    println("Magnitude size: ${magnitude.size}")
    val phase : List<Int> = (0 until magnitude.size).map { (0 until length).random()}

    val epiFundamentalMod = synthesizeSpectrum(
        fundamentalFrequency = 1.0,
        magnitude            = magnitude,
        phase                = phase,
        length               = length
    )//.map { abs(it)}
    val epiFundamental = sineByFrequency(
        amplitude = 1.0,
        frequency = frequency,
        size      = length
    ).zip(epiFundamentalMod) { a,b -> a * b}
    val peak = epiFundamentalMod.max()
    println("Peak: $peak")

    //return epiFundamental

    //TODO: generate a matrix of modulation signals using the epi-fundamental envelope
    val percentages = Ocarina.getPercentages()
    println("Percentages sum (should be 100.0): ${percentages.sum()}")
    val harmonicFormantModDfts:List<List<Double>> = compress(
        baseSignal          = magnitude,
        magnitudePercentage = percentages
    )

    //TODO: Verify the compression output sum is equal to the complement of the input signal

    /*validateMatrix(
        baseSignal = epiFundamental,//magnitude,
        results    = harmonicFormantModDfts
    )*/

    //return epiFundamental

    val underwaveform: MutableList<Double> = MutableList(length) { 0.0}
    for(h in harmonicFormantModDfts.indices) {
        println("Harmonic $h synthesizing...")
        val harmonicNumber = h + 1

        val specialOffset = (peak * (percentages[h]/100.0))

        val formantMod = synthesizeSpectrum(
            fundamentalFrequency = 1.0,
            magnitude            = harmonicFormantModDfts[h],
            phase                = phase,
            length               = length
        ).map { specialOffset - it }

        val formant = sineByFrequency(
            amplitude = 1.0,
            frequency = frequency*harmonicNumber,
            size      = length
        ).zip(formantMod) { a,b -> a * b}

        for(j in formant.indices) {
            underwaveform[j] += formant[j]
        }
    }

    return epiFundamental

    val preFundamental =  epiFundamental.add(underwaveform)

    //TODO: Add the epifundamental formant to the underwaveform return the result. With a fundamental the ratio should be 29 instead of 30

    val underwaveformAverage = underwaveform.map { abs(it)}.average()
    val fundamental = sineByFrequency(
        amplitude = underwaveformAverage * 30,
        frequency = frequency,
        size      = underwaveform.size
    )
    return fundamental.add(preFundamental)
}

fun compress(
    baseSignal: List<Double>,
    magnitudePercentage: List<Double>
) : List<List<Double>> {
    val peak              = baseSignal.map { abs(it)}.max()
    var currentComplement = baseSignal.map { peak - it } //TODO: same except baseSignal[0] = peak, each compression signal should be negated
    val magnitudeSum      = currentComplement.sum()//baseSignal.sum() / (magnitudePercentage[0]/100.0)//currentComplement.sum()
    val output            = mutableListOf<List<Double>>()

    //TODO: How would you seriously do this with functional programming?
    for(h in magnitudePercentage.indices) {
        println("H[$h]: ${magnitudePercentage[h]}...")
        /** The last compression signal should take whatever is remaining but how does this fit into the amplitude scaling?
         * My guess is that if the harmonics used in the magnitudes list really encompass 100% of the magnitude then the remainder
         * will match the magnitude value of magnitude[N -1] identically.
         *
         * **/
        if(h == magnitudePercentage.lastIndex) {
            output.add(currentComplement)
            break
        }

        var currentCapacityTotal  = currentComplement.sum()
        var currentMagnitudeTotal = magnitudeSum * (magnitudePercentage[h] / 100.0)
        //println("currentMagnitudeTotal: $currentMagnitudeTotal")

        val compressionSignal: MutableList<Double> = mutableListOf()
        for(i in currentComplement.indices) {
            currentCapacityTotal -= currentComplement[i]
            val choice = if (i != currentComplement.lastIndex) {
                val min = currentMagnitudeTotal - currentCapacityTotal
                val max = max(currentComplement[i], currentMagnitudeTotal)
                min + (Math.random() * (max - min))
            } else {
                currentMagnitudeTotal
            }
            compressionSignal.add(choice)
            currentMagnitudeTotal -= choice
        }

        output.add(compressionSignal)
        currentComplement = currentComplement.zip(compressionSignal) { a, b -> a - b}
    }
    return output
}
fun compressionSpread(data: List<Double>, baseFrequency: Double, numHarmonics: Int) : List<Double> {
    val fout = FilterUtils.bandpassFIR(
        input                 = data,
        minFreqHz             = 0.1,//(444.0) - (444/2.0),
        maxFreqHz             = baseFrequency + (baseFrequency/2.0),
        impulseResponseLength = Constants.SAMPLE_RATE + 1
    )
    val windows = getWaves(data = fout).map { it.size}
    val dfts: List<MutableList<Double>> = List(windows.size) { (0 until numHarmonics).map { 0.0}.toMutableList()}

    val formants: MutableList<List<Double>> = mutableListOf()
    for(i in 0 until numHarmonics) {
        val harmonicNumber = i + 4
        println("HarmonicNumber: $harmonicNumber...")
        val fout = FilterUtils.bandpassFIR(
            input                 = data,
            minFreqHz             = (harmonicNumber*baseFrequency) - (baseFrequency/2.0),
            maxFreqHz             = (harmonicNumber*baseFrequency) + (baseFrequency/2.0),
            impulseResponseLength = Constants.SAMPLE_RATE + 1
        )
        formants.add(fout)
    }

    //TODO: How to group the values over the given windows
    for(formantNumber in formants.indices) {
        val currentFormant = formants[formantNumber]
        var pos = 0
        for(windowIndex in 0 until windows.size) {
            val w = windows[windowIndex]
            dfts[windowIndex][formantNumber] = currentFormant.subList(pos, pos + w).map {abs(it)}.average()
            pos += w
        }
    }
    val dftsNormalized = dfts.map { normalize(scale = 5000.0, input = it)}
    val spread         = dftsNormalized.map { it.sum()}
    val peaks          = dftsNormalized.map { findPeakPosition(input = it)}

    listData(data = peaks)

    return spread
}

/*fun complexFormant(
    frequency: Double,
    bandwidth: Double,
    length: Int
) : List<Double> {
    val outerFormant = simpleFormant(
        amplitude = 5000.0,
        bandwidth = bandwidth,
        length    = length ,
        booster   = 0.0,
        frequency = frequency
    )

    val innerFormant = carrierSig(centerFrequency = frequency, length = length)

    val preFilterResult = outerFormant.add(
        withAbsAverage(
            average = outerFormant.map {abs(it)}.average() * 4.0,
            input   = innerFormant
        )
    )

    return FilterUtils.bandpassFIR(
        input                 = preFilterResult,
        minFreqHz             = (frequency) - (bandwidth/2.0),
        maxFreqHz             = (frequency) + (bandwidth/2.0),
        impulseResponseLength = Constants.SAMPLE_RATE + 1
    )
}*/

fun extractVibratoSignal(targetFrequency: Double, data: List<Double>) : List<Double> {
    val fout = FilterUtils.bandpassFIR(
        input                 = data,
        minFreqHz             = (targetFrequency) - (444.0/2.0),
        maxFreqHz             = (targetFrequency) + (444.0/2.0),
        impulseResponseLength = Constants.SAMPLE_RATE + 1
    )
    val norm   = normalizeWaves(scale = 5000.0, input = fout)
    val waves  = getWaves(data = norm)
    val sums   = waves.map {it.map{abs(it)}.sum()}
    val re     = relativeEnergy(centerFrequency = targetFrequency, normalizedEnergies = sums)
    val dist   = re.map { (100.0 - it)}
    //val freqs  = re.map { ((it*0.95)/100.0)*targetFrequency}
    val f       = targetFrequency
    val freqs   = dist.map { f - ((it/100.0)*f)}
    return halfPeriodFrequenciesToWaves(frequencies = freqs)
}
fun relativeEnergy(
    centerFrequency: Double,
    normalizedEnergies: List<Double>
): List<Double> {
    val halfPeriod = singleHalfPeriod(amplitude = 5000.0, frequency = centerFrequency).sum()
    println("halfPeriod: ${halfPeriod}, example: ${normalizedEnergies[400]}")
    return normalizedEnergies.map {
        (it / halfPeriod) * 100.0
    }
}
fun carrierSig(
    centerFrequency: Double,
    length         : Int
) : List<Double> {
    val periods  = ((length / (Constants.SAMPLE_RATE / centerFrequency))*2).toInt()

    val m: List<Double> = exponentialFunction(
        amplitude = 5000.0,
        coefficient = -1.0/192.0,//-1.0 / 192,//512,//192.0,
        length      = (periods/2).toInt()
    ).toMutableList()
    //val m        = listOf(0.0) + (0..(periods/2).toInt()).map {5000.0}
    val idft     = inversePolarWithRandomizedPhase(m = m, length = periods)

    val average  = idft.average()
    val centered = idft.map { it - average}

    val desiredAbsAverage = (30.0/888.0)  * centerFrequency // This is used because 30 at 888.0 sounded good during testing

    val scaled = withAbsAverage(average = desiredAbsAverage, input = centered)

    val frequencies =  scaled.map { it + centerFrequency}

    //val std = standardDeviation(data = frequencies)
    //println("STD: $std")

    return halfPeriodFrequenciesToWaves(frequencies)

    //val fsig    = sineByCycles(amplitude = centerFrequency*0.10, cycles = 1.0, size = periods.toInt()).map { it + centerFrequency}

    /*return halfPeriodFrequenciesToWaves(
        frequencies  = sineByCycles(
            amplitude = centerFrequency*0.05,
            cycles    = 0.5,
            size      = periods.toInt()
        ).map { it + centerFrequency }
    )*/

    /*val o       = mutableListOf<Double>()
    var direction = 1.0
    while(o.size < length) {
        val frequency = options[options.indices.random()]
        o.addAll(singleHalfPeriod(amplitude = direction, frequency = frequency))
        direction *= -1
    }
    return o*/
}

/*fun synth(
    numHarmonics    : Int,
    frequency       : Double,
    length          : Int
) : List<Double> {
    val underwaveform: MutableList<Double> = MutableList(length) { 0.0}

    for(h in 2 until numHarmonics) {
        if(h-2 == Ocarina.dft.size) {
            println("Not enough dft values for the desired number of harmonics. Ending early...")
            break
        }
        val currentAmplitude = Ocarina.dft[h - 2]
        val currentFrequency = h*frequency
        println("Current harmonic $h at ${currentFrequency}...")
        if(currentFrequency > (Constants.SAMPLE_RATE/2.0)) {
            println("CurrentFrequency: $currentFrequency, Your Nyquist math is off!!!! ending early")
            break
            //throw RuntimeException("CurrentFrequency: $currentFrequency, Your Nyquist math is off!!!!")
        }

        /*val formant = if(h == 3) {
            sineByFrequency(
                amplitude = 5000.0,
                frequency = currentFrequency,
                size      = length
            )
        } else {
            formant(
                amplitude = 5000.0,
                frequency = currentFrequency,
                bins      = 200,
                length    = length
            )
        }*/
        val formant = complexFormant(
            frequency = currentFrequency,
            bandwidth = frequency,
            length    = length,
        )

        val scaledFormant = withAbsAverage(
            average = currentAmplitude,
            input   = formant
        )

        /** TODO Verify this is working with simple formant equation
         * that operates correctly regardless of whether its sound is authentic **/
        for(j in scaledFormant.indices) {
            if(j == underwaveform.size - 1) {
                break //TODO:
                underwaveform.add(scaledFormant[j])
            } else {
                underwaveform[j] += scaledFormant[j]
            }
        }
    }

    val underwaveformAverage = underwaveform.map { abs(it)}.average()
    val fundamental = sineByFrequency(
        amplitude = underwaveformAverage * 30,
        frequency = frequency,
        size      = underwaveform.size
    )
    return fundamental.add(underwaveform)
}
*/
fun bandwidthAnalysisList(
    data     : List<Double>,
    frequency: Double,
    bandwidth: Double,
    harmonics: Int
) : List<Double> {
    return (1..harmonics).map {
        bandwidthAnalysis(
            data = data,
            frequency = frequency * it,
            bandwidth = bandwidth
        )
    }
}
fun bandwidthAnalysis(
    data     : List<Double>,
    frequency: Double,
    bandwidth: Double
) : Double {
    val spectrum = inharmonicDft(
        x                 = data,
        baseFrequency     = frequency,
        bandwidth         = bandwidth.toInt(),
        frequencyDistance = 1.0
    ).magnitude
    val left  = spectrum.subList(0,spectrum.size/2).reversed()
    val right = spectrum.subList(spectrum.size/2, spectrum.size)
    val dft   = normalize(scale = 5000.0, input = left.add(right))

    /*val peakIndex = findPeakPosition(input = dft)
    val peakPower = pow(dft[peakIndex], 2.0)
    var dropPosition = -1
    for(i in dft.indices) {
        if(i == peakIndex) { continue }
        val power = pow(dft[i], 2.0)
        if(power <= peakPower * 0.1) {
            dropPosition = i
            break
        }
    }*/
    //listData(data = dft)

    println("lower band sum: ${dft.subList(0, 10).sum()/5000.0}, Rest: ${dft.subList(10, dft.size).sum()/5000.0}")

    val bw = dft.sum()/5000.0
    println("f: $frequency, BW: $bw, sum ${dft.sum()}")
    return bw
}

/*fun formant(
    amplitude: Double,
    frequency: Double,
    bins: Int,
    length: Int
) : List<Double> {
    val magnitude: MutableList<Double> = exponentialFunction(
        amplitude = amplitude,
        coefficient = -1.0/192.0,//-1.0 / 192,//512,//192.0,
        length      = bins
    ).toMutableList()

    val initialModulationSignal = synthesizeSpectrum(
        fundamentalFrequency = 1.0,
        magnitude            = magnitude,
        length               = length
    )

    val dc = initialModulationSignal.map { abs(it)}.max() * 1.0
    val modulationSignal = initialModulationSignal.map { it + dc}

    return modulationSignal

    /*val carrierSignal = sineByFrequency(
        amplitude = 1.0,
        frequency = frequency,
        size      = length
    )*/
    val carrierSignal = carrierSig(
        centerFrequency = frequency,
        length          = length
    )

    bandwidthAnalysis(
        data      = carrierSignal,
        frequency = frequency,
        bandwidth = 444.0
    )

    return carrierSignal.zip(modulationSignal) { a, b -> a * b}
}*/


fun filterDft(
    initialHarmonic: Int,
    fundamental : Double,
    lastHarmonic: Int,
    data        : List<Double>
)  {
    val avgs: MutableList<Double> = mutableListOf()
    var harmonicNumber = initialHarmonic
    while(harmonicNumber <= lastHarmonic) {
        val result = FilterUtils.bandpassFIR(
            input                 = data,
            minFreqHz             = harmonicNumber * fundamental - (fundamental / 2.0),
            maxFreqHz             = harmonicNumber * fundamental + (fundamental / 2.0),
            impulseResponseLength = Constants.SAMPLE_RATE + 1
        )
        val amplitude = result.map{abs(it)}.average()
        println(amplitude)
        avgs.add(amplitude)
        harmonicNumber++
    }

    println("scaled")
    val scaled = normalize(scale = 5000.0, input = avgs)
    listData(data = scaled)

    println("list")
    for(s in scaled) {
        println("$s,")
    }
}

fun withAbsAverage(
    average: Double,
    input: List<Double>
) : List<Double> {
    val absAverage = input.map { abs(it)}.average()
    return input.map { it * (average/absAverage)}
}

/**
 *
 * coefficient 128 produces the loud clap I'm often looking for. 96 sounds more like it when the time window is extended
 * and 192 just seems over the top
 *
 * **/
private fun simpleFormant(
    amplitude: Double,
    frequency: Double,
    bandwidth: Double,
    booster: Double,
    coefficient: Double = -1.0/(64.0),//-1.0/128.0,//96.0,
    length:    Int,
    phase: List<Int>
) : List<Double> {
    val magnitude = exponentialHillSignal(
        amplitude   = amplitude,
        coefficient = coefficient,//-1.0/48.0,//-1.0/192.0,
        length      = bandwidth.toInt() + 1
    )

    println("Magnitude.size: ${magnitude.size}")
    println("phase.size: ${phase.size}")

    val outerFormant = synthesizeInharmonicDft(
        magnitude         = magnitude,
        baseFrequency     = frequency,
        bandwidth         = bandwidth,
        frequencyDistance = 1.0,
        phaseOffset       = phase,
        length            = length
    )
    //return outerFormant
    //Default booster is 4
    val innerFormant = sineByFrequency(
        //offset = (0 .. length).random().toDouble(),
        amplitude = booster * outerFormant.map { abs(it)}.average(),//magnitude.sum() * (1/6.0),
        frequency = frequency,
        size      = length
    )
    return innerFormant.add(outerFormant)
}

private fun exponentialHillSignal(
    amplitude:   Double,
    coefficient: Double,
    length:      Int
) : List<Double> {
    val signal = exponentialFunction(
        amplitude   = amplitude,
        coefficient = coefficient,
        length      = length / 2
    )
    return if(length % 2 == 0) {
        signal.reversed() + signal
    } else {
        signal.reversed() + amplitude + signal
    }
}

private fun exponentialFunction(amplitude: Double, coefficient: Double, length: Int) : List<Double> {
    val o: MutableList<Double> = mutableListOf()
    for(i in 0 until length) {
        val value = amplitude* Math.E.pow(coefficient * i.toDouble())
        o.add(value)
    }
    return o
}
