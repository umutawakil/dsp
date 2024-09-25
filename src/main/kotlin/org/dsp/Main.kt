package org.dsp

import org.dsp.analysis.DiscreteFourierTransform
import org.dsp.file.FileUtils
import org.dsp.filters.FilterUtils
import org.dsp.wave.WaveUtils
import java.io.*
import kotlin.math.*
//import java.util.Random
import kotlin.random.Random


fun main() {
    val resourceDirectory     = "/Users/umutawakil/Documents/Git/dsp/src/main/resources"
    val readFile              = File("$resourceDirectory/normalized.wav")//image-only-idft.wav") //normalized-idft-1.wav,normalized.wav,test-sample.wav
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
    val waveLength: MutableList<Double>  = mutableListOf()
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

    for (i in 0 until waveLengths.size step 2) {
        if(i + 1 >= waveLengths.size) { break }
        waveLength.add(waveLengths[i] + waveLengths[i + 1])
    }

    val waves: MutableList<MutableList<Double>> = mutableListOf()
    var l = 0
    for(i in waveLengths.indices) {
        val temp: MutableList<Double> = mutableListOf()
        for(j in 0 until waveLengths[i].toInt()) {
            temp.add(s[l].toDouble())
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
    val full: MutableList<Double> = mutableListOf()

    for (i in waveLengths.indices) {
        if (i % 2 == 0) {
            //println("e($i) v: ${waveLengths[i]}")
            even.add((waveLengths[i] - 52.0).toDouble())

        } else {
            odd.add((waveLengths[i] - 52.0).toDouble())
        }
        full.add(waveLengths[i] - 52.0)
    }

   /** ----------------------------------------------------------------**/

    val newWaves: MutableList<List<Double>> = mutableListOf()
    for(i in 0 until waves.size step 2) {
        if(i + 1 == waves.size) { break }
       newWaves.add(waves[i] + waves[i + 1])
    }
    val o: MutableList<Double> = mutableListOf()
    for(i in newWaves.indices) {
        /*val dft = DiscreteFourierTransform.dftRectangular(x = newWaves[i])
        val timeDomain = DiscreteFourierTransform.inverseDftRectangular(size = 100, R = dft.first, I = dft.second)
        o.addAll(timeDomain)*/

        o.addAll(interpolate(waveform = newWaves[i], delta = 100 - newWaves[i].size))
    }

    val rawWaves: MutableList<List<Double>> = mutableListOf()//getWaves(data = o)
    for(i in 0 until o.size step 50) {
        rawWaves.add(o.subList(i, i + 50))
    }

    val staticWaves: MutableList<List<Double>> = mutableListOf()
    for(i in 0 until rawWaves.size step 2) {
        if(i + 1 == rawWaves.size) { break }
        staticWaves.add(rawWaves[i] + rawWaves[i + 1])
    }

    val waveDiff:MutableList<List<Double>> = mutableListOf()
    for(i in 1 until staticWaves.size) {
        waveDiff.add(staticWaves[i].zip(staticWaves[i - 1]) { a , b -> a - b})
    }

    val genesisWave       = staticWaves[0]
    val transferAmplitude = 100.0
    val pitchedNoiseWaves = pitchedNoiseWaves(amplitude = transferAmplitude, waveLengths = (0..((waveLengths.size/2))).map {100.0})
    val delayLineNoise    = FilterUtils.delayLine(gain = -0.9, delay = 100, x = whiteNoise(range = transferAmplitude, size = 100*400), iterations = 1).map { transferAmplitude *it}

    val dNoise: MutableList<List<Double>> = mutableListOf()
    for(i in 0 until delayLineNoise.size step 100) {
        dNoise.add(delayLineNoise.subList(i, i + 100))
    }

    val integratedWaves: MutableList<List<Double>> = mutableListOf()
    val periodDerivative = waveDiff//.map { normalize(scale= transferAmplitude, input = it)}//dNoise

    /** Integrate the static signal with the period to period derivative signal (This is a key step)
     *
     */
    var tempIntegrationBuff: List<Double> = genesisWave.toList()
    for(i in periodDerivative.indices) {
        tempIntegrationBuff = tempIntegrationBuff.zip(periodDerivative[i]) { a, b -> a + b}
        integratedWaves.add(tempIntegrationBuff)
    }

    val stretchedWaves: MutableList<List<Double>> = mutableListOf()
    var externalWaveLengthSourceCounter = 0
    for (i in integratedWaves.indices) {
        val periodWaves = listOf(
            integratedWaves[i].subList(0,integratedWaves[i].size/2),
            integratedWaves[i].subList(integratedWaves[i].size/2, integratedWaves[i].size)
        )
        if ((periodWaves[0].size != periodWaves[1].size)) {
            throw RuntimeException("Period wave sizes do not match correctly. p0: ${periodWaves[0].size}, p1: ${periodWaves[1].size}")
        }
        if(periodWaves[0].size != 50) {
            listData(input = integratedWaves[i])
            throw RuntimeException("Period($i) -> Period wave size is incorrect: ${periodWaves[0].size}, integratedWave: ${integratedWaves[i].size}")
        }
        if(periodWaves.size != 2) { break }

        if(externalWaveLengthSourceCounter + 1 >= waveLengths.size) { break }
        stretchedWaves.add(
            interpolate(
                waveform = periodWaves[0],
                delta = (waveLengths[externalWaveLengthSourceCounter] - periodWaves[0].size).toInt()
            )
        )
        stretchedWaves.add(
            interpolate(
                waveform = periodWaves[1],
                delta = (waveLengths[externalWaveLengthSourceCounter + 1] - periodWaves[1].size).toInt()
            )
        )
        externalWaveLengthSourceCounter += 2
    }

    /*val harmonicSignal = harmonicSynthesizer(
        fundamentalAmplitude = 5000.0,
        basePeriodLength     = 100,
        baseDft              = dftAvg,
        waveLengths          = waveLengths
    )*/


    writeSamplesToFile(
        fileName = "$resourceDirectory/output.wav",
        input    =  normalize(
                        scale = 5000.0,
                        input = randomSynth(testPeriod = staticWaves[200], waveLength = 100, size = 100*400, vibrato = waveLengths)//stretchedWaves.flatten()//waveDiff.map{it + listOf(2500.0)}.flatten()//pitchedNoiseWaves.flatten()//integratedWaves.flatten()//staticTestSignalWaves.flatten()
                    )
    )
    return
}

fun staticWaveStats(staticWaves: List<List<Double>>) {
    val dftAvg: MutableList<Double> = mutableListOf()
    val numHarmonics = DiscreteFourierTransform.dft(x= staticWaves[0]).size
    for(h in 0 until numHarmonics) {
        var temp = 0.0
        for(i in staticWaves.indices) {
            val dft = DiscreteFourierTransform.dft(x = staticWaves[i])
            temp += dft[h]
        }
        temp /= staticWaves.size
        dftAvg.add(temp)
    }

    val dft = dftAvg//DiscreteFourierTransform.dft(x = staticWaves[200])
    val sum1: Double = dft.subList(2, 11).sum()
    val avgX1 = sum1 / 9.0//(dft.size - 11.0)
    println("Lower Harmonics (Excluding DC and fundamental from division) -> Sum: $sum1, avg: $avgX1")

    val sum2: Double = dft.subList(11, dft.size).sum()
    val avgX2 = sum2 / (dft.size - 9.0)
    println("Higher Harmonics (Excluding DC and fundamental from division) -> Sum: $sum2, avg: $avgX2")

    println("DFT average")
    listData(input = dftAvg)
    println("DFT of static wave 200")
    listData(input = DiscreteFourierTransform.dft(x = staticWaves[200]))
}

//TODO: A dft might be best for this then something else for the sub harmonics
//
fun randomSynth(testPeriod: List<Double>, waveLength: Int, size: Int, vibrato: List<Double>) : List<Double> {
    var o: MutableList<Double> = MutableList(size) { 0.0 }
    val harmonics = listOf(
        0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0
        /*15.0,
        16.0,
        17.0,
        18.0,
        19.0,
        20.0*/
    )

    for(h in harmonics) {
        val currentWaveLength: Int = (waveLength/h).toInt()
        val numPeriods = (size / currentWaveLength).toInt()
        //println("WaveLength: $currentWaveLength, harmonic: $h, numPeriods: $numPeriods")
        for(p in 0 until numPeriods) {
            val amplitude = if(h == 1.0) { 0.0 } else { 100.0 } //1.0
            val startingPosition = p*currentWaveLength
            for(i in 0 until  currentWaveLength) {
                o[i + startingPosition] += amplitude*sin((2*Math.PI*i)/currentWaveLength)
            }
        }
    }

    val periods = o.size / testPeriod.size
    o  = mutableListOf()
    repeat(periods) {
        o.addAll(testPeriod)
    }

    val waves: MutableList<List<Double>> = mutableListOf()
    var c = 0
    val half = waveLength/2
    var direction = 1.0
    while(true) {
        if(c + half - 1 >= o.size) {
            break
        }
        val halfWave = o.subList(c, c + half)
        /*val normalizedUnderWaveform = normalize(scale = 100.0, input = halfWave)
        val fundamentalWaveform = sineByCycles(amplitude = direction*5000.0, cycles = 0.5, size = half)
        waves.add(fundamentalWaveform.zip(normalizedUnderWaveform) { a, b -> a + b})*/
        waves.add(halfWave)
        c += half
        direction *= -1
    }

    val interpolatedWaves: MutableList<List<Double>> = mutableListOf()
    for(w in vibrato.indices) {
        if(w == waves.size) {
            println("Breaking early!!! Vibrato signal and number of waves don't match. w: $w, waves.size: ${waves.size}")
            break
        }
        val delta = (vibrato[w] - waves[w].size).toInt()
        //println("$w, Delta: $delta, waveSize: ${waves[w].size}")
        interpolatedWaves.add(interpolate(waveform = waves[w], delta = delta))
    }
    return interpolatedWaves.flatten()
}

fun harmonicSynthesizer(
    fundamentalAmplitude: Double,
    basePeriodLength: Int,
    baseDft: List<Double>,
    waveLengths: List<Double>
) : List<Double> {
    val o: MutableList<Double> = mutableListOf()
    var dftBuffer: MutableList<Double> = baseDft.toMutableList()

    //listData(input = dftBuffer)

    //return emptyList()
    dftBuffer = normalize(scale = 5000.0, input = dftBuffer).toMutableList()

    for(i in 0 until waveLengths.size step 2) {
        dftBuffer[0] = 0.0
        dftBuffer[1] = 0.0



        //TODO: synthesize a period recursively saving the DFT to reduce los
        val transferSignal = generateRandomDftFromDFT(amplitude = 100.0, rangePercent = 85.0, inputDFT = dftBuffer)
        //listData(input = transferSignal)
        plotSignal(signal = DiscreteFourierTransform.inverseDFTRaw(m = transferSignal))
        println()
       // dftBuffer = normalize(scale = 5000.0, input = dftBuffer.zip(transferSignal) { a, b -> a + b}).toMutableList()
        //dftBuffer = dftBuffer.zip(transferSignal) { a, b -> a + b}.toMutableList()

        dftBuffer[0] = 0.0
        dftBuffer[1] = 0.0//fundamentalAmplitude
        val timeDomain: List<Double> = DiscreteFourierTransform.inverseDftWithLength(m = dftBuffer, length = basePeriodLength)

        //TODO: Interpolate the parts with the wavelength signal

        //val periodLength = waveLengths[i] + waveLengths[i + 1]

        //listData(input = timeDomain + listOf(300.0))
        o.addAll(timeDomain)
    }

    return o
}

private fun generateRandomDftFromDFT(
    amplitude: Double,
    rangePercent: Double,
    inputDFT: List<Double>
) : List<Double> {
    val o: MutableList<Double> = mutableListOf()
    for(i in inputDFT.indices) {
        /*val rangeInt = rangePercent.toInt()
        val percentChange = (-rangeInt..rangeInt).random()/100.0
        val result =  inputDFT[i] + (percentChange*inputDFT[i])
        println("Prev: ${inputDFT[i]}, New: ${result}")
        o.add(result)*/

        o.add((0..1000).random().toDouble())
    }
   // listData(input = o)
    println()
    return normalize(scale = amplitude, input = o)
}

fun whiteNoise(size: Int, mean: Double = 0.0, stdDev: Double = 1.0): List<Double> {
    val random = java.util.Random()
    return List(size) {
        random.nextGaussian() * stdDev + mean
    }
}

fun cosineByCycles(amplitude: Double, cycles: Double, size: Int) : List<Double> {
    return (0 until size).map { (-1*amplitude * cos((2*Math.PI*cycles*(it + size)/size))) + amplitude }
}

fun sineByCycles(amplitude: Double, cycles: Double, size: Int) : List<Double> {
    return (0 until size).map { amplitude * sin((2*Math.PI*cycles*it)/size)}
}

fun interpolate(waveform: List<Double>, delta: Int): List<Double> {
    if (delta == 0) return waveform

    val outputSize = waveform.size + delta
    if (outputSize <= 0) {
        println("Empty List WTF!!!: $outputSize")
        return emptyList()
    }

    val result = mutableListOf<Double>()
    val inputSize = waveform.size.toDouble()

    for (i in 0 until outputSize) {
        val inputIndex = i * (inputSize - 1) / (outputSize - 1)
        val lowerIndex = inputIndex.toInt()
        val upperIndex = minOf(lowerIndex + 1, waveform.lastIndex)

        if (lowerIndex == upperIndex) {
            result.add(waveform[lowerIndex])
        } else {
            val fraction = inputIndex - lowerIndex
            val interpolatedValue = waveform[lowerIndex] * (1 - fraction) + waveform[upperIndex] * fraction
            result.add(interpolatedValue)
        }
    }

    return result
}

/** TODO: Shouldn't amplitude modulation be considered? At low amplitudes that the subharmonic waveforms operate it
 * there indeed seems to be amplitude modulation and it's not perceived.
 */

fun subharmonicSynth(staticSignal: List<Double>, waveLength: Int, size: Int) : List<Double> {
    val noiseWaves = pitchedNoiseWaves(waveLength = waveLength, periods = size/waveLength)
    val o: MutableList<Double> = mutableListOf()

    var temp = normalize(scale = 5000.0, input = staticSignal).toMutableList()
    for(i in noiseWaves.indices) {
        println("N: ${noiseWaves[i].size}, temp: ${temp.size}")
        o.addAll(temp.zip(normalize(scale = 1000.0, input = noiseWaves[i])) { a, b -> a + b }.toMutableList())
        //o.addAll(temp)
        //temp = temp.zip(normalize(scale = 1000.0, input = noiseWaves[i])) { a, b -> a + b }.toMutableList()
    }

    return o
}

fun pitchedNoiseWaves(amplitude: Double, waveLengths: List<Double>) : List<List<Double>> {
    if(waveLengths.size % 2 != 0) { throw RuntimeException("Number of wavelengths supplied must be even. Num waves: ${waveLengths.size}") }

    val o: MutableList<List<Double>> = mutableListOf()
    for(i in 0 until waveLengths.size step 2) {
        if(i + 1 == waveLengths.size) { break }
        val wA = waveLengths[i]
        val wB = waveLengths[i + 1]
        val noiseA = whiteNoise(range = amplitude, size = wA.toInt())//.map { amplitude*(it + 1.0) }
        val noiseB = whiteNoise(range = amplitude, size = wB.toInt()).map { amplitude*(it + 1.0) }
        o.add(noiseA)
        o.add(noiseB.map {-1.0*it})
    }
    return o
}
fun pitchedNoiseWaves(waveLength: Int, periods: Int) : List<List<Double>> {
    val halfPeriod = waveLength/2
    val o: MutableList<List<Double>> = mutableListOf()
    while(o.size < periods) {
        val noise = whiteNoise(range = 1.0, size = halfPeriod)//.map { it + 1.0 }
        o.add(noise)
        o.add(noise.map {-1.0*it})
    }
    return o
}
/*fun pitchedNoise(waveLength: Int, size: Int) : List<Double> {
    val halfPeriod = waveLength/2
    val o: MutableList<Double> = mutableListOf()
    while(o.size < size) {
        o.addAll(whiteNoise(size = halfPeriod).map { it + 1.0 } )
        o.addAll(whiteNoise(size = halfPeriod).map {-1.0*(it + 1.0)})
    }
    return o
}*/

fun mirrorPhaser(waveLengths: List<Double>): List<Double> {
    /** Pretty cool phaser I must say **/
    val o: MutableList<Double> = mutableListOf()
    val sizer = 100*800
    var oi = 0
    var inc = 0
    while(oi < sizer) {
        if(inc + 1 >= waveLengths.size) { break  }
        val wlenA = waveLengths[inc].toInt()
        val wlenB = waveLengths[inc + 1].toInt()

        //TODO: Add the half sine multiplication after confirming all is well

        val pos:List<Double>  = normalize(scale = 5000.0, input = brownNoise(size = wlenA).map{abs(it)})
        val neg: List<Double> = normalize(scale = 5000.0, input = brownNoise(size = wlenB).map { abs(it) *-1})
        o.addAll(pos)
        o.addAll(neg)
        oi += pos.size + neg.size
        inc += 2
    }
    return o
}

fun halfSine(size: Int) : List<Double> {
    return (0 until size).map {
        sin(
            (2.0*Math.PI*it)/(2.0*size)
        )
    }
}

fun whiteNoise(range: Double, size: Int) : List<Double> {
    val o: MutableList<Double> = MutableList(size) {0.0}
    for(i in o.indices) {
        o[i] = (0..range.toInt()).random().toDouble()
    }
    val avg = o.average()

    return o.map { (it - avg) + range}
}

fun brownianNoise(size: Int): List<Int> {
    val noise = mutableListOf<Int>()
    var currentValue = 5 // Start in the middle of the range

    noise.add(currentValue)

    for (i in 1 until size) {
        // Generate a step of -1, 0, or 1
        val step = Random.nextInt(-1, 2)

        // Update the current value and constrain it to the range [0, 10]
        currentValue = (currentValue + step).coerceIn(0, 10)

        noise.add(currentValue)
    }

    return noise
}

fun brownianNoise(size: Int, stepSize: Double, minValue: Double = 0.0, maxValue: Double = 10.0): List<Double> {
    val noise = mutableListOf<Double>()
    var currentValue = (minValue + maxValue) / 2 // Start in the middle of the range

    noise.add(currentValue)

    for (i in 1 until size) {
        // Generate a step that's exactly -stepSize, 0, or +stepSize
        val step = when (Random.nextInt(3)) {
            0 -> -stepSize
            1 -> 0.0
            else -> stepSize
        }

        // Calculate the new value
        var newValue = currentValue + step

        // If the new value is outside the range, clip it to the nearest boundary
        newValue = newValue.coerceIn(minValue, maxValue)

        currentValue = newValue
        noise.add(currentValue)
    }

    return noise
}

fun brownNoise(size: Int): List<Double> {
    val noise = mutableListOf<Double>()
    var lastValue = 0.0
    val scale = sqrt(1.0 / size)

    for (i in 0 until size) {
        // Generate a random value between -1 and 1
        val whiteNoise = Random.nextDouble(-1.0, 1.0)

        // Integrate the white noise
        lastValue += whiteNoise * scale

        // Ensure the values stay within a reasonable range
        lastValue = lastValue.coerceIn(-1.0, 1.0)

        noise.add(lastValue)
    }

    val avg = noise.average()

    return noise.map{it - avg}
}

//TODO: The organic waveLength extraction algorithm may be leaving off the last half period
fun getWaves(data: List<Double>) : List<List<Double>> {
    val o: MutableList<List<Double>> = mutableListOf()
    var direction = if (data[0] >= 0) { 1 } else { -1 }
    var temp: MutableList<Double> = mutableListOf()
    for(i in data.indices) {
        if(data[i] * direction < 0 ) {
        //if(data[i] * direction <= 0 && (i > 0)) {
            direction *= -1
            o.add(temp)
            temp = mutableListOf()
        }
        temp.add(data[i])
    }

    //TODO: Currently you can lose a wave at the end
    /*if(temp.size != 0) {
        println("T: ${temp.size}, o: ${o.size}")
        o.add(temp)
    }*/
    return o
}

//TODO: The organic waveLength extraction algorithm may be leaving off the last half period
fun getPeriods(data: List<Double>) : List<Int> {
    val o: MutableList<Int> = mutableListOf()
    var sign = data[0]
    var currentWaveLength = 0

    for (i in data.indices) {
        val currentSample = data[i]

        if(currentSample * sign >= 0) {
            currentWaveLength++
        } else {
            o.add(currentWaveLength)
            currentWaveLength = 1
        }
        if(currentSample != 0.0) {
            sign = currentSample
        }
    }
    if(currentWaveLength != 1) {
        o.add(currentWaveLength)
    }
    return o
}
fun findPeakPosition(input: List<Double>) : Int {
    var peakIndex = 0
    for(i in input.indices) {
        if(abs(input[peakIndex]) < abs(input[i])) {
            peakIndex = i
        }
    }
    return peakIndex
}
fun findPeak(input: List<Double>) : Double {
    return input[findPeakPosition(input = input)]
}
fun normalize(scale: Double, input: List<Double>) : List<Double> {
    var max = 0.0
    for(i in input.indices) {
        if(abs(max) < abs(input[i])) {
            max = abs(input[i])
        }
    }
    if(max == 0.0) { return  input }
    return input.map { scale * (it/max)}
}
fun listData(input: List<*>) {
    for(i in input.indices) {
        println("$i ${input[i]}")
    }
}
fun createPitchNote(windowSizeInSamples: Int, base: Double) : List<Double> {
    val atomicVibrato = calculateAtomicVibrato(
        scale               = 1.0,
        windowSizeInSamples = windowSizeInSamples,
        base                = base//52.0
    )
    val so: MutableList<Double> = mutableListOf()
    for(i in atomicVibrato.first.indices) {
        if(i >= atomicVibrato.second.size) {break } //TODO: Can they be equalized?
        so.add(atomicVibrato.first[i])
        so.add(atomicVibrato.second[i])
    }
    return so
}

fun generateVibratoSignal(noBase: Boolean, base: Double,size: Int) : List<Double> {
    val INTERNAL_BASE = 52.0
    val paritySignal = generateParitySumSignal(size = size/2)
    val diffSignal   = (0 until size).map { (0..1).random()}
    val output: MutableList<Double> = mutableListOf()

    for(i in paritySignal.indices) {
        if(diffSignal[i] == 0) {
            val even = (paritySignal[i] / 2).toInt()
            val odd  = (paritySignal[i] / 2).toInt()
            output.add(even.toDouble())
            output.add(odd.toDouble())
        } else {
            var diff = 0
            while(diff == 0) {
                diff = (-1..1).random()
            }
            val even = (paritySignal[i] / 2).toInt() + diff
            val odd  = paritySignal[i] - even
            output.add(even.toDouble())
            output.add(odd)
        }
    }

    if(noBase) {
        for (i in output.indices) {
            output[i] =  (output[i] / INTERNAL_BASE)
        }
    } else {
        for (i in output.indices) {
            output[i] = (base * (output[i] / INTERNAL_BASE)) + base
        }
    }

    return output
}

//TODO: Should be able to choose what span of the histogram is desired if the initial envelope is not desired
//TODO: As well as a starting point
fun generateParitySumSignal(size: Int) : List<Double> {
    val histogramMap = mapOf(
        16.0 to 1,
        15.0 to 2,
        14.0 to 1,
        13.0 to 1,
        12.0 to 1,
        11.0 to 3,
        10.0 to 1,
        9.0 to 5,
        8.0 to 20,
        7.0 to 35,
        6.0 to 32,
        5.0 to 123,
        4.0 to 144,
        3.0 to 23,
        2.0 to 4,
        1.0 to 1,
        0.0 to 1
    )
    val histogramList = histogramMap.toList()

    val STEP_SIZE = 1.0
    val options: MutableList<Double> = mutableListOf()
    for(p in histogramList) {
        repeat(p.second) {
            options.add(p.first)
        }
    }
    var tempOptions = options.toMutableList()
    var tempHistogramMap = histogramMap.toMutableMap()

    val o: MutableList<Double> = mutableListOf()
    o.add(histogramList[0].first)
    /*tempOptions.removeAt(index = 0)
    tempHistogramMap[histogramList[0].first] = histogramMap[histogramList[0].first]!! - 1*/
    //var check = false

    while(o.size != size) {
        while(true) {
            val index = (0 until options.size).random()
            val result = options[index]

            if (abs(result - o.last()) <= STEP_SIZE) {
                o.add(result)
                /*if(check) {
                    tempHistogramMap[result] = tempHistogramMap[result]!! - 1
                    if (tempHistogramMap[result] == 0) {
                        tempOptions = tempOptions.subList(index + 1, tempOptions.size)
                    } else {
                        tempOptions.removeAt(index = index)
                    }
                    if (tempOptions.size == 0) {
                        //println("No more options")
                        tempOptions = options.toMutableList()
                        tempHistogramMap = histogramMap.toMutableMap()
                        check = false
                    }
                }*/
                break
            }
        }
    }
    return o
}

fun calculateAtomicVibrato(
    scale: Double,
    windowSizeInSamples: Int,
    base: Double
) : Pair<List<Double>, List<Double>> {
    val histogramsEven: List<Map<Double, Int>> = listOf(
        mapOf(5.0 to 32, 4.0 to 27, 3.0 to 29, 2.0 to 10, 1.0 to 3),
        mapOf(4.0 to 2, 3.0 to 23, 2.0 to 42, 1.0 to 23, 0.0 to 7, -1.0 to 3),
        mapOf(4.0 to 2, 3.0 to 18, 2.0 to 43, 1.0 to 28, 0.0 to 8, -1.0 to 1),
        mapOf(4.0 to 2, 3.0 to 16, 2.0 to 38, 1.0 to 27, 0.0 to 12, -1.0 to 4, -2.0 to 1)
    )

    val histogramsParitySum: List<Map<Double, Int>> = listOf(
        mapOf(10.0 to 10, 9.0 to 6, 8.0 to 19, 7.0 to 33, 6.0 to 26, 5.0 to 7),
        mapOf(6.0 to 8, 5.0 to 40, 4.0 to 48, 3.0 to 4 ),
        mapOf(6.0 to 4, 5.0 to 36, 4.0 to 48, 3.0 to 12),
        mapOf(6.0 to 3, 5.0 to 34, 4.0 to 41, 3.0 to 12, 2.0 to 3, 1.0 to 5, 0.0 to 1)
    )

    println("Calculating even signal pdf...")
    val evenSignal = calculateSignalFromPdf(
        initialValue        = 5.0,
        scale               = scale,
        windowSizeInSamples = windowSizeInSamples,
        base                = base,
        histograms          = histogramsEven
    ).toMutableList()

    var evenWindowSum = 0.0
    var p = 0
    while(evenWindowSum < windowSizeInSamples) {
        evenWindowSum += evenSignal[p] + base
        p++
        println("evenWindowSum: $evenWindowSum, windowSizeInSamples: $windowSizeInSamples")
    }
    println("Even signal number count: $p")

    println("\r\nCalculating parity Sum pdf...")
    val paritySumSignal = calculateSignalFromPdf(
        initialValue        = 10.0,
        scale               = scale,
        windowSizeInSamples = windowSizeInSamples,//(1.5 * windowSizeInSamples).toInt(),
        base                = base,
        histograms          = histogramsParitySum
    )
    /*var sumB = 0.0
    for(t in 0 until 100) {
        sumB += paritySumSignal[t]
    }
    println("SumB: $sumB")*/
    val oddSignal = paritySumSignal.zip(evenSignal) { a, b -> a - b}.toMutableList()

    /** TODO: Might make sense to move this logic into the signal generation code **/
    for (i in evenSignal.indices) {
        evenSignal[i] += base
    }
    for (i in oddSignal.indices) {
        oddSignal[i] += base
    }

    return Pair(first = evenSignal, second = oddSignal)
}

fun calculateSignalFromPdf(
    initialValue: Double,
    scale: Double,
    windowSizeInSamples: Int,
    base: Double,
    histograms: List<Map<Double, Int>>
) : List<Double> {
    val maxDifference          = scale* 1.0 //2.0
    var prev                   = scale*initialValue//10.0//15.0
    val o: MutableList<Double> = mutableListOf()
    for (j in histograms.indices) {
        val histogram = histograms[j]
        val options: MutableList<Double> = mutableListOf()
        for(entry in histogram) {
            repeat(entry.value) {
                options.add(scale*entry.key)
            }
        }
        options.sort() //TODO: Is this really necessary?

        var windowSum = 0.0
        var p = 0
        //while(if(usePeriods){ p < periods } else { windowSum < windowSizeInSamples}) {
        //while(windowSum < windowSizeInSamples) {
        for (i in 0 until windowSizeInSamples) {
            var trials = 0
            if(windowSum > 0) {
            //if(i != 0) {
                var index = (0 until options.size).random()
                while (abs(options[index] - prev) > maxDifference) {
                    index = (0 until options.size).random()
                    trials++
                    if (trials > options.size) {
                        println("trial#:$trials, prev: $prev option: ${options[index]}, diff: ${abs(options[index] - prev)}")
                    }
                }

                //TODO: Cleanup this redundancy so the if statement can be fed to index as a single line
                prev = options[index]
                o.add(prev)
                windowSum += prev + base

            } else  { /** [Important to prevent audible waves]A new window is starting so chose the option closest to the last point in the last window!! **/
                var smallestDiffIndex = 0
                for (oi in options.indices) {
                    if(abs(options[smallestDiffIndex] - prev) > abs(options[oi] - prev)) {
                        smallestDiffIndex = oi
                    }
                }
                prev = options[smallestDiffIndex]
                o.add(prev)
                windowSum += prev + base
            }
            p++
        }
        //println("WindowSum: $windowSum, P: $p")
    }
    return o
}
fun getHistogramStats(start: Int, length: Int, signal: List<Double>) : Map<Double, Int> {
    val histogram: MutableMap<Double, Int> = mutableMapOf()
    for (i in 0 until length) {
        if(i + start >= signal.size) {
            println("Exiting early")
            break
        }
        val value = signal[i + start]
        histogram[value] = histogram.getOrDefault(value, 0) + 1
    }
    val total:Double = histogram.values.sum().toDouble()
    for(k in histogram.keys.sorted()) {
        val value:Int = histogram[k]!!
        println("value: $k, Count: $value,  ${100.0*(value/total)}%")
    }
    //println("Total: ${total}")
    return histogram
}
fun plotSignal(scale: Double = 100.0, signal: List<Double>) {
    val x = normalize(scale = scale, input = signal)
    for(i in 0 until x.size) {
         plot(index = i, value = x[i])
    }
}
fun plot(index: Int, value: Double) {
    print("$index ")
    for(i in 0 until abs(value.toInt())) {
        if(value.toInt() < 0) {
            print("-")
        }
        if(value.toInt() > 0) {
            print("+")
        }
        if(value == 0.0) {
            print("0")
        }
    }
    print(" $value")
    println()
}

fun plotSignals(signalA: List<Double>, signalB: List<Double>) {
    for(i in signalA.indices) {
        if(i == signalB.size) { break }
        plots(index = i, valueA = signalA[i], valueB = signalB[i])
    }
}
fun plots(index: Int, valueA: Double, valueB: Double) {
    print("$index ")
    plotInline(value = valueA)
    print("  ")
    plotInline(value = valueB)
    println()
}
fun plotInline(value: Double) {
    for(i in 0 until abs(value.toInt())) {
        if(value.toInt() < 0) {
            print("-")
        }
        if(value.toInt() > 0) {
            print("+")
        }
        if(value == 0.0) {
            print(" ")
        }
    }
}

fun writeSamplesToFile(fileName: String, input: List<Double>) {
    val samples = input.map {it.toInt().toShort()}
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
fun fourBytesToInt(bytes: ByteArray, pos: Int): Int {
    return (bytes[pos].toInt() and 0xFF) or
            ((bytes[pos + 1].toInt() and 0xFF) shl 8) or
            ((bytes[pos + 2].toInt() and 0xFF) shl 16) or
            ((bytes[pos + 3].toInt() and 0xFF) shl 24)
}

fun bytesToSample(bytes: ByteArray, pos: Int) : Short {
    return ((bytes[pos].toInt() and 0xFF) or ((bytes[pos + 1].toInt() and 0xFF) shl 8)).toShort()
}
