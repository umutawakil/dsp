package org.dsp

import org.dsp.analysis.DiscreteFourierTransform
import org.dsp.analysis.DiscreteFourierTransform.Companion.SpectralBin
import org.dsp.analysis.DiscreteFourierTransform.Companion.dft
import org.dsp.analysis.DiscreteFourierTransform.Companion.inverseSpecialDftPolar
import org.dsp.analysis.DiscreteFourierTransform.Companion.inverseSpecialDftRectangular
import org.dsp.analysis.DiscreteFourierTransform.Companion.specializedDft
import org.dsp.analysis.DiscreteFourierTransform.Companion.synthesizeDft
import org.dsp.analysis.WaveformAnalyzer.Companion.findPeak
import org.dsp.analysis.WaveformAnalyzer.Companion.findPeaks
import org.dsp.analysis.WaveformAnalyzer.Companion.getWaves
import org.dsp.analysis.WaveformAnalyzer.Companion.getWavesP
import org.dsp.analysis.WaveformAnalyzer.Companion.listData
import org.dsp.analysis.WaveformAnalyzer.Companion.plotData
import org.dsp.config.Constants
import org.dsp.filters.FilterUtils
import org.dsp.modulation.WaveformEffect.Companion.normalize
import org.dsp.modulation.WaveformEffect.Companion.normalizeWaves
import org.dsp.signals.SignalGenerator.Companion.cosineByFrequency
import org.dsp.signals.SignalGenerator.Companion.halfSine
import org.dsp.signals.SignalGenerator.Companion.sineByCycles
import org.dsp.signals.SignalGenerator.Companion.sineByFrequency
import org.dsp.signals.SignalGenerator.Companion.sineByWaveLength
import org.dsp.signals.SignalGenerator.Companion.tailoredHalfSine
import org.dsp.tools.add

import org.dsp.wave.WaveUtils
import kotlin.math.*
import kotlin.random.Random
fun main() {
    val resourceDirectory  = "/Users/umutawakil/Documents/Git/dsp/src/main/resources"
    val waveData = WaveUtils.getFileData(
        baseHalfPeriod    = 27.0, //TODO: This is kludgy as hell
        resourceDirectory = resourceDirectory,
        fileName          = "normalized.wav"//"harmonic-2.wav"//"test-normalized.wav"//"harmonic-2.wav"//"harmonic-2.wav"//"normalized.wav"//"harmonic-2.wav"//"normalized.wav"//"filter-out-1.wav"
    )

    val baseFrequency         = 444.0
    val baseWavelength        = (Constants.SAMPLE_RATE/baseFrequency)
    val targetHarmonic        = 4
    val targetFrequency       = targetHarmonic * baseFrequency
    val targetWaveLength      = baseWavelength / targetHarmonic
    val length                = Constants.SAMPLE_RATE//waveData.rawFileData.size//Constants.SAMPLE_RATE + 1
    val impulseResponseLength = length//* 2
    val sBandwidth            = baseFrequency//44.0//baseFrequency//44.0//baseFrequency/4 //(baseFrequency/4) May be all we need or perhaps even (baseFrequency/k) or some function of this.
    val frequencyDistance     = 1.0

    /*val fOut1 = FilterUtils.bandpassFIR(
        input                 = waveData.rawFileData,//.map { it * 1000.0},
        impulseResponseLength = impulseResponseLength,
        minFreqHz             = targetFrequency - (sBandwidth/2),
        maxFreqHz             = targetFrequency + (sBandwidth/2)
    )

    WaveUtils.writeSamplesToFile(
        fileName = "$resourceDirectory/output.wav",
        input    =  normalize(scale = 5000.0, input = fOut1)
    )
    return*/


    val synthAmp = createFormantFrequencyAmplitudeEnvelope(
        amplitude   = 1.0,
        coefficient = -1.0/ 64, //(16*decayFactor),
        length      = sBandwidth.toInt()
    ).toMutableList()

    /*val postFilter = FilterUtils.bandpassFIR(
        input                 = waveData.rawFileData,//.map { it * 1000.0},
        impulseResponseLength = impulseResponseLength + 1,
        minFreqHz             = targetFrequency - (sBandwidth/2),
        maxFreqHz             = targetFrequency + (sBandwidth/2)
    )*/

    //listData(data = synthAmp)
    //return

    val dft = specializedDft(
        x                 = normalize(scale = 5000.0, input = waveData.rawFileData),
        baseFrequency     = targetFrequency,
        sampleRate        = Constants.SAMPLE_RATE,
        bandwidth         = sBandwidth.toInt(),
        frequencyDistance = frequencyDistance
    )

    val nR    = dft.real.filter { it < 0}.sum()
    val pR    = dft.real.filter { it > 0}.sum()
    val peakR = findPeak(input = dft.real)

    val nI    = dft.imaginary.filter{it < 0}.sum()
    val pI    = dft.imaginary.filter{it > 0}.sum()
    val peakI = findPeak(input = dft.imaginary)

    println("nR: $nR, pR: $pR, ratio: ${(nR/pR)*100.0}, PeakR: $peakR, PeakR: ${(peakR/((abs(nR) + abs(pR))) * 100.0)}")
    println("nI: $nI, pR: $pI, ratio: ${(nI/pI)*100.0}, PeakI: $peakI, PeakR: ${(peakI/((abs(nI) + abs(pI))) * 100.0)}")
    //return

    val a = sineByFrequency(amplitude = -15.0,  frequency = 400.0, size = length)
    val b = sineByFrequency(amplitude = 5.0, frequency = 401.0, size = length)
    val c = sineByFrequency(amplitude = 5.0, frequency = 402.0, size = length)
    val d = sineByFrequency(amplitude = 5.0, frequency = 403.0, size = length)
    val e = sineByFrequency(amplitude = 5.0, frequency = 404.0, size = length)
    val f = sineByFrequency(amplitude = 5.0,  frequency = 405.0, size = length)

    val y = a.add(b).add(c).add(d).add(e).add(f)

    val m  = sineByFrequency(amplitude = 1.0, frequency = 0.5, size = length).map { (it) + 1}
    val resa = a.zip(m) { i,j -> i * j}
    val resb = b.zip(m) { i,j -> i * j}
    val resc = c.zip(m) { i,j -> i * j}
    val resd = d.zip(m) { i,j -> i * j}
    val rese = e.zip(m) { i,j -> i * j}

    var nOrgMag = normalize(scale = 5000.0, input = dft.magnitude).toMutableList()
    for(i in 0 until 44) {
        //nOrgMag[(nOrgMag.size/2) - 22 + i] = 0.0
    }
    //nOrgMag[17] = 0.0
    val nSynMag = normalize(scale = 5000.0, input = synthAmp)

    for(i in nOrgMag.indices) {
       // println("($i) O: ${nOrgMag[i]}, S: ${nSynMag[i]}")
    }
    //return

    val orgPeak         = nOrgMag.max()
    val orgMinusPeak    = nSynMag.sum() - orgPeak
    val peakToRestRatio = 100.0 * (orgPeak / orgMinusPeak)

    println("OrgPeak: $orgPeak, OrgMinusPeak: $orgMinusPeak, ratio: $peakToRestRatio")

    val synthPeak            = normalize(scale = 5000.0, input = synthAmp).max()
    val synthMinusPeak       = normalize(scale = 5000.0, input = synthAmp).sum() - synthPeak
    val synthPeakToRestRatio = 100.0 * (synthPeak / synthMinusPeak)

    println("synthPeak: $synthPeak, synthMinusPeak: $synthMinusPeak, synthPeakToRestRatio: $synthPeakToRestRatio")

    val idft = inverseSpecialDftPolar(
        magnitude         = dft.magnitude,
        phaseInPercent    = dft.phaseInPercent,
        baseFrequency     = targetFrequency,
        bandwidth         = sBandwidth.toInt(),
        frequencyDistance = frequencyDistance,
        sampleRate        = Constants.SAMPLE_RATE,
        length            = length
    )

    val wf = getWaves(data = idft)
    println("Waves found: ${wf.size}")
    val nw = wf.map { normalize(scale = 5000.0, input = it)}.flatten()

    /*val newDft = specializedDft(
        x                 = normalize(scale = 5000.0, input = nw),
        baseFrequency     = targetFrequency,
        sampleRate        = Constants.SAMPLE_RATE,
        bandwidth         = sBandwidth.toInt(),
        frequencyDistance = frequencyDistance
    ) */

    val wls = wf.map { it.size }
    val peaks = findPeaks(input = wf.flatten())
    var direction = 1
    val output: MutableList<Double> = mutableListOf()
    for(i in wls.indices) {
        val w = wls[i]
        val amplitude = direction.toDouble() // peaks[i]
        output.addAll(tailoredHalfSine(amplitude = amplitude , size = w))
        //output.addAll(tailoredHalfSine(amplitude = direction.toDouble(), size = (targetWaveLength/2).toInt()))
        direction *= -1
    }

    val postFilter = normalize(scale = 5000.0, input = FilterUtils.bandpassFIR(
        input                 = output,//.map { it * 1000.0},
        impulseResponseLength = impulseResponseLength + 1,
        minFreqHz             = targetFrequency - (sBandwidth/2),
        maxFreqHz             = targetFrequency + (sBandwidth/2)
    )
    )

    val baseSignal = normalize(scale = 2500.0, input = synthesizeDft(
        magnitude         = synthAmp,
        baseFrequency     = targetFrequency,
        bandwidth         = sBandwidth.toInt(),
        frequencyDistance = frequencyDistance,
        sampleRate        = Constants.SAMPLE_RATE,
        length            = length
    )
    )

    WaveUtils.writeSamplesToFile(
        fileName = "$resourceDirectory/output-istvan.wav",
        input    =  normalize(scale = 5000.0, input = postFilter.add(baseSignal))//.add(b).add(resc))
    )
    return
    /*val u: MutableList<Double> = mutableListOf()
    u.addAll(tailoredHalfSine(amplitude = 1.0, size = 13))
    u.addAll(tailoredHalfSine(amplitude = -1.0, size = 14))
    u.addAll(tailoredHalfSine(amplitude = 1.0, size = 13))

    val g = getWaves(data = u)
    listData(data = g.map { it.size})
    println("")
    listData(data = u)
    return */

    val waves2 = getWaves(data = output)
    val ws2 = waves2.map { it.size}
    for(i in ws2.indices) {
        //println("w2: ${ws2[i]}")
        println("($i) w1: ${wls[i]}, w2: ${ws2[i]}")
        if(wls[i] != ws2[i]) {
            throw RuntimeException("It's not working")
        }
    }
    return



    //val baseSound = sineByFrequency(amplitude = 5.0,  frequency = targetFrequency, size = length)

    val sdft = normalizeWaves(input = synthesizeDft(
        magnitude         = synthAmp,//tr,
        baseFrequency     = targetFrequency,
        bandwidth         = sBandwidth.toInt(),
        frequencyDistance = frequencyDistance,
        sampleRate        = Constants.SAMPLE_RATE,
        length            = length
    ), scale = 5000.0)

    WaveUtils.writeSamplesToFile(
        fileName = "$resourceDirectory/output.wav",
        input    =  normalize(scale = 5000.0, input = sdft)//.add(b).add(resc))
    )
    return
    /*var nSynthAmp = synthAmp.map { 1.0 }.toMutableList()
    nSynthAmp = (0 until nSynthAmp.size).map {
        if(it % 2 == 0) {
            nSynthAmp[it]*-1
        } else {
            nSynthAmp[it]
        }
    }.toMutableList()

    for(i in 0 until 44) {
        //nOrgMag[(nOrgMag.size/2) - 22 + i] = 0.0
        nSynthAmp[(nSynthAmp.size/2) - 22 + i] = 0.0
    }

    var innerFormantMagnitude = 0.0
    for(i in 0 until 44) {
        //nOrgMag[(nOrgMag.size/2) - 22 + i] = 0.0
        innerFormantMagnitude += dft.magnitude[(dft.magnitude.size/2) - 22 + i]
    }
    val outerFormantMagnitude = dft.magnitude.sum() - innerFormantMagnitude
    val innerOuterRatio = (innerFormantMagnitude / outerFormantMagnitude)* 100.0
    println("OuterMag: $outerFormantMagnitude, inner: $innerFormantMagnitude, ratio: $innerOuterRatio %")

    var formant = (0..sBandwidth.toInt() - 1).map { 282.1694676478179 / (sBandwidth.toInt() - 44)}.toMutableList()
    for(i in 0 until 44) {
        //nOrgMag[(nOrgMag.size/2) - 22 + i] = 0.0
        formant[(formant.size/2) - 22 + i] = (1.5*113.98468856670495) / 44.0
    }

    val scale = 20
    val mid = (formant.size/2)
    //formant[mid - 2] = -3*scale * formant[mid - 2]
    formant[mid - 1] = -(scale/3) * formant[mid - 1]
    formant[mid]     = 3*(scale/3) * formant[mid]
    formant[mid + 1] = (scale/3) * formant[mid + 1]
    /*formant[mid + 2] = scale * formant[mid + 2] */

    var innerFormantMagnitudeX = 0.0
    for(i in 0 until 44) {
        //nOrgMag[(nOrgMag.size/2) - 22 + i] = 0.0
        innerFormantMagnitudeX += formant[(formant.size/2) - 22 + i]
    }
    val outerFormantMagnitudeX = formant.sum() - innerFormantMagnitudeX
    val innerOuterRatioX       = (innerFormantMagnitudeX / outerFormantMagnitudeX)* 100.0
    println("OuterMagX: $outerFormantMagnitudeX, inner: $innerFormantMagnitudeX, ratio: $innerOuterRatioX %")

    val sdft = synthesizeDft(
        magnitude         = dft.magnitude,//tr,
        baseFrequency     = targetFrequency,
        bandwidth         = sBandwidth.toInt(),
        frequencyDistance = frequencyDistance,
        sampleRate        = Constants.SAMPLE_RATE,
        length            = length
    )

    val averageMagnitude = findPeaks(input = sdft).map { abs(it) }.average()
    val waves = getWaves(data = sdft)
    val newWaves: MutableList<List<Double>> = mutableListOf()
    for(w in waves) {
        val peak = findPeak(input = w)
        val addition = if(peak <0) { averageMagnitude * -1 } else { averageMagnitude }
        newWaves.add(normalize(scale = peak + addition, input = w))
    }
    val out  = getWaves(data = sdft).flatten()
    */

    //val fm = fmSynthesis(fc = 444.0, fm = 0.5, modulationIndex = 20.0, length = length)

    /**
     * Magnitude:
     *  SumR has a net magnitude of 72.128 and SumI 75.1.88.
     *  So a net difference of 2.077% of the total magnitude between them lives in  the imaginary
     *
     *  Amplitude:
     *  sumRP: 33795.63458707377, sumRN: -38333.159536218656, rdiff: 6.290864840175652
     * sumIP: 42086.04605245062, sumIN: -33102.68982113939, idiff: 11.947742074576656
     *
     * On closer inspection above it appears that the negative portion of the real overtakes is positive by 6.29 percent
     * The positive portion of the imaginary overtakes the negative by 11% and result is the amplitude remaining is twice
     * that in magnitude of whats left over in the real but remember that little difference sobrando between them is only 2% of
     * the grand total.
     *
     * So in summary it seems they are both close to being half negative and half positive with their total magnitudes
     * being nearly identical in value. I find it hard to believe a 2% difference can be heard in amplitude given the logorathimic
     * nature of human hearing and potential 18% difference in their "DC sum" or 6% more posoitive than negative in one and 12
     * percent more in another I am also doubtful contributes to the DC balance I'm looking for.
     *
     * The key takeway is that despite the wild multiplicative noise between the two virtually equal magnitude
     * is preserved between them. Perhaps this can be achieved with zero average multiplicative noise on both sides
     * then checking that the positive and negative distributions are the same.
     *
     * The imaginary and real parts are two extremely chaotic signals fluctuating very wildly but what might be the silver lining
     * is that each of the crazy looking signals is basically zero balanced and although also being wildly different
     * from each other they preserve the same net magnitude which is impressive. This could indeed mean a linear curve
     * could be used as the base, and multiplicative noise added to create a real or imaginary part then multiplicative noise
     * added again to get the sibling part.
     *  */

    /*val sumRP = dft.real.map { if(it> 0) {it*1000.0} else { 0.0} }.sum()
    val sumRN = dft.real.map { if(it< 0) {it*1000.0} else { 0.0} }.sum()
    val rdiff = (abs(sumRP + sumRN)/(abs(sumRP) + abs(sumRN))) * 100.0

    val sumIP = dft.imaginary.map { if(it> 0) {it*1000.0} else { 0.0} }.sum()
    val sumIN = dft.imaginary.map { if(it< 0) {it*1000.0} else { 0.0} }.sum()
    val idiff = (abs(sumIP + sumIN)/(abs(sumIP) + abs(sumIN))) * 100.0

    println("sumRP: $sumRP, sumRN: $sumRN, rdiff: $rdiff")
    println("sumIP: $sumIP, sumIN: $sumIN, idiff: $idiff") // SumR has a net magnitude of 72.128 and SumI 75.1.88
    return*/

    //val cdiff = complexDifference(real = dft.magnitude, imaginary = dft.imaginary)
    //val average = cdiff.average()
    //val std     = standardDeviation(data = dft.real.map { 1000.0 * it})
    //val average = dft.real.map { 1000.0 * it}.average()

    //REAL Average: -103.12556702601995, std: 1950.8326615119718, ssd: 4.5381557015333724E8
    //IMAG: Average: 204.16718707525504, std: 2241.4242517265507, ssd: 5.560425835471629E8

    /*val ml = generateRandomNumbers(
        average                 = average,
        standardDeviation       = std,
        potentialMin            = -1000.0,
        potentialMax            = 1000.0,
        length                  = cdiff.size
    )*/

}

//c(t) = Asin(fc*i + Isin(fm*i))

fun fmSynthesis(fc: Double, modulationIndex: Double, fm: Double, length: Int) : List<Double>{
    val o: MutableList<Double> = mutableListOf()
    for(i in 0 until length) {
        val modulator = modulationIndex * sin((2*Math.PI*i) / (Constants.SAMPLE_RATE / fm) )
        val value     = sin(((2*Math.PI*i) / (Constants.SAMPLE_RATE / fc)) + modulator )
        o.add(value)
    }
    return o
}
class InitialRelationship(val frequency: Double,val amplitude: Double,val position: Int)
class Relationship(val frequency: Double,val amplitude: Double, val proportion: Double, val position: Int)

/**
 * TODO: This function greatly assumes a frequency distance of 1 hertz.
 */
fun beatEnvelopTransform(amplitude: List<Double>) : List<List<Relationship>> {
    val initialRelationships: MutableList<MutableList<InitialRelationship>> = mutableListOf()
    for(i in amplitude.indices) {
        //TODO: Loop to the left collecting relationships (probably while loops)
        initialRelationships.add(mutableListOf())
        var distance = 1
        while(i - distance >= 0) {
            initialRelationships[i].add(
                InitialRelationship(
                    frequency = distance.toDouble(),
                    amplitude = amplitude[i - distance],
                    position  = i - distance
                )
            )
            distance++
        }

        //TODO: Loop to the right collecting relationships
        distance = 1
        while(i + distance < amplitude.size) {
            initialRelationships[i].add(
                InitialRelationship(
                    frequency = distance.toDouble(),
                    amplitude = amplitude[i + distance],
                    position  = i + distance
                )
            )
            distance++
        }
    }

    //TODO: Build list of relationships by looping over the list of initial relationships

    val relationships: MutableList<List<Relationship>> = mutableListOf()
    val totals: MutableList<Double> = mutableListOf()
    //TODO: Are we sure we want absolute sum?
    for(ir in initialRelationships) {
        val totalForNode = ir.map { abs(it.amplitude)}.sum() //TODO: Are we sure we want absolute sum?
        totals.add(totalForNode)
        relationships.add(
            ir.map {
                Relationship(
                    frequency  = it.frequency,
                    amplitude  = it.amplitude,
                    position   = it.position,
                    proportion = abs(it.amplitude)/totalForNode //TODO: Do we need to preserve sign?
                )
            }
        )
    }


    //TODO: Loop through the list of relationships to determine whats going into the buffer and whats going to envelopeDFT
    //val envelopeDft: List<Double> // each index represents (hz - 1) = 1

    var totalPossibleContributions = 0.0
    val relationshipStats: MutableMap<Int,Double> = mutableMapOf()
    val usedContributions: MutableSet<String> = mutableSetOf()
    for(i in relationships.indices) {
        val currentNodesRelationships = relationships[i]
        if(abs(totals[i]) < abs(amplitude[i])) {
            println("Less than!!!")
        }

        var contributionTotal = if(abs(totals[i]) < abs(amplitude[i])) {
            abs(totals[i])
        } else {
            abs(amplitude[i])
        }
        //println("ContributionTotal for ${amplitude[i]} is $contributionTotal")

        for(currentRelationship in currentNodesRelationships) {
            val minPos          = min(i, currentRelationship.position)
            val maxPos          = max(i, currentRelationship.position)
            val relationshipKey = "$minPos-$maxPos"
            if(usedContributions.contains(relationshipKey)) { continue} else {
                usedContributions.add(relationshipKey)
            }

            //TODO: we need to determine the polarity of the contribution and also apply proportion
            val polarity = if(currentRelationship.amplitude * amplitude[i] < 0.0) { -1 } else { 1 }
            val contribution = contributionTotal * currentRelationship.proportion * polarity
            //val contribution = min(abs(amplitude[minPos]), abs(amplitude[maxPos])) * polarity//contributionTotal * currentRelationship.proportion * polarity
            //val contribution = min(abs(amplitude[minPos]), abs(amplitude[maxPos])) * currentRelationship.proportion * polarity
            println("ContributionTotal for ${amplitude[i]} at frequency: ${currentRelationship.frequency}, pos: ${currentRelationship.position}, is $contribution")
            relationshipStats[currentRelationship.frequency.toInt()] = relationshipStats.getOrPut(currentRelationship.frequency.toInt()) { 0.0 } + contribution
            totalPossibleContributions += abs(contribution)
        }

        //TODO: Each relationship for each node is turned into a contribution for each "frequency" based on proportion
        // and taking polarity into consideration.

    }
    for(pair in relationshipStats) {
        println("Frequency: ${pair.key}, total: ${pair.value}")
    }
    //TODO: The buffer is the sum all the final relationship sums (the envelope DFT summed) subtracted from the
    // amplitude summed (absolute sum) perhaps

    val estimatedPeak     = amplitude.map { abs(it)}.sum()
    val fluctuationRegion = relationshipStats.values.map { abs(it) }.sum()
    val scaledSum         = fluctuationRegion * sqrt(2.0)
    val estimatedBuffer   = estimatedPeak - fluctuationRegion

    val fluxPercent = 100.0 - (((totalPossibleContributions - fluctuationRegion)/totalPossibleContributions) * 100)
    println("Flux percent: $fluxPercent")
    println("TotalPossibleContributions: $totalPossibleContributions")
    println("Sum: ${fluctuationRegion}, scaledSum: $scaledSum,")
    println("Scaled estimate: ${(scaledSum/estimatedPeak) * 100} % fluctuation")
    println("Buffer: $estimatedBuffer out of $estimatedPeak, or fluctuations in ${(fluctuationRegion/estimatedPeak)* 100} % of signal")

    return relationships
}

fun generateRandomNumbers(
    average: Double,
    standardDeviation: Double,
    potentialMin: Double,
    potentialMax: Double,
    length: Int
): List<Double> {
    val result = mutableListOf<Double>()
    val random = java.util.Random()

    // Ensure that potentialMin and potentialMax allow for a reasonable range
    if (potentialMin > potentialMax) {
        throw IllegalArgumentException("potentialMin must be less than or equal to potentialMax")
    }

    // Generate values using Gaussian distribution
    while (result.size < length) {
        // Generate a candidate value using Gaussian distribution
        val candidate = random.nextGaussian() * standardDeviation + average

        // Check if the candidate is within the potential range
        if (candidate in potentialMin..potentialMax) {
            result.add(candidate)
        }

        // Optional: Break if too many attempts are made to avoid infinite loops
        if (result.size >= length || result.size > 10 * length) {
            break
        }
    }

    // Final filtering to ensure all values are within bounds
    return result.filter { it in potentialMin..potentialMax }
}

/*fun generateMatchingList(
    average: Double,
    sumOfSquaredDifferences: Double,
    standardDeviation: Double,
    potentialMin: Double,
    potentialMax: Double,
    length: Int
): List<Double> {
    val result = mutableListOf<Double>()
    val random = java.util.Random()

    // Calculate variance
    val variance = sumOfSquaredDifferences / (length - 1)

    while (result.size < length) {
        // Generate a candidate value using Gaussian distribution
        val candidate = random.nextGaussian() * standardDeviation + average

        // Check if the candidate is within the potential range
        if (candidate in potentialMin..potentialMax) {
            result.add(candidate)
        }
    }

    // Adjust the list to match the exact average and sum of squared differences
    val currentAverage = result.average()
    val currentSumOfSquaredDifferences = result.sumOf { (it - currentAverage).pow(2) }

    val scaleFactor = sqrt(sumOfSquaredDifferences / currentSumOfSquaredDifferences)
    val adjustedList = result.map { (it - currentAverage) * scaleFactor + average }

    return adjustedList
}*/

fun changeRate(input: List<Double>) : List<Double> {
    return (1..input.size-1).map { ((input[it] - input[it - 1]) / input[it-1]) * 100.0}

    /*val o: MutableList<Double> = mutableListOf()
    for(i in 1 until input.size) {
        o.add(((input[i] - input[i - 1]) / input[i]) * 100.0)
    }
    return o*/
}
fun complexDifference(real: List<Double>, imaginary: List<Double>) : List<Double> {
    return real.zip(imaginary) {r, i ->
        (i/r) * 100.0
    }
}
fun simplePhaseGenerator(step: Double, max: Double, length: Int) : List<Double> {
    val out: MutableList<Double> = mutableListOf()
    var runningTotal             = step
    var currentDirection         = 1

    for(i in 0 until length) {
        out.add(runningTotal)
        if(runningTotal > 0) {
            runningTotal += step * currentDirection
        } else {
            runningTotal += step * 3 * currentDirection
        }

        if(abs(runningTotal) >= max) {
            currentDirection *= -1
        }
    }
    return out
}

/*fun generateRandomNumbers(length: Int, standardDeviation: Double, average: Double): List<Double> {
    val random = Random.Default
    val result = mutableListOf<Double>()

    for (i in 0 until length) {
        val u1 = random.nextDouble()
        val u2 = random.nextDouble()

        val z = sqrt(-2.0 * ln(u1)) * cos(2.0 * PI * u2)
        val value = z * standardDeviation + average

        result.add(value)
    }

    return abs(result)
}*/

fun generatePositiveRandomNumbers(length: Int, standardDeviation: Double, average: Double): List<Double> {
    val random = Random.Default
    val result = mutableListOf<Double>()

    // Calculate the parameters for the truncated normal distribution
    val a = -average / standardDeviation  // Lower bound (in standard normal terms)

    while (result.size < length) {
        val u = random.nextDouble()
        val v = random.nextDouble()

        // Box-Muller transform
        var z = sqrt(-2.0 * ln(u)) * cos(2.0 * PI * v)

        // Truncate the distribution
        if (z > a) {
            // Transform z back to our desired distribution
            val value = z * standardDeviation + average
            result.add(value)
        }
    }

    return result
}
/*fun generatePositiveRandomNumbers(length: Int, standardDeviation: Double, average: Double): List<Double> {
    val random = Random.Default
    val result = mutableListOf<Double>()

    while (result.size < length) {
        val u1 = random.nextDouble()
        val u2 = random.nextDouble()

        val z = sqrt(-2.0 * ln(u1)) * cos(2.0 * PI * u2)
        val value = z * standardDeviation + average

        if (value > 0) {
            result.add(value)
        }
    }

    // Adjust the values to match the desired average and standard deviation
    val actualAverage = result.average()
    val actualStdDev = calculateStandardDeviation(result)

    return result.map { value ->
        ((value - actualAverage) / actualStdDev * standardDeviation + average).coerceAtLeast(0.0)
    }
}

// Helper function to calculate standard deviation
fun calculateStandardDeviation(numbers: List<Double>): Double {
    val mean = numbers.average()
    val variance = numbers.map { (it - mean).pow(2) }.average()
    return sqrt(variance)
}*/

fun symmetricDifference(input: List<Double>) : List<Double> {
    val midPoint = input.size / 2
    val out: MutableList<Double> = mutableListOf()
    for (i in 0 until midPoint) {
        out.add(abs(input[i] - input[input.size - 1 - i]))
    }
    return out
}

fun symmetricDifferenceChange(input: List<Double>) : List<Double> {
    val midPoint = input.size / 2
    val out: MutableList<Double> = mutableListOf()
    for (i in 0 until midPoint) {
        val value = 100 * (abs(input[i] - input[input.size - 1 - i])/input[input.size - 1 - i])
        out.add(value)
    }
    return out
}

fun symmetricSum(input: List<Double>) : List<Double> {
    val midPoint = input.size / 2
    val out: MutableList<Double> = mutableListOf()
    for (i in 0 until midPoint) {
        //val value = sqrt(input[i].pow(2.0) + input[input.size - 1 - i].pow(2.0))
        out.add(input[i] + input[input.size - 1 - i])
        //out.add(value)
    }
    return out
}

//TODO: How does this statistically differ form the SSD? ( I believe for whatever reason SSD places serious emphasis on adjacent jumps)
fun standardDeviation(data: List<Double>): Double {
    // Ensure the list has at least one element to calculate standard deviation
    if (data.isEmpty()) return 0.0

    // Calculate the mean
    val mean = data.sum() / data.size

    // Calculate the sum of squared differences from the mean
    val sumOfSquaredDifferences = data.sumOf { (it - mean).let { diff -> diff * diff } }

    // Calculate and return the standard deviation
    return Math.sqrt(sumOfSquaredDifferences / data.size) // For population standard deviation
    // return Math.sqrt(sumOfSquaredDifferences / (numbers.size - 1)) // For sample standard deviation
}
fun sumOfSquaredDifferences(data: List<Double>): Double {
    // Ensure the list has at least two elements to calculate differences
    if (data.size < 2) return 0.0

    // Compute the sum of squared differences
    return data.zipWithNext { a, b -> (b - a).let { it * it } }.sum()
}

/**
 * These are the stats for the organic signal while the synthetic with bounded noise is still almost double in these roughness measures.
 * It's pretty clear that you need to limit the amount of amplitude pushed into any block of time and 5% of total period envelope is a
 * good base time slot. Large jumps do occur but mostly toward the ends. I'll also add that the new synthetic version is using a bounded noise
 * so the results are vastly improved compared to pure white noise as phase.
 *
 * (Organic)
 * std: 1.9899503959952227, ssd: 90.68670315742034
 *
 * (Synthetic with simple bounded noise for phase)
 * std: 2.2182426329562137, ssd: 148.8759035280646
 *
 * (Synthetic with pure noise for phase)
 * std: 3.0901474808513125, ssd: 234.32443094174167
 *
 * Note that this was all for just the central formants. The transitions will most likely be much more severe on full formants
 *
 * For full formants
 *
 * (Organic)
 * std: 1.6961009697693594, ssd: 32.6877454659895
 *
 * (Full Synthetic with white noise for phase)
 * std: 1.7274920155477012, ssd: 78.43651094007332
 *
 * (Full Synthetic with bounded noise for phase)
 *
 * std: 1.8621174234024198, ssd: 81.5434579659213
 *
 * With harmonics added the difference betweent he bounded noise and unbounded phase noise is pretty minimal
 *
 * -In the organic full analysis exculding the start and finish the transition from any adjacent block of 5% is never more than 1%
 * and in many cases its 0.5 to 0.7% with the last to pushing it into the 1.8 average. The synth is significantly more jumpy with the
 * additional harmonics added and the SSD metric is a much better representation of that as it seems to put more weight on the transitions
 * from point to point where in the case of the standard deviation they arrive at a similar average when they are very different in roughness.
 *
 * -It could be described as a gradual rise to a peak of 8 then it starts to drop off
 *
 */

fun amplitudePhaseDistributionTransform(magnitude: List<Double>, phaseInPercent: List<Double>, divisionSize: Int) {
    val divisions: List<Double>         = (-50..50 step divisionSize).toList().map { it.toDouble()}
    val tempBucket: MutableList<Double> = MutableList(divisions.size) { 0.0 }
    val scaledMagnitude: List<Double>   = normalize(scale = 5000.0, input = magnitude)
    val bag: MutableMap<Double, MutableList<Double>> = mutableMapOf()

    println("------ Scaled Magnitude BEGIN----")
    listData(data = scaledMagnitude)
    println("------ Scaled Magnitude END -----")

    for(k in phaseInPercent.indices) {
        for(d in divisions.indices) {
            val minBound = if(d == 0) { divisions[d] } else { divisions[d - 1] }
            if((minBound <= phaseInPercent[k]) && (phaseInPercent[k] <= divisions[d])) {
                /*if(dft.magnitude[k].toInt() == 1956 ) {
                    println("d($d), CurrentDivision: ${divisions[d]}, minBound: $minBound, phase: ${phase[k]}, division: ${divisions[d]}")
                }*/
                tempBucket[d] += magnitude[k]

                bag.getOrPut(divisions[d]) { mutableListOf() }.add(scaledMagnitude[k]) //TODO: should this be an option
            }
        }
    }

    /** Purely for debugging to ensure the right values are going into the right divisions **/
    for(k in bag.keys.sorted()) {
        println("Phase: $k, SUM: ${bag[k]!!.sum()}, Size: ${bag[k]!!.size}")
        for(p in bag[k]!!) {
            println("Division: $k, Val: $p")
        }
        println()
    }

    println("TEMP BUCKET")
    listData(data = tempBucket)
    println()

    val usableData = tempBucket.subList(1, tempBucket.size) //TODO: This algorithm is sloppy. Must be a better way that doesn't create a useless preceding element in the list.
    val sum        = usableData.sum()
    val dist       = usableData.map { (it/sum)*100 }

    /*** Display the results. Currently, there are no reasons to return the data **/
    listData(data = dist)

    println("SUM: ${dist.sum()}")

    //std --> synth: 3.4, org:  2.821
    //ssd --> synth: 856, org: 285
    val ssd = sumOfSquaredDifferences(data = dist)
    val std = standardDeviation(data = dist)
    println("std: $std, ssd: $ssd")
}

fun ocarinaDft(baseFrequency: Double) : List<FreqPair> {
    val numHarmonics = (Constants.SAMPLE_RATE / baseFrequency).toInt()
    val keyHarmonics: List<FreqPair> = listOf(
        FreqPair(amp = 0.0, freq = baseFrequency*0),
        FreqPair(amp = 168.0, freq = baseFrequency*2),
        FreqPair(amp = 136.0, freq = baseFrequency*3),
        FreqPair(amp = 40.0,  freq = baseFrequency*4),
        FreqPair(amp = 46.0,  freq = baseFrequency*5),
        FreqPair(amp = 53.0,  freq = baseFrequency*6),
        FreqPair(amp = 20.0,  freq = baseFrequency*7),
        FreqPair(amp = 40.0,  freq = baseFrequency*8),
        FreqPair(amp = 7.0,   freq = baseFrequency*9),
        FreqPair(amp = 17.0,  freq = baseFrequency*10)
    )

    val upperHarmonics: MutableList<FreqPair> = mutableListOf()
    var nextUpperHarmonic = keyHarmonics.size + 1
    val initialAmp = 1.0
    for(i in nextUpperHarmonic  until numHarmonics) {
        val currentAmplitude = initialAmp - (initialAmp/(numHarmonics)*i)
        upperHarmonics.add(FreqPair(amp = currentAmplitude, freq = i*baseFrequency))
    }

    return (keyHarmonics + upperHarmonics)//.subList(0, 20)
}
fun compressSignalWithWavelengths(
    input: List<Double>,
    waveLengths: List<Int>
) : List<Double> {
    val targetAmplitude = 5000.0
    val out             = mutableListOf<Double>()
    var direction       = 1
    var pos             = 0

    for (i in waveLengths.indices) {
        val compressedPeriod = compressHalfPeriodHelper(
            targetAmplitude = targetAmplitude*direction,
            input           = input,
            start           = pos,
            end             = pos + waveLengths[i]
        )
        out.addAll(compressedPeriod)
        direction *= -1
        pos += waveLengths[i]
    }
    return out
}

fun compressSignal(
    targetAmplitude: Double,
    baseWavelength: Double,
    input: List<Double>
) : List<Double> {
    val halfWavelength  = (baseWavelength / 2).toInt()

    val out           = mutableListOf<Double>()
    val numPeriods    = input.size / halfWavelength
    var direction     = 1

    for(p in 0 until numPeriods) {
        val compressedPeriod = compressHalfPeriodHelper(
            targetAmplitude = direction*targetAmplitude,
            input           = input,
            start           = halfWavelength*p,
            end             = (halfWavelength*p) + halfWavelength
        )
        out.addAll(compressedPeriod)
        direction *= -1
    }
    return out
}

fun compressHalfPeriodHelper(
    input: List<Double>,
    targetAmplitude: Double,
    start: Int,
    end: Int
) : List<Double> {
    val currentUnderWaveform = input.subList(start, end)
    val iteration1           = compressHalfPeriod(amplitude = targetAmplitude, stepSize = 100.0, input = currentUnderWaveform)
    val iteration2           = compressHalfPeriod(amplitude = targetAmplitude, stepSize = 10.0, input  = iteration1)

    return compressHalfPeriod(amplitude = targetAmplitude, stepSize = 1.0, input   = iteration2)
}
fun compressHalfPeriod(amplitude: Double, stepSize: Double, input: List<Double>) : List<Double> {
    var temp: List<Double>
    var currentAmplitude = 0.0//TODO: Start at zero for those moments its already in range  //direction * range
    var iteration        = 0
    var previousDistance = 0.0

    while(true) {
        //val inputMax = input.max()
        temp = input.zip(halfSine(amplitude = currentAmplitude, size = input.size)) {a,b -> a + b }
        val peak = findPeak(input = temp)
        val currentDistance = abs(peak - amplitude)

        //TODO: Need to check if the currentDistance is larger than the previous meaning an overshoot occurred

        if (currentDistance <= stepSize) {
            break
        }
        if((currentDistance > previousDistance) && (previousDistance > 0.0)) {
            println("Overshoot detected -> previous: $previousDistance, currentDistance: $currentDistance, target: $amplitude, current: $currentAmplitude, step: $stepSize")
            break
        }
        previousDistance = currentDistance

        var dir = 1
        if(amplitude > peak) {
            currentAmplitude += stepSize
        } else if(amplitude < peak ) {
            currentAmplitude -= stepSize
            dir = -1
        } else {
            throw RuntimeException("This should never execute given the check above")
        }

        //if(iteration > 50) {
            //println("i $iteration, t: $amplitude, iMax: $inputMax, d: $currentDistance, cA: $currentAmplitude, p: $peak, range: $stepSize, d: $dir")
        //}
        iteration++

        /** So far this doesn't seem to get hit anymore **/
        if(iteration > 1000) {
            println("iteration: $iteration, target: $amplitude, distance: $currentDistance, currentAmplitude: $currentAmplitude, peak: $peak, range: $stepSize, dir: $dir")
            //listData(data = input)
            //println()
            //listData(data = temp)
            throw RuntimeException("Computation limit reached")
        }
    }
    return temp
}
fun ocarinaSynth(baseFrequency: Double, length: Int) : List<Double> {
    val baseWavelength = Constants.SAMPLE_RATE / baseFrequency
    val numHarmonics    = ((baseWavelength / 2) + 1).toInt()
    //TODO: Need a DFT list of the first 10 harmonics (set fundamental to zero)
    val keyHarmonics: List<FreqPair> = listOf(
        FreqPair(amp = 0.0, freq = baseFrequency*0),
        FreqPair(amp = 168.0, freq = baseFrequency*2),
        FreqPair(amp = 136.0, freq = baseFrequency*3),
        FreqPair(amp = 40.0,  freq = baseFrequency*4),
        FreqPair(amp = 46.0,  freq = baseFrequency*5),
        FreqPair(amp = 53.0,  freq = baseFrequency*6),
        FreqPair(amp = 20.0,  freq = baseFrequency*7),
        FreqPair(amp = 40.0,  freq = baseFrequency*8),
        FreqPair(amp = 7.0,   freq = baseFrequency*9),
        FreqPair(amp = 17.0,  freq = baseFrequency*10)
    )

    val upperHarmonics: MutableList<FreqPair> = mutableListOf()
    var nextUpperHarmonic = keyHarmonics.size + 1
    for(i in nextUpperHarmonic  until numHarmonics) {
        val currentAmplitude = 6.0 - (5.0/(numHarmonics)*i)
        upperHarmonics.add(FreqPair(amp = currentAmplitude, freq = i*baseFrequency))
    }

    println("\r\nkey harmonics:")
    for(i in keyHarmonics.indices) {
        println("A: ${keyHarmonics[i].amp}, F: ${keyHarmonics[i].freq}")
    }
    println("\r\nupper harmonics:")
    for(i in upperHarmonics.indices) {
        println("$i, A: ${upperHarmonics[i].amp}, F: ${upperHarmonics[i].freq}")
    }

    //TODO: Loop through the DFT and generate corresponding formant and add it to the main buffer (called underWaveform)
    val frequencyBinData = (keyHarmonics + upperHarmonics).subList(0, 20)
    //var underWaveform: List<Double> = List(length) {0.0}
    var k = 0

    /**TODO: Is it worth the performance gains to do this procedurally? **/
    val underWaveform: List<Double> = frequencyBinData.fold((0..length).map {0.0}) { acc: List<Double>, currentBin ->
        acc.zip(
            createFormant(
                amplitude     = currentBin.amp,
                coefficient   = -(1.0/16),//-(1.0/32),
                bandWidth     = 444,
                baseFrequency = currentBin.freq,
                length        = length
            )
        ) { a: Double, b: Double ->  a + b }


       /* println("Formant: $k, of ${frequencyBinData.size}")
        k++
        val formant: List<Double> = createFormant(
            amplitude     = currentBin.amp,
            coefficient   = -(1.0/16),//-(1.0/32),
            formantSize   = 444,
            baseFrequency = currentBin.freq,
            length        = length
        )
        return acc.zip(formant) { a: Double, b: Double ->  a + b }*/
    }
    println("Under-waveform: ${underWaveform.size}")

    /*for(currentBin in frequencyBinData) {
        println("Formant: $k, of ${frequencyBinData.size}")
        val currentFormant = formantToSpectrum(
            amplitude     = currentBin.amp,
            coefficient   = -(1.0/16),//-(1.0/32),
            formantSize   = 444,
            baseFrequency = currentBin.freq,
            length        = length,
            iterations    = 1
        )
        underWaveform = underWaveform.zip(currentFormant) { a, b -> a + b }
        k++
    }*/

    val normalizedUnderWaveform = normalize(scale = 600.0, input = underWaveform)
    val uavg = normalizedUnderWaveform.map {abs(it)}.average()
    println("UAVG: $uavg")
    //val vib = Vibrato.genericVibratoSignal(baseHalfWavelength = baseWavelength/2.0, durationInSamples = length)
    //listData(data = vib)
    val fundamental = sineByWaveLength(amplitude = 5000.0, waveLength = baseWavelength, size = length)
   // val result =  fundamental.zip(normalizedUnderWaveform) { a, b -> a + b }

    val result = compressSignal(targetAmplitude = 5000.0, baseWavelength = baseWavelength, input = normalizedUnderWaveform)

    println("Result: ${result.size}")

    return result
}

fun createFormant(
    amplitude: Double,
    coefficient: Double,
    bandWidth: Int,
    baseFrequency: Double,
    length: Int,
    distanceBetweenPoints: Double = 1.0
): List<Double> {
    val spectrumAmplitudeEnvelope = createFormantFrequencyAmplitudeEnvelope(
        amplitude   = amplitude,
        coefficient = coefficient,
        length      = bandWidth
    )
    listData(data = spectrumAmplitudeEnvelope)

    val spectra: MutableList<FreqPair> = mutableListOf()
    for (i in spectrumAmplitudeEnvelope.indices) {
        spectra.add(
            FreqPair(
                amp = spectrumAmplitudeEnvelope[i],
                freq = baseFrequency - (distanceBetweenPoints*bandWidth) + (i*distanceBetweenPoints)
            )
        )
    }
    return synthesizeSpectra(spectra = spectra, length = length)
}

private class PhaseMagnitudeCount(val phase: Double, val currentMagnitude: Double)

//TODO: How can this be done in a functional and immutable way for easier testing?
fun synthesizeFormant(
    amplitude: Double,
    coefficient: Double,
    bandWidth: Double,
    baseFrequency: Double,
    length: Int,
    frequencyDistance: Double
) : List<Double> {
    val magnitudeCurve = createFormantFrequencyAmplitudeEnvelope(
        amplitude   = amplitude,
        coefficient = coefficient,
        length      = bandWidth.toInt()
    )

    val frequencies: List<Double> = generateFrequencies(
        baseFrequency     = baseFrequency,
        bandWidth         = bandWidth,
        frequencyDistance = frequencyDistance
    )
    val sortedMagAndFreq: List<FreqPair> =  magnitudeCurve.zip(frequencies) { a, b -> FreqPair(amp = a, freq = b) }.sortedBy{ it.amp }
    val phases                           = (-50..50 step 5).map{ it.toDouble()}.toList() //TODO: step needs to be 5

    val phaseMagCounts: MutableList<PhaseMagnitudeCount> = phases.map {
        PhaseMagnitudeCount(
            phase            = it,
            currentMagnitude = 0.0
        )
    }.toMutableList()
    val phaseBucket:MutableMap<Double,MutableList<FreqPair>> = phases.fold(initial = mutableMapOf()) {
        acc, curr -> acc.getOrPut(curr) { mutableListOf() }
        acc
    }

    var pos = 0
    while(pos < sortedMagAndFreq.size) {
        val currentPhaseToAugment = phaseMagCounts[0].phase
        phaseBucket[currentPhaseToAugment]!!.add(sortedMagAndFreq[pos])

        phaseMagCounts[0] = PhaseMagnitudeCount(
            phase            = phaseMagCounts[0].phase,
            currentMagnitude = phaseMagCounts[0].currentMagnitude + sortedMagAndFreq[pos].amp
        )
        pos++
        phaseMagCounts.sortBy { it.currentMagnitude }
    }

    val sum  = phaseMagCounts.map { it.currentMagnitude }.sum()
    val dist = phaseMagCounts.map{ (it.currentMagnitude/sum) * 100 }
    //listData(data = dist)

    val ssd = sumOfSquaredDifferences(data = dist)
    val std = standardDeviation(data = dist)
    //println("std: $std, ssd: $ssd")

    val spectralBins: MutableList<SpectralBin> = mutableListOf()
    for(k in phaseBucket.keys) {
        val freqPairs: List<FreqPair> = phaseBucket[k]!!
            spectralBins.addAll(freqPairs.map {
                SpectralBin(magnitude = it.amp, frequency = it.freq, phase = k)
            }
        )
    }

    /*val sortedMagByFreq   = spectralBins.sortedBy { it.frequency }.map { it.magnitude }
    val sortedPhaseByFreq = spectralBins.sortedBy { it.frequency }.map { it.phase }
    amplitudePhaseDistributionTransform(
        magnitude      = sortedMagByFreq,
        phaseInPercent = sortedPhaseByFreq,
        divisionSize   = 5
    )
    */

    return synthesizeSpectralBins(bins = spectralBins, length = length)
}

//TODO: Not sure if this function will ever have value outside of a test
fun balancePhase(
    magnitudeCurve: List<Double>,
    frequencies: List<Double>,
    length: Int
) : List<Double> {
    val sortedMagAndFreq: List<FreqPair> =  magnitudeCurve.zip(frequencies) { a, b -> FreqPair(amp = a, freq = b) }.sortedBy{ it.amp }
    val phases                           = (-50..50 step 5).map{ it.toDouble()}.toList() //TODO: step needs to be 5

    val phaseMagCounts: MutableList<PhaseMagnitudeCount> = phases.map {
        PhaseMagnitudeCount(
            phase            = it,
            currentMagnitude = 0.0
        )
    }.toMutableList()
    val phaseBucket:MutableMap<Double,MutableList<FreqPair>> = phases.fold(initial = mutableMapOf()) {
            acc, curr -> acc.getOrPut(curr) { mutableListOf() }
        acc
    }

    var pos = 0
    while(pos < sortedMagAndFreq.size) {
        val currentPhaseToAugment = phaseMagCounts[0].phase
        phaseBucket[currentPhaseToAugment]!!.add(sortedMagAndFreq[pos])

        phaseMagCounts[0] = PhaseMagnitudeCount(
            phase            = phaseMagCounts[0].phase,
            currentMagnitude = phaseMagCounts[0].currentMagnitude + sortedMagAndFreq[pos].amp
        )
        pos++
        phaseMagCounts.sortBy { it.currentMagnitude }
    }

    val sum  = phaseMagCounts.map { it.currentMagnitude }.sum()
    val dist = phaseMagCounts.map{ (it.currentMagnitude/sum) * 100 }
    //listData(data = dist)

    val ssd = sumOfSquaredDifferences(data = dist)
    val std = standardDeviation(data = dist)
    //println("std: $std, ssd: $ssd")

    val spectralBins: MutableList<SpectralBin> = mutableListOf()
    for(k in phaseBucket.keys) {
        val freqPairs: List<FreqPair> = phaseBucket[k]!!
        spectralBins.addAll(freqPairs.map {
                SpectralBin(magnitude = it.amp, frequency = it.freq, phase = k)
            }
        )
    }

    return synthesizeSpectralBins(bins = spectralBins, length = length)
}

fun generateFrequencies(
    baseFrequency: Double,
    bandWidth: Double,
    frequencyDistance: Double
) : List<Double>{
    val start                            = baseFrequency - (bandWidth/2.0)
    val stop                             = baseFrequency + (bandWidth/2.0)
    var currentFreq                      = start
    val frequencies: MutableList<Double> = mutableListOf()

    frequencies.add(currentFreq)
    while (currentFreq < stop) { //TODO: Should this be <= ?
        currentFreq += frequencyDistance
        frequencies.add(currentFreq)
    }
    return frequencies
}

//TODO: Why length + 1
fun createFormantFrequencyAmplitudeEnvelope(amplitude: Double, coefficient: Double, length: Int) : List<Double> {
    //val signal = interlace(scale = 2.0, input = expF(amplitude = amplitude, coefficient = coefficient, length = length + 1 ))
    //return signal
    val functionalLength = (length / 2) + 1
    val signal = expF(amplitude = amplitude, coefficient = coefficient, length = functionalLength )
    val left   = signal.reversed().subList(0, signal.size - 1).map { it }
    val middle = 0.00//amplitude
    val right  = signal.subList(1, signal.size)

    return left + middle + right
}

/** Coefficient value of -0.25 seems **/
fun expF(amplitude: Double, coefficient: Double,  length: Int) : List<Double> {
    val o: MutableList<Double> = mutableListOf()
    for(i in 0 until length) {
        val value = amplitude* Math.E.pow(coefficient * i.toDouble())
        o.add(value)
    }
    return o
}

class FreqPair(val amp: Double, val freq: Double)


fun synthesizeSpectra(spectra: List<FreqPair>, length: Int) : List<Double> {
    val o: MutableList<Double> = MutableList(length) { 0.0 }
    for(k in spectra.indices) {
        val currentWaveLength = Constants.SAMPLE_RATE/ spectra[k].freq
        val amplitude         = spectra[k].amp
        val ps = (0..length).random()
        //println("currentWaveLength: $currentWaveLength, amp: ${spectra[k].amp}")
        for(i in o.indices) {
            o[i] += amplitude*sin(((2*Math.PI*(i + ps))/currentWaveLength))
        }
    }
    return o
}

fun synthesizeSpectralBins(
    bins: List<SpectralBin>,
    length: Int
) : List<Double> {
    val o: MutableList<Double> = MutableList(length) { 0.0 }
    for(k in bins.indices) {
        val currentWaveLength = Constants.SAMPLE_RATE/ bins[k].frequency
        val amplitude         = bins[k].magnitude
        val currentOffset     = (bins[k].phase/100) * currentWaveLength
        for(i in o.indices) {
            o[i] += amplitude*cos(((2*Math.PI*(i - currentOffset))/currentWaveLength))
        }
    }
    return o
}

fun triangle(centerAmplitude: Double, minAmplitude: Double, length: Int): List<Double> {
    val result = mutableListOf<Double>()
    val midpoint = length / 2
    for (i in 0 until midpoint) {
        val value = minAmplitude + (centerAmplitude - minAmplitude) * (i.toDouble() / midpoint)
        result.add(value)
    }
    if (length % 2 == 0) {
        result.add(centerAmplitude)
    }
    for (i in midpoint until length) {
        val value = centerAmplitude - (centerAmplitude - minAmplitude) * ((i - midpoint).toDouble() / midpoint)
        result.add(value)
    }
    return result
}

