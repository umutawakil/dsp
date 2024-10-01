package org.dsp


import org.dsp.wave.WaveUtils

fun main() {
    val resourceDirectory  = "/Users/umutawakil/Documents/Git/dsp/src/main/resources"
    val waveData = WaveUtils.getFileData(resourceDirectory = resourceDirectory, fileName = "normalized.wave")



   /** ----------------------------------------------------------------**/

    /*writeSamplesToFile(
        fileName = "$resourceDirectory/output.wav",
        input    =  normalize(
                        scale = 5000.0,
                        input = sineE//rand2// randomSynth(waveLength = base, size = base*400, vibrato = generateVibratoSignal(noBase = false, base = base.toDouble(), size = base*400))//stretchedWaves.flatten()//waveDiff.map{it + listOf(2500.0)}.flatten()//pitchedNoiseWaves.flatten()//integratedWaves.flatten()//staticTestSignalWaves.flatten()
                    )
    ) */
}

/** Skeleton code: Kept only because it uses the interpolation method masterfully **/
/*fun randomSynth(waveLength: Int, size: Int, vibrato: List<Double>) : List<Double> {
    val o: MutableList<Double> = mutableListOf()


    /**Normalize the under-waveform and add in the fundamental **/
    val waves: MutableList<List<Double>> = mutableListOf()
    var c = 0
    val half = waveLength/2
    var direction = 1.0
    while(true) {
        if(c + half - 1 >= o.size) {
            break
        }
        val halfWave = o.subList(c, c + half)
        val normalizedUnderWaveform = normalize(scale = 200.0, input = halfWave)
        val fundamentalWaveform = sineByCycles(amplitude = direction*5000.0, cycles = 0.5, size = half)
        waves.add(fundamentalWaveform.zip(normalizedUnderWaveform) { a, b -> a + b})
        //waves.add(halfWave)
        c += half
        direction *= -1
    }

    /** Apply vibrato to the generated signal **/
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
} */

