package org.dsp

import org.dsp.sounds.WhiteNoiseGenerator
import org.dsp.wave.WaveUtils
import java.io.*
fun main() {
    val resourceDirectory = "/Users/umutawakil/Documents/Git/dsp/src/main/resources"
    val file             = File("$resourceDirectory/test.wav")
    file.delete()
    file.createNewFile()
    val os: OutputStream = FileOutputStream(file)
    val bos              = BufferedOutputStream(os)
    val outFile          = DataOutputStream(bos)
    val numberOfSamples  = 44100*5

    /** Create header**/
    WaveUtils.writeWaveHeader(
        outFile         = outFile,
        numberOfSamples = numberOfSamples
    )

    /** Create data and write it here **/
    WhiteNoiseGenerator.writeWhiteNoise(
        outFile    = outFile,
        numSamples = numberOfSamples
    )

    outFile.flush()
    outFile.close()
}