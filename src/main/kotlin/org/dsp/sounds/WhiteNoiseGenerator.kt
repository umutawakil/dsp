package org.dsp.sounds

import org.dsp.file.FileUtils
import java.io.DataOutputStream
import kotlin.math.pow

class WhiteNoiseGenerator {
    companion object {
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
    }
}