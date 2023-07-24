package org.dsp.wave

import java.io.OutputStream

class WaveUtils {
    companion object {
        fun writeWaveHeader(outFile: OutputStream, numberOfSamples: Int) {
            val totalDataLen   = 36
            val channels       = 1
            val sampleRate     = 44100
            val numChannels    = 1
            val bitsPerSample  = 16
            val byteRate       =  sampleRate * numChannels * (bitsPerSample/8)
            val totalAudioLen  = (numberOfSamples * channels * bitsPerSample) / 8

            val header = ByteArray(44)

            header[0] = 'R'.code.toByte() // RIFF/WAVE header

            header[1] = 'I'.code.toByte()
            header[2] = 'F'.code.toByte()
            header[3] = 'F'.code.toByte()
            header[4] = (totalDataLen and 0xff).toByte()
            header[5] = (totalDataLen shr 8 and 0xff).toByte()
            header[6] = (totalDataLen shr 16 and 0xff).toByte()
            header[7] = (totalDataLen shr 24 and 0xff).toByte()
            header[8] = 'W'.code.toByte()
            header[9] = 'A'.code.toByte()
            header[10] = 'V'.code.toByte()
            header[11] = 'E'.code.toByte()
            header[12] = 'f'.code.toByte() // 'fmt ' chunk

            header[13] = 'm'.code.toByte()
            header[14] = 't'.code.toByte()
            header[15] = ' '.code.toByte()
            header[16] = 16 // 4 bytes: size of 'fmt ' chunk

            header[17] = 0
            header[18] = 0
            header[19] = 0
            header[20] = 1 // format = 1

            header[21] = 0
            header[22] = channels.toByte()
            header[23] = 0
            header[24] = (sampleRate and 0xff).toByte()
            header[25] = (sampleRate shr 8 and 0xff).toByte()
            header[26] = (sampleRate shr 16 and 0xff).toByte()
            header[27] = (sampleRate shr 24 and 0xff).toByte()
            header[28] = (byteRate and 0xff).toByte()
            header[29] = (byteRate shr 8 and 0xff).toByte()
            header[30] = (byteRate shr 16 and 0xff).toByte()
            header[31] = (byteRate shr 24 and 0xff).toByte()
            header[32] = (2 * 16 / 8).toByte() // block align

            header[33] = 0
            header[34] = 16 // bits per sample

            header[35] = 0
            header[36] = 'd'.code.toByte()
            header[37] = 'a'.code.toByte()
            header[38] = 't'.code.toByte()
            header[39] = 'a'.code.toByte()
            header[40] = (totalAudioLen and 0xff).toByte()
            header[41] = (totalAudioLen shr 8 and 0xff).toByte()
            header[42] = (totalAudioLen shr 16 and 0xff).toByte()
            header[43] = (totalAudioLen shr 24 and 0xff).toByte()

            outFile.write(header, 0, 44)
        }

    }
}