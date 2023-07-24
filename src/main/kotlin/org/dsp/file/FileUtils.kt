package org.dsp.file

import java.io.DataOutputStream

class FileUtils {
    companion object {
        fun writeShortLE(out: DataOutputStream, value: Short) {
            val byteA = value.toInt() and 0xFF
            val byteB = value.toInt() shr 8 and 0xFF
            out.writeByte(byteA)
            out.writeByte(byteB)
        }
    }
}