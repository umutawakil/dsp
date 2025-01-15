package org.dsp.tools

import java.io.File

class Utils {
    companion object {


    }
}

fun List<Double>.add(other: List<Double>): List<Double> {
    return this.zip(other) { a, b -> a + b }
}

fun writeDoublesToFile(doubles: List<Double>, fileName: String) {
    try {
        val file = File(fileName)
        file.delete()
        file.createNewFile()
        file.bufferedWriter().use { writer ->
            doubles.forEach { number ->
                writer.write(number.toString())
                writer.newLine()
            }
        }
        println("Successfully wrote ${doubles.size} doubles to $fileName")

    } catch (e: Exception) {
        throw RuntimeException("An error occurred while writing to the file: ${e.message}")
    }
}

fun writeStringsToFile(strings: List<String>, fileName: String) {
    try {
        File(fileName).bufferedWriter().use { writer ->
            strings.forEach { str ->
                writer.write(str)
                writer.newLine()
            }
        }
        println("Successfully wrote ${strings.size} strings to $fileName")
    } catch (e: Exception) {
        println("An error occurred while writing to the file: ${e.message}")
    }
}