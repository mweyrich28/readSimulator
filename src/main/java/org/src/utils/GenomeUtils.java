package org.src.utils;

public class GenomeUtils {

    public static final int BASE_A = 65;
    public static final int BASE_T = 84;
    public static final int BASE_G = 71;
    public static final int BASE_C = 67;

    public static String revComplement(String dnaSeq) {
        int len = dnaSeq.length();
        char[] result = new char[len];

        for (int i = 0; i < len; i++) {
            char complement;
            switch (dnaSeq.charAt(len - 1 - i)) {
                case 'A' -> complement = 'T';
                case 'T' -> complement = 'A';
                case 'G' -> complement = 'C';
                case 'C' -> complement = 'G';
                default ->
                        throw new IllegalArgumentException("Invalid DNA sequence character: " + dnaSeq.charAt(len - 1 - i));
            }
            result[i] = complement;
        }
        return new String(result);
    }

    public static byte[] revComplement(byte[] dnaSeq) {
        int len = dnaSeq.length;
        byte[] newSeq = new byte[len];
        for (int i = 0; i < len; i++) {
            byte complement;
            switch (dnaSeq[len - 1 - i]) {
                case BASE_A -> complement = BASE_T;
                case BASE_T -> complement = BASE_A;
                case BASE_G -> complement = BASE_C;
                case BASE_C -> complement = BASE_G;
                default ->
                        throw new IllegalArgumentException("Invalid DNA sequence character: " + dnaSeq[len - 1 - i]);
            }
            newSeq[i] = complement;
        }
        return newSeq;
    }

    public static byte[] subsetArray(byte[] arrToFill, byte[] sourceArr, int start, int end) {
        int pos = 0;
        for (int i = start; i <= end; i++) {
            arrToFill[pos] = sourceArr[i];
            pos++;
        }
        return arrToFill;
    }

    public static byte[] buildFastqEntry(String seqId, byte[] seqBytes, String qId, byte[] qualityBytes, int readId) {
        byte[] seqIdBytes = seqId.getBytes();
        byte[] qIdBytes = qId.getBytes();

        // Calculate total length
        // If readId != 0, we'll add the extra '\n' at the beginning
        int totalLength = seqIdBytes.length + seqBytes.length + 1 + qIdBytes.length + qualityBytes.length;
        if (readId != 0) {
            totalLength += 1; // Add extra space for the newline at the beginning
        }

        byte[] entryData = new byte[totalLength];
        int index = 0;

        if (readId != 0) {
            // If readId != 0, add a newline at the beginning
            entryData[index++] = '\n';
        }

        // Add sequence ID (@0, @1, etc.)
        System.arraycopy(seqIdBytes, 0, entryData, index, seqIdBytes.length);
        index += seqIdBytes.length;

        // Add sequence
        System.arraycopy(seqBytes, 0, entryData, index, seqBytes.length);
        index += seqBytes.length;

        // Newline after sequence
        entryData[index++] = '\n';

        // Add quality identifier (+0, +1, etc.)
        System.arraycopy(qIdBytes, 0, entryData, index, qIdBytes.length);
        index += qIdBytes.length;

        // Add quality scores
        System.arraycopy(qualityBytes, 0, entryData, index, qualityBytes.length);

        return entryData;
    }
}
