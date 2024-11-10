package org.src.utils;

public class GenomeUtils {

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

    private static int getNumericValue(char nucleotide) {
        return switch (nucleotide) {
            case 'A' -> 0;
            case 'C' -> 1;
            case 'G' -> 2;
            case 'T' -> 3;
            default -> throw new IllegalArgumentException("Invalid nucleotide: " + nucleotide);
        };
    }
}
