package readSimulator.utils;

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
}
