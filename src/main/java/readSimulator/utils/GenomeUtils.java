package readSimulator.utils;

public class GenomeUtils {

    public static String revComplement(String dnaSeq) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < dnaSeq.length(); i++) {
            switch (dnaSeq.charAt(i)) {
                case 'A' -> sb.append('T');
                case 'T' -> sb.append('A');
                case 'G' -> sb.append('C');
                case 'C' -> sb.append('G');
            }
        }
        return sb.reverse().toString();
    }
}
