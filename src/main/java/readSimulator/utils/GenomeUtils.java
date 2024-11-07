package readSimulator.utils;

public class GenomeUtils {

    public static String revComplement(String dnaSeq) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < dnaSeq.length(); i++) {
            switch (dnaSeq.charAt(i)) {
                case 'A':
                    sb.append('T');
                    break;
                case 'T':
                    sb.append('A');
                    break;
                case 'G':
                    sb.append('C');
                    break;
                case 'C':
                    sb.append('G');
                    break;
            }
        }
        return sb.reverse().toString();
    }
}
