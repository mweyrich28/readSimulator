package readSimulator.utils;

import java.io.*;
import java.util.ArrayList;

public class LineFilter {
    public static ArrayList<String> readExonLines(File r) throws IOException {
        ArrayList<String> fileLines = new ArrayList<>();
        BufferedReader buff = new BufferedReader(new FileReader(r));
        String line;
        StringBuilder relevantCol = new StringBuilder();

        while((line = buff.readLine()) != null){
            // only read in lines with CDS / exon
            // skip comments
            if (line.charAt(0) == '#') {
                fileLines.add(line);
                continue;
            }

            // save memory by selecting correct col based on tabCount
            int tabCount = 0;
            relevantCol.setLength(0); // reset sb
            for (int i = 0; i < line.length(); i++) {
                if (line.charAt(i) == '\t') {
                    tabCount++;
                }
                else if (tabCount == 2) {
                    relevantCol.append(line.charAt(i));
                }

                if ("exon".contentEquals(relevantCol) || "CDS".contentEquals(relevantCol)) {
                    fileLines.add(line);
                    break;
                }

                if (tabCount == 3) {
                    break;
                }
            }

        }
        return fileLines;
    }

    public static void writeFile(String fileName, ArrayList<String> content) throws IOException {
        BufferedWriter bw = new BufferedWriter(new FileWriter(new File(fileName)));

        for (int i = 0; i < content.size(); i++) {
            bw.write(content.get(i));
            bw.newLine();
        }
        bw.close();
    }
}



