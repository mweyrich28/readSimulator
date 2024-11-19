package org.src;

import org.src.utils.FileUtils;
import org.src.utils.GenomeUtils;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.HashMap;

public class GenomeSequenceExtractor {

    private RandomAccessFile fasta;

    private HashMap<String, ArrayList<Long>> fastaIdx;

    public GenomeSequenceExtractor(String fastaidxPath, String fastaPath) throws IOException {
        this.fastaIdx = new HashMap<>();
        this.fasta = new RandomAccessFile(new File(fastaPath), "r");
        File fastaIdxFile = new File(fastaidxPath);
        initIndex(fastaIdxFile);
    }

    public void initIndex(File fastaIdxFile) throws IOException {
        // simply store coordinates in a HashMap
        ArrayList<String> lines = FileUtils.readLines(fastaIdxFile);
        for (String line: lines) {
            String[] components = line.split("\t");
            ArrayList<Long> values = new ArrayList<>();
            for (int i = 1; i < components.length; i++) {
                values.add(Long.parseLong(components[i]));
            }
            this.fastaIdx.put(components[0], values);
        }
    }
    public String getSequence(String chr, int start, int stop) throws IOException {
        // this shouldn't happen
        if (start > stop) {
            throw new IllegalArgumentException("Start position must be less than or equal to stop position");
        }

        // skip chrs
        if (!fastaIdx.containsKey(chr)) {
            return null;
        }

        // get index
        ArrayList<Long> faIdx = fastaIdx.get(chr);
        long startInFasta = faIdx.get(1);
        long basesPerLine = faIdx.get(2);
        long bytesPerLine = faIdx.get(3);

        long startLine = (start - 1) / basesPerLine;
        long startPos = (start - 1) % basesPerLine;
        long endLine = (stop - 1) / basesPerLine;

        // calc offset
        long offset = startInFasta + (startLine * bytesPerLine) + startPos;

        // determine total bytes to read accurately
        int length = stop - start + 1;
        int linesToRead = (int)(endLine - startLine + 1);
        int bytesToRead = length + linesToRead;

        byte[] buffer = new byte[bytesToRead];

        // seek and read with robust error handling
        this.fasta.seek(offset);
        int totalBytesRead = 0;
        int bytesRemaining = bytesToRead;

        while (bytesRemaining > 0) {
            int read = this.fasta.read(buffer, totalBytesRead, bytesRemaining);
            if (read == -1) break;  // EOF
            totalBytesRead += read;
            bytesRemaining -= read;
        }

        // init sb
        StringBuilder sequence = new StringBuilder(length);
        int basesRead = 0;

        for (int i = 0; i < totalBytesRead && basesRead < length; i++) {
            char currentChar = (char) buffer[i];
            if (currentChar != '\n' && currentChar != '\r') {
                sequence.append(currentChar);
                basesRead++;
            }
        }
        return sequence.toString();
    }

    public String extractRegion(ArrayList<String> regions, String chromosome, char strand, boolean isRw) throws IOException {
        StringBuilder seq = new StringBuilder();
        for (int j = 0; j < regions.size(); j++) {
            String[] test = regions.get(j).split("-"); //
            int start = Integer.parseInt(test[0]);
            int end = Integer.parseInt(test[1]);
            String subSeq = this.getSequence(chromosome, start, end - 1);
            seq.append(subSeq);
        }
        if (strand == '+' && isRw || strand == '-' && !isRw) {
            return GenomeUtils.revComplement(seq.toString());
        }
        else {
            return seq.toString();
        }
    }
}
