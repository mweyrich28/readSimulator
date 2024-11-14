package org.src;

import org.src.utils.FileUtils;

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

    public byte[] getSequence(String chr, int start, int stop) throws IOException {
        // skip non existent chromosomes
        if (!fastaIdx.containsKey(chr)) {
            return null;
        }

        // get idx data
        ArrayList<Long> faIdx = fastaIdx.get(chr);
        long startInFasta = faIdx.get(1);
        long basesPerLine = faIdx.get(2);
        long bytesPerLine = faIdx.get(3);

        // get file coordinates
        long startLine = start / basesPerLine;
        long startPos = (start - 1) % basesPerLine;

        // cal offset in bytes
        long offset = startInFasta + (startLine * bytesPerLine) + startPos;

        // how many nucleotides do I need to read
        int length = stop - start + 1;

        // how many lines do we expect
        int expectedLines = (int) ((startPos + length) / basesPerLine);
        if ((startPos + length) % basesPerLine == 0) expectedLines--;

        // total bytes to read including newlines
        int bytesToRead = length + expectedLines;

        // jump to correct pos
        this.fasta.seek(offset);

        byte[] buffer = new byte[bytesToRead];
        int totalBytesRead = 0;

        // read all bytes
        while (totalBytesRead < bytesToRead) {
            int read = this.fasta.read(buffer, totalBytesRead, bytesToRead - totalBytesRead);
            if (read == -1) break;  // EOF
            totalBytesRead += read;
        }

        // build sequence, skipping newlines
        // StringBuilder sequence = new StringBuilder(length);
        byte[] byteSequence = new byte[length];
        int basesRead = 0;

        for (int i = 0; i < totalBytesRead && basesRead < length; i++) {
            if (buffer[i] != '\n' && buffer[i] != '\r' && buffer[i] != 0) {
                byteSequence[basesRead] = buffer[i];
                basesRead++;
            }
        }

        return byteSequence;
    }
}
