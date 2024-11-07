package readSimulator;

import readSimulator.utils.FileUtils;
import readSimulator.utils.GenomeUtils;

import java.io.File;
import java.io.FileNotFoundException;
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

    public String getSequence1(String chr, int start, int stop, boolean revStrand) throws IOException {
        // skip chrs that don't exist in idx
        if (!fastaIdx.containsKey(chr)) {
            return null;
        }

        long startInFasta = fastaIdx.get(chr).get(2);
        this.fasta.seek(startInFasta + start);

        // how far do i read
        int length = stop - start;
        StringBuilder sequence = new StringBuilder();

        int bytesRead = 0;
        while (bytesRead < length) {
            String line = this.fasta.readLine();
            // EOF
            if (line == null) break;

            // just in case
            // line = line.trim();

            // how much further
            int remaining = length - bytesRead;
            if (line.length() > remaining) {
                sequence.append(line.substring(0, remaining));
                break;
            } else {
                sequence.append(line);
                bytesRead += line.length();
            }
        }
        if (revStrand) {
            // TODO: check if this is valid
            return sequence.reverse().toString();
        } else {
            return sequence.toString();
        }
    }
    public String getSequence2(String chr, int start, int stop, boolean revStrand) throws IOException {
        // Ensure the chromosome exists in the index
        if (!fastaIdx.containsKey(chr)) {
            return null;
        }

        // Get FASTA index information for the chromosome
        ArrayList<Long> faIdx = fastaIdx.get(chr);
        long startInFasta = faIdx.getFirst();
        int lineWidth = 60;
        int bytesPerLine = 61;

        // Calculate the number of lines before the start position (zero-based)
        long lineOffset = start / lineWidth;
        long positionInLine = start % lineWidth;

        // Calculate the offset in the FASTA file
        long offset = startInFasta + lineOffset * bytesPerLine + positionInLine;

        // Seek to the calculated offset in the FASTA file
        this.fasta.seek(offset);

        // Calculate the length to read
        int length = stop - start + 1;
        StringBuilder sequence = new StringBuilder();
        int bytesRead = 0;

        while (bytesRead < length) {
            String line = this.fasta.readLine();
            System.out.println(line);
            if (line == null) break; // End of file

            int remaining = length - bytesRead;
            if (line.length() > remaining) {
                sequence.append(line, 0, remaining);
                break;
            } else {
                sequence.append(line);
                bytesRead += line.length();
            }
        }

        // Reverse complement if needed
        if (revStrand) {
            return GenomeUtils.revComplement(sequence.toString());
        } else {
            return sequence.toString();
        }
    }
    public String getSequence4(String chr, int start, int stop, boolean revStrand) throws IOException {
        // Ensure the chromosome exists in the index
        if (!fastaIdx.containsKey(chr)) {
            return null;
        }

        // Get FASTA index information for the chromosome
        ArrayList<Long> faIdx = fastaIdx.get(chr);
        long startInFasta = faIdx.getFirst();
        final int BASES_PER_LINE = 60;
        final int BYTES_PER_LINE = 61; // 60 bases + 1 newline character

        // Calculate the starting position
        long lineOffset = start / BASES_PER_LINE;
        long positionInLine = start % BASES_PER_LINE;
        long offset = startInFasta + (lineOffset * BYTES_PER_LINE) + positionInLine;

        // Calculate the length to read
        int length = stop - start + 1;

        // Seek to the calculated offset in the FASTA file
        this.fasta.seek(offset);

        // Read the sequence
        StringBuilder sequence = new StringBuilder(length);
        int basesRead = 0;
        byte[] buffer = new byte[8192]; // Use buffered reading for better performance

        while (basesRead < length) {
            int bytesToRead = Math.min(buffer.length, (length - basesRead) + ((length - basesRead) / BASES_PER_LINE) + 1);
            int bytesRead = this.fasta.read(buffer, 0, bytesToRead);

            if (bytesRead == -1) break; // End of file

            // Process the bytes read
            for (int i = 0; i < bytesRead && basesRead < length; i++) {
                // Skip newline characters
                if (buffer[i] != '\n' && buffer[i] != '\r') {
                    sequence.append((char) buffer[i]);
                    basesRead++;
                }
            }
        }

        String result = sequence.toString().toUpperCase();

        // Reverse complement if requested
        if (revStrand) {
            result = GenomeUtils.revComplement(result);
        }

        return result;
    }

    public String getSequence(String chr, int start, int stop, boolean revStrand) throws IOException {
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
        long offset = startInFasta;
        offset += (startLine * bytesPerLine);
        offset += startPos;

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
        StringBuilder sequence = new StringBuilder(length);
        int basesRead = 0;

        for (int i = 0; i < totalBytesRead && basesRead < length; i++) {
            if (buffer[i] != '\n' && buffer[i] != '\r') {
                sequence.append((char) buffer[i]);
                basesRead++;
            }
        }

        String result = sequence.toString().toUpperCase();

        if (revStrand) {
            result = GenomeUtils.revComplement(result);
        }

        return result;
    }
    // Method to get a sequence from the FASTA file
    public String getSequenceaa(String chr, int start, int end, boolean revStrand) throws IOException {
        // Check if the chromosome exists in the index map
        if (!this.fastaIdx.containsKey(chr)) {
           return null;
        }

        // Get the file offset for the start of the chromosome
        long chrOffset = fastaIdx.get(chr).get(1);

        // Calculate the offset to seek to within the chromosome
        long seekPosition = chrOffset + start;

        // Seek to the start position
        fasta.seek(seekPosition);

        // Read the sequence of the specified length
        int length = end - start + 1;
        byte[] buffer = new byte[length];
        fasta.read(buffer, 0, length);

        // Convert the buffer to a String and return it
        if (revStrand) {
           return GenomeUtils.revComplement(new String(buffer));
        } else {
            return new String(buffer);
        }
    }

    // Close method to release file resources
    public void close() throws IOException {
        if (fasta != null) {
            fasta.close();
        }
    }
}
