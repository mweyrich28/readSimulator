package readSimulator;

import org.apache.commons.math3.distribution.NormalDistribution;
import readSimulator.utils.FileUtils;
import readSimulator.utils.GenomeUtils;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

public class ReadSimulator {

    private static final char[] NUCLEOTIDES = {'A', 'T', 'C', 'G'};
    private final Random random = new Random();
    private Genome genome;
    private HashMap<String, HashMap<String, Integer>> readCounts = new HashMap<>();

    public ReadSimulator(int length,
                         int frlength,
                         int SD,
                         double mutRate,
                         String gtfPath,
                         String readCountsPath,
                         String fastaPath,
                         String idxPath,
                         String od) throws IOException
    {

        this.initReadCounts(readCountsPath);
        this.genome = new Genome(idxPath, fastaPath);
        this.genome.readGTF(gtfPath);
        this.genome.initTargetGeneSeqs(readCounts);
        this.generateReads(length, frlength, SD, mutRate);
    }

    public void initReadCounts(String readCountsPath) throws IOException {
        ArrayList<String> lines = FileUtils.readLines(new File(readCountsPath));

        for (String line: lines) {
            // skip header
            if (line.startsWith("gene")) {
                continue;
            }

            String[] components = line.split("\t");
            HashMap<String, Integer> innerMap = new HashMap<>();
            if (readCounts.containsKey(components[0])) {
                readCounts.get(components[0]).put(components[1], Integer.parseInt(components[2]));
            } else {
                readCounts.put(components[0], innerMap);
                innerMap.put(components[1], Integer.parseInt(components[2]));
            }
        }
    }


    public void generateReads(int length, int frlength, int SD, double mutRate) {
        int readId = 0;
        NormalDistribution normalDist = new NormalDistribution(frlength, SD);

        for (String geneKey : readCounts.keySet()) {

            for (String transcriptKey : readCounts.get(geneKey).keySet()) {
                System.out.println("Generating Reads in " + geneKey + " " + transcriptKey + " " + readId);

                Transcript transcript = this.genome.getGenes().get(geneKey).getTranscriptMap().get(transcriptKey);
                String transcriptSeq = transcript.getTranscriptSeq();
                int sampleAmount = readCounts.get(geneKey).get(transcriptKey);

                for (int i = 0; i < sampleAmount; i++) {

                    int fragmentLength;
                    // get frlength
                    do {
                        fragmentLength = (int) Math.round(normalDist.sample());
                    } while (fragmentLength < length || fragmentLength > transcriptSeq.length());

                    // calc max pos start pos
                    int maxStartPos = transcriptSeq.length() - fragmentLength;
                    int randomStartPos = random.nextInt(maxStartPos + 1);

                    // get fragment
                    String fragment = transcriptSeq.substring(randomStartPos, randomStartPos + fragmentLength);

                    // get fw and rw reads
                    String fwSeqRead = fragment.substring(0, length);
                    String rwSeqRead = GenomeUtils.revComplement(fragment.substring(fragment.length() - length, fragment.length()));
                    int fwStart = randomStartPos;
                    int fwEnd = randomStartPos + length;
                    int rwStart = randomStartPos + fragmentLength - length;
                    int rwEnd = randomStartPos + fragmentLength;
                    Read fwRead = new Read(fwSeqRead, fwStart, fwEnd, readId, false);
                    Read rwRead = new Read(rwSeqRead, rwStart, rwEnd, readId, true);
                    mutateRead(fwRead, mutRate);
                    mutateRead(rwRead, mutRate);
                    transcript.addReads(fwRead, rwRead);

                    // System.out.println(transcriptKey);
                    // System.out.println(fwSeqRead);
                    // System.out.println(rwSeqRead);
                    // System.out.println(transcriptSeq);
                    // System.out.println(fragment);
                    // System.out.println();

                    readId++;
                }

            }

        }
    }

    public void mutateRead(Read read, double mutRate) {
        StringBuilder mutatedSeq = new StringBuilder();
        String orgSeq = read.getReadSeq();

        for (int i = 0; i < orgSeq.length(); i++) {
            char originalNucleotide = orgSeq.charAt(i);
            char newNucleotide = originalNucleotide;

            // Check if a mutation should occur
            if (random.nextDouble() < mutRate) {
                // Choose a new nucleotide that is not the same as the original
                do {
                    newNucleotide = NUCLEOTIDES[random.nextInt(NUCLEOTIDES.length)];
                } while (newNucleotide == originalNucleotide);
            }

            // Append the chosen nucleotide to mutatedSeq
            mutatedSeq.append(newNucleotide);
        }

        // Set the mutated sequence in the Read object
        read.setReadSeq(mutatedSeq.toString());
    }
}
