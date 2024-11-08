package readSimulator;

import org.apache.commons.math3.distribution.NormalDistribution;
import readSimulator.utils.FileUtils;
import readSimulator.utils.GenomeUtils;

import java.io.File;
import java.io.IOException;
import java.util.*;

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
        this.genome.readGTF(gtfPath, readCounts);
        this.genome.initTargetGeneSeqs(readCounts);
        this.generateReads(length, frlength, SD, mutRate);
        System.out.println();
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
                    // String fragment = transcriptSeq.substring(randomStartPos, randomStartPos + fragmentLength);

                    // get fw and rw reads
                    String fwSeqRead = transcriptSeq.substring(randomStartPos, randomStartPos + length);
                    String rwSeqRead = GenomeUtils.revComplement(transcriptSeq.substring(randomStartPos + fragmentLength - length, randomStartPos + fragmentLength));
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
                    // System.out.println(fwRead.getReadSeq());
                    // System.out.println(rwSeqRead);
                    // System.out.println(rwRead.getReadSeq());
                    // System.out.println(transcriptSeq);
                    // System.out.println(transcriptSeq.substring(randomStartPos, randomStartPos + fragmentLength));
                    // System.out.println();

                    readId++;
                }

            }

        }
    }

    public void mutateRead(Read read, double mutRate) {
        char[] seqArr = read.getReadSeq();
        int seqLength = seqArr.length;

        // Calculate the expected number of mutations
        int numMutations = (int) Math.round(mutRate * seqLength);

        // Randomly select unique positions to mutate
        HashSet<Integer> mutationPositions = new HashSet<>();
        while (mutationPositions.size() < numMutations) {
            mutationPositions.add(random.nextInt(seqLength));
        }

        // Perform mutations at the selected positions
        for (int pos : mutationPositions) {
            char originalNucleotide = seqArr[pos];
            char newNucleotide;

            // Choose a new nucleotide that is different from the original
            do {
                newNucleotide = NUCLEOTIDES[random.nextInt(NUCLEOTIDES.length)];
            } while (newNucleotide == originalNucleotide);

            seqArr[pos] = newNucleotide; // Mutate the character at the position
            read.addMutPos(pos);
        }
    }

}
