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
        this.generateReads(length, frlength, SD);
        this.mutateReads(mutRate);
    }

    public void mutateReads(double mutRate) {
        for (String geneKey : readCounts.keySet()) {
            for (String transcriptKey : readCounts.get(geneKey).keySet()) {
            }
        }

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


    public void generateReads(int length, int frlength, int SD) {
        int readId = 0;
        NormalDistribution normalDist = new NormalDistribution(frlength, SD);
        Random random = new Random();

        for (String geneKey : readCounts.keySet()) {

            for (String transcriptKey : readCounts.get(geneKey).keySet()) {

                String transcriptSeq = this.genome.getGenes().get(geneKey).getTranscriptMap().get(transcriptKey).getTranscriptSeq();
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
                    System.out.println(transcriptKey);
                    System.out.println(fwSeqRead);
                    System.out.println(rwSeqRead);
                    System.out.println(transcriptSeq);
                    System.out.println(fragment);
                    System.out.println();
                }
            }
        }
    }
}
