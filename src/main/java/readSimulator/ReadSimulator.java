package readSimulator;

import org.apache.commons.math3.distribution.NormalDistribution;
import readSimulator.utils.FileUtils;
import readSimulator.utils.GenomeUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
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
        // this.generateReads(length, frlength, SD, mutRate);
        this.generateReads(length, frlength, SD, mutRate, 100, od);
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


    public void generateReads(int length, int frlength, int SD, double mutRate) throws IOException {
        int readId = 0;
        NormalDistribution normalDist = new NormalDistribution(frlength, SD);
        StringBuilder mutSeqBuilder = new StringBuilder();
        StringBuilder summaryBuilder = new StringBuilder();
        ArrayList<Read> fwReads = new ArrayList<>();
        ArrayList<Read> rwReads = new ArrayList<>();
        ArrayList<String> summaryList = new ArrayList<>();

        for (String geneKey : readCounts.keySet()) {
            Gene currGene = genome.getGenes().get(geneKey);

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
                    int fwEnd = randomStartPos + length - 1;
                    int rwStart = randomStartPos + fragmentLength - length;
                    int rwEnd = randomStartPos + fragmentLength - 1;
                    Read fwRead = new Read(fwSeqRead, fwStart, fwEnd, readId, false);
                    Read rwRead = new Read(rwSeqRead, rwStart, rwEnd, readId, true);
                    mutateRead(fwRead, mutRate, mutSeqBuilder);
                    mutateRead(rwRead, mutRate, mutSeqBuilder);
                    ArrayList<String> fwRegVec = getGenomicRegion(fwRead, transcript, currGene.getStrand());
                    ArrayList<String> rwRegVec = getGenomicRegion(rwRead, transcript, currGene.getStrand());

                    summaryBuilder.append(readId + "\t");
                    summaryBuilder.append(currGene.getChr() + "\t");
                    summaryBuilder.append(currGene.getGeneId() + "\t");
                    summaryBuilder.append(transcript.getTranscriptId() + "\t");
                    summaryBuilder.append(String.join("|", fwRegVec) + "\t");
                    summaryBuilder.append(String.join("|", rwRegVec) + "\t");
                    summaryBuilder.append(fwStart+ "-" + fwEnd + "\t");
                    summaryBuilder.append(rwStart+ "-" + rwEnd + "\t");
                    summaryBuilder.append(String.join(",", fwRead.getMutPos()));
                    summaryList.add(summaryBuilder.toString());
                    summaryBuilder.setLength(0);

                    fwReads.add(fwRead);
                    rwReads.add(rwRead);

                    // maybe don't need that
                    // transcript.addReads(fwRead, rwRead);

                    // if(currGene.getStrand() == '+') {
                    //     System.out.println("ID " + transcriptKey + " " + currGene.getStrand());
                    //     System.out.println("S  " + transcriptSeq);
                    //     System.out.println("FR " + transcriptSeq.substring(randomStartPos, randomStartPos + fragmentLength));
                    //     System.out.println("FW " + fwSeqRead);
                    //     System.out.println(fwRead.getStartInTranscript() + "-" + fwRead.getStopInTranscript());
                    //     ArrayList<String> reg = getGenomicRegion(fwRead, transcript, currGene.getStrand());
                    //     System.out.println(reg);
                    //     System.out.println("FM " + fwRead.getReadSeq());
                    //     System.out.print("LO " );
                    //     for (int j = 0; j < reg.size(); j++) {
                    //         String[] test = reg.get(j).split("-");
                    //         if (test.length == 1) {
                    //             int val = Integer.parseInt(test[0]);
                    //             System.out.print(genome.getGSE().getSequence(genome.getGenes().get(geneKey).getChr(), val, val));
                    //         }
                    //         else {
                    //             int start = Integer.parseInt(test[0]);
                    //             int end = Integer.parseInt(test[1]);
                    //             System.out.print(genome.getGSE().getSequence(genome.getGenes().get(geneKey).getChr(), start, end));
                    //         }
                    //     }
                    //     System.out.println();
                    //     System.out.println(rwRead.getStartInTranscript() + "-" + rwRead.getStopInTranscript());
                    //     reg = getGenomicRegion(rwRead, transcript, currGene.getStrand());
                    //     System.out.println(reg);
                    //     System.out.println("RM " + rwRead.getReadSeq());
                    //     System.out.print("LO " );
                    //     for (int j = 0; j < reg.size(); j++) {
                    //         String[] test = reg.get(j).split("-");
                    //         if (test.length == 1) {
                    //             int val = Integer.parseInt(test[0]);
                    //             System.out.print(genome.getGSE().getSequence(genome.getGenes().get(geneKey).getChr(), val, val));
                    //         }
                    //         else {
                    //             int start = Integer.parseInt(test[0]);
                    //             int end = Integer.parseInt(test[1]);
                    //             System.out.print(GenomeUtils.revComplement(genome.getGSE().getSequence(genome.getGenes().get(geneKey).getChr(), start, end)));
                    //         }
                    //     }
                    //     System.out.println("\n--------------------------------------------------------------------------------------------------------------------------------");
                    // } else {
                    //     System.out.println("ID " + transcriptKey + " " + currGene.getStrand());
                    //     System.out.println("S  " + transcriptSeq);
                    //     System.out.println("FR " + transcriptSeq.substring(randomStartPos, randomStartPos + fragmentLength));
                    //     System.out.println("FW " + fwSeqRead);
                    //     System.out.println(fwRead.getStartInTranscript() + "-" + fwRead.getStopInTranscript());
                    //     ArrayList<String> reg = getGenomicRegion(fwRead, transcript, currGene.getStrand());
                    //     System.out.println(reg);
                    //     System.out.println("FM " + fwRead.getReadSeq());
                    //     System.out.print("LO " );
                    //     for (int j = 0; j < reg.size(); j++) {
                    //         String[] test = reg.get(j).split("-");
                    //         if (test.length == 1) {
                    //             int val = Integer.parseInt(test[0]);
                    //             System.out.print(GenomeUtils.revComplement(genome.getGSE().getSequence(genome.getGenes().get(geneKey).getChr(), val, val)));
                    //         }
                    //         else {
                    //             int start = Integer.parseInt(test[0]);
                    //             int end = Integer.parseInt(test[1]);
                    //             System.out.print(GenomeUtils.revComplement(genome.getGSE().getSequence(genome.getGenes().get(geneKey).getChr(), start, end)));
                    //         }
                    //     }
                    //     System.out.println();
                    //     System.out.println(rwRead.getStartInTranscript() + "-" + rwRead.getStopInTranscript());
                    //     reg = getGenomicRegion(rwRead, transcript, currGene.getStrand());
                    //     System.out.println(reg);
                    //     System.out.println("RM " + rwRead.getReadSeq());
                    //     System.out.print("LO " );
                    //     for (int j = 0; j < reg.size(); j++) {
                    //         String[] test = reg.get(j).split("-");
                    //         if (test.length == 1) {
                    //             int val = Integer.parseInt(test[0]);
                    //             System.out.print(GenomeUtils.revComplement(genome.getGSE().getSequence(genome.getGenes().get(geneKey).getChr(), val, val)));
                    //         }
                    //         else {
                    //             int start = Integer.parseInt(test[0]);
                    //             int end = Integer.parseInt(test[1]);
                    //             System.out.print(genome.getGSE().getSequence(genome.getGenes().get(geneKey).getChr(), start, end));
                    //         }
                    //     }
                    //     System.out.println("\n--------------------------------------------------------------------------------------------------------------------------------");

                    // }


                    readId++;
                }
            }

        }
    }
    public void generateReads(int length, int frlength, int SD, double mutRate, int batchSize, String od) throws IOException {

        // create od if not existent
        File outputDir = new File(od);
        if (!outputDir.exists()) {
            outputDir.mkdirs();
        }

        int readId = 0;

        NormalDistribution normalDist = new NormalDistribution(frlength, SD);

        // init string builders
        StringBuilder summaryBuilder = new StringBuilder();
        StringBuilder fwFastqBuilder = new StringBuilder();
        StringBuilder rwFastqBuilder = new StringBuilder();
        StringBuilder mutateSeqBuilder = new StringBuilder();

        // for each file init buff writer
        BufferedWriter summaryWriter = new BufferedWriter(new FileWriter(od + File.separator + "read.mappinginfo"));
        BufferedWriter fwFastqWriter = new BufferedWriter(new FileWriter(od + File.separator + "fw.fastq"));
        BufferedWriter rwFastqWriter = new BufferedWriter(new FileWriter(od + File.separator + "rw.fastq"));
        summaryWriter.write("readid\tchr\tgene\ttranscript\tfw_regvec\trw_regvec\tt_fw_regvec\tt_rw_regvec\tfw_mut\trw_mut");




        for (String geneKey : readCounts.keySet()) {
            Gene currGene = genome.getGenes().get(geneKey);

            for (String transcriptKey : readCounts.get(geneKey).keySet()) {
                Transcript transcript = this.genome.getGenes().get(geneKey).getTranscriptMap().get(transcriptKey);
                String transcriptSeq = transcript.getTranscriptSeq();
                int sampleAmount = readCounts.get(geneKey).get(transcriptKey);

                for (int i = 0; i < sampleAmount; i++) {
                    int fragmentLength;
                    do {
                        fragmentLength = (int) Math.round(normalDist.sample());
                    } while (fragmentLength < length || fragmentLength > transcriptSeq.length());

                    int maxStartPos = transcriptSeq.length() - fragmentLength;
                    int randomStartPos = random.nextInt(maxStartPos + 1);

                    String fwSeqRead = transcriptSeq.substring(randomStartPos, randomStartPos + length);
                    String rwSeqRead = GenomeUtils.revComplement(transcriptSeq.substring(randomStartPos + fragmentLength - length, randomStartPos + fragmentLength));
                    int fwStart = randomStartPos;
                    int fwEnd = randomStartPos + length - 1;
                    int rwStart = randomStartPos + fragmentLength - length;
                    int rwEnd = randomStartPos + fragmentLength - 1;

                    Read fwRead = new Read(fwSeqRead, fwStart, fwEnd, readId, false);
                    Read rwRead = new Read(rwSeqRead, rwStart, rwEnd, readId, true);

                    mutateRead(fwRead, mutRate, mutateSeqBuilder);
                    mutateRead(rwRead, mutRate, mutateSeqBuilder);

                    ArrayList<String> fwRegVec = getGenomicRegion(fwRead, transcript, currGene.getStrand());
                    ArrayList<String> rwRegVec = getGenomicRegion(rwRead, transcript, currGene.getStrand());

                    // concat summary col
                    summaryBuilder.append("\n").append(readId).append("\t")
                            .append(currGene.getChr()).append("\t")
                            .append(currGene.getGeneId()).append("\t")
                            .append(transcript.getTranscriptId()).append("\t")
                            .append(String.join("|", fwRegVec)).append("\t")
                            .append(String.join("|", rwRegVec)).append("\t")
                            .append(fwStart).append("-").append(fwEnd).append("\t")
                            .append(rwStart).append("-").append(rwEnd).append("\t")
                            .append(String.join(",", fwRead.getMutPos()));

                    // format read seqs in fastq files
                    String quality = "I".repeat(fwSeqRead.length());
                    // avoid line breaks at end of file
                    if (readId == 0) {
                        fwFastqBuilder.append("@").append(readId).append("\n").append(fwSeqRead).append("\n").append("+").append(readId).append("\n").append(quality);
                        rwFastqBuilder.append("@").append(readId).append("\n").append(rwSeqRead).append("\n").append("+").append(readId).append("\n").append(quality);
                    }
                    else {
                        fwFastqBuilder.append("\n").append("@").append(readId).append("\n").append(fwSeqRead).append("\n").append("+").append(readId).append("\n").append(quality);
                        rwFastqBuilder.append("\n").append("@").append(readId).append("\n").append(rwSeqRead).append("\n").append("+").append(readId).append("\n").append(quality);

                    }

                    readId++;

                    // batch writing
                    if (readId % batchSize == 0) {
                        summaryWriter.write(summaryBuilder.toString());
                        fwFastqWriter.write(fwFastqBuilder.toString());
                        rwFastqWriter.write(rwFastqBuilder.toString());

                        // flush sbs
                        summaryBuilder.setLength(0);
                        fwFastqBuilder.setLength(0);
                        rwFastqBuilder.setLength(0);
                    }
                }
            }
        }

        if (!summaryBuilder.isEmpty()) {
            summaryWriter.write(summaryBuilder.toString());
        }
        if (!fwFastqBuilder.isEmpty()) {
            fwFastqWriter.write(fwFastqBuilder.toString());
        }
        if (!rwFastqBuilder.isEmpty()) {
            rwFastqWriter.write(rwFastqBuilder.toString());
        }

        summaryWriter.flush();
        fwFastqWriter.flush();
        rwFastqWriter.flush();
    }

    public ArrayList<String> getGenomicRegion(Read read, Transcript transcript, char strand) {
        ArrayList<Exon> exons = transcript.getExonList();
        ArrayList<String> genomicRegions = new ArrayList<>();

        int readStart = read.getStartInTranscript();
        int readEnd = read.getStopInTranscript();

        int transcriptPosition = 0; // Tracks position in the transcript (1-based)

        for (Exon exon : exons) {
            int exonStart = exon.getGenomicStart();
            int exonEnd = exon.getGenomicEnd();
            int exonLength = exonEnd - exonStart + 1;

            // Check if the read overlaps with this exon in transcript coordinates
            if (transcriptPosition + exonLength - 1 >= readStart && transcriptPosition <= readEnd) {
                // Calculate overlap within transcript coordinates
                int overlapStartInTranscript = Math.max(readStart, transcriptPosition);
                int overlapEndInTranscript = Math.min(readEnd, transcriptPosition + exonLength - 1);

                // Map overlap to genomic coordinates based on strand
                int overlapStartInGenomic;
                int overlapEndInGenomic;

                if (strand == '+') {
                    // Forward strand: Calculate genomic positions from exon start
                    overlapStartInGenomic = exonStart + (overlapStartInTranscript - transcriptPosition);
                    overlapEndInGenomic = exonStart + (overlapEndInTranscript - transcriptPosition);
                } else {
                    // Reverse strand: Calculate genomic positions from exon end (reversed order)
                    overlapStartInGenomic = exonEnd - (overlapStartInTranscript - transcriptPosition);
                    overlapEndInGenomic = exonEnd - (overlapEndInTranscript - transcriptPosition);
                }

                // Format and add the genomic region
                if (overlapStartInGenomic == overlapEndInGenomic) {
                    genomicRegions.add(String.valueOf(overlapStartInGenomic));
                } else {
                    // On the reverse strand, ensure we output coordinates in decreasing order
                    if (strand == '+') {
                        genomicRegions.add(overlapStartInGenomic + "-" + overlapEndInGenomic);
                    } else {
                        genomicRegions.add(overlapEndInGenomic + "-" + overlapStartInGenomic);
                    }
                }
            }

            // Advance transcript position for the next exon
            transcriptPosition += exonLength;
        }

        if (read.isRw() && genomicRegions.size() > 1) {
            Collections.reverse(genomicRegions);
        }
        return genomicRegions;
    }


    public void mutateRead(Read read, double mutRate, StringBuilder sb) {
        String seq = read.getReadSeq();
        sb.setLength(0);
        sb.append(seq);
        int seqLength = seq.length();

        // Calculate the expected number of mutations
        int numMutations = (int) Math.round(mutRate * seqLength);

        // Randomly select unique positions to mutate
        HashSet<Integer> mutationPositions = new HashSet<>();
        while (mutationPositions.size() < numMutations) {
            mutationPositions.add(random.nextInt(seqLength));
        }

        // Perform mutations at the selected positions
        for (int pos : mutationPositions) {
            char originalNucleotide = seq.charAt(pos);
            char newNucleotide;

            // Choose a new nucleotide that is different from the original
            do {
                newNucleotide = NUCLEOTIDES[random.nextInt(NUCLEOTIDES.length)];
            } while (newNucleotide == originalNucleotide);

            sb.setCharAt(pos ,newNucleotide); // Mutate the character at the position
            read.addMutPos(pos);
        }
        read.setReadSeq(sb.toString());
    }
}
