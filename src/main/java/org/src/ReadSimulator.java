package org.src;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.src.utils.FileUtils;
import org.src.utils.GenomeUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class ReadSimulator {

    private static final char[] NUCLEOTIDES = {'A', 'T', 'C', 'G'};
    private static boolean debug = false;
    private final SplittableRandom splittableRandom = new SplittableRandom();
    private final Genome genome;
    private final HashMap<String, HashMap<String, Integer>> readCounts = new HashMap<>();

    public ReadSimulator(int length, int frlength, int SD, double mutRate, String gtfPath, String readCountsPath, String fastaPath, String idxPath, String od, boolean debug, String transcriptomePath) throws IOException {
        ReadSimulator.debug = debug;
        this.initReadCounts(readCountsPath);
        this.genome = new Genome(idxPath, fastaPath);
        this.genome.readGTF(gtfPath, readCounts); // O(n) n = lines
        this.genome.initTargetGeneSeqs(readCounts);
        if (debug && transcriptomePath != null) {
            genome.validateTranscripts("/home/malte/projects/gobi/readSimulator/inputFiles/Homo_sapiens.GRCh37.75.cdna.all.fa", readCounts);
        }
        this.generateReads(length, frlength, SD, mutRate, od);
    }

    public void initReadCounts(String readCountsPath) throws IOException {
        ArrayList<String> lines = FileUtils.readLines(new File(readCountsPath));

        for (String line : lines) {
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

    public void generateReads(int length, int frlength, int SD, double mutRate, String od) throws IOException {

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

        // for each file init buff writer
        BufferedWriter summaryWriter = new BufferedWriter(new FileWriter(od + File.separator + "read.mappinginfo"));
        BufferedWriter fwFastqWriter = new BufferedWriter(new FileWriter(od + File.separator + "fw.fastq"));
        BufferedWriter rwFastqWriter = new BufferedWriter(new FileWriter(od + File.separator + "rw.fastq"));
        summaryWriter.write("readid\tchr\tgene\ttranscript\tt_fw_regvec\tt_rw_regvec\tfw_regvec\trw_regvec\tfw_mut\trw_mut");


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
                    int randomStartPos = splittableRandom.nextInt(maxStartPos + 1);

                    String fwSeqRead = transcriptSeq.substring(randomStartPos, randomStartPos + length);
                    String rwSeqRead = GenomeUtils.revComplement(transcriptSeq.substring(randomStartPos + fragmentLength - length, randomStartPos + fragmentLength));

                    // organizing coordinates
                    int fwStart = randomStartPos;
                    int fwEnd = randomStartPos + length - 1;
                    int rwStart = randomStartPos + fragmentLength - length;
                    int rwEnd = randomStartPos + fragmentLength - 1;

                    // TODO: I don't need these objects maybe
                    Read fwRead = new Read(fwSeqRead, fwStart, fwEnd, readId, false);
                    Read rwRead = new Read(rwSeqRead, rwStart, rwEnd, readId, true);

                    if (debug) {
                        ArrayList<String> reg = getGenomicRegion(fwRead, transcript, currGene.getStrand());
                        boolean correctFw = fwRead.getReadSeq().equals(this.genome.getGSE().extractRegion(reg, currGene.getChr(), currGene.getStrand(), false));

                        reg = getGenomicRegion(rwRead, transcript, currGene.getStrand());
                        boolean correctRw = rwRead.getReadSeq().equals(this.genome.getGSE().extractRegion(reg, currGene.getChr(), currGene.getStrand(), true));

                        if (!(correctRw && correctFw)) {
                            System.out.println("DEBUG:");
                            reg = getGenomicRegion(fwRead, transcript, currGene.getStrand());
                            System.out.println("FWR " + reg);
                            System.out.println("FWS " + fwSeqRead);
                            System.out.println("EXT " + this.genome.getGSE().extractRegion(reg, currGene.getChr(), currGene.getStrand(), false));

                            reg = getGenomicRegion(rwRead, transcript, currGene.getStrand());
                            System.out.println("RWR " + reg);
                            System.out.println("RWS " + rwSeqRead);
                            System.out.println("EXT " + this.genome.getGSE().extractRegion(reg, currGene.getChr(), currGene.getStrand(), true));
                            throw new RuntimeException("Read " + readId + " does not match with its GenomicReadVector.");
                        }
                    }

                    mutateRead(fwRead, mutRate);
                    mutateRead(rwRead, mutRate);

                    ArrayList<String> fwRegVec = getGenomicRegion(fwRead, transcript, currGene.getStrand());
                    ArrayList<String> rwRegVec = getGenomicRegion(rwRead, transcript, currGene.getStrand());


                    // concat summary col
                    summaryBuilder.append("\n").append(readId).append("\t")
                            .append(currGene.getChr()).append("\t")
                            .append(currGene.getGeneId()).append("\t")
                            .append(transcript.getTranscriptId()).append("\t")
                            .append(fwStart).append("-").append(fwEnd + 1).append("\t") // + 1 because 1 based
                            .append(rwStart).append("-").append(rwEnd + 1).append("\t") // + 1 because 1 based
                            .append(String.join("|", fwRegVec)).append("\t")
                            .append(String.join("|", rwRegVec)).append("\t")
                            .append(String.join(",", fwRead.getMutPos())).append("\t")
                            .append(String.join(",", rwRead.getMutPos()));

                    // format read seqs in fastq files
                    String quality = "I".repeat(fwSeqRead.length());
                    // avoid line breaks at end of file
                    if (readId == 0) {
                        fwFastqBuilder.append("@").append(readId).append("\n")
                                .append(fwRead.getReadSeq()).append("\n")
                                .append("+").append(readId).append("\n")
                                .append(quality);
                        rwFastqBuilder.append("@").append(readId).append("\n")
                                .append(rwRead.getReadSeq()).append("\n")
                                .append("+").append(readId).append("\n")
                                .append(quality);
                    } else {
                        fwFastqBuilder.append("\n")
                                .append("@").append(readId).append("\n")
                                .append(fwRead.getReadSeq()).append("\n").append("+").append(readId).append("\n")
                                .append(quality);
                        rwFastqBuilder.append("\n")
                                .append("@").append(readId).append("\n")
                                .append(rwRead.getReadSeq()).append("\n").append("+").append(readId).append("\n")
                                .append(quality);

                    }

                    readId++;

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

        summaryWriter.close();
        fwFastqWriter.close();
        rwFastqWriter.close();
    }

    public ArrayList<String> getGenomicRegion(Read read, Transcript transcript, char strand) {
        ArrayList<Exon> exons = transcript.getExonList();
        ArrayList<String> genomicRegions = new ArrayList<>();

        int readStart = read.getStartInTranscript();
        int readEnd = read.getStopInTranscript();

        // tracks position in the transcript (1-based)
        int transcriptPosition = 0;

        for (Exon exon : exons) {
            int exonStart = exon.getGenomicStart();
            int exonEnd = exon.getGenomicEnd();
            int exonLength = exonEnd - exonStart + 1;

            // check if the read overlaps with this exon in transcript coordinates
            if (transcriptPosition + exonLength - 1 >= readStart && transcriptPosition <= readEnd) {
                // calculate overlap within transcript coordinates
                int overlapStartInTranscript = Math.max(readStart, transcriptPosition);
                int overlapEndInTranscript = Math.min(readEnd, transcriptPosition + exonLength - 1);

                // map overlap to genomic coordinates based on strand
                int overlapStartInGenomic;
                int overlapEndInGenomic;

                if (strand == '+') {
                    // +: calculate genomic positions from exon start
                    overlapStartInGenomic = exonStart + (overlapStartInTranscript - transcriptPosition);
                    overlapEndInGenomic = exonStart + (overlapEndInTranscript - transcriptPosition);
                } else {
                    // -: calculate genomic positions from exon end (reversed order)
                    overlapStartInGenomic = exonEnd - (overlapStartInTranscript - transcriptPosition);
                    overlapEndInGenomic = exonEnd - (overlapEndInTranscript - transcriptPosition);
                }

                // add genomic region
                if (overlapStartInGenomic == overlapEndInGenomic) {
                    genomicRegions.add(overlapStartInGenomic + "-" + (overlapStartInGenomic + 1));
                } else {
                    if (strand == '+') {
                        genomicRegions.add(overlapStartInGenomic + "-" + (overlapEndInGenomic + 1));
                    } else {
                        // reverse strand:coordinates need to be in decreasing order
                        genomicRegions.add(overlapEndInGenomic + "-" + (overlapStartInGenomic + 1));
                    }
                }
            }

            // walk exon dist in transcript
            transcriptPosition += exonLength;
        }

        if (strand == '-' && genomicRegions.size() > 1) {
            // since exons are stored in rev order for transcript of rev strand
            // we have to reverse our order if length > 1
            Collections.reverse(genomicRegions);
        }

        return genomicRegions;
    }

    public void mutateRead(Read read, double mutRate) {
        String seq = read.getReadSeq();
        StringBuilder sb = new StringBuilder(seq);
        int seqLength = seq.length();

        // convert to prob
        double mutProb = mutRate / 100;

        // Check each position for mutation based on probability
        for (int i = 0; i < seqLength; i++) {
            if (splittableRandom.nextDouble() < mutProb) {
                char originalNucleotide = seq.charAt(i);
                char newNucleotide;
                do {
                    newNucleotide = NUCLEOTIDES[splittableRandom.nextInt(NUCLEOTIDES.length)];
                } while (newNucleotide == originalNucleotide);
                sb.setCharAt(i, newNucleotide);
                read.addMutPos(i);
            }
        }
        read.setReadSeq(sb.toString());
    }
}
