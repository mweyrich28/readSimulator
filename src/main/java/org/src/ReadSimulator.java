package org.src;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.src.utils.FileUtils;
import org.src.utils.GenomeUtils;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class ReadSimulator {

    private static final char[] NUCLEOTIDES = {'A', 'T', 'C', 'G'};
    private static final byte[] BYTE_NUCLEOTIDES = {65, 84, 71, 67};

    private final SplittableRandom splittableRandom = new SplittableRandom();
    private final Genome genome;
    private final HashMap<String, HashMap<String, Integer>> readCounts = new HashMap<>();

    public ReadSimulator(int length, int frlength, int SD, double mutRate, String gtfPath, String readCountsPath, String fastaPath, String idxPath, String od) throws IOException {
        this.initReadCounts(readCountsPath);
        this.genome = new Genome(idxPath, fastaPath);
        this.genome.readGTF(gtfPath, readCounts);
        this.genome.initTargetGeneSeqs(readCounts); // check
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
        StringBuilder idBuilder = new StringBuilder();
        StringBuilder qidBuilder = new StringBuilder();

        // for each file init buff writer
        BufferedWriter summaryWriter = new BufferedWriter(new FileWriter(od + File.separator + "read.mappinginfo"));
        BufferedOutputStream fwOutputStream = new BufferedOutputStream(new FileOutputStream(od + File.separator + "fw.fastq"));
        BufferedOutputStream rwOutputStream = new BufferedOutputStream(new FileOutputStream(od + File.separator + "rw.fastq"));
        summaryWriter.write("readid\tchr\tgene\ttranscript\tt_fw_regvec\tt_rw_regvec\tfw_regvec\trw_regvec\tfw_mut\trw_mut");

        // pre format quality
        String quality = "I".repeat(length);
        byte[] qualityBytes = quality.getBytes();

        // allocate read arrays
        byte[] rwByteRead = new byte[length];
        byte[] fwByteRead = new byte[length];
        for (String geneKey : readCounts.keySet()) {

            Gene currGene = genome.getGenes().get(geneKey);

            for (String transcriptKey : readCounts.get(geneKey).keySet()) {

                Transcript transcript = this.genome.getGenes().get(geneKey).getTranscriptMap().get(transcriptKey);
                byte[] byteTranscriptSeq = transcript.getByteTranscriptSeq();
                int sampleAmount = readCounts.get(geneKey).get(transcriptKey);

                for (int i = 0; i < sampleAmount; i++) {

                    int fragmentLength;

                    do {
                        fragmentLength = (int) Math.round(normalDist.sample());
                    } while (fragmentLength < length || fragmentLength > byteTranscriptSeq.length);

                    int maxStartPos = byteTranscriptSeq.length - fragmentLength;
                    int randomStartPos = splittableRandom.nextInt(maxStartPos + 1);

                    // organizing coordinates
                    int fwStart = randomStartPos;
                    int fwEnd = randomStartPos + length - 1;
                    int rwStart = randomStartPos + fragmentLength - length;
                    int rwEnd = randomStartPos + fragmentLength - 1;


                    GenomeUtils.subsetArray(fwByteRead, byteTranscriptSeq, fwStart, fwEnd);

                    GenomeUtils.subsetArray(rwByteRead, byteTranscriptSeq, rwStart, rwEnd);
                    rwByteRead = GenomeUtils.revComplement(rwByteRead);


                    // TODO: I don't need these objects maybe
                    Read fwRead = new Read(fwByteRead, fwStart, fwEnd, readId, false);
                    Read rwRead = new Read(rwByteRead, rwStart, rwEnd, readId, true);

                    mutateByteRead(fwRead, mutRate);
                    mutateByteRead(rwRead, mutRate);

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

                    String seqId = idBuilder.append("@").append(readId).append("\n").toString();
                    String qId = qidBuilder.append("+").append(readId).append("\n").toString();

                    byte[] fwSeqData = GenomeUtils.buildFastqEntry(seqId, fwRead.getReadSeqBytes(), qId, qualityBytes, readId);
                    byte[] rwSeqData = GenomeUtils.buildFastqEntry(seqId, rwRead.getReadSeqBytes(), qId, qualityBytes, readId);

                    fwOutputStream.write(fwSeqData);
                    rwOutputStream.write(rwSeqData);

                    // reset sb
                    idBuilder.setLength(0);
                    qidBuilder.setLength(0);

                    readId++;

                    summaryWriter.write(summaryBuilder.toString());

                    // reset
                    summaryBuilder.setLength(0);
                }
            }
        }

        if (!summaryBuilder.isEmpty()) {
            summaryWriter.write(summaryBuilder.toString());
        }
        summaryWriter.flush();
        summaryWriter.close();
        rwOutputStream.flush();
        fwOutputStream.flush();
        rwOutputStream.close();
        fwOutputStream.close();
    }

    public ArrayList<String> getGenomicRegion(Read read, Transcript transcript, char strand) {
        ArrayList<Exon> exons = transcript.getExonList();
        ArrayList<String> genomicRegions = new ArrayList<>();

        int readStart = read.getStartInTranscript();
        int readEnd = read.getStopInTranscript();

        // tracks position in the transcript (1-based)
        int transcriptPosition = 0;

        for (int i = 0; i < exons.size(); i++) {
            Exon exon = exons.get(i);
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

    public void mutateRead(Read read, double mutRate, StringBuilder sb) {
        String seq = read.getReadSeq();
        sb.setLength(0);
        sb.append(seq);
        int seqLength = seq.length();

        // convert to prob
        double mutProb = mutRate / 100;

        // Check each position for mutation based on probability
        ArrayList<Integer> mutationPositions = new ArrayList<>();
        for (int i = 0; i < seqLength; i++) {
            if (splittableRandom.nextDouble() < mutProb) {
                mutationPositions.add(i);
            }
        }

        // mutate selected positions
        for (int pos : mutationPositions) {
            char originalNucleotide = seq.charAt(pos);
            char newNucleotide;
            do {
                newNucleotide = NUCLEOTIDES[splittableRandom.nextInt(NUCLEOTIDES.length)];
            } while (newNucleotide == originalNucleotide);
            sb.setCharAt(pos, newNucleotide);
            read.addMutPos(pos);
        }

        read.setReadSeq(sb.toString());
    }

    public void mutateByteRead(Read read, double mutRate) {
        byte[] seq = read.getReadSeqBytes();
        int seqLength = seq.length;

        // convert to prob
        double mutProb = mutRate / 100;
        // Check each position for mutation based on probability
        for (int i = 0; i < seqLength; i++) {
            if (splittableRandom.nextDouble() < mutProb) {
                byte originalNucleotide = seq[i];
                byte newNucleotide;
                do {
                    newNucleotide = BYTE_NUCLEOTIDES[splittableRandom.nextInt(BYTE_NUCLEOTIDES.length)];
                } while (newNucleotide == originalNucleotide);
                seq[i] = newNucleotide;
                read.addMutPos(i);
            }
        }
    }
}
