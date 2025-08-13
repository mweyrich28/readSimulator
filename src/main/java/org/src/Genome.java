package org.src;

import org.src.utils.FileUtils;
import org.src.utils.GenomeUtils;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.SplittableRandom;

public class Genome {
    private final HashMap<String, Gene> genes;
    private GenomeSequenceExtractor GSE;
    private final SplittableRandom splittableRandom = new SplittableRandom(42);
    private static final char[] NUCLEOTIDES = {'A', 'T', 'C', 'G'};

    public Genome() {
        this.genes = new HashMap<>();
    }

    public Genome(String fastaIdx, String fasta) throws IOException {
        this.genes = new HashMap<>();
        this.GSE = new GenomeSequenceExtractor(fastaIdx, fasta);

    }

    public HashMap<String, Gene> getGenes() {
        return genes;
    }

    public void initTargetGeneSeqs(HashMap<String , HashMap<String, Integer>> readCounts, double mutRate) throws IOException {
        for (String geneKey : readCounts.keySet()) {
            Gene gene = this.genes.get(geneKey);
            if (gene == null) {
                continue;
            }
            String chr = gene.getChr();
            int start = gene.getStart();
            int end = gene.getEnd();
            String seq = GSE.getSequence(chr, start, end);
            if (seq == null) {
                continue;
            }

            gene.setSequence(seq);
            if (mutRate != 0.0) {
                simulateMutations(gene, mutRate);
            }

            // for all transcripts of gene generate exons and trans seq
            for (String transcriptKey : readCounts.get(geneKey).keySet()) {
                Transcript transcript = gene.getTranscriptMap().get(transcriptKey);
                if (transcript == null){
                    continue;
                }
                StringBuilder transcriptSeq = new StringBuilder();

                for (int i = 0; i < transcript.getExonList().size(); i++) {
                    Exon exon = null;
                    if (gene.getStrand() == '-') {
                        exon = transcript.getExonList().get(transcript.getExonList().size() - i - 1);
                    } else {
                        exon = transcript.getExonList().get(i);
                    }

                    int relStart = exon.getGenomicStart() - gene.getStart();
                    int relEnd = exon.getGenomicEnd() - gene.getStart() + 1;

                    transcriptSeq.append(gene.getSeq(), relStart, relEnd);
                }
                if (gene.getStrand() == '-') {
                    transcript.setTranscriptSeq(GenomeUtils.revComplement(transcriptSeq.toString()));
                } else {
                    transcript.setTranscriptSeq(transcriptSeq.toString());
                }
            }
        }
    }

    public void printSeq(String geneId) {
        System.out.println(this.genes.get(geneId).getSeq());
    }

    public void readGTF(String pathToGtf, HashMap<String, HashMap<String, Integer>> readCounts) throws IOException {
        // sanity check vars
        Gene lastGene = null;
        int exonCounter = 0;

        BufferedReader buff = new BufferedReader(new FileReader(pathToGtf));
        String line;

        while((line = buff.readLine()) != null) {
            // skipp all lines that don t contain a relevant gene id
            if (!FileUtils.filterLine(line, readCounts)) {
               continue;
            }

            // extract main components (line split by \t)
            String[] mainComponents = line.split("\t");
            // split attributes again at ";"
            String[] attributeEntries = mainComponents[mainComponents.length - 1].split(";");

            // get newGeneId of current line
            String newGeneId = FileUtils.parseGTFAttributes(attributeEntries, "gene_id");
            if (newGeneId.length() != 15) {
                // gene name contains a version at suffix
                newGeneId= newGeneId.split("\\.")[0];
            }

            // check if we hit a new gene
            if (mainComponents[2].equals("gene")) {
                // update gene and continue with next gtf line
                int geneStart = Integer.parseInt(mainComponents[3]);
                int geneEnd = Integer.parseInt(mainComponents[4]);
//                String geneName = FileUtils.parseGTFAttributes(attributeEntries, "gene_name");
                String geneName = FileUtils.parseGTFAttributes(attributeEntries, "gene_id");
                if (geneName.length() != 15) {
                    // gene name contains a version at suffix
                    geneName = geneName.split("\\.")[0];
                }
                String chr = mainComponents[0].replace("chr", "");
                String annot = FileUtils.parseGTFAttributes(attributeEntries, "gene_type");
                char strand = mainComponents[6].charAt(0);
                lastGene = new Gene(newGeneId, geneStart, geneEnd, geneName, chr, strand, annot);
                genes.put(lastGene.getGeneId(), lastGene);

                continue;
            }


            // did we hit a new transcript
            if (mainComponents[2].equals("transcript")) {

                // only add cds to current transcript
                String transcriptId = FileUtils.parseGTFAttributes(attributeEntries, "transcript_id");
                if (transcriptId.length() != 15) {
                    // gene name contains a version at suffix
                    transcriptId= transcriptId.split("\\.")[0];
                }

                // add gene to genome
                genes.put(lastGene.getGeneId(), lastGene);

                // add new transcript to current gene
                int transcriptStart = Integer.parseInt(mainComponents[3]);
                int transcriptStop = Integer.parseInt(mainComponents[4]);
                Transcript transcript = new Transcript(transcriptId, mainComponents[2], transcriptStart, transcriptStop);
                lastGene.addTranscript(transcript);

                // reset
                exonCounter = 0;

            }
            // add exon to last transcript
            else if (mainComponents[2].equals("exon")) {
                int start = Integer.parseInt(mainComponents[3]);
                int end = Integer.parseInt(mainComponents[4]);
                lastGene.getLastTranscript().addExon(
                        start,
                        end,
                        exonCounter
                );
                exonCounter++;
            }
        }
    }

    public void validateTranscripts(String transcriptomePath, HashMap<String, HashMap<String, Integer>> readCounts) throws IOException {
       HashMap<String, StringBuilder> transcriptome = new HashMap<>();

        BufferedReader buff = new BufferedReader(new FileReader(transcriptomePath));
        String line;
        // init relevant transcript seqs
        String transcriptID = null;
        String geneID = null;
        boolean relevant = false;
        while((line = buff.readLine()) != null) {
           if(line.startsWith(">")) {
               String[] components = line.split(" ");
               transcriptID = components[0].substring(1);
               geneID = components[3].substring(5);
               if (readCounts.containsKey(geneID) && readCounts.get(geneID).containsKey(transcriptID)) {
                   transcriptome.put(transcriptID, new StringBuilder());
                   relevant = true;
                   continue;
               } else {
                   relevant = false;
                   continue;
               }
           }

           if (relevant) {
               transcriptome.get(transcriptID).append(line.trim());
           }
        }

        int count = 0;
        for (String geneId: readCounts.keySet()) {
            for (String transcriptId : readCounts.get(geneId).keySet()) {
                String expected = transcriptome.get(transcriptId).toString();
                String actual = this.genes.get(geneId).getTranscriptMap().get(transcriptId).getTranscriptSeq();
                count++;
                if (!actual.equals(expected)) {
                    throw new RuntimeException("Transcript " + transcriptId + " of gene " + geneId + " does not correspond to reference in transcriptome.\nREF: " +
                            expected + "\nGOT: " + actual );
                }
            }
        }
        System.out.println("DEBUG: ALL " + count + " TRANSCRIPTS HAVE BEEN SUCCESSFULLY VALIDATED VIA PROVIDED TRANSCRIPTOME.");
    }

    public GenomeSequenceExtractor getGSE() {
        return GSE;
    }

    public void simulateMutations(Gene gene, double mutRate) {
        String seq = gene.getSeq();
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
                gene.addMutation(i);
            }
        }
        gene.setSequence(sb.toString());
    }
}
