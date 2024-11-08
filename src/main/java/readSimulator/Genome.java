package readSimulator;

import readSimulator.utils.FileUtils;
import readSimulator.utils.GenomeUtils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class Genome {
    private final HashMap<String, Gene> genes;
    private GenomeSequenceExtractor GSE;

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

    public void initTargetGeneSeqs(HashMap<String , HashMap<String, Integer>> readCounts) throws IOException {
        for (String geneKey : readCounts.keySet()) {
            Gene gene = this.genes.get(geneKey);
            String chr = gene.getChr();
            int start = gene.getStart();
            int end = gene.getEnd();
            String seq = GSE.getSequence(chr, start, end);
            gene.setSequence(seq);

            // for all transcripts of gene generate exons and trans seq
            for (String transcriptKey : readCounts.get(geneKey).keySet()) {
                Transcript transcript = gene.getTranscriptMap().get(transcriptKey);
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



    public void readGTF(String pathToGtf, HashMap<String, HashMap<String, Integer>> readCounts) throws IOException {
        // get lines of gtf

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

            // check if we hit a new gene
            if (mainComponents[2].equals("gene")) {
                // update gene and continue with next gtf line
                int geneStart = Integer.parseInt(mainComponents[3]);
                int geneEnd = Integer.parseInt(mainComponents[4]);
                String geneName = FileUtils.parseGTFAttributes(attributeEntries, "gene_name");
                String chr = mainComponents[0];
                char strand = mainComponents[6].charAt(0);
                lastGene = new Gene(newGeneId, geneStart, geneEnd, geneName, chr, strand);
                genes.put(lastGene.getGeneId(), lastGene);

                continue;
            }


            // did we hit a new transcript
            if (mainComponents[2].equals("transcript")) {

                // only add cds to current transcript
                String transcriptId = FileUtils.parseGTFAttributes(attributeEntries, "transcript_id");

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
            // if we don't have "transcript" in mainComp[2], we are either in CDS
            // or exon of last transcript
            else {
                // add exon to last transcript
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

    public GenomeSequenceExtractor getGSE() {
        return GSE;
    }
}
