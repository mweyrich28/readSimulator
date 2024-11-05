package readSimulator;

import org.src.utils.FileUtils;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

public class Genome {
    private String name;
    private String version;
    private ArrayList<Gene> proteinCodingGenes;

    public Genome() {
        this.proteinCodingGenes = new ArrayList<>();
    }

    public Genome(String name, String version) {
        this.name = name;
        this.version = version;
        this.proteinCodingGenes = new ArrayList<>();
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getVersion() {
        return version;
    }

    public void setVersion(String version) {
        this.version = version;
    }

    public ArrayList<Gene> getProteinCodingGenes() {
        return proteinCodingGenes;
    }

    public ArrayList<String> generateESSE() {
        ArrayList<String> events = new ArrayList<>();
        StringBuilder sb = new StringBuilder();
        sb.append("id\t").append("symbol\t").append("chr\t").append("strand\t").append("nprots\t").append("ntrans\t").append("SV\t").append("WT\t").append("SV_prots\t").append("WT_prots\t").append("min_skipped_exon\t").append("max_skipped_exon\t").append("min_skipped_bases\t").append("max_skipped_bases");

        // add colNames
        events.add(sb.toString());

        for (Gene gene : proteinCodingGenes) {
            if (gene.getStrand() == '-') {
               gene.invertTranscripts();
            }
            gene.generateIntrons();
            events.addAll(gene.getEvents());
        }

        return events;
    }

    // TODO: improve for memory
    public String parseAttributes(String[] attributeEntries, String attributeName) {
        for (int i = 0; i < attributeEntries.length; i++) {
            String trimmedEntry = attributeEntries[i].trim();

            int posSpace = trimmedEntry.indexOf(' ');
            String attributeKey = trimmedEntry.substring(0, posSpace);

            String attributeVal;

            if (trimmedEntry.endsWith("\"")) {
                attributeVal = trimmedEntry.substring(posSpace + 2, trimmedEntry.length() - 1);
            }
            // there are entries like exon_number, which are not encapsulated in quotes
            else {
                attributeVal = trimmedEntry.substring(posSpace + 1, trimmedEntry.length() - 1);
            }

            if (attributeKey.equals(attributeName)) {
                return attributeVal;
            }
        }
        return null;
    }

    public void readGTFCDS(String pathToGtf) throws IOException {
        // get lines of gtf
        ArrayList<String> lines = FileUtils.readExonLines(new File(pathToGtf));

        // sanity check vars
        Gene lastGeneId = null;
        String lastTranscriptId = null;
        int cdsCounter = 0;

        for (int i = 0; i < lines.size() - 1; i++) {
            String currLine = lines.get(i);
            // skip potential comments
            if (currLine.startsWith("#")) {
                continue;
            }

            // extract main components (line split by \t)
            String[] mainComponents = currLine.split("\t");
            // split attributes again at ";"
            String[] attributeEntries = mainComponents[mainComponents.length - 1].split(";");

            // get newGeneId of current line
            String newGeneId = parseAttributes(attributeEntries, "gene_id");

            // check if we hit a new gene
            if (lastGeneId == null || !newGeneId.equals(lastGeneId.getGeneId())) {
                // update gene and continue with next gtf line
                int geneStart = Integer.parseInt(mainComponents[3]);
                int geneEnd = Integer.parseInt(mainComponents[4]);
                String geneName = parseAttributes(attributeEntries, "gene_name");
                String chr = mainComponents[0];
                char strand = mainComponents[6].charAt(0);
                lastGeneId = new Gene(newGeneId, geneStart, geneEnd, geneName, chr, strand);

                continue;
            }

            // only add cds to current transcript
            String transcriptId = parseAttributes(attributeEntries,"transcript_id");
            if (mainComponents[2].equals("CDS")) {
                String cdsIdKey = "protein_id";
                String cdsId = parseAttributes(attributeEntries, cdsIdKey);
                if (cdsId == null) {
                    // fall back to ccds id
                    cdsId = parseAttributes(attributeEntries, "ccdsid");
                }
                if (cdsId == null) {
                    // fall back to empty id
                    cdsId = "NaN";
                }
                // check if we are in a new transcript
                if(lastGeneId.getTranscripts().isEmpty()) { // if gene transcripts are empty, just add new transcript

                    // add gene to genome (based on if it is p coding or not)
                    this.proteinCodingGenes.add(lastGeneId);

                    cdsCounter = 0; // reset cdsCounter

                    // add new transcript to current gene
                    Transcript transcript = new Transcript(transcriptId, mainComponents[2]);
                    lastGeneId.addTranscript(transcript);

                    // add cds to current transcript
                    lastGeneId.getLastTranscript().addCds(
                            cdsId,
                            Integer.parseInt(mainComponents[3]),
                            Integer.parseInt(mainComponents[4]),
                            cdsCounter
                    );
                    cdsCounter++;
                }
                // else check if we are still in the same transcript
                else if (transcriptId.equals(lastGeneId.getLastTranscript().getTranscriptId())) {
                    lastGeneId.getLastTranscript().addCds(
                            cdsId,
                            Integer.parseInt(mainComponents[3]),
                            Integer.parseInt(mainComponents[4]),
                            cdsCounter
                    );
                    cdsCounter++;
                }
                // else we add a new transcript
                else {
                    cdsCounter = 0; // reset cdsCounter
                    // add new transcript to current gene
                    lastGeneId.addTranscript(new Transcript(transcriptId, mainComponents[2]));
                    // add cds to current transcript
                    lastGeneId.getLastTranscript().addCds(
                            cdsId,
                            Integer.parseInt(mainComponents[3]),
                            Integer.parseInt(mainComponents[4]),
                            cdsCounter
                    );
                    cdsCounter++;
                }
            }

            // here we increment the ntrans count, we currently don't need to save exons.
            // by only having a count, we save space
            if (mainComponents[2].equals("CDS") || mainComponents[2].equals("exon")) {
                if (!(transcriptId.equals(lastTranscriptId))) {
                    lastGeneId.incTrans();
                    lastTranscriptId = transcriptId;
                }
            }
        }
    }
}
