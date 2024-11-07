package readSimulator;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class Transcript {
    private final String transcriptId;
    private final int start;
    private final int stop;
    private final ArrayList<Exon> exonList;
    private final String transcriptType;

    private final ArrayList<String> fwReads;
    private final ArrayList<String> rwReads;
    private String transcriptSeq; // patched together using its exons

    public Transcript(String transcriptId, String transcriptType, int start, int stop) {
        this.transcriptId = transcriptId;
        this.transcriptType = transcriptType;
        this.start = start;
        this.stop = stop;
        this.exonList = new ArrayList<>();
        this.fwReads = new ArrayList<>();
        this.rwReads = new ArrayList<>();
    }

    public void addExon(String exonid, int start, int end, int pos) {
        Exon cds = new Exon(exonid, start, end, pos);
        exonList.add(cds);
    }

    public String getTranscriptId() {
        return transcriptId;
    }

    public ArrayList<Exon> getExonList() {
        return this.exonList;
    }

    public void reversCdsList() {
        Collections.reverse(this.exonList);
    }

    public void setTranscriptSeq(String transcriptSeq) {
        this.transcriptSeq = transcriptSeq;
    }

    public String getTranscriptSeq() {
        return this.transcriptSeq;
    }
}
