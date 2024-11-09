package org.src;

import java.util.ArrayList;
import java.util.Collections;

public class Transcript {
    private final String transcriptId;
    private final int start;
    private final int stop;
    private final ArrayList<Exon> exonList;
    private final String transcriptType;

    private String transcriptSeq; // patched together using its exons

    public Transcript(String transcriptId, String transcriptType, int start, int stop) {
        this.transcriptId = transcriptId;
        this.transcriptType = transcriptType;
        this.start = start;
        this.stop = stop;
        this.exonList = new ArrayList<>();
    }

    public void addExon(int start, int end, int pos) {
        Exon exon = new Exon(start, end, pos, end-start + 1);
        exonList.add(exon);
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
