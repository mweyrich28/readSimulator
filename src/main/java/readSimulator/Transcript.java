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

    private final ArrayList<Read> fwReads;
    private final ArrayList<Read> rwReads;
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

    public void addExon(int start, int end, int pos) {
        Exon cds = new Exon(start, end, pos, end-start + 1);
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

    public ArrayList<Read> getFwReads() {
        return fwReads;
    }

    public ArrayList<Read> getRwReads() {
        return rwReads;
    }

    public void addReads(Read fw, Read rw) {
        this.fwReads.add(fw);
        this.rwReads.add(rw);
    }
}
