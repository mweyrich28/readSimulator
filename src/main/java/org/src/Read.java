package org.src;

import java.util.ArrayList;

public class Read {
    private boolean isRw;
    private String readSeq;
    private byte[] readSeqBytes;
    private int id;
    private int startInTranscript;
    private int stopInTranscript;
    private int length;
    private ArrayList<String> mutPos = new ArrayList<>();

    public Read(String seq, int start, int end, int id, boolean isRw) {
        this.isRw = isRw;
        this.readSeq = seq;
        this.id = id;
        this.startInTranscript = start;
        this.stopInTranscript = end;
        this.length = stopInTranscript - startInTranscript;
    }

    public Read(byte[] seq, int start, int end, int id, boolean isRw) {
        this.isRw = isRw;
        this.readSeqBytes = seq;
        this.id = id;
        this.startInTranscript = start;
        this.stopInTranscript = end;
        this.length = stopInTranscript - startInTranscript;
    }
    public String getReadSeq() {
        return readSeq;
    }

    public void setReadSeq(String readSeq) {
        this.readSeq = readSeq;
    }

    public void addMutPos(int pos) {
        this.mutPos.add(Integer.toString(pos));
    }

    public ArrayList<String> getMutPos() {
        return mutPos;
    }

    public int getStartInTranscript() {
        return startInTranscript;
    }

    public int getStopInTranscript() {
        return stopInTranscript;
    }

    public int getLength() {
        return length;
    }

    public boolean isRw() {
        return isRw;
    }

    public byte[] getReadSeqBytes() {
        return readSeqBytes;
    }

    public void setReadSeqBytes(byte[] readSeqBytes) {
        this.readSeqBytes = readSeqBytes;
    }
}
