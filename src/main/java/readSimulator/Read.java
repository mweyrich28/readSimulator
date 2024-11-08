package readSimulator;

import java.util.ArrayList;

public class Read {
    private boolean isRw;
    private char[] readSeq;
    private int id;
    private int startInFragment;
    private int stopInFragment;
    private String mutPos = "";

    public Read(String seq, int start, int end, int id, boolean isRw) {
        this.isRw = isRw;
        this.readSeq = initSeq(seq);
        this.id = id;
        this.startInFragment = start;
        this.stopInFragment = end;
    }

    public char[] initSeq(String seq) {
        char[] seqArr = new char[seq.length()];
        for (int i = 0; i < seq.length(); i++) {
           seqArr[i] = seq.charAt(i);
        }
        return seqArr;
    }

    public char[] getReadSeq() {
        return readSeq;
    }

    public void addMutPos(int pos) {
        this.mutPos += pos;
    }
}
