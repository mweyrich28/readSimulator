package readSimulator;

public class Read {
    private boolean isRw;
    private String readSeq;
    private int id;
    private int startInFragment;
    private int stopInFragment;

    public Read(String seq, int start, int end, int id, boolean isRw) {
        this.isRw = isRw;
        this.readSeq = seq;
        this.id = id;
        this.startInFragment = start;
        this.stopInFragment = end;
    }

    public String getReadSeq() {
        return readSeq;
    }

    public void setReadSeq(String readSeq) {
        this.readSeq = readSeq;
    }
}
