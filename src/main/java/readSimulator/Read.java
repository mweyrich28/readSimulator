package readSimulator;

public class Read {
    private boolean isRw;
    private String readSeq;
    private int id;
    private int start;
    private int end;

    public Read(String seq, int start, int end, int id, boolean isRw) {
        this.isRw = isRw;
        this.readSeq = seq;
        this.id = id;
        this.start = start;
        this.end = end;
    }
}
