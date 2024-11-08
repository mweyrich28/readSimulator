package readSimulator;

import java.util.Objects;

public class Exon {
    private int length;
    private final int genomicStart;
    private final int genomicEnd;
    private int pos;

    public Exon(int start, int end, int pos, int length) {
        this.length = length;
        this.genomicStart = start;
        this.genomicEnd = end;
        this.pos = pos;
    }

    public int getGenomicStart() {
        return genomicStart;
    }

    public int getGenomicEnd() {
        return genomicEnd;
    }

    public int getPos() {
        return pos;
    }

    public void setPos(int pos) {
        this.pos = pos;
    }

    @Override
    public String toString() {
        return this.genomicStart + ":" + this.genomicEnd + " " + "[" + this.pos +"] " + this.length;
    }

    public int getLength() {
        return length;
    }
}
