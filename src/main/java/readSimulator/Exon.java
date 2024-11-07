package readSimulator;

import java.util.Objects;

public class Exon {
    private final String id;
    private final int genomicStart;
    private final int genomicEnd;
    private int relStart;
    private int relEnd;
    private int pos;

    public Exon(String id, int start, int end, int pos) {
        this.id = id;
        this.genomicStart = start;
        this.genomicEnd = end;
        this.pos = pos;
    }

    public String getId() {
        return id;
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
    public boolean equals(Object o) {
        if (this == o) return true;  // Check if they are the same object
        if (o == null || getClass() != o.getClass()) return false;  // Check if the other object is a Cds instance

        Exon cds = (Exon) o;
        return genomicStart == cds.genomicStart &&
                genomicEnd == cds.genomicEnd &&
                Objects.equals(id, cds.id);
    }

    @Override
    public int hashCode() {
        return Objects.hash(id, genomicStart, genomicEnd);
    }

    @Override
    public String toString() {
        return this.id + " " + this.genomicStart + ":" + this.genomicEnd + " " + "[" + this.pos +"]";
    }

    public void setRelEnd(int relEnd) {
        this.relEnd = relEnd;
    }

    public int getRelStart() {
        return relStart;
    }

    public int getRelEnd() {
        return relEnd;
    }

    public void setRelStart(int relStart) {
        this.relStart = relStart;
    }
}
