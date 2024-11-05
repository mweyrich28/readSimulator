package readSimulator;

import java.util.Objects;

public class CodingDnaSequence {
    private final String id;
    private final int start;
    private final int end;
    private int pos;

    public CodingDnaSequence(String id, int start, int end, int pos) {
        this.id = id;
        this.start = start;
        this.end = end;
        this.pos = pos;
    }

    public String getId() {
        return id;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
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

        CodingDnaSequence cds = (CodingDnaSequence) o;
        return start == cds.start &&
                end == cds.end &&
                Objects.equals(id, cds.id);
    }

    @Override
    public int hashCode() {
        return Objects.hash(id, start, end);
    }

    @Override
    public String toString() {
        return this.id + " " + this.start + ":" + this.end + " " + "[" + this.pos +"]";
    }
}
