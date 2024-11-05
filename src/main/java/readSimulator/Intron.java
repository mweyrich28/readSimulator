package readSimulator;

public class Intron {
    private final int start;
    private final int end;

    public Intron(int start, int end) {
        this.start = start;
        this.end = end;
    }

    public int getEnd() {
        return end;
    }

    public int getStart() {
        return start;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;  // Same object reference
        if (obj == null || getClass() != obj.getClass()) return false;

        Intron intron = (Intron) obj;
        return start == intron.start && end == intron.end;
    }

    @Override
    public int hashCode() {
        return 31 * start + end;
    }
}
