package readSimulator;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

public class Transcript {
    private final String transcriptId;
    private final HashMap<Integer, CodingDnaSequence> cdsEndIndices;
    private final HashMap<Integer, CodingDnaSequence> cdsStartIndices;
    private final HashMap<String , CodingDnaSequence> cdsIdMap;
    private final ArrayList<CodingDnaSequence> cdsList;
    private final String transcriptType;

    public Transcript(String transcriptId, String transcriptType) {
        this.transcriptId = transcriptId;
        this.transcriptType = transcriptType;
        this.cdsEndIndices = new HashMap<>();
        this.cdsStartIndices = new HashMap<>();
        this.cdsIdMap = new HashMap<>();
        this.cdsList = new ArrayList<>();
    }

    public void addCds(String cdsid, int start, int end, int pos) {
        CodingDnaSequence cds = new CodingDnaSequence(cdsid, start, end, pos);

        // easily check if transcript has a cds ending at I_s / starting at I_e
        cdsEndIndices.put(end, cds);
        cdsStartIndices.put(start, cds);
        cdsIdMap.put(cds.getId(), cds);
        cdsList.add(cds);
    }

    public String getTranscriptId() {
        return transcriptId;
    }

    public ArrayList<CodingDnaSequence> getCdsList() {
        return this.cdsList;
    }

    public HashMap<Integer, CodingDnaSequence> getCdsEndIndices() {
        return cdsEndIndices;
    }

    public HashMap<Integer, CodingDnaSequence> getCdsStartIndices() {
        return cdsStartIndices;
    }

    public void reversCdsList() {
        Collections.reverse(this.cdsList);
    }

    public HashMap<String, CodingDnaSequence> getCdsIdMap() {
        return cdsIdMap;
    }
}
