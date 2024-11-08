package readSimulator;

import java.io.IOException;

public class Main {
    public static void main(String[] args) throws IOException {
        // add argparser
        // ArgumentParser parser = ArgumentParsers.newFor("ExonSkipping").build().defaultHelp(true).description("Usage:\n\t-gtf <path-to-gtf>\n\t-out <path-to-out-tsv>");
        // try {
            // parser.addArgument("-gtf").required(true).help("Path to Gene Transfer Format File.");
            // parser.addArgument("-o").required(true).help("Specify Output File Name.");
            // Namespace ns = parser.parseArgs(args);
            // String gtfPath = ns.getString("gtf");
            // String outPath = ns.getString("o");

        // }
        // print usage entry if not all required args were provided
        // catch (ArgumentParserException e) {
        //     parser.printHelp();
        // }
        int length = 75; // read length
        int frlength = 20;
        int SD = 80;
        double mutRate = 1.0 / 100; // 0.01
        // String gtfPath = "./inputFiles/test.GRCh37.75.gtf";
        // String readCountsPath = "./inputFiles/simul.readcons";
        String gtfPath = "./inputFiles/Homo_sapiens.GRCh37.75.gtf";
        String readCountsPath = "inputFiles/readcounts.simulation";
        String fastaPath = "./inputFiles/Homo_sapiens.GRCh37.75.dna.toplevel.fa";
        String idxPath = "./inputFiles/Homo_sapiens.GRCh37.75.dna.toplevel.fa.fai";
        String od = "output";
        ReadSimulator r = new ReadSimulator(length , frlength , SD , mutRate , gtfPath , readCountsPath , fastaPath , idxPath , od);
        // Genome genome = new Genome(idxPath, fastaPath);
        // String seq = genome.getGSE().getSequence("14",22907539, 22907546);
    }
}
