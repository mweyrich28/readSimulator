package readSimulator;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

import java.io.IOException;

public class Main {
    public static void main(String[] args) throws IOException {
        // add argparser
        ArgumentParser parser = ArgumentParsers.newFor("ExonSkipping").build().defaultHelp(true).description("Usage:\n\t-gtf <path-to-gtf>\n\t-out <path-to-out-tsv>");
        try {
            parser.addArgument("-gtf").required(true).help("Path to Gene Transfer Format File.");
            parser.addArgument("-od").required(true).help("Specify Output File Name.");
            parser.addArgument("-length").required(true).help("Specify Read Length.").type(Integer.class);;
            parser.addArgument("-frlength").required(true).help("Specify Mean Fragment Length.").type(Integer.class);;
            parser.addArgument("-SD").required(true).help("Specify Standard Error for Fragment Length.").type(Integer.class);;
            parser.addArgument("-mutationrate").required(true).help("Specify Mutationrate for Reads. 1.0 corresponds to 1%.").type(Double.class);;
            parser.addArgument("-fasta").required(true).help("Path to Fasta.");
            parser.addArgument("-fidx").required(true).help("Path to Fasta Index. This should correspond to your provided Fasta.");
            parser.addArgument("-readcounts").required(true).help("A TSV containing Entries of Transcripts to be simulated x Amount.");

            Namespace ns = parser.parseArgs(args);
            String gtfPath = ns.getString("gtf");
            String od = ns.getString("od");
            int length = ns.getInt("length");
            int frlength = ns.getInt("frlength");
            int SD = ns.getInt("SD");
            double mutRate = ns.getDouble("mutationrate");
            String fastaPath = ns.getString("fasta");
            String idxPath = ns.getString("fidx");
            String readCountsPath = ns.getString("readcounts");

            ReadSimulator r = new ReadSimulator(length , frlength , SD , mutRate , gtfPath, readCountsPath , fastaPath , idxPath , od);

        }
        // print usage entry if not all required args were provided
        catch (ArgumentParserException e) {
            parser.printHelp();
        }


        // int length = 75; // read length
        // int frlength = 20;
        // int SD = 80;
        // double mutRate = 1.0 / 100; // 0.01
        // String readCountsPath = "./inputFiles/simul.readcons";
        // String gtfPath = "./inputFiles/Homo_sapiens.GRCh37.75.gtf";
        // String fastaPath = "./inputFiles/Homo_sapiens.GRCh37.75.dna.toplevel.fa";
        // String idxPath = "./inputFiles/Homo_sapiens.GRCh37.75.dna.toplevel.fa.fai";
        // String od = "output";
        // ReadSimulator r = new ReadSimulator(length , frlength , SD , mutRate , gtfPath , readCountsPath , fastaPath , idxPath , od);

        // Genome genome = new Genome(idxPath, fastaPath);
        // String seq = genome.getGSE().getSequence("14",22907539, 22907546);
    }
}
