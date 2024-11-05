package readSimulator;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

import org.apache.commons.math3.distribution.NormalDistribution;
import readSimulator.utils.LineFilter;

import java.io.IOException;
import java.util.ArrayList;

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
    }
}
