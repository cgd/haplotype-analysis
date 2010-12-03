/*
 * Copyright (c) 2010 The Jackson Laboratory
 * 
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software.  If not, see <http://www.gnu.org/licenses/>.
 */
package org.jax.haplotype.analysis;

import java.io.BufferedOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jax.geneticutil.data.BasePairInterval;
import org.jax.geneticutil.data.PartitionedIntervalSet;
import org.jax.haplotype.analysis.jaxbfactory.JaxbHaplotypeAssociationTestFactory;
import org.jax.haplotype.analysis.jaxbfactory.JaxbPhylogenyAssociationTestFactory;
import org.jax.haplotype.jaxbgenerated.CommaSeparatedHaplotypeAssociationTestOutputType;
import org.jax.haplotype.jaxbgenerated.HaplotypeAssociationExperimentDesign;
import org.jax.haplotype.jaxbgenerated.HaplotypeAssociationTestOutputType;
import org.jax.haplotype.jaxbgenerated.HaplotypeAssociationTestType;
import org.jax.haplotype.jaxbgenerated.PhylogenyAssociationTestOutputType;
import org.jax.haplotype.jaxbgenerated.PhylogenyAssociationTestType;
import org.jax.haplotype.jaxbgenerated.TabDelimitedPhylogenyAssociationTestOutputType;
import org.jax.haplotype.phylogeny.data.PhylogenyTestResult;
import org.jax.haplotype.phylogeny.data.PhylogenyTreeNode;
import org.jax.util.datastructure.SetUtilities;
import org.jax.util.io.CharacterDelimitedParser;

/**
 * Main command line entry point for doing HAM analysis
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class HaplotypeAssociationMappingMain
{
    /**
     * our logger
     */
    private static final Logger LOG = Logger.getLogger(
            HaplotypeAssociationMappingMain.class.getName());
    
    private String analysisDesignFile = null;
    
    /**
     * Constructor
     */
    public HaplotypeAssociationMappingMain()
    {
    }
    
    /**
     * Getter for the design file to perform analysis on
     * @return the analysisDesignFile
     */
    public String getAnalysisDesignFile()
    {
        return this.analysisDesignFile;
    }
    
    /**
     * Setter for the design file to perform analysis on
     * @param analysisDesignFile the analysisDesignFile to set
     */
    public void setAnalysisDesignFile(String analysisDesignFile)
    {
        this.analysisDesignFile = analysisDesignFile;
    }
    
    /**
     * Perform the HAM tests
     */
    private void performAnalysis()
    {
        try
        {
            JAXBContext jaxbContext = JAXBContext.newInstance(
                    HaplotypeAssociationExperimentDesign.class.getPackage().getName());
            Unmarshaller unmarshaller = jaxbContext.createUnmarshaller();
            
            HaplotypeAssociationExperimentDesign haplotypeAssociationExperimentDesign =
                (HaplotypeAssociationExperimentDesign)unmarshaller.unmarshal(
                    new FileInputStream(this.getAnalysisDesignFile()));
            
            // take care of writing our haplotype test results to the outputs
            List<HaplotypeAssociationTestOutputType> haplotypeOutputList =
                haplotypeAssociationExperimentDesign.getHaplotypeAssociationTestOutput();
            this.generateHaplotypeAssociationTestOutputs(haplotypeOutputList);
            
            // take care of writing the phylogeny test results to the outputs
            List<PhylogenyAssociationTestOutputType> phylogenyTestOutputs =
                haplotypeAssociationExperimentDesign.getPhylogenyAssociationTestOutput();
            this.generatePhylogenyAssociationTestOutputs(phylogenyTestOutputs);
        }
        catch(JAXBException ex)
        {
            LOG.log(Level.SEVERE,
                    "failed to parse experiment design file",
                    ex);
        }
        catch(IOException ex)
        {
            LOG.log(Level.SEVERE,
                    "file IO error",
                    ex);
        }
    }
    
    /**
     * Take care of writing the phylogeny test results data sources to the given
     * outputs
     * @param phylogenyTestOutputs
     *          the output list
     */
    private void generatePhylogenyAssociationTestOutputs(
            List<PhylogenyAssociationTestOutputType> phylogenyTestOutputs)
    {
        try
        {
            Map<String, Map<Integer, List<PhylogenyTestResult>>> phylogenyTestResultsCache =
                new HashMap<String, Map<Integer, List<PhylogenyTestResult>>>();
            
            for(PhylogenyAssociationTestOutputType currPhyloTestOut: phylogenyTestOutputs)
            {
                if(currPhyloTestOut instanceof TabDelimitedPhylogenyAssociationTestOutputType)
                {
                    TabDelimitedPhylogenyAssociationTestOutputType currTabDelimPhyloTestOut =
                        (TabDelimitedPhylogenyAssociationTestOutputType)currPhyloTestOut;
                    
                    PrintStream fileOut = new PrintStream(
                            new BufferedOutputStream(new FileOutputStream(
                                    currTabDelimPhyloTestOut.getFileLocation())));
                    
                    PhylogenyAssociationTestType phylogenyAssociationTestEntity =
                        (PhylogenyAssociationTestType)currPhyloTestOut.getPhylogenyAssociationTestId();
                    Map<Integer, List<PhylogenyTestResult>> matchingResults = phylogenyTestResultsCache.get(
                            phylogenyAssociationTestEntity.getId());
                    if(matchingResults == null)
                    {
                        PhylogenyAssociationTest phyloTest =
                            JaxbPhylogenyAssociationTestFactory.getPhylogenyAssociationTest(
                                    phylogenyAssociationTestEntity);
                        matchingResults = phyloTest.getTestResults();
                        phylogenyTestResultsCache.put(
                                phylogenyAssociationTestEntity.getId(),
                                matchingResults);
                    }
                    
                    System.out.println("# matching results: " + matchingResults.size());
                    for(Entry<Integer, List<PhylogenyTestResult>> chromoTestEntry:
                        matchingResults.entrySet())
                    {
                        for(PhylogenyTestResult phylogenyTestResult: chromoTestEntry.getValue())
                        {
                            BasePairInterval interval =
                                phylogenyTestResult.getPhylogenyInterval().getInterval();
                            PhylogenyTreeNode phylogeny =
                                phylogenyTestResult.getPhylogenyInterval().getPhylogeny();
                            fileOut.print(
                                    interval.getStartInBasePairs() + "\t" +
                                    interval.getExtentInBasePairs() + "\t" +
                                    phylogeny.toNewickFormat());
                            fileOut.print('\t');
                            fileOut.print(phylogenyTestResult.getPValue());
                            fileOut.println();
                        }
                    }
                    fileOut.close();
                }
                else
                {
                    throw new IllegalStateException(
                            "unknown output type: " +
                            currPhyloTestOut.getClass().getName());
                }
            }
        }
        catch(IOException ex)
        {
            LOG.log(Level.SEVERE,
                    "file IO error",
                    ex);
        }
    }

    /**
     * Take care of writing the haplotype test results data sources to the given
     * outputs
     * @param haplotypeOutputList
     *          the output list
     */
    private void generateHaplotypeAssociationTestOutputs(
            List<HaplotypeAssociationTestOutputType> haplotypeOutputList)
    {
        try
        {
            // create an initially empty cache for test results
            Map<String, HaplotypeEquivalenceClassTestResult[]> haplotypeTestResultsCache =
                new HashMap<String, HaplotypeEquivalenceClassTestResult[]>();
            
            for(HaplotypeAssociationTestOutputType currHapTestOut: haplotypeOutputList)
            {
                if(currHapTestOut instanceof CommaSeparatedHaplotypeAssociationTestOutputType)
                {
                    CommaSeparatedHaplotypeAssociationTestOutputType currCommaSeparatedHapTestOut =
                        (CommaSeparatedHaplotypeAssociationTestOutputType)currHapTestOut;
                    
                    PrintStream fileOut = new PrintStream(
                            new BufferedOutputStream(new FileOutputStream(
                            currCommaSeparatedHapTestOut.getFileLocation())));
                    
                    fileOut.println("pvalue,chromosome,block_start_bp,block_extent_bp,block_end_bp,block_middle_bp,equiv_class_extent_bp,strains");
                    
                    HaplotypeAssociationTestType haplotypeAssociationTestEntity =
                        (HaplotypeAssociationTestType)currHapTestOut.getHaplotypeAssociationTestId();
                    HaplotypeEquivalenceClassTestResult[] matchingResults = haplotypeTestResultsCache.get(
                            haplotypeAssociationTestEntity.getId());
                    if(matchingResults == null)
                    {
                        HaplotypeAssociationTest hapAssocTest = JaxbHaplotypeAssociationTestFactory.getHaplotypeAssociationTest(
                                haplotypeAssociationTestEntity);
                        matchingResults = hapAssocTest.getEquivalenceClassTestResults();
                        haplotypeTestResultsCache.put(
                                haplotypeAssociationTestEntity.getId(),
                                matchingResults);
                    }
                    
                    for(int i = 0; i < matchingResults.length; i++)
                    {
                        PartitionedIntervalSet currHaplotypeEquivClass =
                            matchingResults[i].getHaplotypeEquivalenceClass();
                        BitSet currStrainsBits =
                            currHaplotypeEquivClass.getStrainBitSet();
                        
                        for(BasePairInterval currSnpBlock:
                            currHaplotypeEquivClass.getSnpIntervals())
                        {
                            fileOut.print(matchingResults[i].getPValue());
                            fileOut.print(CharacterDelimitedParser.DEFAULT_DELIMITER_CHAR);
                            
                            fileOut.print(currSnpBlock.getChromosomeNumber());
                            fileOut.print(CharacterDelimitedParser.DEFAULT_DELIMITER_CHAR);
                            
                            fileOut.print(currSnpBlock.getStartInBasePairs());
                            fileOut.print(CharacterDelimitedParser.DEFAULT_DELIMITER_CHAR);
                            
                            fileOut.print(currSnpBlock.getExtentInBasePairs());
                            fileOut.print(CharacterDelimitedParser.DEFAULT_DELIMITER_CHAR);
    
                            fileOut.print(currSnpBlock.getEndInBasePairs());
                            fileOut.print(CharacterDelimitedParser.DEFAULT_DELIMITER_CHAR);
                            
                            fileOut.print(
                                    currSnpBlock.getStartInBasePairs() +
                                    (currSnpBlock.getExtentInBasePairs() / 2.0));
                            fileOut.print(CharacterDelimitedParser.DEFAULT_DELIMITER_CHAR);
                            
                            fileOut.print(currHaplotypeEquivClass.getCumulativeExtentInBasePairs());
                            fileOut.print(CharacterDelimitedParser.DEFAULT_DELIMITER_CHAR);
    
                            fileOut.print(SetUtilities.bitSetToBinaryString(
                                    currStrainsBits));
                            fileOut.println();
                        }
                    }
                    
                    fileOut.close();
                }
                else
                {
                    throw new IllegalStateException(
                            "unknown output type: " +
                            currHapTestOut.getClass().getName());
                }
            }
        }
        catch(IOException ex)
        {
            LOG.log(Level.SEVERE,
                    "file IO error",
                    ex);
        }
    }

    /**
     * the main application function for haplotype association mapping
     * @param args
     *          command line arguments
     */
    public static void main(String[] args)
    {
        // Deal with the options.
        CommandLineParser parser = new GnuParser();
        Options options = new Options();
        CommandLine commandLine = null;
        
        final Option helpOption;
        {
            helpOption = new Option("help", "Print this help and exit");
            helpOption.setRequired(false);
            options.addOption(helpOption);
        }
        
        final Option designfileOption;
        {
            designfileOption = new Option("designfile", "the analysis design file");
            designfileOption.setRequired(false);
            designfileOption.setArgs(1);
            designfileOption.setArgName("file name");
            options.addOption(designfileOption);
        }
        
        try
        {
            commandLine = parser.parse(options, args);
            
            // See if we just need to print the help options.
            if(commandLine.hasOption(helpOption.getOpt()))
            {
                HelpFormatter helpFormatter = new HelpFormatter();
                helpFormatter.printHelp("ham-analysis", options);
            }
            else
            {
                if(commandLine.hasOption(designfileOption.getOpt()))
                {
                    String fileName =
                        commandLine.getOptionValue(designfileOption.getOpt());
                    
                    HaplotypeAssociationMappingMain mainInstance =
                        new HaplotypeAssociationMappingMain();
                    mainInstance.setAnalysisDesignFile(fileName);
                    mainInstance.performAnalysis();
                }
                else
                {
                    System.out.println("initialization failed");
                    HelpFormatter helpFormatter = new HelpFormatter();
                    helpFormatter.printHelp("ham-analysis", options);
                }
            }
        }
        catch(ParseException ex)
        {
            LOG.log(Level.SEVERE,
                    "initialization failed",
                    ex);
            HelpFormatter helpFormatter = new HelpFormatter();
            helpFormatter.printHelp("ham-analysis", options);
        }
    }
}
