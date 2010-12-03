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

package org.jax.ham.analysis;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;

import org.jax.geneticutil.data.PartitionedInterval;
import org.jax.geneticutil.data.PartitionedIntervalSet;
import org.jax.haplotype.analysis.MPDIndividualStrainPhenotypeParser;
import org.jax.haplotype.analysis.SexFilter;
import org.jax.haplotype.analysis.StrainBinaryPartitionSignificanceTester;
import org.jax.haplotype.data.ChromosomeDataSource;
import org.jax.haplotype.data.CommaSeparatedChromosomeDataSource;
import org.jax.haplotype.inference.HaplotypeEquivalenceClassCreator;
import org.jax.haplotype.inference.IntervalScanningHaplotypeEstimator;
import org.jax.haplotype.io.GenotypeParser;
import org.jax.util.math.StatisticUtilities;
import org.junit.Test;


/**
 * A unit test for the significance tester.
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class HaplotypeSignificanceTesterTest
{
    /**
     * 
     * @throws FileNotFoundException
     * @throws IOException
     */
    @Test
    public void testSignificanceTest() throws FileNotFoundException, IOException
    {
//        MPDIndividualStrainPhenotypeParser phenoParser =
//            new MPDIndividualStrainPhenotypeParser();
//        Map<String, Double> meanPhenotypes = phenoParser.parseMeanPhenotypeFromStream(
//                new FileInputStream("data/mpd-diastolic-bp.csv"));
//        for(Entry<String, Double> currPhenoEntry: meanPhenotypes.entrySet())
//        {
//            System.out.println("Strain:     " + currPhenoEntry.getKey());
//            System.out.println("Mean Pheno: " + currPhenoEntry.getValue());
//        }
        int minNumSnps = 3;
        int minNumChromo = 4;
        IntervalScanningHaplotypeEstimator slidingWindowEstimator;
        GenotypeParser parser = new GenotypeParser();
        List<PartitionedInterval> allHaplotypeBlocks =
            new ArrayList<PartitionedInterval>();
        int cumulativeChromosomeEquivClasses = 0;
        MPDIndividualStrainPhenotypeParser phenoParser =
            new MPDIndividualStrainPhenotypeParser();
        Set<String> phenoStrainNames = phenoParser.parseAvailableStrainNames(
                new FileInputStream("data/phenotype/peters1-rbc.csv"));
        Set<String> genoStrainNames = new HashSet<String>();
        for(int i = 1; i <= 20; i++)
        {
            String numOrX;
            if(i == 20)
            {
                numOrX = "X";
            }
            else
            {
                numOrX = Integer.toString(i);
            }
            
            String fileName = "PBHMM2_SNP_chr" + numOrX + "_genoData.csv";
            File file = new File("data");
            file = new File(file, "genotype");
            file = new File(file, "old-genodata");
            file = new File(file, fileName);
            ChromosomeDataSource chromosomeDataSource = new CommaSeparatedChromosomeDataSource(
                    file,
                    i);
            genoStrainNames = chromosomeDataSource.getAvailableStrains();
//            Set<StrainChromosome> genotype = parser.parseGenotypeFromStream(
//                    new FileInputStream(file),
//                    true,
//                    phenoStrainNames);
//            System.out.println("chromosomes parsed: " + genotype.size());
//            for(StrainChromosome currChromosome: genotype)
//            {
//                genoStrainNames.add(currChromosome.getStrainName());
//            }
            
            String[] sortedStrainNames = genoStrainNames.toArray(new String[0]);
            Arrays.sort(sortedStrainNames);
            slidingWindowEstimator =
                new IntervalScanningHaplotypeEstimator(
                        minNumSnps,
                        minNumChromo);
            List<PartitionedInterval> estimatedHaplotypes = slidingWindowEstimator.estimateHaplotypeBlocks(
                    chromosomeDataSource.getSdpInputStream(sortedStrainNames),
                    chromosomeDataSource.getSnpPositionInputStream());
            Collections.sort(estimatedHaplotypes);
            
            System.out.println(
                    "Number of haplotype blocks for chromosome " + numOrX + ": " +
                    estimatedHaplotypes.size());
            List<PartitionedIntervalSet> currEquivalenceClasses =
                HaplotypeEquivalenceClassCreator.createEquivalenceClassesFromBlocks(
                        estimatedHaplotypes);
            System.out.println(
                    "Number of unique strain groups for chromosome " + numOrX +
                    ": " + currEquivalenceClasses.size());
            cumulativeChromosomeEquivClasses += currEquivalenceClasses.size();
            
            allHaplotypeBlocks.addAll(estimatedHaplotypes);
        }
        
        List<PartitionedIntervalSet> allEquivalenceClasses =
            HaplotypeEquivalenceClassCreator.createEquivalenceClassesFromBlocks(
                allHaplotypeBlocks);
        System.out.println(
                "total number of equivalence classes: " + allEquivalenceClasses.size());
        System.out.println(
                "cumulative chromosome only equivalence classes: " +
                cumulativeChromosomeEquivClasses);
        
        Set<String> genoPhenoStrainSet = new HashSet<String>(genoStrainNames);
        genoPhenoStrainSet.retainAll(phenoStrainNames);
        String[] sortedGenoPhenoStrains = genoPhenoStrainSet.toArray(new String[0]);
        Arrays.sort(sortedGenoPhenoStrains);
        
        Set<String> phenotypeNames = phenoParser.parseAvailablePhenotypes(
                new FileInputStream("data/phenotype/peters1-rbc.csv"));
        for(String currPhenoName: phenotypeNames)
        {
            System.out.println("phenotype name: " + currPhenoName);
            
            Map<String, List<Double>> strainToPhenoMap =
                phenoParser.parsePhenotypesFromStream(
                        currPhenoName,
                        new FileInputStream("data/phenotype/peters1-rbc.csv"),
                        SexFilter.AGNOSTIC,
                        genoStrainNames);
            
//            List<Double> randomPhenotypeValueBucket = new ArrayList<Double>();
            System.out.println("Num strains: " + strainToPhenoMap.size());
            for(Entry<String, List<Double>> currEntry: strainToPhenoMap.entrySet())
            {
                System.out.print(currEntry.getKey() + " Values:");
                for(double currPheno: currEntry.getValue())
                {
                    System.out.print(" " + currPheno);
                }
                System.out.println();
//                randomPhenotypeValueBucket.addAll(currEntry.getValue());
            }
            
            double[][] orderedPhenotypeResponses = new double[sortedGenoPhenoStrains.length][];
            for(int i = 0; i < sortedGenoPhenoStrains.length; i++)
            {
                String currStrain = sortedGenoPhenoStrains[i];
                List<Double> currResponseList = strainToPhenoMap.get(currStrain);
                double[] currResponseArray = new double[currResponseList.size()];
                for(int j = 0; j < currResponseList.size(); j++)
                {
                    currResponseArray[j] = currResponseList.get(j);
                }
                
                orderedPhenotypeResponses[i] = currResponseArray;
            }
            
            StrainBinaryPartitionSignificanceTester significanceTester =
                new StrainBinaryPartitionSignificanceTester();
            
            double[] pValues = significanceTester.tTestMultipleResponseSignificance(
                    allEquivalenceClasses,
                    orderedPhenotypeResponses);
            Arrays.sort(pValues);
            System.out.println(
                    "min p-value: " +
                    pValues[0]);
            
            double[] normalizedPValues = significanceTester.normalizedTestMultipleResponseSignificance(
                    allEquivalenceClasses,
                    orderedPhenotypeResponses);
            Arrays.sort(normalizedPValues);
            System.out.println(
                    "normalized min p-value: " +
                    normalizedPValues[0]);
            
            System.out.print("running permutations: ");
//            int permCount = 10000;
            int permCount = 200;
            double[] permutationMinPValues = new double[permCount];
            double[] normalizedPermutationMinPValues = new double[permCount];
//            Random rand = new Random();
//            for(int i = 0; i < permCount; i++)
//            {
//                System.out.print("*");
//                
//                List<List<Double>> randomPhenotypes = new ArrayList<List<Double>>(
//                        strainToPhenoMap.values());
//                StatisticUtilities.shuffle(randomPhenotypes, rand);
//                List<String> strainNames = new ArrayList<String>(strainToPhenoMap.keySet());
//                Map<String, List<Double>> randomStrainToPhenoMap =
//                    new HashMap<String, List<Double>>();
//                for(int j = 0; j < randomPhenotypes.size(); j++)
//                {
//                    randomStrainToPhenoMap.put(
//                            strainNames.get(j),
//                            randomPhenotypes.get(j));
//                }
//                
//                double[] currPermSignificanceValues = significanceTester.testMultipleResponseSignificance(
//                        allEquivalenceClasses,
//                        randomStrainToPhenoMap);
//                Arrays.sort(currPermSignificanceValues);
//                permutationMinPValues[i] =
//                    currPermSignificanceValues[0];
//                
//                double[] normalizedCurrPermSigValues = significanceTester.normalizedTestMultipleResponseSignificance(
//                        allEquivalenceClasses,
//                        randomStrainToPhenoMap);
//                Arrays.sort(normalizedCurrPermSigValues);
//                normalizedPermutationMinPValues[i] =
//                    normalizedCurrPermSigValues[0];
//            }
            System.out.println();
            
            Arrays.sort(permutationMinPValues);
            Arrays.sort(normalizedPermutationMinPValues);
//            System.out.println("Permutation p-values:");
//            for(double currPermVal: permutationMinPValues)
//            {
//                System.out.println(" " + currPermVal);
//            }
            
            double pValuePValue = StatisticUtilities.calculatePValueForDataPoint(
                    pValues[0],
                    permutationMinPValues);
            double normalizedPValuePValue = StatisticUtilities.calculatePValueForDataPoint(
                    normalizedPValues[0],
                    normalizedPermutationMinPValues);
            System.out.println("p-value p-value: " + pValuePValue);
            System.out.println("normalized p-value p-value: " + normalizedPValuePValue);
        }
    }
}
