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

import java.io.Serializable;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jax.geneticutil.data.PartitionedInterval;
import org.jax.geneticutil.data.PartitionedIntervalSet;

/**
 * Interface used by haplotype significance tests
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class HaplotypeAssociationTest implements Serializable
{
    /**
     * every {@link java.io.Serializable} is supposed to have one of these
     */
    private static final long serialVersionUID = 1260525589152046482L;
    
    private final String name;
    private final HaplotypeDataSource haplotypeDataSource;
    private final PhenotypeDataSource phenotypeDataSource;
    
    /**
     * Constructor
     * @param name
     *          the name to use for this test (should be unique among tests)
     * @param haplotypeDataSource
     *          the haplotype data source to use
     * @param phenotypeDataSource
     *          the phenotype data source to use
     */
    public HaplotypeAssociationTest(
            String name,
            HaplotypeDataSource haplotypeDataSource,
            PhenotypeDataSource phenotypeDataSource)
    {
        this.name = name;
        this.haplotypeDataSource = haplotypeDataSource;
        this.phenotypeDataSource = phenotypeDataSource;
    }
    
    /**
     * Getter for the name of this test
     * @return the name
     */
    public String getName()
    {
        return this.name;
    }
    
    /**
     * The haplotype data source
     * @return the haplotypeDataSource
     */
    public HaplotypeDataSource getHaplotypeDataSource()
    {
        return this.haplotypeDataSource;
    }
    
    /**
     * The phenotype data source
     * @return the phenotypeDataSource
     */
    public PhenotypeDataSource getPhenotypeDataSource()
    {
        return this.phenotypeDataSource;
    }
    
    /**
     * get the strain names that are common between the
     * {@link #getPhenotypeDataSource()} and {@link #getHaplotypeDataSource()}
     * @return
     *          the common strain set
     */
    public Set<String> getCommonStrains()
    {
        Set<String> haploStrains =
            this.getHaplotypeDataSource().getAvailableStrains();
        Set<String> phyloStrains =
            this.getPhenotypeDataSource().getPhenotypeData().keySet();
        
        if(haploStrains.equals(phyloStrains))
        {
            // they're equal... no need to perform the intersection operation
            return haploStrains;
        }
        else
        {
            Set<String> intersectionStrains = new HashSet<String>(haploStrains);
            intersectionStrains.retainAll(phyloStrains);
            return intersectionStrains;
        }
    }

    /**
     * Perform the significance test
     * @return
     *          the results of the test
     */
    public HaplotypeEquivalenceClassTestResult[] getEquivalenceClassTestResults()
    {
        StrainBinaryPartitionSignificanceTester haplotypeSignificanceTester =
            new StrainBinaryPartitionSignificanceTester();
        
        // find the common set of strains and sort them
        Set<String> haplotypeStrains =
            this.haplotypeDataSource.getAvailableStrains();
        Map<String, List<Double>> phenotypeDataMap =
            this.phenotypeDataSource.getPhenotypeData();
        int originalPhenoStrainCount = phenotypeDataMap.size();
        phenotypeDataMap.keySet().retainAll(haplotypeStrains);
        
        String[] sortedCommonStrains = phenotypeDataMap.keySet().toArray(
                new String[phenotypeDataMap.size()]);
        Arrays.sort(sortedCommonStrains);
        
        System.out.println("# of haplotype strains: " + haplotypeStrains.size());
        System.out.println("# of phenotype strains: " + originalPhenoStrainCount);
        System.out.println("# of strains in common: " + sortedCommonStrains.length);
        
        List<PartitionedIntervalSet> haplotypeDataList =
            this.haplotypeDataSource.getHaplotypeEquivalenceClassData(
                    phenotypeDataMap.keySet());
        
        double[][] orderedPhenotypData =
            new double[phenotypeDataMap.size()][];
        for(int strainIndex = 0;
            strainIndex < sortedCommonStrains.length;
            strainIndex++)
        {
            List<Double> currPhenoResponseList = phenotypeDataMap.get(
                    sortedCommonStrains[strainIndex]);
            double[] newPhenoResponses =
                new double[currPhenoResponseList.size()];
            for(int responseIndex = 0;
                responseIndex < currPhenoResponseList.size();
                responseIndex++)
            {
                newPhenoResponses[responseIndex] =
                    currPhenoResponseList.get(responseIndex);
            }
            
            orderedPhenotypData[strainIndex] = newPhenoResponses;
        }
        
        double[] pValueResults =
            haplotypeSignificanceTester.tTestMultipleResponseSignificance(
                haplotypeDataList,
                orderedPhenotypData);
        
        System.out.println(
                "completed processing " + pValueResults.length +
                " equivalence classes for significance test: " +
                this.name);
        
        if(haplotypeDataList.size() == pValueResults.length)
        {
            HaplotypeEquivalenceClassTestResult[] newTestResults =
                new HaplotypeEquivalenceClassTestResult[pValueResults.length];
            Iterator<PartitionedIntervalSet> haplotypeIter =
                haplotypeDataList.iterator();
            for(int i = 0; i < pValueResults.length; i++)
            {
                PartitionedIntervalSet currHapEquivClass =
                    haplotypeIter.next();
                newTestResults[i] = new HaplotypeEquivalenceClassTestResult(
                        currHapEquivClass,
                        pValueResults[i]);
            }
            
            return newTestResults;
        }
        else
        {
            throw new IllegalStateException(
                    "haplotype data size doesn't match significance " +
                    "test size: " + haplotypeDataList.size() + " vs " +
                    pValueResults.length);
        }
    }
    
    /**
     * Perform the significance test on haplotypes
     * @param chromosomeNumber
     *          the chromosome number to get results for
     * @return
     *          the results of the test
     */
    public HaplotypeBlockTestResult[] getHaplotypeTestResults(
            int chromosomeNumber)
    {
        StrainBinaryPartitionSignificanceTester haplotypeSignificanceTester =
            new StrainBinaryPartitionSignificanceTester();
        
        // find the common set of strains and sort them
        Set<String> haplotypeStrains =
            this.haplotypeDataSource.getAvailableStrains();
        Map<String, List<Double>> phenotypeDataMap =
            this.phenotypeDataSource.getPhenotypeData();
        int originalPhenoStrainCount = phenotypeDataMap.size();
        phenotypeDataMap.keySet().retainAll(haplotypeStrains);
        
        String[] sortedCommonStrains = phenotypeDataMap.keySet().toArray(
                new String[phenotypeDataMap.size()]);
        Arrays.sort(sortedCommonStrains);
        
        System.out.println("# of haplotype strains: " + haplotypeStrains.size());
        System.out.println("# of phenotype strains: " + originalPhenoStrainCount);
        System.out.println("# of strains in common: " + sortedCommonStrains.length);
        
        List<PartitionedInterval> haplotypeDataList = this.haplotypeDataSource.getHaplotypeData(
                Collections.singleton(chromosomeNumber),
                phenotypeDataMap.keySet());
        
        double[][] orderedPhenotypData =
            new double[phenotypeDataMap.size()][];
        for(int strainIndex = 0;
            strainIndex < sortedCommonStrains.length;
            strainIndex++)
        {
            List<Double> currPhenoResponseList = phenotypeDataMap.get(
                    sortedCommonStrains[strainIndex]);
            double[] newPhenoResponses =
                new double[currPhenoResponseList.size()];
            for(int responseIndex = 0;
                responseIndex < currPhenoResponseList.size();
                responseIndex++)
            {
                newPhenoResponses[responseIndex] =
                    currPhenoResponseList.get(responseIndex);
            }
            
            orderedPhenotypData[strainIndex] = newPhenoResponses;
        }
        
        double[] pValueResults =
            haplotypeSignificanceTester.tTestMultipleResponseSignificance(
                haplotypeDataList,
                orderedPhenotypData);
        
        System.out.println(
                "completed processing " + pValueResults.length +
                " haplotype blocks for significance test: " +
                this.name);
        
        if(haplotypeDataList.size() == pValueResults.length)
        {
            HaplotypeBlockTestResult[] newTestResults =
                new HaplotypeBlockTestResult[pValueResults.length];
            
            Iterator<PartitionedInterval> haplotypeIter =
                haplotypeDataList.iterator();
            for(int i = 0; i < pValueResults.length; i++)
            {
                PartitionedInterval nextHaploBlock = haplotypeIter.next();
                newTestResults[i] = new HaplotypeBlockTestResult(
                        nextHaploBlock,
                        pValueResults[i]);
            }
            
            return newTestResults;
        }
        else
        {
            throw new IllegalStateException(
                    "haplotype data size doesn't match significance " +
                    "test size: " + haplotypeDataList.size() + " vs " +
                    pValueResults.length);
        }
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public String toString()
    {
        return this.getName();
    }
}
