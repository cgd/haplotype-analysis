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

import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jax.geneticutil.data.MultiPartitionedInterval;
import org.jax.haplotype.data.MultiGroupHaplotypeDataSource;

/**
 * A test class that uses the sliding window algorithm
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class MultiGroupHaplotypeAssociationTest implements MultiHaplotypeBlockTest
{
    /**
     * every {@link java.io.Serializable} is supposed to have one of these
     */
    private static final long serialVersionUID = -6951713978946827064L;
    
    private final String name;
    private final MultiGroupHaplotypeDataSource haplotypeDataSource;
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
    public MultiGroupHaplotypeAssociationTest(
            String name,
            MultiGroupHaplotypeDataSource haplotypeDataSource,
            PhenotypeDataSource phenotypeDataSource)
    {
        this.name = name;
        this.haplotypeDataSource = haplotypeDataSource;
        this.phenotypeDataSource = phenotypeDataSource;
    }
    
    /**
     * Getter for the haplotype data source
     * @return the data source
     */
    public MultiGroupHaplotypeDataSource getHaplotypeDataSource()
    {
        return this.haplotypeDataSource;
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
     * {@inheritDoc}
     */
    public int[] getAvailableChromosomes()
    {
        return this.haplotypeDataSource.getAvailableChromosomes();
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
     * Perform the significance tests
     * @param chromosomeNumber
     *          the chromosome number to get results for
     * @return
     *          the results of the test
     */
    public MultiHaplotypeBlockTestResult[] getTestResults(
            int chromosomeNumber)
    {
        // find the common set of strains and sort them
        Set<String> haploStrains =
            this.haplotypeDataSource.getAvailableStrains();
        Map<String, List<Double>> phenotypeDataMap =
            this.phenotypeDataSource.getPhenotypeData();
        int originalPhenoStrainCount = phenotypeDataMap.size();
        phenotypeDataMap.keySet().retainAll(haploStrains);
        
        String[] sortedCommonStrains = phenotypeDataMap.keySet().toArray(
                new String[phenotypeDataMap.size()]);
        Arrays.sort(sortedCommonStrains);
        
        System.out.println("# of genotype strains:  " + haploStrains.size());
        System.out.println("# of phenotype strains: " + originalPhenoStrainCount);
        System.out.println("# of strains in common: " + sortedCommonStrains.length);
        
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
        
        StrainMultiPartitionSignificanceTester tester =
            new StrainMultiPartitionSignificanceTester();
        try
        {
            List<MultiPartitionedInterval> haplotypeDataList =
                this.haplotypeDataSource.getHaplotypeData(
                        chromosomeNumber,
                        phenotypeDataMap.keySet());
            double[] pValueResults =
                tester.fTestMultipleResponseSignificance(
                        haplotypeDataList,
                        orderedPhenotypData);
            
            System.out.println(
                    "completed processing " + pValueResults.length +
                    " haplotype blocks for significance test: " +
                    this.name);
            
            if(haplotypeDataList.size() == pValueResults.length)
            {
                MultiHaplotypeBlockTestResult[] newTestResults =
                    new MultiHaplotypeBlockTestResult[pValueResults.length];
                
                Iterator<MultiPartitionedInterval> haplotypeIter =
                    haplotypeDataList.iterator();
                for(int i = 0; i < pValueResults.length; i++)
                {
                    MultiPartitionedInterval nextHaploBlock = haplotypeIter.next();
                    newTestResults[i] = new MultiHaplotypeBlockTestResult(
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
        catch(RuntimeException ex)
        {
            // allow runtime exceptions to pass through
            throw ex;
        }
        catch(Exception ex)
        {
            // this may be a bad thing to do. we're repackaging the exception
            // and throwing it as a runtime
            throw new RuntimeException(ex);
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
