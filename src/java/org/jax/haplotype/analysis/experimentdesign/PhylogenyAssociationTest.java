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

package org.jax.haplotype.analysis.experimentdesign;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;

import org.jax.haplotype.analysis.PhylogenySignificanceTester;
import org.jax.haplotype.phylogeny.data.PhylogenyInterval;
import org.jax.haplotype.phylogeny.data.PhylogenyTestResult;
import org.jax.haplotype.phylogeny.data.PhylogenyTreeEdgeWithRealValue;
import org.jax.haplotype.phylogeny.data.PhylogenyTreeNode;

/**
 * A test description for phylogeny association. To get the results use
 * {@link #getTestResults()}.
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class PhylogenyAssociationTest implements Serializable
{
    /**
     * every {@link java.io.Serializable} is supposed to have one of these
     */
    private static final long serialVersionUID = -7921255440767326768L;

    private static final PhylogenySignificanceTester PHYLOGENY_TESTER =
        new PhylogenySignificanceTester();
    
    private final String name;
    
    private final PhylogenyDataSource phylogenyDataSource;
    
    private final PhenotypeDataSource phenotypeDataSource;

    /**
     * Constructor
     * @param name
     *          the name of the test (should be unique among tests)
     * @param phylogenyDataSource
     *          the data source for phylogenies
     * @param phenotypeDataSource
     *          the data source for phenotypes
     */
    public PhylogenyAssociationTest(
            String name,
            PhylogenyDataSource phylogenyDataSource,
            PhenotypeDataSource phenotypeDataSource)
    {
        this.name = name;
        this.phylogenyDataSource = phylogenyDataSource;
        this.phenotypeDataSource = phenotypeDataSource;
    }
    
    /**
     * Get the phylogeny tester that is used by this test
     * @return the phylogenyTester
     *          the tester
     */
    public PhylogenySignificanceTester getPhylogenyTester()
    {
        return PHYLOGENY_TESTER;
    }
    
    /**
     * Getter for the test name
     * @return the testName
     */
    public String getName()
    {
        return this.name;
    }
    
    /**
     * get the strain names that are common between the
     * {@link #getPhenotypeDataSource()} and {@link #getPhylogenyDataSource()}
     * @return
     *          the common strain set
     */
    public Set<String> getCommonStrains()
    {
        Set<String> haploStrains =
            this.getPhylogenyDataSource().getAvailableStrains();
        Set<String> phyloStrains =
            this.getPhenotypeDataSource().getPhenotypeData().keySet();
        
        if(haploStrains.equals(phyloStrains))
        {
            // they're equal... no need to 
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
     * Getter for the phenotype data used
     * @return the phenotypeDataSource
     */
    public PhenotypeDataSource getPhenotypeDataSource()
    {
        return this.phenotypeDataSource;
    }
    
    /**
     * Getter for the phylogeny data used
     * @return the phylogenyDataSource
     */
    public PhylogenyDataSource getPhylogenyDataSource()
    {
        return this.phylogenyDataSource;
    }
    
    /**
     * Get an array of test results
     * @return
     *          the test results
     */
    public Map<Integer, List<PhylogenyTestResult>> getTestResults()
    {
        PhylogenySignificanceTester phyloTester =
            new PhylogenySignificanceTester();
        
        // find the common set of strains
        Set<String> phyloStrains =
            this.phylogenyDataSource.getAvailableStrains();
        Map<String, List<Double>> phenotypeDataMap =
            this.phenotypeDataSource.getPhenotypeData();
        int originalPhenoStrainCount = phenotypeDataMap.size();
        phenotypeDataMap.keySet().retainAll(phyloStrains);
        
        System.out.println("# of phylogeny strains: " + phyloStrains.size());
        System.out.println("# of phenotype strains: " + originalPhenoStrainCount);
        System.out.println("# of strains in common: " + phenotypeDataMap.size());
        
        Map<Integer, List<PhylogenyInterval>> phyloData = this.phylogenyDataSource.getPhylogenyData(
                phenotypeDataMap.keySet());
        Map<Integer, List<PhylogenyTestResult>> testResults =
            new HashMap<Integer, List<PhylogenyTestResult>>(phyloData.size());
        for(Entry<Integer, List<PhylogenyInterval>> phyloEntry: phyloData.entrySet())
        {
            List<PhylogenyTestResult> currTestResultList =
                new ArrayList<PhylogenyTestResult>(phyloEntry.getValue().size());
            testResults.put(phyloEntry.getKey(), currTestResultList);
            for(PhylogenyInterval currPhyloInterval: phyloEntry.getValue())
            {
                PhylogenyTreeNode phylogenyWithPValue = phyloTester.testMultipleResponseSignificance(
                        currPhyloInterval.getPhylogeny(),
                        phenotypeDataMap);
                PhylogenyTreeEdgeWithRealValue smallestEdge = PhylogenyTreeEdgeWithRealValue.getEdgeWithMininumValue(
                        phylogenyWithPValue);
                
                currTestResultList.add(new PhylogenyTestResult(
                        currPhyloInterval,
                        smallestEdge == null ? 1.0 : smallestEdge.getRealValue()));
            }
        }
        
        return testResults;
    }

    /**
     * Get an array of test results
     * @param chromosomeNumber
     *          the chromosome number
     * @return
     *          the test results for the given chrosome number
     */
    public List<PhylogenyTestResult> getTestResults(int chromosomeNumber)
    {
        // find the common set of strains
        Set<String> phyloStrains =
            this.phylogenyDataSource.getAvailableStrains();
        Map<String, List<Double>> phenotypeDataMap =
            this.phenotypeDataSource.getPhenotypeData();
        int originalPhenoStrainCount = phenotypeDataMap.size();
        phenotypeDataMap.keySet().retainAll(phyloStrains);
        
        System.out.println("# of phylogeny strains: " + phyloStrains.size());
        System.out.println("# of phenotype strains: " + originalPhenoStrainCount);
        System.out.println("# of strains in common: " + phenotypeDataMap.size());
        
        Map<Integer, List<PhylogenyInterval>> phyloData = this.phylogenyDataSource.getPhylogenyData(
                phenotypeDataMap.keySet(),
                Collections.singleton(chromosomeNumber));
        Map<Integer, List<PhylogenyTestResult>> testResults =
            new HashMap<Integer, List<PhylogenyTestResult>>(phyloData.size());
        for(Entry<Integer, List<PhylogenyInterval>> phyloEntry: phyloData.entrySet())
        {
            List<PhylogenyTestResult> currTestResultList =
                new ArrayList<PhylogenyTestResult>(phyloEntry.getValue().size());
            testResults.put(phyloEntry.getKey(), currTestResultList);
            for(PhylogenyInterval currPhyloInterval: phyloEntry.getValue())
            {
                PhylogenyTreeNode phylogenyWithPValue =
                    PHYLOGENY_TESTER.testMultipleResponseSignificance(
                            currPhyloInterval.getPhylogeny(),
                            phenotypeDataMap);
                PhylogenyTreeEdgeWithRealValue smallestEdge =
                    PhylogenyTreeEdgeWithRealValue.getEdgeWithMininumValue(
                            phylogenyWithPValue);
                
                currTestResultList.add(new PhylogenyTestResult(
                        currPhyloInterval,
                        smallestEdge == null ? 1.0 : smallestEdge.getRealValue()));
            }
        }
        
        return testResults.get(chromosomeNumber);
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
