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

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * Pretty much the same as {@link HaplotypeAssociationTest} except that
 * the test results and common strains are cached.
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class CachingHaplotypeAssociationTest extends HaplotypeAssociationTest
{
    private HaplotypeEquivalenceClassTestResult[] equivClassTestResults = null;
    
    private final Map<Integer, HaplotypeBlockTestResult[]> chromoToHaploBlockTestResultsMap =
        new HashMap<Integer, HaplotypeBlockTestResult[]>();
    
    private Set<String> commonStrains = null;
    
    /**
     * Constructor that allows you to create a caching data source from a
     * non-caching data source
     * @param testToCache
     *          the association test that we should cache
     */
    public CachingHaplotypeAssociationTest(
            HaplotypeAssociationTest testToCache)
    {
        this(testToCache.getName(),
             testToCache.getHaplotypeDataSource(),
             testToCache.getPhenotypeDataSource());
    }
    
    /**
     * Constructor
     * @param name
     *          the name of this test
     * @param haplotypeDataSource
     *          the haplotype data source
     * @param phenotypeDataSource
     *          the phenotype data source
     */
    public CachingHaplotypeAssociationTest(
            String name,
            HaplotypeDataSource haplotypeDataSource,
            PhenotypeDataSource phenotypeDataSource)
    {
        super(name, haplotypeDataSource, phenotypeDataSource);
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public synchronized HaplotypeEquivalenceClassTestResult[] getEquivalenceClassTestResults()
    {
        if(this.equivClassTestResults == null)
        {
            this.equivClassTestResults = super.getEquivalenceClassTestResults();
        }
        
        return this.equivClassTestResults;
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public synchronized HaplotypeBlockTestResult[] getHaplotypeTestResults(
            int chromosomeNumber)
    {
        HaplotypeBlockTestResult[] results =
            this.chromoToHaploBlockTestResultsMap.get(chromosomeNumber);
        if(results == null)
        {
            results = super.getHaplotypeTestResults(chromosomeNumber);
            this.chromoToHaploBlockTestResultsMap.put(
                    chromosomeNumber,
                    results);
        }
        
        return results;
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public synchronized Set<String> getCommonStrains()
    {
        if(this.commonStrains == null)
        {
            this.commonStrains = super.getCommonStrains();
        }
        
        return this.commonStrains;
    }
}
