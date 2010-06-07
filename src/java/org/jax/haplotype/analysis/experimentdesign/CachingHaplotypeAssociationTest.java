/*
 * Copyright (c) 2008 The Jackson Laboratory
 *
 * Permission is hereby granted, free of charge, to any person obtaining  a copy
 * of this software and associated documentation files (the  "Software"), to
 * deal in the Software without restriction, including  without limitation the
 * rights to use, copy, modify, merge, publish,  distribute, sublicense, and/or
 * sell copies of the Software, and to  permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be  included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,  EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF  MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY  CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,  TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE  SOFTWARE OR THE
 * USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.jax.haplotype.analysis.experimentdesign;

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
