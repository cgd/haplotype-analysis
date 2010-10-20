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

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.jax.geneticutil.data.PartitionedInterval;
import org.jax.geneticutil.data.PartitionedIntervalSet;
import org.jax.haplotype.data.ChromosomeDataSource;
import org.jax.haplotype.data.GenomeDataSource;
import org.jax.haplotype.inference.HaplotypeEquivalenceClassCreator;
import org.jax.haplotype.inference.HaplotypeEstimator;
import org.jax.haplotype.inference.IntervalScanningHaplotypeEstimator;
import org.jax.haplotype.io.SdpInputStream;

/**
 * A haplotype data source that infers haplotype structure from genotype
 * data.
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class GenotypeInferenceHaplotypeDataSource implements HaplotypeDataSource
{
    /**
     * every {@link java.io.Serializable} is supposed to have one of these
     */
    private static final long serialVersionUID = -3937759676101800546L;

    private static final Logger LOG = Logger.getLogger(
            GenotypeInferenceHaplotypeDataSource.class.getName());
    
    private final GenomeDataSource genomeDataSource;
    
    private final int minimumSnpExtent;
    
    private final int minimumStrainGroupSize;

    private final Set<String> persistentStrainsToAcceptFilter;

    private final String name;
    
    /**
     * Constructor
     * @param name
     *          the name of this data source (can be null)
     * @param genomeDataSource
     *          the genome data source to use
     * @param persistentStrainsToAcceptFilter 
     *          the strains we should accept
     * @param minimumSnpExtent
     *          the minimum extent in SNPs use before calling something
     *          a haplotype block
     * @param minimumStrainGroupSize
     *          the minimum # of strains that have to be in a block
     */
    public GenotypeInferenceHaplotypeDataSource(
            String name,
            GenomeDataSource genomeDataSource,
            Set<String> persistentStrainsToAcceptFilter,
            int minimumSnpExtent,
            int minimumStrainGroupSize)
    {
        this.name = name;
        this.genomeDataSource = genomeDataSource;
        this.persistentStrainsToAcceptFilter = persistentStrainsToAcceptFilter;
        this.minimumSnpExtent = minimumSnpExtent;
        this.minimumStrainGroupSize = minimumStrainGroupSize;
    }
    
    /**
     * {@inheritDoc}
     */
    public String getName()
    {
        return this.name;
    }

    /**
     * {@inheritDoc}
     */
    public Set<String> getAvailableStrains()
    {
        Set<String> availableStrains = this.genomeDataSource.getAvailableStrains();
        
        if(this.persistentStrainsToAcceptFilter != null)
        {
            availableStrains.retainAll(this.persistentStrainsToAcceptFilter);
        }
        
        return availableStrains;
    }
    
    /**
     * {@inheritDoc}
     */
    public int[] getAvailableChromosomes()
    {
        return this.genomeDataSource.getAvailableChromosomes();
    }
    
    /**
     * {@inheritDoc}
     */
    public List<PartitionedInterval> getHaplotypeData(
            Set<Integer> chromosomesToAccept,
            Set<String> strainsToAcceptFilter)
    {
        try
        {
            if(strainsToAcceptFilter == null)
            {
                strainsToAcceptFilter = this.getAvailableStrains();
            }
            else
            {
                strainsToAcceptFilter = new HashSet<String>(strainsToAcceptFilter);
                strainsToAcceptFilter.retainAll(this.getAvailableStrains());
            }
            
            String[] strainsToAcceptArray = strainsToAcceptFilter.toArray(
                    new String[strainsToAcceptFilter.size()]);
            Arrays.sort(strainsToAcceptArray);
            
            HaplotypeEstimator haplotypeEstimator =
                new IntervalScanningHaplotypeEstimator(
                        this.minimumSnpExtent,
                        this.minimumStrainGroupSize);
            
            List<PartitionedInterval> haplotypeBlocks =
                new ArrayList<PartitionedInterval>();
            Map<Integer, ? extends ChromosomeDataSource> chromosomeDataSources =
                this.genomeDataSource.getChromosomeDataSources();
            List<Integer> chromosomeNumbers = new ArrayList<Integer>(
                    chromosomeDataSources.keySet());
            Collections.sort(chromosomeNumbers);
            
            // apply chromosome filtering
            if(chromosomesToAccept != null)
            {
                chromosomeNumbers.retainAll(chromosomesToAccept);
            }
            
            for(Integer chromosomeNumber: chromosomeNumbers)
            {
                ChromosomeDataSource currChromoDataSource =
                    chromosomeDataSources.get(chromosomeNumber);
                SdpInputStream sdpInputStream = currChromoDataSource.getSdpInputStream(
                        strainsToAcceptArray);
                List<PartitionedInterval> currHaplotypeData =
                    haplotypeEstimator.estimateHaplotypeBlocks(
                            sdpInputStream,
                            currChromoDataSource.getSnpPositionInputStream());
                haplotypeBlocks.addAll(currHaplotypeData);
            }
            
            return haplotypeBlocks;
        }
        catch(IOException ex)
        {
            LOG.log(Level.SEVERE,
                    "Passing IOException on as a runtime exception",
                    ex);
            throw new RuntimeException(ex);
        }
    }

    /**
     * {@inheritDoc}
     */
    public List<PartitionedIntervalSet> getHaplotypeEquivalenceClassData(
            Set<String> strainsToAcceptFilter)
    {
        List<PartitionedInterval> haplotypeBlocks = this.getHaplotypeData(
                null,
                strainsToAcceptFilter);
        
        return HaplotypeEquivalenceClassCreator.createEquivalenceClassesFromBlocks(
                haplotypeBlocks);
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
