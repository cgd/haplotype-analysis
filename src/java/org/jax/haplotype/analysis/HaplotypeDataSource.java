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
import java.util.List;
import java.util.Set;

import org.jax.geneticutil.data.PartitionedInterval;
import org.jax.geneticutil.data.PartitionedIntervalSet;

/**
 * The interface for haplotype data sources
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public interface HaplotypeDataSource extends Serializable
{
    /**
     * Getter for the name of this haplotype data source
     * @return
     *          the name of the haplotype data source (can be null)
     */
    public String getName();
    
    /**
     * Get the equivalence class data
     * @param strainsToAccept
     *          the strains we should consider when reading in the data. use
     *          null if you want to use all available strains
     * @return
     *          the equivalence classes. the bit ordering used will be
     *          the sorted ordering of the strains that make up the intersection
     *          of {@code strainsToAccept} and {@link #getAvailableStrains()}
     */
    public List<PartitionedIntervalSet> getHaplotypeEquivalenceClassData(
            Set<String> strainsToAccept);
    
    /**
     * Get the haplotype blocks
     * @param chromosomesToAccept
     *          the chromosomes to allow through
     * @param strainsToAccept
     *          the strains to allow through
     * @return
     *          the list of haplotype blocks given the strain and chromosome
     *          filters
     */
    public List<PartitionedInterval> getHaplotypeData(
            Set<Integer> chromosomesToAccept,
            Set<String> strainsToAccept);
    
    /**
     * Get the available strains in the data
     * @return
     *          the strains
     */
    public Set<String> getAvailableStrains();
    
    /**
     * Getter for the available chromosome numbers
     * @return
     *          the chromosomes
     */
    public int[] getAvailableChromosomes();
}
