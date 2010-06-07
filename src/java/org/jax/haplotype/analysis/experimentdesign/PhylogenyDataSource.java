/*
 * Copyright (c) 2008 The Jackson Laboratory
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package org.jax.haplotype.analysis.experimentdesign;

import java.io.Serializable;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jax.haplotype.phylogeny.data.PhylogenyInterval;

/**
 * A data source for reading phylogeny data
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public interface PhylogenyDataSource extends Serializable
{
    /**
     * Getter for the name of this phylogeny data source
     * @return
     *          the name (can be null)
     */
    public String getName();
    
    /**
     * Get the data
     * @param strainsToAccept
     *          the strains we should consider when reading in the data, null
     *          means all available strains
     * @return
     *          the phylogeny intervals
     */
    public Map<Integer, List<PhylogenyInterval>> getPhylogenyData(
            Set<String> strainsToAccept);
    
    /**
     * Get the data
     * @param strainsToAccept
     *          the strains we should consider when reading in the data
     * @param chromosomesToAccept
     *          the chromosomes we should consider when reading the data
     * @return
     *          the phylogeny intervals
     */
    public Map<Integer, List<PhylogenyInterval>> getPhylogenyData(
            Set<String> strainsToAccept,
            Set<Integer> chromosomesToAccept);
    
    /**
     * Get the available strains in the data
     * @return
     *          the strains
     */
    public Set<String> getAvailableStrains();

    /**
     * Getter for the available chromosomes
     * @return
     *          the chromosomes
     */
    public int[] getAvailableChromosomes();
}
