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

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.jax.geneticutil.data.BasePairInterval;
import org.jax.geneticutil.data.IndexedSnpInterval;
import org.jax.haplotype.data.ChromosomeDataSource;
import org.jax.haplotype.data.CommaSeparatedChromosomeDataSource;
import org.jax.haplotype.data.GenomeDataSource;
import org.jax.haplotype.io.StreamDirection;
import org.jax.haplotype.phylogeny.data.NoValidPhylogenyException;
import org.jax.haplotype.phylogeny.data.PhylogenyInterval;
import org.jax.haplotype.phylogeny.data.PhylogenyTreeNode;
import org.jax.haplotype.phylogeny.inference.IntervalScanner;
import org.jax.haplotype.phylogeny.inference.PhylogenyScanner;

/**
 * A {@link PhylogenyDataSource} that generates the phylogenies by inferring
 * them from the given {@link GenomeDataSource} using the MAX-K algorithm
 * which was developed by the Computational Biology Group at UNC
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class GenotypeInferencePhylogenyDataSource implements PhylogenyDataSource
{
    /**
     * every {@link java.io.Serializable}
     */
    private static final long serialVersionUID = -6931330567617849010L;

    /**
     * our logger
     */
    private static final Logger LOG = Logger.getLogger(
            GenotypeInferencePhylogenyDataSource.class.getName());
    
    private static final IntervalScanner INTERVAL_SCANNER = new IntervalScanner();
    
    private static final PhylogenyScanner PHYLOGENY_SCANNER = new PhylogenyScanner();
    
    private final GenomeDataSource genomeDataSource;

    private final String name;
    
    /**
     * Constructor
     * @param name
     *          the name of this data source (can be null)
     * @param genomeDataSource
     *          the genome data. note that the chromosome data sources in
     *          the genome data need to be of type
     *          {@link CommaSeparatedChromosomeDataSource} or we'll end up
     *          throwing an {@link IllegalStateException} when
     *          {@link #getPhylogenyData(Set)} is called
     */
    public GenotypeInferencePhylogenyDataSource(
            String name,
            GenomeDataSource genomeDataSource)
    {
        this.name = name;
        this.genomeDataSource = genomeDataSource;
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
        return this.genomeDataSource.getAvailableStrains();
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
    public Map<Integer, List<PhylogenyInterval>> getPhylogenyData(
            Set<String> strainsToAccept)
    {
        return this.getPhylogenyData(
                strainsToAccept,
                null);
    }

    /**
     * {@inheritDoc}
     */
    public Map<Integer, List<PhylogenyInterval>> getPhylogenyData(
            Set<String> strainsToAccept,
            Set<Integer> chromosomesToAccept)
    {
        try
        {
            if(strainsToAccept == null)
            {
                strainsToAccept = this.getAvailableStrains();
            }
            
            Map<Integer, ? extends ChromosomeDataSource> chromoDataSources =
                this.genomeDataSource.getChromosomeDataSources();
            List<Integer> sortedChromosomeNumbers = new ArrayList<Integer>(
                    chromoDataSources.keySet());
            if(chromosomesToAccept != null)
            {
                sortedChromosomeNumbers.retainAll(chromosomesToAccept);
            }
            Collections.sort(sortedChromosomeNumbers);
            
            String[] sortedStrains = strainsToAccept.toArray(new String[0]);
            Arrays.sort(sortedStrains);
            
            Map<Integer, List<PhylogenyInterval>> phylogenyData =
                new HashMap<Integer, List<PhylogenyInterval>>(
                        sortedChromosomeNumbers.size());
            for(Integer currChromoNumber: sortedChromosomeNumbers)
            {
                ChromosomeDataSource currChromoDataSource =
                    chromoDataSources.get(currChromoNumber);
                List<IndexedSnpInterval> indexedMaxKIntervals =
                    INTERVAL_SCANNER.maxKScan(
                            currChromoDataSource.getSdpInputStream(sortedStrains),
                            currChromoDataSource.getSdpInputStream(StreamDirection.REVERSE, sortedStrains),
                            currChromoDataSource.getSdpInputStream(sortedStrains));
                int intervalCount = indexedMaxKIntervals.size();
                List<PhylogenyTreeNode> phyloTrees =
                    PHYLOGENY_SCANNER.inferPerfectPhylogenies(
                            currChromoDataSource.getSdpInputStream(sortedStrains),
                            indexedMaxKIntervals);
                List<BasePairInterval> maxKIntervals = INTERVAL_SCANNER.toOrderedPhysicalIntervals(
                        indexedMaxKIntervals,
                        currChromoDataSource.getSnpPositionInputStream());
                List<PhylogenyInterval> phyloIntervals =
                    new ArrayList<PhylogenyInterval>(maxKIntervals.size());
                for(int i = 0; i < intervalCount; i++)
                {
                    phyloIntervals.add(new PhylogenyInterval(
                            phyloTrees.get(i),
                            maxKIntervals.get(i)));
                }
                
                System.out.println("phylo intervals: " + phyloIntervals.size());
                phylogenyData.put(currChromoNumber, phyloIntervals);
            }
            
            return phylogenyData;
        }
        catch(IOException ex)
        {
            LOG.log(Level.SEVERE,
                    "Failed to build phylogeny intervals",
                    ex);
            return null;
        }
        catch(NoValidPhylogenyException ex)
        {
            LOG.log(Level.SEVERE,
                    "Failed to infer phylogeny from intervals",
                    ex);
            return null;
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
