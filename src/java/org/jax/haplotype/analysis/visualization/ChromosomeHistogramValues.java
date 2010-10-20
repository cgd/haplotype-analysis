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

package org.jax.haplotype.analysis.visualization;

import java.util.List;

import org.jax.geneticutil.data.RealValuedBasePairInterval;

/**
 * Holds the values needed to describe a chromosome histogram
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class ChromosomeHistogramValues
{
    private final int chromosomeNumber;
    
    private final List<? extends RealValuedBasePairInterval> intervals;
    
    private final long startInBasePairs;
    
    private final long extentInBasePairs;
    
    private final HighlightedSnpInterval visualInterval;

    /**
     * Constructor
     * @param chromosomeNumber
     *          the chromosome number
     * @param intervals
     *          the intervals to use
     * @param startInBasePairs
     *          where should we start the graph?
     * @param extentInBasePairs
     *          how far should the graph extend
     * @param visualInterval
     *          the visual interval to use
     */
    public ChromosomeHistogramValues(
            final int chromosomeNumber,
            final List<? extends RealValuedBasePairInterval> intervals,
            final long startInBasePairs,
            final long extentInBasePairs,
            final HighlightedSnpInterval visualInterval)
    {
        this.chromosomeNumber = chromosomeNumber;
        this.intervals = intervals;
        this.startInBasePairs = startInBasePairs;
        this.extentInBasePairs = extentInBasePairs;
        this.visualInterval = visualInterval;
    }
    
    /**
     * Getter for the chromosome number for this histogram
     * @return the chromosome number
     */
    public int getChromosomeNumber()
    {
        return this.chromosomeNumber;
    }
    
    /**
     * Getter for the real-valued intervals for this histogram
     * @return the intervals
     */
    public List<? extends RealValuedBasePairInterval> getIntervals()
    {
        return this.intervals;
    }
    
    /**
     * Getter for the starting position in base pairs
     * @return the start in base pairs
     */
    public long getStartInBasePairs()
    {
        return this.startInBasePairs;
    }
    
    /**
     * Getter for the extent in base pairs
     * @return the extent in base pairs
     */
    public long getExtentInBasePairs()
    {
        return this.extentInBasePairs;
    }
    
    /**
     * Getter for the visual interval
     * @return the visual interval
     */
    public HighlightedSnpInterval getVisualInterval()
    {
        return this.visualInterval;
    }
}
