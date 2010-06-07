/*
 * Copyright (c) 2009 The Jackson Laboratory
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
