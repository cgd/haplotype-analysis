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

package org.jax.haplotype.analysis;

import java.util.List;

import org.jax.geneticutil.data.BasePairInterval;
import org.jax.geneticutil.data.SimpleBasePairInterval;

/**
 * A class for calculating overlap between SNP block lists.
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class SnpBlockOverlapCalculator
{
    private final long calculationLimit;

    /**
     * Constructor for applying no effective constraint on the input size.
     */
    public SnpBlockOverlapCalculator()
    {
        this(Long.MAX_VALUE);
    }
    
    /**
     * Constructor
     * @param calculationLimit
     *          this is the maximum number of calculations we allow where
     *          the number of calculations is
     *          <code>snpBlocksToTest.size() * snpBlocksToTestAgainst.size()</code>
     */
    public SnpBlockOverlapCalculator(final long calculationLimit)
    {
        this.calculationLimit = calculationLimit;
    }
    
    /**
     * Calculate overlap for the given lists of SNP blocks
     * TODO This function should run faster than the current O(m*n) algorithm
     * @param snpBlocksToTest
     *          the blocks that we're testing for overlap.
     *          These should already be sorted with
     *          {@link SimpleBasePairInterval#SNP_INTERVAL_COMPARATOR}
     * @param snpBlocksToTestAgainst
     *          the blocks that we're testing against
     *          These should already be sorted with
     *          {@link SimpleBasePairInterval#SNP_INTERVAL_COMPARATOR}
     * @return
     *          the 0-1 cumulative overlap of the blocks to test vs. the
     *          blocks to test against. this array will always be the
     *          same size as the blocks to test
     */
    public double[] calculateSnpBlockOverlap(
            List<BasePairInterval> snpBlocksToTest,
            List<BasePairInterval> snpBlocksToTestAgainst)
    {
        // calculation limit check
        {
            final int numCalculations =
                snpBlocksToTest.size() * snpBlocksToTestAgainst.size();
            if(numCalculations > this.calculationLimit)
            {
                throw new IllegalArgumentException(
                        "The input data size exceeds the current limit " +
                        "(# rows in file 1) * (# rows in file 2) must be " +
                        "less than or equal to " + this.calculationLimit +
                        " but the current value of " + numCalculations +
                        " is exceeds this limit.");
            }
        }
        
        double[] overlapResults = new double[snpBlocksToTest.size()];
        
        // loop through the blocks to test
        for(int blockToTestIndex = 0;
            blockToTestIndex < overlapResults.length;
            blockToTestIndex++)
        {
            // initialize some data
            double currentOverlap = 0.0;
            BasePairInterval blockToTest =
                snpBlocksToTest.get(blockToTestIndex);
            long blockToTestExtent = blockToTest.getExtentInBasePairs();
            
            if(blockToTestExtent > 0)
            {
                for(BasePairInterval blockToTestAgainst: snpBlocksToTestAgainst)
                {
                    long overlap = blockToTest.getOverlapInBasePairs(
                            blockToTestAgainst);
                    if(overlap >= 1)
                    {
                        currentOverlap +=
                            overlap / (double)blockToTestExtent;
                    }
                }
            }
            
            // set the accumulated overlap for the block to test
            overlapResults[blockToTestIndex] = currentOverlap;
        }
        
        return overlapResults;
    }
}
