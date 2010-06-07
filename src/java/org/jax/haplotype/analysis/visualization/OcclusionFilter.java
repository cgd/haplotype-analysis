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

package org.jax.haplotype.analysis.visualization;

import java.util.ArrayList;
import java.util.List;

import org.jax.geneticutil.data.RealValuedBasePairInterval;

/**
 * A class for performing occlusion filtering on real valued intervals
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class OcclusionFilter
{
    /**
     * Removes occluded intervals where we define occluded as an interval that
     * is fully contained by another interval whose value is greater
     * @param <T>
     *          the type
     * @see RealValuedBasePairInterval#getRealValue()
     * @see RealValuedBasePairInterval#contains(org.jax.geneticutil.data.BasePairInterval)
     * @param realValuedIntervals
     *          the values to filter. this list must be presorted using
     *          {@link RealValuedBasePairInterval}'s comparison method
     * @param invertOcclusionComparison
     *          setting this to true basically acts as if you had
     *          changed every {@link RealValuedBasePairInterval#getRealValue()}
     *          to 1/{@link RealValuedBasePairInterval#getRealValue()}
     * @return
     *          the list sans any occluded intervals
     */
    public static <T extends RealValuedBasePairInterval> List<T> filterOutOccludedIntervals(
            List<T> realValuedIntervals,
            boolean invertOcclusionComparison)
    {
        boolean[] occludedBools = new boolean[realValuedIntervals.size()];
        int unoccludedCount = occludedBools.length;
        for(int i = 0; i < occludedBools.length; i++)
        {
            RealValuedBasePairInterval currInterval =
                realValuedIntervals.get(i);
            
            // look behind you for occluded intervals. there can only be
            // occlusion behind you if the starting point is the same because
            // we are assuming that the input is presorted
            for(int j = i - 1; j >= 0; j--)
            {
                RealValuedBasePairInterval behindMeInternal =
                    realValuedIntervals.get(j);
                if(behindMeInternal.getStartInBasePairs() == currInterval.getStartInBasePairs())
                {
                    if(!occludedBools[j] &&
                       occludes(currInterval, behindMeInternal, invertOcclusionComparison))
                    {
                        // the interval is occluded by this interval
                        occludedBools[j] = true;
                        unoccludedCount--;
                    }
                }
                else
                {
                    // don't need to look behind me anymore
                    break;
                }
            }
            
            // look in front of you
            for(int j = i + 1; j < occludedBools.length; j++)
            {
                RealValuedBasePairInterval inFrontOfMeInternal =
                    realValuedIntervals.get(j);
                if(currInterval.intersects(inFrontOfMeInternal))
                {
                    if(!occludedBools[j] &&
                       occludes(currInterval, inFrontOfMeInternal, invertOcclusionComparison))
                    {
                        occludedBools[j] = true;
                        unoccludedCount--;
                    }
                }
                else
                {
                    // we're definately not going to occlude anything else
                    // in front
                    break;
                }
            }
        }
        
        List<T> unoccludedIntervals = new ArrayList<T>(unoccludedCount);
        for(int j = 0; j < occludedBools.length; j++)
        {
            if(!occludedBools[j])
            {
                unoccludedIntervals.add(realValuedIntervals.get(j));
            }
        }
        
        return unoccludedIntervals;
    }
    
    private static boolean occludes(
            RealValuedBasePairInterval interval1,
            RealValuedBasePairInterval interval2,
            boolean invertOcclusionComparison)
    {
        if(invertOcclusionComparison)
        {
            return interval1.getRealValue() < interval2.getRealValue() &&
                   interval1.contains(interval1);
        }
        else
        {
            return interval1.getRealValue() > interval2.getRealValue() &&
                   interval1.contains(interval1);
        }
    }
}
