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
