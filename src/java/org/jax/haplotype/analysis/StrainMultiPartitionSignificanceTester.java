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

package org.jax.haplotype.analysis;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math.MathException;
import org.apache.commons.math.stat.inference.OneWayAnova;
import org.apache.commons.math.stat.inference.OneWayAnovaImpl;
import org.jax.geneticutil.data.MultiGroupStrainPartition;
import org.jax.util.datastructure.SequenceUtilities;
import org.jax.util.math.StatisticUtilities;

/**
 * Performs an F-test on multiple partitions
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class StrainMultiPartitionSignificanceTester
{
    private final OneWayAnova oneWayAnova = new OneWayAnovaImpl();
    
    /**
     * Test the significance of the given responses
     * @param strainPartitions
     *          the partitions to test
     * @param strainMultipleResponses
     *          the responses... strain names are keys and responses
     *          are values
     * @return
     *          test p-values. this will be an array as long as the
     *          input partitions
     */
    public double[] fTestMultipleResponseSignificance(
            List<? extends MultiGroupStrainPartition> strainPartitions,
            double[][] strainMultipleResponses)
    {
        // for right now just average the responses. Later on we can
        // do something smarter
        double[] strainResponses = new double[strainMultipleResponses.length];
        for(int i = 0; i < strainMultipleResponses.length; i++)
        {
            strainResponses[i] = StatisticUtilities.calculateMean(
                    strainMultipleResponses[i]);
        }
        
        return this.fTestSingleResponseSignificance(
                strainPartitions,
                strainResponses);
    }
    
    /**
     * Test the significance of the given responses
     * @param strainPartitions
     *          the partitions to test
     * @param strainResponses
     *          the responses
     * @return
     *          test p-values. this will be an array as long as the
     *          input partitions
     */
    public double[] fTestSingleResponseSignificance(
            List<? extends MultiGroupStrainPartition> strainPartitions,
            double[] strainResponses)
    {
        double[] significanceValues =
            new double[strainPartitions.size()];
        for(int currPartitionIndex = 0;
            currPartitionIndex < strainPartitions.size();
            currPartitionIndex++)
        {
            MultiGroupStrainPartition currPartition =
                strainPartitions.get(currPartitionIndex);
            
            // segregate the responses
            short[] strainGroups = currPartition.getStrainGroups();
            if(strainGroups.length != strainResponses.length)
            {
                throw new IllegalArgumentException(
                        "Expect the strain group array length to match the strain " +
                        "response count but " + strainGroups.length + " and " +
                        strainResponses.length + " do not match");
            }
            
            Map<Short, List<Double>> segregatedResponseLists =
                new HashMap<Short, List<Double>>();
            for(int i = 0; i < strainGroups.length; i++)
            {
                List<Double> groupResponses =
                    segregatedResponseLists.get(strainGroups[i]);
                if(groupResponses == null)
                {
                    groupResponses = new ArrayList<Double>();
                    segregatedResponseLists.put(strainGroups[i], groupResponses);
                }
                
                groupResponses.add(strainResponses[i]);
            }
            
            List<double[]> segregatedResponses = new ArrayList<double[]>(
                    segregatedResponseLists.size());
            for(List<Double> groupResponses: segregatedResponseLists.values())
            {
                if(groupResponses.size() >= 2)
                {
                    segregatedResponses.add(SequenceUtilities.toDoubleArray(
                            groupResponses));
                }
            }
            
            // perform F-test on segregated responses
            if(segregatedResponses.size() < 2)
            {
                significanceValues[currPartitionIndex] = 1.0;
            }
            else
            {
                try
                {
                    significanceValues[currPartitionIndex] =
                        this.oneWayAnova.anovaPValue(segregatedResponses);
                }
                catch(MathException ex)
                {
                    throw new IllegalStateException(ex);
                }
            }
        }
        
        return significanceValues;
    }
}
