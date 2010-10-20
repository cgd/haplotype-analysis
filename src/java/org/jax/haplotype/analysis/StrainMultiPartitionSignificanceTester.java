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
