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

import java.util.BitSet;
import java.util.List;

import org.apache.commons.math.MathException;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math.stat.inference.TTest;
import org.apache.commons.math.stat.inference.TTestImpl;
import org.jax.geneticutil.data.BinaryStrainPartition;
import org.jax.geneticutil.data.PartitionedIntervalSet;
import org.jax.util.math.StatisticUtilities;

/**
 * Class that tests a bunch of strain partitions for their significance
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class StrainBinaryPartitionSignificanceTester
{
    /**
     * T-Test for significance
     */
    private final TTest tTester = new TTestImpl();
    
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
    public double[] tTestMultipleResponseSignificance(
            List<? extends BinaryStrainPartition> strainPartitions,
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
        
        return this.tTestSingleResponseSignificance(
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
    public double[] tTestSingleResponseSignificance(
            List<? extends BinaryStrainPartition> strainPartitions,
            double[] strainResponses)
    {
        double[] significanceValues =
            new double[strainPartitions.size()];
        for(int currPartitionIndex = 0;
            currPartitionIndex < strainPartitions.size();
            currPartitionIndex++)
        {
            BinaryStrainPartition currPartition =
                strainPartitions.get(currPartitionIndex);
            
            // segregate the responses
            BitSet currStrainBitSet = currPartition.getStrainBitSet();
            int currStrainCount = currStrainBitSet.cardinality();
            double[] insidePartitionResponses =
                new double[currStrainCount];
            double[] outsidePartitionResponses =
                new double[strainResponses.length - currStrainCount];
            
            int currInsidePartitionCursor = 0;
            int currOutsidePartitionCursor = 0;
            for(int responseIndex = 0;
                responseIndex < strainResponses.length;
                responseIndex++)
            {
                if(currStrainBitSet.get(responseIndex))
                {
                    // this response is in the partition
                    insidePartitionResponses[currInsidePartitionCursor] =
                        strainResponses[responseIndex];
                    currInsidePartitionCursor++;
                }
                else
                {
                    // this response is outside the partition
                    outsidePartitionResponses[currOutsidePartitionCursor] =
                        strainResponses[responseIndex];
                    currOutsidePartitionCursor++;
                }
            }
            
            assert currInsidePartitionCursor == insidePartitionResponses.length;
            assert currOutsidePartitionCursor == outsidePartitionResponses.length;
            
            // perform t-test on segregated responses
            if(insidePartitionResponses.length <= 2 ||
               outsidePartitionResponses.length <= 2)
            {
                significanceValues[currPartitionIndex] = 1.0;
            }
            else
            {
                DescriptiveStatistics insidePartitionResponseSummary =
                    new DescriptiveStatistics();
                for(double currInsideResponseValue: insidePartitionResponses)
                {
                    insidePartitionResponseSummary.addValue(currInsideResponseValue);
                }
                
                DescriptiveStatistics outsidePartitionResponseSummary =
                    new DescriptiveStatistics();
                for(double currOutsideResponseValue: outsidePartitionResponses)
                {
                    outsidePartitionResponseSummary.addValue(currOutsideResponseValue);
                }
                
                try
                {
                    double pValue = this.tTester.tTest(
                            insidePartitionResponseSummary,
                            outsidePartitionResponseSummary);
                    
                    significanceValues[currPartitionIndex] = pValue;
                }
                catch(MathException ex)
                {
                    throw new IllegalStateException(ex);
                }
            }
        }
        
        return significanceValues;
    }
    
    /**
     * Test the significance of the given responses
     * @param genomicPartitions
     *          the partitions to test
     * @param strainMultipleResponses
     *          multiple responses. major index is for strains minor index is
     *          for response
     * @return
     *          significance values. this will be an array as long as the
     *          input partitions
     */
    public double[] normalizedTestMultipleResponseSignificance(
            List<? extends PartitionedIntervalSet> genomicPartitions,
            double[][] strainMultipleResponses)
    {
        double[] strainResponses = new double[strainMultipleResponses.length];
        for(int i = 0; i < strainMultipleResponses.length; i++)
        {
            strainResponses[i] = StatisticUtilities.calculateMean(
                    strainMultipleResponses[i]);
        }
        
        return this.normalizedTestSingleResponseSignificance(
                genomicPartitions,
                strainResponses);
    }
    
    /**
     * Test the significance of the given responses
     * @param genomicPartitions
     *          the partitions to test
     * @param strainResponses
     *          the responses
     * @return
     *          significance values. this will be an array as long as the
     *          input partitions
     */
    public double[] normalizedTestSingleResponseSignificance(
            List<? extends PartitionedIntervalSet> genomicPartitions,
            double[] strainResponses)
    {
        double[] sinificanceValues =
            new double[genomicPartitions.size()];
        for(int currPartitionIndex = 0;
            currPartitionIndex < genomicPartitions.size();
            currPartitionIndex++)
        {
            BinaryStrainPartition currPartition =
                genomicPartitions.get(currPartitionIndex);
            
            // segregate the responses
            BitSet currPartitionStrainBitSet = currPartition.getStrainBitSet();
            int currPartitionChromoCount = currPartitionStrainBitSet.cardinality();
            double[] insidePartitionResponses =
                new double[currPartitionChromoCount];
            double[] outsidePartitionResponses =
                new double[strainResponses.length - currPartitionChromoCount];
            
            int currInsidePartitionCursor = 0;
            int currOutsidePartitionCursor = 0;
            for(int responseIndex = 0;
                responseIndex < strainResponses.length;
                responseIndex++)
            {
                if(currPartitionStrainBitSet.get(responseIndex))
                {
                    // this response is in the partition
                    insidePartitionResponses[currInsidePartitionCursor] =
                        strainResponses[responseIndex];
                    currInsidePartitionCursor++;
                }
                else
                {
                    // this response is outside the partition
                    outsidePartitionResponses[currOutsidePartitionCursor] =
                        strainResponses[responseIndex];
                    currOutsidePartitionCursor++;
                }
            }
            
            assert currInsidePartitionCursor == insidePartitionResponses.length;
            assert currOutsidePartitionCursor == outsidePartitionResponses.length;
            
            // perform t-test on segregated responses
            if(insidePartitionResponses.length <= 2 ||
               outsidePartitionResponses.length <= 2)
            {
                sinificanceValues[currPartitionIndex] = 1.0;
            }
            else
            {
                DescriptiveStatistics insidePartitionResponseSummary =
                    new DescriptiveStatistics();
                for(double currInsideResponseValue: insidePartitionResponses)
                {
                    insidePartitionResponseSummary.addValue(currInsideResponseValue);
                }
                
                DescriptiveStatistics outsidePartitionResponseSummary =
                    new DescriptiveStatistics();
                for(double currOutsideResponseValue: outsidePartitionResponses)
                {
                    outsidePartitionResponseSummary.addValue(currOutsideResponseValue);
                }
                
                try
                {
                    double pValue = this.tTester.tTest(
                            insidePartitionResponseSummary,
                            outsidePartitionResponseSummary);
                    
                    // reduce the pValue relative to cumulative extent
                    pValue /=
                        genomicPartitions.get(currPartitionIndex).getCumulativeExtentInBasePairs();
                    
                    sinificanceValues[currPartitionIndex] = pValue;
                }
                catch(MathException ex)
                {
                    throw new IllegalStateException(ex);
                }
            }
        }
        
        return sinificanceValues;
    }
}
