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

import org.jax.geneticutil.data.PartitionedIntervalSet;

/**
 * A test result
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class HaplotypeEquivalenceClassTestResult
{
    private final PartitionedIntervalSet haplotypeEquivalenceClass;
    
    private final double pValue;

    /**
     * Constructor
     * @param haplotypeEquivalenceClass
     *          the equivalence class that was tested
     * @param pValue
     *          the p-value returned by the test
     */
    public HaplotypeEquivalenceClassTestResult(
            PartitionedIntervalSet haplotypeEquivalenceClass,
            double pValue)
    {
        this.haplotypeEquivalenceClass = haplotypeEquivalenceClass;
        this.pValue = pValue;
    }
    
    /**
     * Getter for the equivalence class that was tested
     * @return
     *          the haplotype equivalence class
     */
    public PartitionedIntervalSet getHaplotypeEquivalenceClass()
    {
        return this.haplotypeEquivalenceClass;
    }
    
    /**
     * Getter for the p-value of the test
     * @return
     *          the p-value
     */
    public double getPValue()
    {
        return this.pValue;
    }
}
