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

package org.jax.haplotype.analysis.experimentdesign;

import org.jax.geneticutil.data.CompositeRealValuedBasePairInterval;
import org.jax.geneticutil.data.MultiGroupStrainPartition;
import org.jax.geneticutil.data.MultiPartitionedInterval;

/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class MultiHaplotypeBlockTestResult
extends CompositeRealValuedBasePairInterval
implements MultiGroupStrainPartition
{
    /**
     * every {@link java.io.Serializable} is supposed to have one of these
     */
    private static final long serialVersionUID = -360759604320672055L;

    /**
     * Constructor
     * @param interval
     *          the interval
     * @param pValue
     *          the p-value
     */
    public MultiHaplotypeBlockTestResult(
            MultiPartitionedInterval interval,
            double pValue)
    {
        super(interval, pValue);
    }
    
    /**
     * Just returns {@link #getRealValue()}
     * @return
     *          the p-value for this test
     */
    public double getPValue()
    {
        return this.getRealValue();
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public MultiPartitionedInterval getDelegateInterval()
    {
        return (MultiPartitionedInterval)super.getDelegateInterval();
    }
    
    /**
     * {@inheritDoc}
     */
    public short[] getStrainGroups()
    {
        return this.getDelegateInterval().getStrainGroups();
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public String toString()
    {
        return "p-value=" + this.getPValue() + ", " +
               "interval=" + this.getDelegateInterval();
    }
}
