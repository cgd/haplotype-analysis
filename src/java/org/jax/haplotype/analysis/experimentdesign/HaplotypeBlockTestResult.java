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

import java.util.BitSet;

import org.jax.geneticutil.data.BinaryStrainPartition;
import org.jax.geneticutil.data.CompositeRealValuedBasePairInterval;
import org.jax.geneticutil.data.PartitionedInterval;

/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class HaplotypeBlockTestResult
extends CompositeRealValuedBasePairInterval
implements BinaryStrainPartition
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
    public HaplotypeBlockTestResult(PartitionedInterval interval, double pValue)
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
    public PartitionedInterval getDelegateInterval()
    {
        return (PartitionedInterval)super.getDelegateInterval();
    }
    
    /**
     * {@inheritDoc}
     */
    public BitSet getStrainBitSet()
    {
        return this.getDelegateInterval().getStrainBitSet();
    }
}
