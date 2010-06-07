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
