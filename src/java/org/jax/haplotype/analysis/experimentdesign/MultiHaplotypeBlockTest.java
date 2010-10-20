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

import java.io.Serializable;
import java.util.Set;

/**
 * An interface for representing a significance test of the given haplotype
 * data with the given phenotype data
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public interface MultiHaplotypeBlockTest extends Serializable
{
    /**
     * Getter for the name of this test
     * @return the name
     */
    public String getName();
    
    /**
     * The phenotype data source
     * @return the phenotypeDataSource
     */
    public PhenotypeDataSource getPhenotypeDataSource();
    
    /**
     * get the strain names that are common between the composite data sources
     * @return
     *          the common strain set
     */
    public Set<String> getCommonStrains();
    
    /**
     * Perform the significance tests
     * @param chromosomeNumber
     *          the chromosome number to get results for
     * @return
     *          the results of the test
     */
    public MultiHaplotypeBlockTestResult[] getTestResults(
            int chromosomeNumber);

    /**
     * Getter for the chromosomes available in this test
     * @return
     *          the chromosome numbers
     */
    public int[] getAvailableChromosomes();
}
