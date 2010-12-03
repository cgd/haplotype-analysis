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

package org.jax.haplotype.analysis.jaxbfactory;

import java.io.IOException;

import org.jax.haplotype.analysis.HaplotypeAssociationTest;
import org.jax.haplotype.analysis.HaplotypeDataSource;
import org.jax.haplotype.analysis.PhenotypeDataSource;
import org.jax.haplotype.jaxbgenerated.HaplotypeAssociationTestType;
import org.jax.haplotype.jaxbgenerated.HaplotypeDataSourceType;
import org.jax.haplotype.jaxbgenerated.PhenotypeDataSourceType;

/**
 * Class responsible for generating haplotype association tests from the
 * given JAXB description.
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class JaxbHaplotypeAssociationTestFactory
{
    /**
     * Get the haplotype association test for the given JAX-B
     * descriptions
     * @param jaxbHaplotypeAssociationTest
     *          the JAX-B description of the haplotype association test
     * @return
     *          the matching test
     * @throws IOException
     *          if we get an exception reading/writing data
     */
    public static HaplotypeAssociationTest getHaplotypeAssociationTest(
            HaplotypeAssociationTestType jaxbHaplotypeAssociationTest)
            throws IOException
    {
        HaplotypeDataSourceType jaxbHaplotypeDataSource =
            (HaplotypeDataSourceType)jaxbHaplotypeAssociationTest.getHaplotypeDataSourceId();
        PhenotypeDataSourceType jaxbPhenotypeDataSource =
            (PhenotypeDataSourceType)jaxbHaplotypeAssociationTest.getPhenotypeDataSourceId();
        
        HaplotypeDataSource haplotypeDataSource =
            JaxbHaplotypeDataSourceFactory.getHaplotypeDataSource(
                    jaxbHaplotypeDataSource);
        PhenotypeDataSource phenotypeDataSource =
            JaxbPhenotypeDataSourceFactory.getPhenotypeDataSource(
                    jaxbPhenotypeDataSource);
        
        HaplotypeAssociationTest haplotypeAssociationTest =
            new HaplotypeAssociationTest(
                    jaxbHaplotypeAssociationTest.getId(),
                    haplotypeDataSource,
                    phenotypeDataSource);
        
        return haplotypeAssociationTest;
    }
}
