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

import org.jax.haplotype.analysis.experimentdesign.PhenotypeDataSource;
import org.jax.haplotype.analysis.experimentdesign.PhylogenyAssociationTest;
import org.jax.haplotype.analysis.experimentdesign.PhylogenyDataSource;
import org.jax.haplotype.jaxbgenerated.PhenotypeDataSourceType;
import org.jax.haplotype.jaxbgenerated.PhylogenyAssociationTestType;
import org.jax.haplotype.jaxbgenerated.PhylogenyDataSourceType;

/**
 * A factory for creating phylogeny association tests
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class JaxbPhylogenyAssociationTestFactory
{
    /**
     * Factory method for creating a java association test from the given
     * XML definition
     * @param phylogenyAssociationTestElement
     *          the XML definition for the test
     * @return
     *          the test
     * @throws IOException
     *          if we run into trouble doing I/O
     */
    public static PhylogenyAssociationTest getPhylogenyAssociationTest(
            PhylogenyAssociationTestType phylogenyAssociationTestElement)
            throws IOException
    {
        PhylogenyDataSourceType phylogenyDataSourceElement =
            (PhylogenyDataSourceType)phylogenyAssociationTestElement.getPhylogenyDataSourceId();
        PhylogenyDataSource phylogenyDataSource =
            JaxbPhylogenyDataSourceFactory.getPhylogenyDataSource(
                    phylogenyDataSourceElement);
        PhenotypeDataSourceType phenoDataSourceElement =
            (PhenotypeDataSourceType)phylogenyAssociationTestElement.getPhenotypeDataSourceId();
        PhenotypeDataSource phenotypeDataSource =
            JaxbPhenotypeDataSourceFactory.getPhenotypeDataSource(
                    phenoDataSourceElement);
        
        return new PhylogenyAssociationTest(
                phylogenyAssociationTestElement.getId(),
                phylogenyDataSource,
                phenotypeDataSource);
    }
}
