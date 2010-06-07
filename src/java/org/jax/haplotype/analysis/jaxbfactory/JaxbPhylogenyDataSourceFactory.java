/*
 * Copyright (c) 2008 The Jackson Laboratory
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package org.jax.haplotype.analysis.jaxbfactory;

import org.jax.haplotype.analysis.experimentdesign.GenotypeInferencePhylogenyDataSource;
import org.jax.haplotype.analysis.experimentdesign.PhylogenyDataSource;
import org.jax.haplotype.data.GenomeDataSource;
import org.jax.haplotype.data.JaxbGenomeDataSourceFactory;
import org.jax.haplotype.jaxbgenerated.GenomeDataSourceType;
import org.jax.haplotype.jaxbgenerated.GenotypeInferencePhylogenyDataSourceType;
import org.jax.haplotype.jaxbgenerated.PhylogenyDataSourceType;

/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class JaxbPhylogenyDataSourceFactory
{
    /**
     * Build a phylogeny data source from the given JAXB definition
     * @param phylogenyDataSourceElement
     *          the JAXB XML element
     * @return
     *          the real data source object
     */
    public static PhylogenyDataSource getPhylogenyDataSource(PhylogenyDataSourceType phylogenyDataSourceElement)
    {
        if(phylogenyDataSourceElement instanceof GenotypeInferencePhylogenyDataSourceType)
        {
            return getPhylogenyDataSource(
                (GenotypeInferencePhylogenyDataSourceType)phylogenyDataSourceElement);
        }
        else
        {
            throw new IllegalArgumentException(
                    "don't know how to create data source for: " +
                    phylogenyDataSourceElement.getClass().getName());
        }
    }
    
    /**
     * Build a genotype inference phylogeny data source object
     * @param genoInferencePhyloDataSourceElement
     *          the JAXB element that acts as a definition for the object
     *          that this function builds
     * @return
     *          the phylogeny data source
     */
    public static GenotypeInferencePhylogenyDataSource getPhylogenyDataSource(
            GenotypeInferencePhylogenyDataSourceType genoInferencePhyloDataSourceElement)
    {
        GenomeDataSourceType genomeDataSourceElement = (GenomeDataSourceType)genoInferencePhyloDataSourceElement.getGenomeDataSourceId();
        GenomeDataSource genomeDataSource = JaxbGenomeDataSourceFactory.getGenomeDataSource(
                genomeDataSourceElement);
        return new GenotypeInferencePhylogenyDataSource(
                genoInferencePhyloDataSourceElement.getName(),
                genomeDataSource);
    }
}
