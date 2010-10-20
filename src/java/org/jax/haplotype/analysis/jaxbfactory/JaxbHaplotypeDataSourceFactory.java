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
import java.util.HashSet;
import java.util.Set;

import org.jax.haplotype.analysis.experimentdesign.GenotypeInferenceHaplotypeDataSource;
import org.jax.haplotype.analysis.experimentdesign.HaplotypeDataSource;
import org.jax.haplotype.data.GenomeDataSource;
import org.jax.haplotype.data.JaxbGenomeDataSourceFactory;
import org.jax.haplotype.jaxbgenerated.GenomeDataSourceType;
import org.jax.haplotype.jaxbgenerated.GenotypeInferenceHaplotypeDataSourceType;
import org.jax.haplotype.jaxbgenerated.HaplotypeDataSourceType;
import org.jax.haplotype.jaxbgenerated.StrainType;

/**
 * A factory that takes in JAXB definitions for a haplotype data source
 * and returns the actual data source
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class JaxbHaplotypeDataSourceFactory
{
    /**
     * Get the haplotype data source for the given JAXB parameter
     * @param jaxbHaplotypeDataSource
     *          the JAXB definition of the haplotype data source
     * @return
     *          the actual haplotype data source
     * @throws IOException
     */
    public static HaplotypeDataSource getHaplotypeDataSource(
            HaplotypeDataSourceType jaxbHaplotypeDataSource)
    throws IOException
    {
        if(jaxbHaplotypeDataSource instanceof GenotypeInferenceHaplotypeDataSourceType)
        {
            return getHaplotypeDataSource(
                    (GenotypeInferenceHaplotypeDataSourceType)jaxbHaplotypeDataSource);
        }
        else
        {
            throw new RuntimeException(
                    "Unknown haplotype data source type: " +
                    jaxbHaplotypeDataSource.getClass().getName());
        }
    }
    
    /**
     * Get the haplotype data source for the given JAXB parameter
     * @param jaxbHaplotypeDataSource
     *          the JAXB definition of the haplotype data source
     * @return
     *          the actual haplotype data source
     * @throws IOException
     */
    public static HaplotypeDataSource getHaplotypeDataSource(
            GenotypeInferenceHaplotypeDataSourceType jaxbHaplotypeDataSource)
    throws IOException
    {
        Set<String> strainNamesToAccept = null;
        if(!jaxbHaplotypeDataSource.getStrainToAcceptFilter().isEmpty())
        {
            strainNamesToAccept = new HashSet<String>(
                    jaxbHaplotypeDataSource.getStrainToAcceptFilter().size());
            for(StrainType strainType: jaxbHaplotypeDataSource.getStrainToAcceptFilter())
            {
                strainNamesToAccept.add(strainType.getStrainName());
            }
        }
        
        GenomeDataSourceType genomeDataSourceElement =
            (GenomeDataSourceType)jaxbHaplotypeDataSource.getGenomeDataSourceId();
        GenomeDataSource genomeDataSource = JaxbGenomeDataSourceFactory.getGenomeDataSource(
                genomeDataSourceElement);
        
        return new GenotypeInferenceHaplotypeDataSource(
                jaxbHaplotypeDataSource.getName(),
                genomeDataSource,
                strainNamesToAccept,
                (int)jaxbHaplotypeDataSource.getMinimumSnpExtent(),
                (int)jaxbHaplotypeDataSource.getMinimumStrainGroupSize());
    }
}
