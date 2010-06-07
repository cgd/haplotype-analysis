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

import java.io.File;
import java.io.IOException;

import org.jax.haplotype.analysis.experimentdesign.MPDIndividualPhenotypeDataSource;
import org.jax.haplotype.analysis.experimentdesign.PhenotypeDataSource;
import org.jax.haplotype.analysis.experimentdesign.SexFilter;
import org.jax.haplotype.jaxbgenerated.MpdIndividualPhenotypeDataSourceType;
import org.jax.haplotype.jaxbgenerated.PhenotypeDataSourceType;
import org.jax.haplotype.jaxbgenerated.SexConstraintType;

/**
 * A factory class for getting a "native" phenotype data source out of a
 * JAXB type
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class JaxbPhenotypeDataSourceFactory
{
    /**
     * Construct a "native" phenotype data source from the given JAXB type
     * @param jaxbPhenotypeDataSource
     *          the JAXB definition for the phenotype data source
     * @return
     *          the "native" phenotype data source
     * @throws IOException
     *          if we catch an IO exception while constructing the data source
     */
    public static PhenotypeDataSource getPhenotypeDataSource(
            PhenotypeDataSourceType jaxbPhenotypeDataSource) throws IOException
    {
        if(jaxbPhenotypeDataSource instanceof MpdIndividualPhenotypeDataSourceType)
        {
            MpdIndividualPhenotypeDataSourceType commaSeperatedDataSource =
                (MpdIndividualPhenotypeDataSourceType)jaxbPhenotypeDataSource;
            SexFilter sexFilter = jaxbSexConstraintToNativeSexFilter(
                    jaxbPhenotypeDataSource.getSexConstraint());
            
            return new MPDIndividualPhenotypeDataSource(
                    commaSeperatedDataSource.getName(),
                    new File(commaSeperatedDataSource.getFileLocation()),
                    commaSeperatedDataSource.getPhenotype(),
                    null,
                    sexFilter);
        }
        else
        {
            throw new RuntimeException(
                    "unknown phenotype data source type: ");
        }
    }
    
    /**
     * A convenience function for going between the jaxb type and native
     * java type for sex constrainst
     * @param jaxbSexConstraint
     *          the jaxb type
     * @return
     *          the java type
     */
    public static SexFilter jaxbSexConstraintToNativeSexFilter(
            SexConstraintType jaxbSexConstraint)
    {
        switch(jaxbSexConstraint)
        {
            case ALLOW_MALE: return SexFilter.ALLOW_MALE;
            
            case ALLOW_FEMALE: return SexFilter.ALLOW_FEMALE;
            
            case SEX_AGNOSTIC: return SexFilter.AGNOSTIC;
            
            default: return null;
        }
    }
}
